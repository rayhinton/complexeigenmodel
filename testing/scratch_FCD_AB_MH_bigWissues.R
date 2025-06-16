# MH proposals for A, B

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
# source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_data-covar-dense.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

# Check if package is installed, install if not
if (!require("truncdist", quietly = TRUE)) {
    # Set up personal library
    user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
    if (!dir.exists(user_lib)) {
        dir.create(user_lib, recursive = TRUE)
    }
    .libPaths(c(user_lib, .libPaths()))
    
    # Install package
    install.packages("truncdist", 
                     repos = "https://cran.rstudio.com/",
                     dependencies = TRUE)
    library(truncdist)
}

library(foreach)
library(doParallel)

# proposal functions ------------------------------------------------------

ab_minmaxs <- function(x, width_scale = 1, extend_hilo = TRUE) { 
    P <- length(x)
    
    mins <- x[2:(P-1)] - (x[2:(P-1)] - x[3:P])/2 * width_scale
    maxs <- x[2:(P-1)] + (x[1:(P-2)] - x[2:(P-1)])/2 * width_scale
    
    if (extend_hilo) {
        # special extended ranges for the highest and lowest values
        maxs[1] <- x[2] + (x[1] - x[2]) * width_scale
        mins[P-2] <- x[P-1] - (x[P-1] - x[P]) * width_scale
    }
    
    return(list(mins = mins, maxs = maxs))
}

# parameters --------------------------------------------------------------

doPlot <- FALSE
extendRanges <- FALSE

Ukseed <- 21052025
# Ukseed <- 22052025

P <- 8
d <- 4
K <- 4
al_be_upper <- 2
al_be_lower <- 1
bing_its <- 5000
gs <- 200

S_its <- 200
S_print <- 10
MCits <- 1e5
# should the foreach loop output "inorder"? can potentially speed up
inorder <- TRUE
batch_count <- 10
inParallel <- TRUE

# true parameter values, as an override of the randomness
# w_0 <- 4

# w prior and proposal parameters
eta0_w <- 2
tau02_w <- 1/1000
wp_sd <- 1

# how should initialization be done? 
# - FALSE means the parameters are initialized to their true values
initialRandom <- TRUE

# should parameters be sampled or kept constant?
alphaConstant <- TRUE
betaConstant <- TRUE
wConstant <- FALSE

# FCD by MC integration functions ----------------------------------------

### functions to simulate from complex Stiefel manifold
# return n independent samples from the standard complex normal distribution
rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

# draw matrix distributed uniformly on V_d^P(F_betaf), where 
# - betaf = 1 --> Real
# - betaf = 2 --> Complex
runif_stiefel <- function(P, d, betaf) {
    # real
    if (betaf == 1) {
        X1 <- matrix(rnorm(P*d), ncol = d)
        # complex
    } else if (betaf == 2) {
        X1 <- matrix(rscnorm(P*d), ncol = d)        
    }
    
    # take QR decomposition - not guaranteed to be unique due to numerical methods
    X1qr <- qr(X1)
    QX1 <- qr.Q(X1qr)
    # extract sign of the diagonal elements of R
    D <- sign(Re(diag(qr.R(X1qr))))
    # transform by Q = QS, where S = diag(D). this transformation guarantees Q is unique.
    
    # Qfinal <- t(t(QX1) * D) # faster than matrix multiplication, for larger matrices
    Qfinal <- QX1 %*% diag(D)
    
    return(Qfinal)
}

# calculate LogSumExp
logsumexp <- function(x) {
    xstar <- max(x)
    lse <- xstar + log( sum(exp( x - xstar )) )
    return(lse)
}

# Approximate 0F0 by MC integration
hgf0F0MC <- function(x, y, betaf, MCits = 1e4, logscale = FALSE) {
    P <- length(x)
    d <- length(y)
    
    # compute the log integrand
    logintgd <- rep(NA, MCits)
    for (i in 1:MCits) {
        U <- runif_stiefel(P, d, betaf = betaf)
        logintgd[i] <- Re(crossprod(x, (Conj(U) * U) %*% y))
    }
    
    # return result on linear or log scales
    if (logscale == FALSE) {
        return(mean(exp(logintgd)))
    } else {
        # use LogSumExp trick for log scale approximation
        return(logsumexp(logintgd) - log(MCits))
    }
}

# generate parameters ------------------------------------------------------

# set.seed(8052025) 
set.seed(13052025)

# V matrix parameter
# V_0
V_0 <- runitary(P, P)

V_0

# w scalar
# w_0
# w_0 <- rgamma(1, 1, 1)
if (!exists("w_0")) {
    w_0 <- rgamma(1, eta0_w/2, rate = tau02_w)
} else {
    rgamma(1, eta0_w/2, rate = tau02_w)
}

# "previous" alpha, beta samples
alpha_0 <- c(al_be_upper,
             runif(P-2, al_be_lower, al_be_upper) |> sort(decreasing = TRUE),
             al_be_lower)
beta_0 <- c(al_be_upper,
            runif(d-2, al_be_lower, al_be_upper) |> sort(decreasing = TRUE),
            al_be_lower)
alpha_0
beta_0
w_0

A_0 <- diag(sqrt(w_0) * alpha_0)
B_0 <- diag(sqrt(w_0) * beta_0)
G_0 <- V_0 %*% A_0 %*% t(Conj(V_0))
diag(G_0) <- Re(diag(G_0))

G_0

# generate Uk data --------------------------------------------------------

# U_k matrix parameters
# U_k_0
U_ks <- array(NA, c(P, d, K))
set.seed(Ukseed)
Ukinit <- runitary(P, d)
for (i in 1:bing_its) {
    if (i %% 500 == 0) {
        print(paste0("i = ", i))
    }
    # for rectangular Bingham matrices
    Ukinit <- rcmb(Ukinit, G_0, B_0)
    # for square Bingham matrices
    # Ukinit <- rcBingUP_gibbs(Ukinit, G_0, B_0, Imtol = 1e4*.Machine$double.eps)
}
for (s in 1:(K*100)) {
    # for rectangular Bingham matrices
    Ukinit <- rcmb(Ukinit, G_0, B_0)
    # for square Bingham matrices
    # Ukinit <- rcBingUP_gibbs(Ukinit, G_0, B_0, Imtol = 1e4*.Machine$double.eps)
    if (s %% 100 == 0) {
        print(paste0("k = ", s/100))
        U_ks[, , s/100] <- Ukinit
    }
}

print(U_ks[, , K])

# load("966112.RData")

# set up alpha and beta samples array -------------------------------------

# this can be used inside to the function, to derive w based on the convention
# ws <- (As[1] / al_be_upper)^2
# alphas <- As / sqrt(ws)
# betas <- Bs / sqrt(ws)

# calculate matrix M, sum of function of V and U_k matrices (constant for all A, B)
M <- matrix(0, nrow = P, ncol = d)
for (k in 1:K) {
    M <- M + Re( Conj( t(Conj(V_0))%*%U_ks[, , k] ) * t(Conj(V_0))%*%U_ks[, , k] )
}

# initialize
alpha_S <- matrix(NA, P, S_its)
alpha_S[c(1, P), ] <- c(al_be_upper, al_be_lower)
beta_S <- matrix(NA, d, S_its)
beta_S[c(1, d), ] <- c(al_be_upper, al_be_lower)
w_S <- rep(NA, S_its)

# initialize the first sample
if (initialRandom) {
    set.seed(17052025)
    alpha_S[2:(P-1), 1] <- sort(runif(P-2, 1, 2), decreasing = TRUE)
    beta_S[2:(d-1), 1] <- sort(runif(d-2, 1, 2), decreasing = TRUE)
    w_S[1] <- rgamma(1, eta0_w/2, rate = tau02_w)
} else {
    alpha_S[, 1] <- alpha_0
    beta_S[, 1] <- beta_0
    w_S[1] <- w_0
}

alpha_S[, 1:5]
beta_S[, 1:5]
w_S[1]

# alpha FCD sampling ------------------------------------------------------

if (inParallel) {
    parallel::detectCores()
    
    cluster <- makeCluster(batch_count)
    registerDoParallel(cluster)
} else {
    `%dopar%` <- foreach::`%do%`
}

# one step
# s <- 2
alphas <- ifelse(rep(alphaConstant, P), alpha_0, alpha_S[, 1])
betas <- ifelse(rep(betaConstant, d), beta_0, beta_S[, 1])
ws <- ifelse(wConstant, w_0, w_S[1])
accCounts <- rep(0, 3)

(log0F0s <- hgf0F0MC(sqrt(ws) * alphas, sqrt(ws) * betas, 
                    betaf = 2, MCits = MCits, logscale = TRUE))
print(log0F0s)

# for (s in 2:S_its) {

    if (s %% S_print == 0) {
        print(paste0("s = ", s))
        print(paste0("Acc. rates: " ,
                     paste0(signif(accCounts/(s-1), 3), collapse = ", ")))
        # print(alpha_S[, s-1])
    }
    
    ##### 
    # alpha proposal
    #####
    
    # independent Uniform proposals
    # for each alphaj, where j = (2, P-1)
    # need a minimum, maximum for Uniform distribution
    # minimum = alphaj - widthj
    # maximum = alphaj + widthj
    
    # maybe could use pmin, to calculate these minimums simultaneously, i.e. element by element?
    # widths <- pmin(alphas[1:(P-2)] - alphas[2:(P-1)], 
    #                alphas[2:(P-1)] - alphas[3:P]) / 2 * 0.95
    # mins <- alphas[2:(P-1)] - widths
    # maxs <- alphas[2:(P-1)] + widths
    
    ### moving uniform proposals 
    # s_minmaxs <- ab_minmaxs(alphas, extend_hilo = extendRanges)
    # # proposal
    # alphap <- c(2, runif(P-2, s_minmaxs$mins, s_minmaxs$maxs), 1)
    # # bounds based on the proposal
    # p_minmaxs <- ab_minmaxs(alphap, extend_hilo = extendRanges)
    # boundsCondition <- all(p_minmaxs$mins <= alphas[2:(P-1)] & 
    #                            alphas[2:(P-1)] <= p_minmaxs$maxs)
    ### end of moving uniform proposals
    
    if (alphaConstant) {
        alphas <- alpha_0
        alpha_S[, s] <- alphas
    } else {
        
        ### uniform order statistics
        alphap <- c(2, sort(runif(P-2, 1, 2), decreasing = TRUE), 1)
        boundsCondition <- TRUE
        ### end of uniform order statistics
        
        if ( boundsCondition ) {
            ### moving uniform proposal densities
            # logJ_p <- sum(log( 1/(s_minmaxs$maxs - s_minmaxs$mins) ))
            # logJ_s <- sum(log( 1/(p_minmaxs$maxs - p_minmaxs$mins) ))
            ### uniform order statistics proposal densities
            logJ_p <- 0
            logJ_s <- 0
            
            # proceed with the rest of sampling
            
            # logr, log acceptance ratio
            # log0F0p <- hgf0F0MC(sqrt(ws) * alphap, sqrt(ws) * betas, 
            #                     betaf = 2, MCits = MCits, logscale = TRUE)
            
            batch_means <- vector("list", batch_count)
            batch_means <- foreach(i = 1:batch_count, .inorder = inorder) %dopar% {
                hgf0F0MC(sqrt(ws) * alphap, sqrt(ws) * betas, 
                         betaf = 2, MCits = MCits/batch_count, logscale = TRUE)
            }
            log0F0p <- logsumexp(unlist(batch_means)) - log(batch_count)
            
            traceterm <- t(alphap - alphas) %*% (ws * M %*% betas)
            
            logr <- -K*( log0F0p - log0F0s ) + traceterm + logJ_s - logJ_p
            
            # accept or reject
            if (log(runif(1)) <= logr) {
                # accept
                alphas <- alphap
                alpha_S[, s] <- alphas
                log0F0s <- log0F0p
                accCounts[1] <- accCounts[1] + 1
            } else {
                alpha_S[, s] <- alphas
            }
            
        } else {
            # reject outright: do not do the rest of sampling
            # set the previous value to the current value
            alpha_S[, s] <- alphas
        } # end of Boundary Condition for alpha
    } # end of constant or proposal condition
    
    ##### 
    # beta proposal
    #####
    
    if (betaConstant) {
        betas <- beta_0
        beta_S[, s] <- betas
    } else {
        
        ### uniform order statistics
        betap <- c(2, sort(runif(d-2, 1, 2), decreasing = TRUE), 1)
        boundsCondition <- TRUE
        ### end of uniform order statistics
        
        if ( boundsCondition ) {
            ### moving uniform proposal densities
            # logJ_p <- sum(log( 1/(s_minmaxs$maxs - s_minmaxs$mins) ))
            # logJ_s <- sum(log( 1/(p_minmaxs$maxs - p_minmaxs$mins) ))
            ### uniform order statistics proposal densities
            logJ_p <- 0
            logJ_s <- 0
            
            # proceed with the rest of sampling
            
            # logr, log acceptance ratio
            # log0F0p <- hgf0F0MC(sqrt(ws) * alphap, sqrt(ws) * betas, 
            #                     betaf = 2, MCits = MCits, logscale = TRUE)
            
            batch_means <- vector("list", batch_count)
            batch_means <- foreach(i = 1:batch_count, .inorder = inorder) %dopar% {
                hgf0F0MC(sqrt(ws) * alphas, sqrt(ws) * betap, 
                         betaf = 2, MCits = MCits/batch_count, logscale = TRUE)
            }
            log0F0p <- logsumexp(unlist(batch_means)) - log(batch_count)
            
            # traceterm <- t(alphap - alphas) %*% (ws * M %*% betas)
            traceterm <- ws * t(alphas) %*% M %*% (betap - betas)
            
            logr <- -K*( log0F0p - log0F0s ) + traceterm + logJ_s - logJ_p
            
            # accept or reject
            if (log(runif(1)) <= logr) {
                # accept
                betas <- betap
                beta_S[, s] <- betas
                log0F0s <- log0F0p
                accCounts[2] <- accCounts[2] + 1
            } else {
                beta_S[, s] <- betas
            }
            
        } else {
            # reject outright: do not do the rest of sampling
            # set the previous value to the current value
            beta_S[, s] <- betas
        } # end of Boundary Condition for beta
    } # end of constant or proposal condition
    
    #####
    # w proposal
    #####
    
    if (wConstant) {
        ws <- w_0
        w_S[s] <- ws
    } else {
        
        # make a proposal
        (wp <- truncdist::rtrunc(1, "norm", a = 0, b = Inf, mean = ws, sd = wp_sd))
        
        # calculate new 0F0 value
        batch_means <- vector("list", batch_count)
        batch_means <- foreach(i = 1:batch_count, .inorder = inorder) %dopar% {
            hgf0F0MC(sqrt(wp) * alphas, sqrt(wp) * betas, 
                     betaf = 2, MCits = MCits/batch_count, logscale = TRUE)
        }
        (log0F0p <- logsumexp(unlist(batch_means)) - log(batch_count))
        
        # calculate the rest of terms in acceptance ratio
        traceterm <- (wp - ws) * (t(alphas) %*% M %*% betas)
        
        logprior_w <- dgamma(wp, eta0_w/2, rate = tau02_w, log = TRUE) - 
            dgamma(ws, eta0_w/2, rate = tau02_w, log = TRUE)
        logprop_w <- 
            truncdist::dtrunc(ws, "norm", 0, Inf, mean = wp, sd = wp_sd, log = TRUE) -
            truncdist::dtrunc(wp, "norm", 0, Inf, mean = ws, sd = wp_sd, log = TRUE)
        
        # calculate acceptance ratio
        (logr <- -K*(log0F0p - log0F0s) + traceterm + logprior_w + logprop_w)
        
        # accept or reject
        if (log(runif(1)) <= logr) {
            # accept
            ws <- wp
            w_S[s] <- ws
            log0F0s <- log0F0p
            accCounts[3] <- accCounts[3] + 1
        } else {
            # reject
            w_S[s] <- ws
        } # end of Boundary Condition for w
    } # end of constant or proposal condition
    
# } # end of s loop

if (inParallel) {
    stopCluster(cluster)
}

accCounts/(S_its-1)

gibbsKeep <- floor(S_its/2):S_its
### investigate Alphas
alpha_0
plot(alpha_S[2, ])
abline(h = alpha_0[2], col = "red")

apply(alpha_S[2:(P-1), gibbsKeep], 1, 
      quantile, probs = c(0.025, 0.5, 0.975)) |> 
    t() |>
    cbind("true" = alpha_0[2:(P-1)])

apply(alpha_S[2:(P-1), gibbsKeep], 1,
      mean) |> 
    cbind("true" = alpha_0[2:(P-1)])

### investigate Betas
beta_0
plot(beta_S[2, ])
abline(h = beta_0[2], col = "red")

apply(beta_S[2:(d-1), gibbsKeep], 1, 
      quantile, probs = c(0.025, 0.5, 0.975)) |> 
    t() |>
    cbind("true" = beta_0[2:(d-1)])

apply(beta_S[2:(d-1), gibbsKeep], 1,
      mean) |> 
    cbind("true" = beta_0[2:(d-1)])

### investigate ws
w_0
plot(w_S)
abline(h = w_0, col = "red")

print(summary(w_S[gibbsKeep]))

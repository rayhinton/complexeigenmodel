# MH proposals for A, B

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
# source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_data-covar-dense.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

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

# Beta proposal investigation ---------------------------------------------

be_mode <- .5
# concentration must be greater than 2
be_conc <- 5

be_a <- be_mode * (be_conc - 2) + 1
be_b <- (1 - be_mode) * (be_conc - 2) + 1

curve(dbeta(x, shape1 = be_a, shape2 = be_b), 
      from = 0, to = 1, n = 201)

# if I simply scale drawn Beta r.v.s, does that match with the correct density?

y_min <- 1.2
y_max <- 1.8

ys <- rbeta(1000, be_a, be_b) * (y_max - y_min) + y_min

plot(density(ys))
curve(dbeta((x - y_min)/(y_max - y_min), 
            shape1 = be_a, shape2 = be_b) / (y_max - y_min), 
      from = y_min, to = y_max, n = 201, add = TRUE, col = "red")

P <- 8

set.seed(3)
alphas <- c(2, sort(runif(P-2, 1, 2), decreasing = TRUE), 1)

ls <- c(alphas[2:(P-2)] - 0.5 * (alphas[2:(P-2)] - alphas[3:(P-1)]), 1)
us <- c(2, alphas[3:(P-1)] + 0.5 * (alphas[2:(P-2)] - alphas[3:(P-1)]))

plot(alphas, rep(1, P), ylim = c(0, 5))
abline(v = ls, lty = 1:(P-2), col = 1:(P-2))
abline(v = us, lty = 1:(P-2), col = 1:(P-2))

modeX <- (alphas[2:(P-1)] - ls) / (us - ls)
a_pars <- modeX * (be_conc - 2) + 1
b_pars <- (1 - modeX) * (be_conc - 2) + 1
for (j in 1:(P-2)) {
    print(j)
    curve(dbeta((x - ls[j])/(us[j] - ls[j]), 
                shape1 = a_pars[j], shape2 = b_pars[j]) / (us[j] - ls[j]), 
          from = ls[j], to = us[j], n = 201, add = TRUE, col = j)
}

set.seed(4)
alpha_p <- c(2, rbeta(P-2, a_pars, b_pars) * (us - ls) + ls, 1)
# alpha_p <- c(2, c(1.85, 1.62, 1.55, 1.4, 1.28, 1.01), 1)

points(alpha_p, rep(1.5, P), pch = 3)
ls_p <- c(alpha_p[2:(P-2)] - 0.5 * (alpha_p[2:(P-2)] - alpha_p[3:(P-1)]), 1)
us_p <- c(2, alpha_p[3:(P-1)] + 0.5 * (alpha_p[2:(P-2)] - alpha_p[3:(P-1)]))

ls_p <= alphas[2:(P-1)] & alphas[2:(P-1)] <= us_p


# mixture of Betas --------------------------------------------------------



# finding different intervals, as a function ------------------------------

P <- 8

set.seed(3)
alphas <- c(2, sort(runif(P-2, 1, 2), decreasing = TRUE), 1)

s_minmaxs <- ab_minmaxs(alphas)
cbind(s_minmaxs$mins, s_minmaxs$maxs, alphas[2:(P-1)])

# proposal
# set.seed(5)
set.seed(8) # 8 results in all proposed values being in the bounds
alphap <- c(2, runif(P-2, s_minmaxs$mins, s_minmaxs$maxs), 1)
# J_p: p proposal density, given s previous
(J_p <- prod(1/(s_minmaxs$maxs - s_minmaxs$mins)))
(logJ_p <- sum(log( 1/(s_minmaxs$maxs - s_minmaxs$mins) )))

p_minmaxs <- ab_minmaxs(alphap)
# J_s: s previous density, given p proposal
# equal to 0, if one of the below are FALSE
p_minmaxs$mins <= alphas[2:(P-1)] & alphas[2:(P-1)] <= p_minmaxs$maxs

cbind(p_minmaxs$mins, p_minmaxs$maxs, alphas[2:(P-1)])

if ( all(p_minmaxs$mins <= alphas[2:(P-1)] & alphas[2:(P-1)] <= p_minmaxs$maxs) ) {
	# equal to this, if all of the above are TRUE
	(J_s <- prod(1/(p_minmaxs$maxs - p_minmaxs$mins)))
    
    (logJ_s <- sum(log( 1/(p_minmaxs$maxs - p_minmaxs$mins) )))
    
    # proceed with the rest of sampling
} else {
	(J_s <- 0)
    
    # do not do the rest of sampling
    # set the previous value to the current value
}

# do some plots to show the problem

plot(alphas, rep(1, P))
abline(v = c(s_minmaxs$mins, s_minmaxs$maxs))
points(alphap, rep(1.1, P), pch = 2, col = 2)
abline(v = c(p_minmaxs$mins, p_minmaxs$maxs), col = 2)


# w MH proposals with truncated Normal ------------------------------------

ws <- 4
wp_sd <- 2

wp <- truncdist::rtrunc(1, "norm", a = 0, b = Inf, mean = ws, sd = wp_sd)

curve(truncdist::dtrunc(x, "norm", 0, Inf, mean = ws, sd = wp_sd),
      from = 0, to = ws + 4*wp_sd, xlab = "w", ylab = "density")
abline(v = wp, lty = 2)

curve(truncdist::dtrunc(x, "norm", 0, Inf, mean = wp, sd = wp_sd),
      from = 0, to = wp + 4*wp_sd, add = TRUE, col = "red")
abline(v = ws, lty = 2, col = "red")

truncdist::dtrunc(ws, "norm", 0, Inf, mean = wp, sd = wp_sd)
truncdist::dtrunc(ws, "norm", 0, Inf, mean = wp, sd = wp_sd)
truncdist::dtrunc(wp, "norm", 0, Inf, mean = ws, sd = wp_sd)

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

S_its <- 2000
MCits <- 1e5
inorder <- TRUE
batch_count <- 10

parallel::detectCores()

cluster <- makeCluster(batch_count)
registerDoParallel(cluster)

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

# w scalar
# w_0
w_0 <- rgamma(1, 1, 1)

# "previous" alpha, beta samples
alpha_0 <- c(al_be_upper,
             runif(P-2, al_be_lower, al_be_upper) |> sort(decreasing = TRUE),
             al_be_lower)
beta_0 <- c(al_be_upper,
            runif(d-2, al_be_lower, al_be_upper) |> sort(decreasing = TRUE),
            al_be_lower)
alpha_0
beta_0

A_0 <- diag(sqrt(w_0) * alpha_0)
B_0 <- diag(sqrt(w_0) * beta_0)
G_0 <- V_0 %*% A_0 %*% t(Conj(V_0))

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

# set up alpha samples array ----------------------------------------------

# this can be used inside to the function, to derive w based on the convention
# ws <- (As[1] / al_be_upper)^2
# alphas <- As / sqrt(ws)
# betas <- Bs / sqrt(ws)
ws <- w_0
betas <- beta_0

# calculate matrix M, sum of function of V and U_k matrices (constant for all A, B)
M <- matrix(0, nrow = P, ncol = d)
for (k in 1:K) {
    M <- M + Re( Conj( t(Conj(V_0))%*%U_ks[, , k] ) * t(Conj(V_0))%*%U_ks[, , k] )
}

# this vector stays constant, when sampling alpha
Mbeta <- M %*% betas

# initialize
alpha_S <- matrix(NA, P, S_its)
alpha_S[c(1, P), ] <- c(al_be_upper, al_be_lower)
alpha_S[, 1:5]

# initialize the first sample
# set.seed(8032025)
set.seed(17052025)
alpha_S[2:(P-1), 1] <- sort(runif(P-2, 1, 2), decreasing = TRUE)
# alpha_S[, 1] <- seq(al_be_upper, al_be_lower, length.out = P)
# TODO unrealistic starting values, but testing
# alpha_S[, 1] <- alpha_0
alpha_S[, 1:5]

# alpha FCD sampling ------------------------------------------------------

log0F0s <- hgf0F0MC(sqrt(ws) * alpha_S[, 1], sqrt(ws) * betas, 
                    betaf = 2, MCits = MCits, logscale = TRUE)

# one step
s <- 2
accCount <- 0
for (s in 2:S_its) {
    if (s %% 100 == 0) {
        print(paste0("s = ", s))
    }
    
    alphas <- alpha_S[, s-1]
    
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
    
    s_minmaxs <- ab_minmaxs(alphas, extend_hilo = extendRanges)
    # proposal
    alphap <- c(2, runif(P-2, s_minmaxs$mins, s_minmaxs$maxs), 1)
    # bounds based on the proposal
    p_minmaxs <- ab_minmaxs(alphap, extend_hilo = extendRanges)
    
    if ( all(p_minmaxs$mins <= alphas[2:(P-1)] & alphas[2:(P-1)] <= p_minmaxs$maxs) ) {
        logJ_p <- sum(log( 1/(s_minmaxs$maxs - s_minmaxs$mins) ))
        logJ_s <- sum(log( 1/(p_minmaxs$maxs - p_minmaxs$mins) ))
        
        # proceed with the rest of sampling
        
        # logr, log acceptance ratio
        # log0F0p <- hgf0F0MC(sqrt(ws) * alphap, sqrt(ws) * betas, 
        #                     betaf = 2, MCits = MCits, logscale = TRUE)
        
        batch_means <- vector("list", batch_count)
        batch_means <- foreach(i = 1:batch_count, .inorder = FALSE) %dopar% {
            hgf0F0MC(sqrt(ws) * alphap, sqrt(ws) * betas, 
                     betaf = 2, MCits = MCits/batch_count, logscale = TRUE)
        }
        log0F0p <- logsumexp(unlist(batch_means)) - log(batch_count)
        
        traceterm <- t(alphap - alphas) %*% (ws * Mbeta)
        
        logr <- -K*( log0F0p - log0F0s ) + traceterm + logJ_s - logJ_p
        
        # accept or reject
        if (log(runif(1)) <= logr) {
            # accept
            alphas <- alphap
            alpha_S[, s] <- alphas
            log0F0s <- log0F0p
            accCount <- accCount + 1
        } else {
            alpha_S[, s] <- alphas
        }
        
    } else {
        # reject outright: do not do the rest of sampling
        # set the previous value to the current value
        alpha_S[, s] <- alphas
    }
    
} # end of s loop
    
stopCluster(cluster)

accCount/(S_its-1)

alpha_0
plot(alpha_S[2, ])
alpha_S[, s]

gibbsKeep <- floor(S_its/2):S_its
# gibbsKeep <- 1:1000

apply(alpha_S[2:(P-1), gibbsKeep], 1, 
      quantile, probs = c(0.025, 0.5, 0.975)) |> 
    t() |>
    cbind("true" = alpha_0[2:(P-1)])

apply(alpha_S[2:(P-1), gibbsKeep], 1,
      mean) |> 
    cbind("true" = alpha_0[2:(P-1)])

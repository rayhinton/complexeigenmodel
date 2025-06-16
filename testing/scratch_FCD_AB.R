# FCD sampler function for A, B matrices

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
# source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_data-covar-dense.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")


# parameters --------------------------------------------------------------

P <- 8
d <- 8
K <- 4
al_be_upper <- 2
al_be_lower <- 1
bing_its <- 5000
S_its <- 10000
gs <- 200
MCits <- 1e2

# original AB FCD function based on approximation -------------------------

AB_gibbs_densCovar <- function(data_list, param_list, gs = 201, al_be_lower = 1,
                               al_be_upper = 2,
							   eta0 = 2, tau02 = 1/1000) {
	P <- param_list$P
	d <- param_list$d
	K <- param_list$K
	Vs <- param_list$Vs
	U_ks <- param_list$U_ks
	As <- param_list$As
	Bs <- param_list$Bs
	
	# all i < j
	firstpairs <- combn(1:P, 2) |> t()
	# such that i < d + 1, equivalent to i <= d
	# TODO the choice of pairs should depend on the smallest value in A, B matrices
	# - if alpha, beta range from 2 to 1: use i <= d
	# - if alpha, beta range from 1 to 0: use i <= d - 1
	Ppairs <- firstpairs[firstpairs[, 1] <= d, ]
	
	# calculate matrix M, sum of function of V and U_k matrices
	M <- matrix(0, nrow = P, ncol = d)
	for (k in 1:K) {
	    M <- M + Conj( t(Conj(Vs))%*%U_ks[, , k] ) * t(Conj(Vs))%*%U_ks[, , k]
	}
	# M should be exactly real, but it comes out complex with 0 imaginary values
	M <- Re(M)
	
	ws <- (param_list$As[1] / al_be_upper)^2
	
	alphas <- As / sqrt(ws)
	betas <- Bs / sqrt(ws)
	
	##########
	# sample alpha
	##########
	# this vector stays constant, when sampling alpha=
	Mbeta <- M %*% betas
	
	# TODO sample a_ind randomly, so we do it in a different order each time
	for (a_ind in sample(2:(P-1))) {

		# for the a index: select all the rows where the 1st or 2nd column is a_ind
		# a_ind <- 3
		# these are the combined columns
		# Ppairs[Ppairs[, 1] == a_ind | Ppairs[, 2] == a_ind, ]
		# these have the indices I actually need
		others <- Ppairs[Ppairs[, 1] == a_ind, 2]
		some <- Ppairs[Ppairs[, 2] == a_ind, 1]

		# for each value in grid xs = seq(0, 1, .1) (or some other sequence),
		# calculate the product of (a[some] - x) and (x - a[others]). the
		# sampled alpha value should be between the surrounding alpha values.
		# since the alphas are sorted decreasing, that means it should be
		# between the alpha values with a_ind+1 and a_ind-1, minimum and maximum
		# respectively. However, since the prior is that the alpha values are
		# order statistics from a Uniform(1, 2) distribution, and thus they are
		# equal with probability 0, and further, since the approximation we are
		# using requires that the alpha values are not equal, we should not give
		# positive probability to the values possibly being equal in the
		# posterior. Therefore, we should not use a sequence starting and ending
		# exactly with the respective min and max alpha values. Instead, we will
		# begin and end the sequences on a small offset from the exact
		# surrounding alpha values.

		# xs <- seq(0, 1, .1)
		# create the min and max values for the sequence. Give them an offset, so that we do not possibly sample equal alpha values.
		min_al <- alphas[a_ind + 1]
		max_al <- alphas[a_ind - 1]
		gint <- (max_al - min_al)/(gs-1)
		min_al <- min_al + gint/2
		max_al <- max_al - gint/2

		xs <- seq(min_al, max_al, length.out = gs)

		# TODO check that this sum of logs is actually calculating the right log density
		aprod1 <- sapply(xs, function(x) {
		    sum(log((alphas[some] - x))) + sum(log(x - alphas[others]))
		    })

		# aprod2
		aprod2 <- ws * xs * Mbeta[a_ind]
		# calculate the proportional density, with a possible additional factor
		proplogdens <- K * aprod1 + aprod2
		if (2 <= a_ind & a_ind <= d) {
		    # optional: aprod3 (only needed if 2 <= a_ind <= d)
		    aprod3 <- -K * ws * alphas[a_ind] * betas[a_ind]
		    proplogdens <- proplogdens + aprod3
		}

		# draw one value from xs (i.e. an alpha value), using the logweights
		alphas[a_ind] <- sample_gumbel(xs, 1, proplogdens)
	} # end of alphas sample for loop
	
	
	##########
	# sample beta
	##########
	
	alphaM <- t(alphas) %*% M

	for (b_ind in sample(2:(d-1))) {
        
        # select all the rows where the 1st or 2nd column is b_ind
        # indexes that are after b_ind (i.e. 2nd)
        b2s <- Ppairs[Ppairs[, 1] == b_ind & Ppairs[, 2] <= d, 2]
        # indexes that are before b_ind (i.e. 1st)
        b1s <- Ppairs[Ppairs[, 2] == b_ind, 1]
        
        # create a grid of beta values, depending on the next highest and lowest
        # betas values. first, determine the upper and lower limits of the
        # grid.
        min_bl <- betas[b_ind + 1]
        max_bl <- betas[b_ind - 1]
        gint <- (max_bl - min_bl)/(gs-1)
        min_bl <- min_bl + gint/2
        max_bl <- max_bl - gint/2
        # make the grid
        xs <- seq(min_bl, max_bl, length.out = gs)
        
        # calculate terms in the log density
        bprod1 <- sapply(xs, function(x) {
            sum(log((betas[b1s] - x))) + sum(log(x - betas[b2s]))
            })
        bprod2 <- (P-d)*log(betas[b_ind])
        bprod3 <- xs * (-K * ws * alphas[b_ind] + alphaM[b_ind])
        
        # calculate the log density
        proplogdens_b <- K*bprod1 + K*bprod2 + bprod3
        
        # draw one value from xs (i.e. beta value), using the logweights
        betas[b_ind] <- sample_gumbel(xs, 1, proplogdens_b)
    }
	
	##########
	# sample w
	##########
	
	fstar <- choose(P, 2) - choose(P-d, 2)
	betatilde <- c(betas, rep(0, P-d))

	w_shape <- K*fstar + eta0/2
	w_rate <- K * (t(alphas) %*% betatilde) - t(alphas) %*% M %*% betas + tau02

	ws1 <- rgamma(1, shape = w_shape, rate = w_rate)
	
	return(list(As = sqrt(ws1) * alphas, 
	            Bs = sqrt(ws1) * betas))
}

# FCD by MC integration functions ----------------------------------------

### functions to simulate from complex Stiefel manifold
# return n independent samples from the standard complex normal distribution
rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

# draw matrix distributed uniformly on V_d^P(C)
rcstiefel <- function(P, d) {
    # matrix of independent standard complex Normals
    X1 <- matrix(rscnorm(P*d), ncol = d)
    
    # take QR decomposition - not guaranteed to be unique due to numerical methods
    X1qr <- qr(X1)
    QX1 <- qr.Q(X1qr)
    # extract sign of the diagonal elements of R
    D <- sign(Re(diag(qr.R(X1qr))))
    # transform by Q = QS, where S = diag(D)
    
    # Qfinal <- t(t(QX1) * D) # faster than matrix multiplication, for larger matrices
    Qfinal <- QX1 %*% diag(D)
    
    return(Qfinal)
}

# exp(matrixStats::logSumExp(log(unlist(intgd))) - log(M))

logsumexp <- function(x) {
    xstar <- max(x)
    lse <- xstar + log( sum(exp( x - xstar )) )
    return(lse)
}

hgf0F0MC <- function(x, y, MCits = 1e4, logscale = FALSE) {
	P <- length(x)
	d <- length(y)
	# intgd <- vector("list", MCits)
	logintgd <- rep(NA, MCits)

	for (i in 1:MCits) {
	# intgd <- foreach(i = 1:M) %dopar% {
	    U <- rcstiefel(P, d)
	    # intgd[i] <- exp(Re(t(x) %*% (Conj(U) * U) %*% y))
	    logintgd[i] <- Re(crossprod(x, (Conj(U) * U) %*% y))
	}
    
	if (logscale == FALSE) {
    	# return(c(Reduce("+", intgd)/MCits))
	    return(mean(exp(logintgd)))
	} else {
	    return(logsumexp(logintgd) - log(MCits))
	}
}


# Nasuda Laplace approximation functions ----------------------------------

# define a function for the beta-multivariate gamma function
mvGamma <- function(c, m, betaf, logscale = FALSE) {
    term1 <- .25*m*(m-1)*betaf*log(pi)
    term2 <- sum(lgamma(c - ((1:m)-1)/2 * betaf))
    
    if (logscale) {
        return(term1 + term2)
    } else {
        return(exp(term1 + term2))
    }
}

omega_mbeta <- function(m, betaf, logscale = FALSE) {
    numer <- m*log(2) + (m^2 * betaf/2)*log(pi)
    denom <- mvGamma(m*betaf/2, m, betaf, logscale = TRUE)
    
    # return on log or linear scale
    if (logscale) {
        return(numer - denom)
    } else {
        return(exp(numer - denom))
    }
}

Omega_ms <- function(ms, betaf, logscale = FALSE) {
    m <- sum(ms)
    
    # calculate the log product over ms values
    omegasum <- 0
    for (j in 1:length(ms)) {
        omegasum <- omegasum + omega_mbeta(ms[j], betaf, logscale = TRUE)
    }
    # divide by this term
    omegam <- omega_mbeta(m, betaf, logscale = TRUE)
    
    # return on log or linear scale
    if (logscale) {
        return(omegasum - omegam)
    } else {
        return(exp(omegasum - omegam))
    }
}

logLA_Nasu_0s <- function(as, bs, betaf) {
    m <- length(as)
    d <- length(bs)
    stopifnot(m >= d)
    
    ppairs <- t(combn(m, 2))
    ii <- ppairs[ppairs[, 1] <= d, 1]
    jj <- ppairs[ppairs[, 1] <= d, 2]
    ss <- length(ii)
    # with this adjustment, the same function handles whether as and bs are the
    # same or different lengths.
    if (m == d) {
        ms <- rep(1, m)
    } else {
        ms <- c(rep(1, d), m-d)
    }
    
    bs0 <- c(bs, rep(0, m-d))
    
    logJAB <- sum(log(2/betaf*(as[ii] - as[jj])*(bs0[ii] - bs0[jj])))
    
    result <- (betaf/2*ss)*log(2/betaf*pi) +
        Omega_ms(ms, betaf, logscale = TRUE) +
        -(betaf/2)*logJAB +
        sum(as*bs0)
    
    return(result)
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
set.seed(21052025)
Ukinit <- runitary(P, d)
for (i in 1:bing_its) {
    if (i %% 500 == 0) {
        print(paste0("i = ", i))
    }
    # for rectangular Bingham matrices
    # Ukinit <- rcmb(Ukinit, G_0, B_0)
    # for square Bingham matrices
    Ukinit <- rcBingUP_gibbs(Ukinit, G_0, B_0, Imtol = 1e4*.Machine$double.eps)
}
for (s in 1:(K*100)) {
    # for rectangular Bingham matrices
    # Ukinit <- rcmb(Ukinit, G_0, B_0)
    # for square Bingham matrices
    Ukinit <- rcBingUP_gibbs(Ukinit, G_0, B_0, Imtol = 1e4*.Machine$double.eps)
    if (s %% 100 == 0) {
        print(paste0("k = ", s/100))
        U_ks[, , s/100] <- Ukinit
    }
}

# test the sampler function ----------------------------------------------

# P <- param_list$P
# d <- param_list$d
# K <- param_list$K
# Vs <- param_list$Vs
# U_ks <- param_list$U_ks
# As <- param_list$As
# Bs <- param_list$Bs

# Vs <- runitary(P, P)
# 
# As <- seq(2, 1, length.out = P)
# Bs <- seq(2, 1, length.out = d)
# 
# # TODO properly generate Uk - this is just a placeholder to make sure the function runs
# U_ks <- array(NA, c(P, d, K))
# for (k in 1:K) {
# 	U_ks[, , k] <- runitary(P, d)
# }
# 
# 
# # TODO placeholder data - the AB FCD does not even require this data
# data_list <- list()
# 
# param_list <- list(
# 	P = P,
# 	d = d,
# 	K = K,
# 	Vs = Vs,
# 	U_ks = U_ks,
# 	As = As,
# 	Bs = Bs
# )
# 
# # test running the function
# AB_gibbs_densCovar(data_list, param_list)

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


# alpha sampling with MC integral for 0F0 ---------------------------------

library(foreach)
library(doParallel)

parallel::detectCores()

cluster <- makeCluster(14)
registerDoParallel(cluster)

print(Sys.time())
for (s in 2:S_its) {
    
    if (s %% 10 == 0) {
        print(paste0("s = ", s, ": ", Sys.time()))
    }
        
    alphas <- alpha_S[, s-1]
    # for (a_ind in 2:(P-1)) {
    for (a_ind in sample(2:(P-1))) {
        
        min_al <- alphas[a_ind + 1]
        max_al <- alphas[a_ind - 1]
        gint <- (max_al - min_al)/(gs-1)
        min_al <- min_al + gint/2
        max_al <- max_al - gint/2
        
        xs <- seq(min_al, max_al, length.out = gs)
        
        alphatmp <- alphas
        # hgfterm <- rep(NA, gs)
        hgfterm <- vector("list", gs)
        
        # for (j in 1:gs) {
        hgfterm <- foreach(j = 1:gs) %dopar% {
        
        # loop over the grid: j in 1:gs
        # compute the 0F0 approximation with one value of alpha
            alphatmp[a_ind] <- xs[j]
            
            # one value of 0F0
            hgfterm[j] <- hgf0F0MC(sqrt(ws)*alphatmp, sqrt(ws)*betas, MCits)
        }
        hgfterm <- unlist(hgfterm)
        # calculate log density
        # -K * log(hgf terms) +
        # w * xs * Mbeta[a_ind]
        logdens <- -K * log(hgfterm) + ws * xs * Mbeta[a_ind]
         
        # sample from the log density
        # alphas[a_ind] <- sample_gumbel(xs, 1, logdens)
        alphas[a_ind] <- sample_gumbel(xs, 1, 
                                       # log(lowess(xs, exp(logdens - max(logdens)))$y)
                                       # lowess(xs, logdens)$y
                                       lowess(xs, logdens - max(logdens))$y
                                       # logdens
                                       )
    }
    alpha_S[, s] <- alphas
}

# Save all objects with a custom filename including timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save.image(file = paste0("test-FCD-AB_", timestamp, ".RData"))
load("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/test-FCD-AB_20250514_201259.RData")

stopCluster(cl = cluster)

# investigate the (log) density from last iteration
log(lowess(xs, exp(logdens - max(logdens)))$y)

lowess(xs, exp(logdens - max(logdens)))$y
lowess(xs, exp(logdens))$y

plot(xs, (logdens - max(logdens)))
lines(lowess(xs, (logdens - max(logdens))))

plot(xs, exp(logdens))
lines(lowess(xs, exp(logdens)))
plot(xs, logdens)
lines(lowess(xs, logdens))

# investigate posterior intervals and sample mixing
gibbsKeep <- floor(S_its/2):S_its

apply(alpha_S[2:(P-1), gibbsKeep], 1, 
      quantile, probs = c(0.025, 0.5, 0.975)) |> 
    t() |>
    cbind("true" = alpha_0[2:(P-1)])

apply(alpha_S[2:(P-1), gibbsKeep], 1,
      mean) |> 
    cbind("true" = alpha_0[2:(P-1)])

plot(alpha_S[2, gibbsKeep], type = "l", main = "alpha2")
plot(alpha_S[3, gibbsKeep], type = "l", main = "alpha3")
plot(alpha_S[6, gibbsKeep], type = "l", main = "alpha6")
plot(alpha_S[7, gibbsKeep], type = "l", main = "alpha7")

plot(density(alpha_S[2, gibbsKeep]))
abline(v = alpha_0[2], col = "red")

plot(density(alpha_S[7, gibbsKeep]))
abline(v = alpha_0[7], col = "red")




# Metropolis-Hastings in Gibbs step for w ---------------------------------

# variance of the proposal distribution
sigma2_wMH <- 50

# w prior parameters
eta0 <- 2
tau02 <- 1/1000

w_its <- 100
w_S <- rep(NA, w_its)
w_S[1] <- ws
hgf_ws <- hgf0F0MC(sqrt(ws)*alphas, sqrt(ws)*betas, MCits = 1e4, logscale = TRUE)

# begin loop
s <- 2
w_acc <- 0
for (s in 2:w_its) {
    
    if (s %% 10 == 0) {
        print(paste0("s = ", s))
    }
    
    # sample w_p from proposal distribution
    # proposal distribution: prior on w
    a_prop <- ws^2 / sigma2_wMH
    l_prop <- ws / sigma2_wMH
    w_p <- rgamma(1, a_prop, rate = l_prop)
    
    # calculate acceptance ratio, r
    # r, first term: 0F0 ratios
    
    # 0F0 w_p = 0F0( sqrt(w_p)*alphas, sqrt(w_p)*betas)
    # 0F0 w_s = 0F0( sqrt(ws)*alphas, sqrt(ws)*betas)
    
    hgf_wp <- hgf0F0MC(sqrt(w_p)*alphas, sqrt(w_p)*betas, MCits = 1e4, logscale = TRUE)
    # hgf_ws <- hgf0F0MC(sqrt(ws)*alphas, sqrt(ws)*betas, MCits = 1e4, logscale = TRUE)
    
    # r, second term: etr ratios
    # M = (sum of V^H %*% U_k hadamard products)
    # C = t(alphas) %*% M %*% betas
    # trace = C*(w_p - ws)
    
    Cs <- t(alphas) %*% M %*% betas
    
    # log acceptance ratio
    logr <- -K*(hgf_wp - hgf_ws) + Cs*(w_p - ws) +
        # non-symmetric proposal portion
        (eta0/2 - a_prop)*log(w_p/ws) + (l_prop - tau02)*(w_p - ws)
    
    # accept
    if (log(runif(1)) <= logr) {
        ws <- w_p
        w_S[s] <- w_p
        w_acc <- w_acc + 1
        
        hgf_ws <- hgf_wp
        
    } else {
        w_S[s] <- ws
    }
}

summary(w_S)
plot(w_S)
w_acc/w_its

w_right <- 2.3
w_try <- 1.63
sigma2_try <- 50
curve(dgamma(x, w_try^2 / sigma2_try, rate = w_try/sigma2_try),
      from = 0, to = max(c(100, 1.5*w_try)))
abline(v = w_right, col = "red")

# alpha FCD by Laplace approximation, round 2 -----------------------------

for (s in 2:S_its) {
    
    if (s %% 100 == 0) {
        print(paste0("s = ", s, ": ", Sys.time()))
    }
    
    alphas <- alpha_S[, s-1]
    # for (a_ind in 2:(P-1)) {
    for (a_ind in sample(2:(P-1))) {
        
        min_al <- alphas[a_ind + 1]
        max_al <- alphas[a_ind - 1]
        gint <- (max_al - min_al)/(gs-1)
        min_al <- min_al + gint/2
        max_al <- max_al - gint/2
        
        xs <- seq(min_al, max_al, length.out = gs)
        
        alphatmp <- alphas
        # hgfterm <- rep(NA, gs)
        hgfterm <- vector("list", gs)
        
        for (j in 1:gs) {
        # hgfterm <- foreach(j = 1:gs) %dopar% {
            
            # loop over the grid: j in 1:gs
            # compute the 0F0 approximation with one value of alpha
            alphatmp[a_ind] <- xs[j]
            
            # one value of 0F0
            # hgfterm[j] <- hgf0F0MC(sqrt(ws)*alphatmp, sqrt(ws)*betas, MCits)
            hgfterm[j] <- logLA_Nasu_0s(sqrt(ws)*alphatmp, sqrt(ws)*betas, 2)
        }
        hgfterm <- unlist(hgfterm)
        # calculate log density
        # -K * log(hgf terms) +
        # w * xs * Mbeta[a_ind]
        logdens <- -K * hgfterm + ws * xs * Mbeta[a_ind]
        
        # sample from the log density
        alphas[a_ind] <- sample_gumbel(xs, 1, logdens)
    }
    alpha_S[, s] <- alphas
}

gibbsKeep <- floor(S_its/2):S_its
# gibbsKeep <- 1:1000

apply(alpha_S[2:(P-1), gibbsKeep], 1, 
      quantile, probs = c(0.025, 0.5, 0.975)) |> 
    t() |>
    cbind("true" = alpha_0[2:(P-1)])

apply(alpha_S[2:(P-1), gibbsKeep], 1,
      mean) |> 
    cbind("true" = alpha_0[2:(P-1)])

plot(alpha_S[2, gibbsKeep], type = "l", main = "alpha2")
plot(alpha_S[3, gibbsKeep], type = "l", main = "alpha3")
plot(alpha_S[6, gibbsKeep], type = "l", main = "alpha6")
plot(alpha_S[7, gibbsKeep], type = "l", main = "alpha7")

plot(density(alpha_S[2, gibbsKeep]))
abline(v = alpha_0[2])


# save.image("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_FCD_AB_redoFCDLaplace.Rdata.RData")
load("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_FCD_AB_redoFCDLaplace.RData")

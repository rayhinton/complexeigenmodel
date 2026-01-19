source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Uk_CMACG.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcomplex_wishart.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/matrix_distances.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/first-time-series/geoSS.R")

closest_semiorth <- function(U) {
 
    A <- Re(U)
    svd_result <- svd(A)
    V <- svd_result$u %*% t(svd_result$v)
    
    return(V)
}

sdf_ar2 <- function(omega, phis, sigma2 = 1) {
    phi1 <- phis[1]
    phi2 <- phis[2]
    
    denom <- 1 + phi1^2 + phi2^2 + 2*phi1*(phi2 - 1) * cos(2*pi*omega) - 
        2*phi2*cos(4*pi*omega)
    
    return(sigma2 / denom)
}

# unnormalized PDF
unnorm_logPDF <- function(x, tau2, ajk, mu, N, logscale = FALSE) {
    logf1 <- -N * log(1 + x) - ajk/(1 + x)
    logf2 <-  -(.5/tau2) * (x - mu)^2
    
    result <- ifelse(x <= 0, -Inf, logf1 + logf2)
    
    if (logscale) {
        return(result)
    } else {
        return(exp(result))
    }
}

logd_Uk <- function(Us, lw, invSigmas, sigma_k2, result_Lambdas, data_list_w, 
                    num_freqs, k) {
    
    # invSigmas <- solve(Sigmas)
    
    # calculate traces and other terms, for the MH acceptance ratio
    # for (l in 1:num_freqs) {
    Skl <- data_list_w[[lw]][[k]]
    lambdakl <- result_Lambdas[, k, lw]
    Omegakl <- diag(lambdakl / (1 + lambdakl))
    
    # trace <- Re(sum(diag( Us %*% Omegakl %*% t(Conj(Us)) %*% Skl )))
    traces <- Re(sum(Conj(Us) * (Skl %*% Us %*% Omegakl)))
    
    # } # end of summing over frequencies
    
    logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*% invSigmas %*% Us )))
    
    logdens <- traces/sigma_k2 - P*logdets
    
    return(logdens)
}

# Generate a smooth positive random function using splines
generate_smooth_positive_function <- function(N = 100, n_knots = 7, 
                                              log_mean = 0, log_sd = 1) {
    # Define evaluation points
    x <- 1:N
    
    # Sample knot locations (including boundaries)
    knot_locs <- sort(c(1, sample(2:(N-1), n_knots - 2, replace = FALSE), N))
    
    # Sample log-values at knots
    log_y_knots <- rnorm(n_knots, mean = log_mean, sd = log_sd)
    
    # Fit smooth spline through knots
    spline_fit <- splinefun(knot_locs, log_y_knots, method = "natural")
    
    # Evaluate at all points and exponentiate
    log_f <- spline_fit(x)
    f <- exp(log_f)
    
    return(list(x = x, f = f, knot_locs = knot_locs, log_y_knots = log_y_knots))
}

generate_VAR1_coef <- function(P, max_eigenvalue = 0.9) {
    # Generate random matrix and scale to control eigenvalues
    A1 <- matrix(rnorm(P * P), P, P)
    
    # Eigen decomposition
    eig <- eigen(A1)
    
    # Scale eigenvalues to be inside unit circle
    eig$values <- eig$values / max(Mod(eig$values)) * max_eigenvalue
    
    # Reconstruct matrix with scaled eigenvalues
    A1 <- Re(eig$vectors %*% diag(eig$values) %*% solve(eig$vectors))
    
    return(A1)
}

generate_AR1_covariance <- function(P, sigma2 = 1, rho = 0.5) {
    # Check validity
    if (abs(rho) >= 1) stop("rho must be in (-1, 1)")
    
    # Create AR(1) covariance matrix
    Sigma <- sigma2 * rho^abs(outer(1:P, 1:P, "-"))
    
    return(Sigma)
}

logdet <- function(X) {
    return( log(EigenR::Eigen_det(X)) )
}

# setup -------------------------------------------------------------------

# dataseed <- 21092025
# dataseed <- 22092025
# dataseed <- 10102025
# dataseed <- 17102025
dataseed <- 28112025

# parseed <- 314
parseed <- 3141
# parseed <- 963456789

P <- 4
d <- 2
K <- 2
Tt <- 1024 # length of time series
LL <- round(sqrt(Tt))
# LL <- round(5)

# options include: 1RW, 2RWPN
Lambda_prior <- "2RWPN"

gibbsIts <- 5000
burnin <- 0.5
gibbsPrint <- 100

num_freqs <- Tt/2 - 1

# time series generation parameters
# number of knots in the function that generates Lambda curves
n_knots <- 4

# Geodesic slice sampling parameters
w <- 10
m <- 2

### Ukl MH tuning parameters
# tau_Uk <- rep(.1, K)
tau_Ukl <- array(0.1, c(K, num_freqs))
num_tau_check <- 20
show_tau_tune_summ <- TRUE
doCayleyZeros <- FALSE
CayleyZeroProb <- 0.5

### Sigmal MH tuning parameters
n_Sig <- rep(50, num_freqs)

# grid of frequencies to calculate over
omegaw <- seq(1/Tt, by = 1/Tt, length.out = num_freqs)

# hyperparameters to the prior for tau2
tau2_a <- 1
tau2_b <- 1

### adapting MH tuning parameters
# what s index is the maximum burnin iteration?
burninS <- floor(gibbsIts * burnin)
# how often should the adaptation happen?
tau_numin <- floor(burninS / (num_tau_check))
# at which s iterations should the adaptation happen?
tau_s_check <- seq(tau_numin, burninS, tau_numin)

# adaptive version
# tau_s_check <- seq(50, gibbsIts, 50)

# generate true parameters ------------------------------------------------

# set.seed(parseed)
set.seed(112)
{
VARpars <- array(NA, c(P, P, K))
# for (k in 1:K) {
    # VARpars[, , k] <- generate_VAR1_coef(P, 0.8)
# }
VARpars[, , 1:K] <- generate_VAR1_coef(P, 0.8)
noiseSigma <- generate_AR1_covariance(P, sigma2 = 1, rho = 0.5)

# sigmakl02, scale parameter
# sigmakl02 <- matrix(rgamma(K*(num_freqs + 2), 1, 1), 
#                     K, num_freqs + 2)
# this makes it constant across all frequencies for each k
sigmakl02 <- matrix(rgamma(K, 1, 1),
                    K, Tt)

fkTR <- array(NA, c(P, P, K, Tt))
U_kl0 <- array(NA, c(P, d, K, Tt))
Lambdakl0 <- array(NA, c(d, K, Tt))

for (k in 1:K) {
    for (t in 1:Tt) {
        Hz <- solve( diag(P) - exp(-1i*2*pi * t/Tt) * VARpars[, , k] )
        fomega <- Hz %*% noiseSigma %*% t(Conj(Hz))
        
        f_evd <- eigen(fomega)
        
        thisU <- f_evd$vectors[, 1:d]
        thisLambda <- f_evd$values[1:d]
        
        fkTR[, , k, t] <- sigmakl02[k, t] * ( thisU %*% diag(thisLambda) %*% 
                                          t(Conj(thisU)) + diag(P) )
        U_kl0[, , k, t] <- thisU
        Lambdakl0[, k, t] <- thisLambda
    }
}

# check Lambdakl0 ---------------------------------------------------------

k <- 1
plot(Lambdakl0[1, k, 1:num_freqs], type = "l", ylab = "lambda",
     main = paste0("k = ", k), ylim = c(0, max(Lambdakl0[, k, ])))
lines(Lambdakl0[2, k, 1:num_freqs], col = 2)

k <- 2
plot(Lambdakl0[1, k, 1:num_freqs], type = "l", ylab = "lambda",
     main = paste0("k = ", k), ylim = c(0, max(Lambdakl0[, k, ])))
lines(Lambdakl0[2, k, 1:num_freqs], col = 2)

# calculate the 2nd differences -------------------------------------------

diffs2nd <- array(NA, c(d, K, Tt))

diffs2nd[, , 2:(Tt-1)] <- Lambdakl0[, , 1:(Tt-2)] - 2*Lambdakl0[, , 2:(Tt-1)] +
    Lambdakl0[, , 3:Tt]
    
Lambdakl0[1, 1, 1] - 2*Lambdakl0[1, 1, 2] + Lambdakl0[1, 1, 3]
diffs2nd[1, 1, 2]

apply(diffs2nd[, , 2:(Tt-1)]^2, c(1, 2), sum)

plot(diffs2nd[1, 1, 2:(num_freqs-1)]^2, type = "l")
plot(diffs2nd[2, 1, 2:(num_freqs-1)]^2, type = "l")

apply(diffs2nd[, , 2:(num_freqs-1)]^2, c(1, 2), max) / 
    apply(diffs2nd[, , 2:(num_freqs-1)]^2, c(1, 2), median)

}


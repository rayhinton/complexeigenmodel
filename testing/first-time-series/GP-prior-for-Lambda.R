source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Lambdak.R")

# Functions ---------------------------------------------------------------

sdf_ar2 <- function(omega, phis, sigma2 = 1) {
    phi1 <- phis[1]
    phi2 <- phis[2]
    
    denom <- 1 + phi1^2 + phi2^2 + 2*phi1*(phi2 - 1) * cos(2*pi*omega) - 
        2*phi2*cos(4*pi*omega)
    
    return(sigma2 / denom)
}

# Parameters and setup ----------------------------------------------------

# dataseed <- 21092025
# dataseed <- 22092025
dataseed <- 10102025
parseed <- 314

P <- 4
d <- 2
K <- 2
Tt <- 1024 # length of time series
LL <- round(sqrt(Tt))

num_freqs <- 384

omegaw <- seq(1/Tt, by = 1/Tt, length.out = num_freqs)

gibbsIts <- 10000

tau2_Lambda <- 0.5
# hyperparameters to the prior for tau2
tau2_a <- 1
tau2_b <- 1

set.seed(parseed)
sigmak02 <- c(5, 10)
Sigma0 <- rFTCW(diag(P), P+1, P)

U_k0 <- array(NA, c(P, d, K))
for (k in 1:K) {
    U <- qr.Q(qr(matrix(rnorm(P*d), ncol = d)))
    U_k0[, , k] <- U
}

# AR(2) parameters
# one pair of SDFs starts out separated by 1
phis <- c(1, -.25, -.75, -.5,
            .75, -.65, -1, -.05) |>
    array(c(2, d, K))

# one pair of SDFs starts out separated by 0.1
# phis <- c(1, -.25, -.75, -.5,
#             .15, -.85, -1, -.05) |>
#     array(c(2, d, K))

Lambdak0_w <- array(NA, c(d, K, num_freqs))

for (w in 1:num_freqs) {
    for (k in 1:K) {
        for (j in 1:d) {
            Lambdak0_w[j, k, w] <- sdf_ar2(omegaw[w], phis[, j, k])
        }
    }
}

plot(Lambdak0_w[1, 2, ], type = "l")
lines(Lambdak0_w[2, 2, ], type = "l", col = "red")

# Lambdak0 <- Lambdak0_w[, 2, ]

# simulate time series ----------------------------------------------------
TS_data <- list()
SDMests <- list()

len_freq <- Tt/2

Utp <- waveslim::sine.taper(Tt, LL)

set.seed(dataseed)
for (k in 1:K) {
    Xm <- matrix(NA, Tt, d)
    
    for (j in 1:d) {
        Xm[, j] <- arima.sim(list(ar = phis[, j, k]), Tt)
    }
    X <- ts(Xm, start = c(1, 1), frequency = 1)
    
    Ymat <- sqrt(sigmak02[k]) * (U_k0[, , k] %*% t(X) + 
                                     matrix(rnorm(P*Tt), ncol = Tt))
    Yts <- ts(t(Ymat), start = c(1, 1), freq = 1)
    
    TS_data[[k]] <- Yts
    
    # estimate SDMs with multitaper
    Y_tp <- apply(Utp, MARGIN = 2, function(u) u*Yts, simplify = FALSE)
    F_tp_list <- lapply(Y_tp, FUN = function(Y) astsa::mvspec(Y,plot = FALSE))
    F_tp1 <- array(0, c(P, P, len_freq))
    
    for (ell in 1:len_freq) {
        for(j in 1:length(F_tp_list)){
            F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
        }
        F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
    }
    
    SDMests[[k]] <- F_tp1 * Tt
}

# Data (exact SDMs, for now) ----------------------------------------------

# generate the K observed P-dim. time series

# Shat, complex pd Hermitian matrices, PxP x Nw
Shat <- array(NA, c(P, P, num_freqs))

for (w in 1:num_freqs){
    # TODO only sampling for k = 1, for now
    Shat[, , w] <- SDMests[[1]][, , w]
}

# Elliptical slice sampling, setup ----------------------------------------

sigma2ESS <- 5
lESS <- 0.04

Kappa <- matrix(NA, num_freqs - 1, num_freqs - 1)

for (i in 1:(num_freqs-1)) {
    for (j in 1:(num_freqs-1)) {
        Kappa[i, j] <- sigma2ESS * exp( -(omegaw[i] - omegaw[j])^2 / (2*lESS^2) )
    }
}

Kappa <- Kappa + diag(1e-8, num_freqs - 1)

CholKappa <- chol(Kappa)

microbenchmark::microbenchmark(
    CholKappa <- chol(Kappa),
    times = 500
)

# to test, need

# Nw is the amount of frequencies
# LL is the amount of tapers
# k is the index of the observation

# all of this needs to be done over all k in 1:K, but I can ignore that for now

# Shat, complex pd Hermitian matrices, PxP x Nw
# Uk, complex semi-unitary matrices, Pxd x Nw
# sigmak2, just one scalar for now, 1
# Lambdak, positive real diagonal matrices, dx Nw

# j is the entry of Lambda that we are sampling
j <- 1
# only sampling for k  = 1
sigmak2 <- sigmak02[1]

# initialize

Lambdak_s <- array(NA, c(d, num_freqs-1, gibbsIts))
Lambdak_s[, , 1] <- Lambdak0_w[, 1, 2:num_freqs]

# ESS, sampling -----------------------------------------------------------

onesum <- 0
Aw <- rep(NA, num_freqs)
Bw <- rep(NA, num_freqs-1)

for (w in 1:num_freqs) {
    # TODO assuming j = 1
    # change this to instead select i =/= j
    for (i in 2:d) {
        onesum <- onesum + -LL * log(Lambdak0_w[i, 1, w] + 1)
        onesum <- onesum +
            -LL/sigmak2 * Lambdak0_w[i, 1, w] + 1 / (Lambdak0_w[i, 1, w] + 1) *
            Re(t(Conj(U_k0[, i, 1])) %*% Shat[, , w] %*% U_k0[, i, 1])
    }
    
    Aw[w] <- -LL/sigmak2 * Re(t(Conj(U_k0[, j, 1])) %*% Shat[, , w] %*% U_k0[, j, 1])   
}

Aw <- Aw[2:num_freqs]

for (s in 2:gibbsIts) {
    # draw nu from prior N(mu, Kappa)
    nuZ <- rnorm(num_freqs - 1)
    nu <- t(CholKappa) %*% nuZ
    
    # set previous Lambda
    Lamthe <- Lambdak_s[j, , s-1]
    
    Bw <- -LL * log(Lamthe + 1) - Lamthe / (Lamthe + 1) * Aw
    sumtwo <- sum(Bw)
    logdatadens <- onesum + sumtwo
    
    logy <- log(runif(1)) + logdatadens
    
    # Draw initial proposal, defining a bracket:
    theta_init <- runif(1, 0, 2*pi)
    theta_min <- theta_init - 2*pi
    theta_max <- theta_init
    
    # New Lambda proposal is based on the drawn theta value:
    (Lamnew <- Lamthe*cos(theta_init) + nu*sin(theta_init))
    
    while(any(Lamnew <= 0)) {
        theta_init <- runif(1, 0, 2*pi)
        theta_min <- theta_init - 2*pi
        theta_max <- theta_init
        
        # New Lambda proposal is based on the drawn theta value:
        Lamnew <- Lamthe*cos(theta_init) + nu*sin(theta_init)
    }
    
    Lamnew
    
    # Decide whether to accept the new Lambda value:
    
    # calculate log data density at proposed Lambda
    Bw <- -LL * log(Lamnew + 1) - Lamnew / (Lamnew + 1) * Aw
    sumtwo <- sum(Bw)
    (logpropdens <- onesum + sumtwo)
    
    # Reject, shrink theta interval and draw new theta
    while (logpropdens <= logy) {
        if (theta_init < 0) {
            theta_min <- theta_init
        } else {
            theta_max <- theta_init
        }
        
        theta_init <- runif(1, theta_min, theta_max)
        
        # New Lambda proposal is based on the drawn theta value:
        Lamnew <- Lamthe*cos(theta_init) + nu*sin(theta_init)
        
        while(any(Lamnew <= 0) ) {
            theta_init <- runif(1, theta_min, theta_max)
            # New Lambda proposal is based on the drawn theta value:
            Lamnew <- Lamthe*cos(theta_init) + nu*sin(theta_init)
        }
        
        # calculate log data density at proposed Lambda
        Bw <- -LL * log(Lamnew + 1) - Lamnew / (Lamnew + 1) * Aw
        sumtwo <- sum(Bw)
        (logpropdens <- onesum + sumtwo)
    }
    
    # Accept
    Lambdak_s[1, , s] <- Lamnew
}

# true Lambda curve
# Lambdak0_w[1, 1, ]

gibbsPostBurn <- (gibbsIts/2):gibbsIts

dim(Lambdak_s)

apply(Lambdak_s[1, , ], c(1), mean)
apply(Lambdak_s[1, , ], c(1), quantile, probs = 0.025)
apply(Lambdak_s[1, , ], c(1), quantile, probs = 0.975)
mean(Lambdak_s[1, 1, ])
quantile(Lambdak_s[1, 1, ], probs = c(0.025, .975))

plot(apply(Lambdak_s[1, , gibbsPostBurn], c(1), mean), type = "l")
lines(Lambdak0_w[1, 1, ], col = "red")
lines(apply(Lambdak_s[1, , gibbsPostBurn], c(1), quantile, probs = 0.025), lty = 2)
lines(apply(Lambdak_s[1, , gibbsPostBurn], c(1), quantile, probs = 0.975), lty = 2)


# other Cholesky calculations ---------------------------------------------

# Efficient square root for stationary covariance on a regular 1D grid
# using circulant embedding

# Function to compute circulant embedding square root
stationary_sqrt <- function(grid_points, cov_function, ...) {
    n <- length(grid_points)
    
    # Compute covariances for first row (stationary => function of distance)
    distances <- abs(grid_points - grid_points[1])
    cov_vec <- cov_function(distances, ...)
    
    # Extend to circulant matrix (embed in size 2n)
    # Need: [cov_vec, 0, reversed cov_vec[-1]]
    circulant_vec <- c(cov_vec, 0, rev(cov_vec[-1]))
    
    # Eigenvalues via FFT (eigenvalues of circulant matrix)
    eigenvalues <- Re(fft(circulant_vec))
    
    # Check for numerical issues (eigenvalues should be non-negative)
    if(any(eigenvalues < -1e-10)) {
        warning("Negative eigenvalues detected - covariance may not be valid")
    }
    eigenvalues <- pmax(eigenvalues, 0)  # Numerical safety
    
    # Square root of eigenvalues
    sqrt_eigenvalues <- sqrt(eigenvalues)
    
    # Return function that applies the square root
    function(z) {
        # z should be standard normal vector of length n
        # Extend z with zeros
        z_extended <- c(z, rep(0, n))
        
        # Apply: R*z = IFFT(sqrt(lambda) * FFT(z))
        fft_z <- fft(z_extended)
        result <- fft(sqrt_eigenvalues * fft_z, inverse = TRUE) / length(z_extended)
        
        # Return first n elements (real part)
        Re(result[1:n])
    }
}

# Gaussian/Squared Exponential covariance
gaussian_cov <- function(d, sigma = 1, lengthscale = 1) {
    sigma^2 * exp(-(d / lengthscale)^2)
}

# Example usage
set.seed(123)
n <- 510
grid <- seq(0, 1, length.out = n)

# Get the square root operator
sqrt_K <- stationary_sqrt(grid, gaussian_cov, sigma = 2, lengthscale = 0.01)

# Generate samples
z <- rnorm(n)
sample <- sqrt_K(z)

# Visualize
plot(grid, sample, type = 'l', 
     main = "GP Sample via FFT-based Square Root",
     xlab = "x", ylab = "f(x)")

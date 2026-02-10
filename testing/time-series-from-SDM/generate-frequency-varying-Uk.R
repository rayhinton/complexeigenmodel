# generate time series with frequency-varying Uk 

source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("functions/matrix_distances.R")

library(splines)
library(expm)

# functions ---------------------------------------------------------------

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

library(expm)

generate_smooth_Uk <- 
    function(P, d, K, Tt, U_star_method, 
             n_basis = 3, scale_base = 0.5, scale_k = 0.1) {
    
    Uk <- array(NA_complex_, c(P, d, K, Tt))
    
    # Frequencies from 0 to π
    n_freq <- Tt/2 + 1
    omega <- 2 * pi * (0:(Tt/2)) / Tt
    
    # Random real orthonormal starting point
    U_0 <- qr.Q(qr(matrix(rnorm(P * d), P, d)))
    
    # Coefficient matrices for smooth base path
    C_R <- array(0, c(P, P, n_basis))  # skew-symmetric
    C_S <- array(0, c(P, P, n_basis))  # symmetric
    
    for (j in 1:n_basis) {
        temp <- matrix(rnorm(P^2), P, P) * scale_base
        C_R[,,j] <- (temp - t(temp)) / 2
        temp <- matrix(rnorm(P^2), P, P) * scale_base
        C_S[,,j] <- (temp + t(temp)) / 2
    }
    
    # Compute base path U*(ω)
    U_star <- array(NA_complex_, c(P, d, n_freq))
    
    if (U_star_method == "frombasis") {
        # Case 1: generate a smooth base path of U
        for (idx in 1:n_freq) {
            om <- omega[idx]
            A <- matrix(0 + 0i, P, P)
            for (j in 1:n_basis) {
                A <- A + cos(j * om) * C_R[,,j]
                A <- A + 1i * sin(j * om) * C_S[,,j]  # vanishes at 0 and π
            }
            U_star[,,idx] <- expm(A) %*% U_0
        }
    } else if (U_star_method == "VAR1") {
        # Case 2: use the U from a VAR(1) model
        VARpars <- generate_VAR1_coef(P, 0.8)
        noiseSigma <- generate_AR1_covariance(P, sigma2 = 1, rho = 0.5)
        
        for (idx in 1:n_freq) {
            om <- omega[idx]
            Hz <- solve( diag(P) - exp(-1i*om) * VARpars)
            fomega <- Hz %*% noiseSigma %*% t(Conj(Hz))
            f_evd <- eigen(fomega)
            
            U_star[, , idx] <- f_evd$vectors[, 1:d]
        }
    }
    
    # Series-specific rotations (real skew-symmetric preserves real-ness at boundaries)
    for (k in 1:K) {
        temp <- matrix(rnorm(P^2), P, P) * scale_k
        A_k <- (temp - t(temp)) / 2
        Q_k <- expm(A_k)
        
        for (idx in 1:n_freq) {
            U_k_omega <- Q_k %*% U_star[,,idx]
            
            if (idx == 1) {
                Uk[,,k,Tt] <- Re(U_k_omega)
            } else if (idx == n_freq) {
                Uk[,,k,Tt/2] <- Re(U_k_omega)
            } else {
                t_idx <- idx - 1
                Uk[,,k,t_idx] <- U_k_omega
                Uk[,,k,Tt - t_idx] <- Conj(U_k_omega)
            }
        }
    }
    
    return(list(Uk = Uk, U_star = U_star))
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

plot_sdm_smoothness <- function(sdm_array, n_diffs = 3, freq = NULL) {
    # sdm_array: P x P x T array of matrices evaluated at T frequencies
    # n_diffs:   number of finite difference orders to compute and plot
    # freq:      optional frequency grid (length T); defaults to equally spaced on [0, pi]
    
    T_len <- dim(sdm_array)[3]
    if (is.null(freq)) freq <- seq(0, pi, length.out = T_len)
    
    # Compute Frobenius norm at each frequency
    frob <- apply(sdm_array, 3, function(mat) norm(mat, type = "F"))
    
    diffs <- list(frob)
    f_grids <- list(freq)
    
    for (d in seq_len(n_diffs)) {
        prev <- diffs[[d]]
        h <- diff(f_grids[[d]])
        diffs[[d + 1]] <- diff(prev) / h
        f_grids[[d + 1]] <- f_grids[[d]][-1]
    }
    
    old.par <- par(no.readonly = TRUE)
    
    n_plots <- n_diffs + 1
    if (n_plots <= 4) {
        par(mfrow = c(n_plots, 1), mar = c(3, 4, 2, 1), oma = c(2, 0, 2, 0))
    } else {
        n_rows <- ceiling(n_plots / 2)
        par(mfrow = c(n_rows, 2), mar = c(3, 4, 2, 1), oma = c(2, 0, 2, 0))
    }
    
    for (d in seq_len(n_plots)) {
        label <- if (d == 1) "||f(ω)||_F" else paste0("Δ^", d - 1, " ||f||_F")
        plot(f_grids[[d]], diffs[[d]], type = "l",
             xlab = "", ylab = label,
             main = if (d == 1) "Frobenius norm of SDM" else "")
    }
    mtext("Frequency", side = 1, outer = TRUE)
    
    par(old.par)
}

# setup -------------------------------------------------------------------

Tt <- 1024
P <- 4
d <- 2
K <- 2

scale_k <- 0.1

len_freq <- Tt/2
LL <- round(sqrt(Tt))

Utp <- waveslim::sine.taper(Tt, LL)

# for the random Lambda curves
n_knots <- 4

set.seed(30102025)

sigmak2 <- rgamma(K, 1, 1)

# generate Uks ------------------------------------------------------------

# generate_smooth_Uk <- function(P, d, K, Tt, n_basis = 3, scale_base = 0.5, scale_k = 0.1) {

Ukall <- generate_smooth_Uk(P, d, K, Tt, U_star_method = "VAR1")
Uk <- Ukall$Uk
dim(Uk)

# Uk <- array(NA, c(P, d, K, Tt))
# 
# for (t in 1:(Tt/2 - 1)) {
#     for (k in 1:K) {
#         Uk[, , k, t] <- rCMACG(P, d, diag(P))
#         Uk[, , k, Tt - t] <- Conj(Uk[, , k, t])
#     }
# }
# 
# for (k in 1:K) {
#     Uk[, , k, Tt/2] <- rMACG(P, d, diag(P))
#     Uk[, , k, Tt] <- rMACG(P, d, diag(P))
# }


# investigate smoothness of U_star and Uk ---------------------------------

U_star_ds <- rep(NA, dim(Ukall$U_star)[3]-1)
Uk_ds <- matrix(NA, K, length(U_star_ds))

for (i in 2:length(U_star_ds)) {
    U_star_ds[i-1] <- frob_dist(Ukall$U_star[, , i-1], Ukall$U_star[, , i])
    
    for (k in 1:K) {
        Uk_ds[k, i-1] <- frob_dist(Uk[, , k, i-1], Uk[, , k, i])
    }
}

plot(U_star_ds, type = "l")
plot(Uk_ds[1, ], type = "l")

# look at entrywise Ustar - it should be smooth
plot(Re(Ukall$U_star[3, 2, ]), type = "l")

# generate Lambdas --------------------------------------------------------

Lambdak <- array(NA, c(d, K, Tt))

for (k in 1:K) {
    for(j in 1:d) {
        Lambdak[j, k, c(Tt, 1:(Tt/2))] <- 
            generate_smooth_positive_function(Tt/2 + 1, n_knots = n_knots, 
                                              log_mean = 1, log_sd = 0.25)$f
        Lambdak[j, k, (Tt/2 + 1) : (Tt - 1)] <- Lambdak[j, k, (Tt/2 - 1) : 1]
    }
    
    Lambdak[, k, ] <- Lambdak[order(Lambdak[, k, 1], decreasing = TRUE), k, ]
}

plot(Lambdak[1, 1, ], type = "l", 
     ylim = c(min(Lambdak[, 1, ]), max(Lambdak[, 1, ])))
lines(Lambdak[2, 1, ], col = "red")

plot(Lambdak[1, 2, ], type = "l",
     ylim = c(min(Lambdak[, 2, ]), max(Lambdak[, 2, ])))
lines(Lambdak[2, 2, ], col = "red")

# calculate true SDMs -----------------------------------------------------

Sl <- array(NA, c(P, P, K, Tt))

for (k in 1:K) {
    for (t in 1:Tt) {
        Sl[, , k, t] <- sigmak2[k] * 
            (Uk[, , k, t] %*% diag(Lambdak[, k, t]) %*% t(Conj(Uk[, , k, t])) +
                  diag(P))
        diag(Sl[, , k, t]) <- Re(diag(Sl[, , k, t]))
    }
}

# plot the diagonal entries for K = 1
plot(Re(Sl[1, 1, 1, ]), type = "l")
lines(Re(Sl[2, 2, 1, ]), col = 2)
lines(Re(Sl[3, 3, 1, ]), col = 3)
lines(Re(Sl[4, 4, 1, ]), col = 4)

# plot the diagonal entries for K = 2
plot(Re(Sl[1, 1, 2, ]), type = "l")
lines(Re(Sl[2, 2, 2, ]), col = 2)
lines(Re(Sl[3, 3, 2, ]), col = 3)
lines(Re(Sl[4, 4, 2, ]), col = 4)

# calculate distances between the SDMs ------------------------------------

SDM_ds <- matrix(NA, K, Tt)

for (t in 2:Tt) {
    for (k in 1:K) {
        SDM_ds[k, t-1] <- frob_dist(Sl[, , k, t-1], Sl[, , k, t])
    }
}

plot(SDM_ds[2, 1:(Tt/2)], type = "l")

# calculate Cholesky decompositions ---------------------------------------

Rfs <- array(NA, c(P, P, K, Tt))

for (k in 1:K) {
    for (t in 1:Tt) {
        # need lower triangular part, take t(Conj())
        # AFAICT, this is the unique Chol. with positive real diagonals
        Rfs[, , k, t] <- t(Conj(EigenR::Eigen_chol(Sl[, , k, t])))
    }
}

# generate random vectors -------------------------------------------------

Zs <- array(NA, c(P, K, Tt))

for (k in 1:K) {
    Zs[, k, 1:(Tt/2 - 1)] <- rscnorm(P * (Tt/2 - 1)) / sqrt(Tt)
    Zs[, k, c(Tt/2, Tt)] <- rnorm(P * 2) / sqrt(Tt)
    
    Zs[, k, (Tt/2 + 1):(Tt - 1)] <- Conj(Zs[, k, (Tt/2 - 1):1])
}


# generate the time series ------------------------------------------------

Yts <- array(0, c(P, K, Tt))

for (k in 1:K) {
    for (t in 1:Tt) {
        for (l in 1:Tt) {
            Yts[, k, t] <- Yts[, k, t] + 
                Rfs[, , k, l] %*% Zs[, k, l] * exp(2 * pi * 1i * l/Tt * t)
        }
    }
}

# verify the imaginary parts are small (i.e. to numerical precision)
max(abs(Im(Yts)))

# set Y to the real components, only
Yts <- Re(Yts)

plot(Yts[1, 1, ], type = "l")
lines(Yts[2, 1, ], col = 2)
lines(Yts[3, 1, ], col = 3)
lines(Yts[4, 1, ], col = 4)

plot(Yts[1, 2, ], type = "l")
lines(Yts[2, 2, ], col = 2)
lines(Yts[3, 2, ], col = 3)
lines(Yts[4, 2, ], col = 4)


# convert one time series to a matrix -------------------------------------

Yts[, 1, 1:10] |> t()

# SDM estimation, previous code -------------------------------------------

SDMests <- list()
TS_data <- list()

# Yts <- ts(t(Ymat), start = c(1, 1), freq = 1)

for (k in 1:K) {
    thisYt <- ts(t(Yts[, k, ]))
    TS_data[[k]] <- thisYt
    
    Y_tp <- apply(Utp, MARGIN = 2, function(u) u*thisYt, simplify = FALSE)
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

# compare estimates and true ----------------------------------------------

k <- 2
j <- 4
plot(Re(SDMests[[k]][j, j, 1:(len_freq - 1)]), type = "l", 
     ylab = "spectral density")
lines(Re(Sl[j, j, k, 1:(Tt/2 - 1)]), lty = 2)

# a test of distances -----------------------------------------------------

for (scale_k in c(0.001, 0.01, 0.1, 1, 10)) {
    temp <- matrix(rnorm(P^2), P, P) * scale_k
    A_k <- (temp - t(temp)) / 2
    Q_k <- expm::expm(A_k)
    cat("scale_k =", scale_k, ": ||Q_k - I||_F =", norm(Q_k - diag(P), "F"), "\n")
}

# After generating Uk, check directly:
cat("Direct Frobenius between U_1 and U_2 at freq 1:", 
    norm(U_kl0[,,1,2] - U_kl0[,,2,2], "F"), "\n")

# And check a few matrices visually
print(U_kl0[,,1,2])
print(U_kl0[,,2,2])


# alternative: VAR(1) as a base path --------------------------------------

temp <- matrix(rnorm(P^2), P, P) * scale_k
A_k <- (temp - t(temp)) / 2
Q_k <- expm(A_k)

par_VAR1 <- generate_VAR1_coef(P, 0.8)
U_VAR1 <- array(NA, c(P, d, Tt))
Uk_mats <- array(NA, c(P, d, Tt))
trunc_SDM_VAR1 <- array(NA, c(P, P, Tt))
Sk_mats <- array(NA, c(P, P, Tt))
orig_SDM_VAR1 <- array(NA, c(P, P, Tt))

true_evals <- matrix(NA, P, Tt)

U_VAR1_ds <- rep(NA, dim(U_VAR1)[3]-1)
Uk_VAR_ds <- matrix(NA, K, length(U_VAR1_ds))
trunc_VAR1_ds <- rep(NA, dim(U_VAR1)[3]-1)
orig_VAR1_ds <- rep(NA, dim(U_VAR1)[3]-1)
Sk_ds <- rep(NA, dim(U_VAR1)[3]-1)

noiseSigma <- generate_AR1_covariance(P, sigma2 = 1, rho = 0.5)

for (t in 1:Tt) {
    Hz <- solve( diag(P) - exp(-1i*2*pi * t/Tt) * par_VAR1 )
    fomega <- Hz %*% noiseSigma %*% t(Conj(Hz))
    
    f_evd <- eigen(fomega)
    
    thisU <- f_evd$vectors[, 1:d]
    thisLambda <- f_evd$values[1:d]
    
    Uk <- Q_k %*% thisU
    Uk_mats[, , t] <- Uk
    
    trunc_SDM_VAR1[, , t] <- thisU %*% diag(thisLambda) %*% t(Conj(thisU)) + 
        diag(P)
    U_VAR1[, , t] <- thisU
    orig_SDM_VAR1[, , t] <- fomega
    
    Sk_mats[, , t] <- Uk %*% diag(thisLambda) %*% t(Conj(Uk)) + diag(P)
    
    true_evals[, t] <- f_evd$values
    # Lambdakl0[, k, t] <- thisLambda
}


for (i in 2:length(U_VAR1_ds)) {
    U_VAR1_ds[i-1] <- evec_Frob_stat(U_VAR1[, , i-1], U_VAR1[, , i])
    orig_VAR1_ds[i-1] <- frob_dist(orig_SDM_VAR1[, , i-1], orig_SDM_VAR1[, , i])
    
    trunc_VAR1_ds[i-1] <- frob_dist(trunc_SDM_VAR1[, , i-1], 
                                    trunc_SDM_VAR1[, , i])
    
    # for (k in 1:K) {
        # Uk_ds[k, i-1] <- frob_dist(Uk[, , k, i-1], Uk[, , k, i])
    # }
}

plot(U_VAR1_ds, type = "l")
# plot(Uk_ds[1, ], type = "l")

plot(orig_VAR1_ds, type = "l")
plot(trunc_VAR1_ds, type = "l")

plot(true_evals[1, ], type = "l", ylim = c(0, max(true_evals)))
lines(true_evals[2, ], col = 2)
lines(true_evals[3, ], col = 3)
lines(true_evals[4, ], col = 4)

plot(true_evals[3, ], type = "l")

plot(Im(orig_SDM_VAR1[1, 2, ]), type = "l")

plot_sdm_smoothness(orig_SDM_VAR1, n_diffs = 7)
plot_sdm_smoothness(trunc_SDM_VAR1, n_diffs = 7)

plot_sdm_smoothness(Sk_mats, n_diffs = 7)

# from the MCMC sampler script
plot_sdm_smoothness(fkTR[, , 1, 1:(Tt/2-1)], n_diffs = 7)

plot_sdm_smoothness(Sl[, , 1, 1:(Tt/2-1)], n_diffs = 7)

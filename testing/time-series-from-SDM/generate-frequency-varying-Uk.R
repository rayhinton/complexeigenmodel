# generate time series with frequency-varying Uk 

source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")

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

generate_smooth_Uk <- function(P, d, K, Tt, n_basis = 3, scale_base = 0.5, scale_k = 0.1) {
    
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
    
    for (idx in 1:n_freq) {
        om <- omega[idx]
        A <- matrix(0 + 0i, P, P)
        for (j in 1:n_basis) {
            A <- A + cos(j * om) * C_R[,,j]
            A <- A + 1i * sin(j * om) * C_S[,,j]  # vanishes at 0 and π
        }
        U_star[,,idx] <- expm(A) %*% U_0
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
    
    return(Uk)
}

# setup -------------------------------------------------------------------

Tt <- 1024
P <- 4
d <- 2
K <- 2

len_freq <- Tt/2
LL <- round(sqrt(Tt))

Utp <- waveslim::sine.taper(Tt, LL)

# for the random Lambda curves
n_knots <- 4

set.seed(30102025)

sigmak2 <- rgamma(K, 1, 1)

# generate Uks ------------------------------------------------------------

# generate_smooth_Uk <- function(P, d, K, Tt, n_basis = 3, scale_base = 0.5, scale_k = 0.1) {

Uk <- generate_smooth_Uk(P, d, K, Tt)
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


# generate time series with frequency-varying Uk 

source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")

library(splines)

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

# setup -------------------------------------------------------------------

Tt <- 512
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

Uk <- array(NA, c(P, d, K, Tt))

for (t in 1:(Tt/2 - 1)) {
    for (k in 1:K) {
        Uk[, , k, t] <- rCMACG(P, d, diag(P))
        Uk[, , k, Tt - t] <- Conj(Uk[, , k, t])
    }
}

for (k in 1:K) {
    Uk[, , k, Tt/2] <- rMACG(P, d, diag(P))
    Uk[, , k, Tt] <- rMACG(P, d, diag(P))
}

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
plot(Re(Sl[2, 2, 1, ]), type = "l")
plot(Re(Sl[3, 3, 1, ]), type = "l")
plot(Re(Sl[4, 4, 1, ]), type = "l")

# plot the diagonal entries for K = 2
plot(Re(Sl[1, 1, 2, ]), type = "l")
plot(Re(Sl[2, 2, 2, ]), type = "l")
plot(Re(Sl[3, 3, 2, ]), type = "l")
plot(Re(Sl[4, 4, 2, ]), type = "l")

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


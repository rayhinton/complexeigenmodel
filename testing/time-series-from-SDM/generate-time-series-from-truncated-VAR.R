source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")

library(splines)

# functions ---------------------------------------------------------------

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

# setup -------------------------------------------------------------------

P <- 4
d <- 2
K <- 1

Tt <- 1024

len_freq <- Tt/2
LL <- round(sqrt(Tt))

Utp <- waveslim::sine.taper(Tt, LL)

# testing -----------------------------------------------------------------

# Example usage
A1 <- generate_VAR1_coef(3, max_eigenvalue = 0.8)
max(Mod(eigen(A1)$values))  # Check: should be ≤ 0.8

# generate ----------------------------------------------------------------

sigmakl2 <- rep(2, Tt)

# sigmakl2 <- rep(NA, Tt)
# sigmakl2[c(1:(Tt/2), Tt)] <- rgamma(Tt/2 + 1, 1, 1)
# sigmakl2[(Tt-1) : (Tt/2+1)] <- sigmakl2[1:(Tt/2 - 1)]

# generate a valid matrix, A1, of VAR(1) coefficients
A1 <- generate_VAR1_coef(P, max_eigenvalue = 0.8)
# generate noiseSigma, the covariance matrix for the noise
# noiseSigma <- diag(P)
noiseSigma <- generate_AR1_covariance(P, sigma2 = 2, rho = 0.5)

fTR <- array(NA, c(P, P, Tt))
Ul <- array(NA, c(P, d, Tt))
Lambdal <- matrix(NA, d, Tt)

# at each frequency, t/Tt, calculate quantities:
for (t in 1:Tt) {
    Hz <- solve( diag(P) - exp(-1i*2*pi * t/Tt) * A1 )
    fomega <- Hz %*% noiseSigma %*% t(Conj(Hz))
    
    f_evd <- eigen(fomega)
    
    thisU <- f_evd$vectors[, 1:d]
    thisLambda <- f_evd$values[1:d]
    
    fTR[, , t] <- sigmakl2[t] * ( thisU %*% diag(thisLambda) %*% t(Conj(thisU)) + diag(P) )
    Ul[, , t] <- thisU
    Lambdal[, t] <- thisLambda
    
}

plot(Lambdal[1, ], type = "l", ylim = c(0, max(Lambdal)))
lines(Lambdal[2, ], col = 2)

Ul[, , 1]
Ul[, , 2]
Ul[, , 3]



# calculate Sigmal at each t ----------------------------------------------

# Sigmal <- array(NA, c(P, P, Tt/2 - 1))
# 
# for(t in 1: (Tt/2 - 1)) {
#     avg_Uk <- apply(Ul[, , t])
# }

# avg_Uk0 <- apply(U_k0, c(1, 2), mean)
# avg_Uk0 <- avg_Uk0 %*% 
#     solve( EigenR::Eigen_sqrt( t(Conj(avg_Uk0)) %*% avg_Uk0 ) )
# 
# avg_Uk0_perp <- (qr(avg_Uk0, complete = TRUE) |> 
#                      qr.Q(complete = TRUE))[, (d+1):P]
# 
# V_Sigma0 <- cbind(avg_Uk0, avg_Uk0_perp)
# Sigma0 <- V_Sigma0 %*% diag(P*(P:1) / sum(P:1)) %*% t(Conj(V_Sigma0))

# calculate Cholesky decompositions ---------------------------------------

Rfs <- array(NA, c(P, P, K, Tt))

for (k in 1:K) {
    for (t in 1:Tt) {
        # need lower triangular part, take t(Conj())
        # AFAICT, this is the unique Chol. with positive real diagonals
        Rfs[, , k, t] <- t(Conj(EigenR::Eigen_chol(fTR[, , t])))
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

# plot(Yts[1, 2, ], type = "l")
# lines(Yts[2, 2, ], col = 2)
# lines(Yts[3, 2, ], col = 3)
# lines(Yts[4, 2, ], col = 4)

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

k <- 1
j <- 2
plot(Re(SDMests[[k]][j, j, 1:(len_freq - 1)]), type = "l", 
     ylab = "spectral density")
lines(Re(fTR[j, j, 1:(Tt/2 - 1)]), lty = 2)

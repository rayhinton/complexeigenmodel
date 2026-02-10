# generate by rotating one VAR(1) path

source("testing/time-series-from-SDM/model-simulation-parameters.R")

source("functions/utility.R")
source("functions/FCD_Uk_CMACG.R")
source("functions/rcomplex_wishart.R")
source("functions/matrix_distances.R")
source("testing/first-time-series/geoSS.R")

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
K <- 2
Tt <- 1024

scale_k <- 0.1

# base path ---------------------------------------------------------------

VAR_pars <- generate_VAR1_coef(P, 0.9)
noiseSigma <- generate_AR1_covariance(P, sigma2 = 1, rho = 0.5)

Ukl0 <- array(NA, c(P, d, K, Tt))
base_VAR1_SDM <- array(NA, c(P, P, Tt))
base_trunc_SDM <- array(NA, c(P, P, Tt))
Skl <- array(NA, c(P, P, K, Tt))
Q_k <- array(NA, c(P, P, K))

for (k in 1:K) {
    temp <- matrix(rnorm(P^2), P, P) * scale_k
    A_k <- (temp - t(temp)) / 2
    Q_k[, , k] <- expm::expm(A_k)
}

for (t in 1:Tt) {
    Hz <- solve( diag(P) - exp(-1i*2*pi * t/Tt) * VAR_pars )
    fomega <- Hz %*% noiseSigma %*% t(Conj(Hz))
    f_evd <- eigen(fomega)
    
    thisU <- f_evd$vectors[, 1:d]
    thisLambda <- f_evd$values[1:d]
    
    base_VAR1_SDM[, , t] <- fomega
    base_trunc_SDM[, , t] <- thisU %*% diag(thisLambda) %*% t(Conj(thisU)) +
        diag(P)
    
    for (k in 1:K) {
        Ukl0[, , k, t] <- Q_k[, , k] %*% thisU
        Skl[, , k, t] <- Ukl0[, , k, t] %*% diag(thisLambda) %*% 
            t(Conj(Ukl0[, , k, t])) + diag(P)
    }
}

# for (k in 1:K) {
#     for (t in 1:Tt) {
#         Hz <- solve( diag(P) - exp(-1i*2*pi * t/Tt) * VARpars[, , k] )
#         fomega <- Hz %*% noiseSigma %*% t(Conj(Hz))
#         
#         f_evd <- eigen(fomega)
#         
#         thisU <- f_evd$vectors[, 1:d]
#         thisLambda <- f_evd$values[1:d]
#         
#         fkTR[, , k, t] <- sigmak02[k] * ( thisU %*% diag(thisLambda) %*% 
#                                               t(Conj(thisU)) + diag(P) )
#         U_kl0[, , k, t] <- thisU
#         Lambdakl0[, k, t] <- thisLambda
#     }
# }


# analyze SDMs ------------------------------------------------------------

plot_sdm_smoothness(base_VAR1_SDM, n_diffs = 7)
plot_sdm_smoothness(base_trunc_SDM, n_diffs = 7)
plot_sdm_smoothness(Skl[, , 1, ], n_diffs = 7)

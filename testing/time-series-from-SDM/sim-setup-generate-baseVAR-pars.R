# sim-setup

# sigmak02 <- rgamma(K, 1, 1)
# 
# fkTR <- array(NA, c(P, P, K, Tt))
# U_kl0 <- array(NA, c(P, d, K, Tt))
# Lambdakl0 <- array(NA, c(d, K, Tt))

VAR_pars <- generate_VAR1_coef(P, 0.9)
noiseSigma <- generate_AR1_covariance(P, sigma2 = 1, rho = 0.5)

U_kl0 <- array(NA, c(P, d, K, Tt))
base_VAR1_SDM <- array(NA, c(P, P, Tt))
base_trunc_SDM <- array(NA, c(P, P, Tt))
fkTR <- array(NA, c(P, P, K, Tt))
Q_k <- array(NA, c(P, P, K))

for (k in 1:K) {
    temp <- matrix(rnorm(P^2), P, P) * U_k_scale_k
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
        U_kl0[, , k, t] <- Q_k[, , k] %*% thisU
        fkTR[, , k, t] <- 
            sigmak02[k] * (U_kl0[, , k, t] %*% diag(thisLambda) %*% 
                               t(Conj(U_kl0[, , k, t])) + diag(P))
        
        Lambdakl0[, k, t] <- thisLambda
    }
}
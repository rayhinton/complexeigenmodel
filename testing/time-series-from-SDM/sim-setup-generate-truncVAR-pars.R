# generate truncated VAR parameters

VARpars <- array(NA, c(P, P, K))
if (all_same_VAR_pars) { 
    VARpars[, , 1:K] <- generate_VAR1_coef(P, 0.8)
} else {
    for (k in 1:K) {
        VARpars[, , k] <- generate_VAR1_coef(P, 0.8)
    }
}

noiseSigma <- generate_AR1_covariance(P, sigma2 = 1, rho = 0.5)

for (k in 1:K) {
    for (t in 1:Tt) {
        Hz <- solve( diag(P) - exp(-1i*2*pi * t/Tt) * VARpars[, , k] )
        fomega <- Hz %*% noiseSigma %*% t(Conj(Hz))
        
        f_evd <- eigen(fomega)
        
        thisU <- f_evd$vectors[, 1:d]
        thisLambda <- f_evd$values[1:d]
        
        fkTR[, , k, t] <- sigmak02[k] * ( thisU %*% diag(thisLambda) %*% 
                                              t(Conj(thisU)) + diag(P) )
        U_kl0[, , k, t] <- thisU
        Lambdakl0[, k, t] <- thisLambda
    }
}
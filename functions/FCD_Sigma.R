
# functions ---------------------------------------------------------------

logdet <- function(X) {
    return( log(EigenR::Eigen_det(X)) )
}

# FCD sampling function ---------------------------------------------------

Sigma_gibbs_densCovar <- 
    function(data_list, param_list, n_Sig, useEigenR = FALSE, 
                byCholesky = FALSE) {
    # extract stuff from param_list
    P <- param_list$P
    d <- param_list$d
    K <- param_list$K
    Sigmas <- param_list$Sigmas
    U_ks <- param_list$U_ks
    
    invSigmas <- solve(Sigmas)
    
    # make proposal
    Sigmap <- rFTCW(Sigmas, n_Sig, P, useEigenR, byCholesky)
    invSigmap <- solve(Sigmap)
        
    # sumlog_p
    sumlog_p <- 0
    sumlog_s <- 0
    for (k in 1:K) {
        sumlog_p <- sumlog_p + 
            logdet( t(Conj(U_ks[, , k])) %*% invSigmap %*% U_ks[, , k] )
        sumlog_s <- sumlog_s + 
            logdet( t(Conj(U_ks[, , k])) %*% invSigmas %*% U_ks[, , k] )
    }
    
    logdens_num <- -d*K* logdet(Sigmap) - P * sumlog_p - 
        n_Sig * logdet(Sigmap) +
        (n_Sig - P) * logdet(Sigmas) - 
        P*n_Sig * log( Re(sum(diag( invSigmap %*% Sigmas))) ) 
        
    logdens_den <- -d*K* logdet(Sigmas) - P * sumlog_s - 
        n_Sig * logdet(Sigmas) +
        (n_Sig - P) * logdet(Sigmap) - 
        P*n_Sig * log( Re(sum(diag( invSigmas %*% Sigmap ))) ) 
    
    logr <- Re(logdens_num - logdens_den)
    
    if (log(runif(1)) <= logr) {
        Sigmas <- Sigmap
        accCount <- TRUE
    } else {
        accCount <- FALSE
    }
    
    return(list(Sigmas = Sigmas, accCount_Sigma = accCount))
}

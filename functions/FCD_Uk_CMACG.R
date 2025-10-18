# sigmak2_gibbs_densCovar <- function(data_list, param_list) {

# functions ---------------------------------------------------------------

# generate a proposal via a small rotation parameterized with Cayley transform
Uk_MH_Cayley <- function(Us, tau_U, doCayleyZeros, CayleyZeroProb) {
    P <- nrow(Us)
    
    Sp <- matrix(0, P, P)
    Sp[upper.tri(Sp)] <- tau_U*rscnorm(choose(P, 2))

    if (doCayleyZeros) {
        Sp[upper.tri(Sp)] <- 
            Sp[upper.tri(Sp)] * rbinom(choose(P, 2), 1, 1 - CayleyZeroProb)
    }

    for(j in 1:(P-1)) {
        for (i in (j+1):P) {
            Sp[i, j] <- -Conj(Sp[j, i])
        }
    }
    Rp <- (diag(P) - Sp) %*% solve(diag(P) + Sp)
    Up <- Rp %*% Us
    
    return(Up)
}

# FCD sampling function ---------------------------------------------------

Uk_CMACG_gibbs_densCovar <- function(data_list, param_list, tau_Uk,
    doCayleyZeros = FALSE, CayleyZeroProb = 0.5) {

    # extract stuff from param_list
    P <- param_list$P
    d <- param_list$d
    K <- param_list$K


    if (length(tau_Uk) == 1) {
        tau_Uk <- rep(tau_Uk, K)
    }

    Sigmas <- param_list$Sigmas
    invSigmas <- solve(Sigmas)

    accCount <- rep(TRUE, K)
    newU_ks <- array(NA, c(P, d, K))

    # for each k
    for (k in 1:K) {
        Y_k <- data_list[[k]]
        Lambda_k <- diag(param_list$Lambda_ks[, k])
        Omega_k <- solve( solve(Lambda_k) + diag(d) )
        sigma_k2 <- param_list$sigma_k2s[k]
        
        Us <- param_list$U_ks[, , k]
        
        # propose U by small rotation of Us
        Up <- Uk_MH_Cayley(Us, tau_Uk[k], doCayleyZeros, CayleyZeroProb)
        
        # calculate terms based on priors
        tracep <- Re(sum(diag( Up %*% Omega_k %*% t(Conj(Up)) %*% Y_k )))
        logdetp <- log(Re(EigenR::Eigen_det( t(Conj(Up)) %*% invSigmas %*% Up )))
        traces <- Re(sum(diag( Us %*% Omega_k %*% t(Conj(Us)) %*% Y_k )))
        logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*% invSigmas %*% Us )))
        
        # calculate acceptance ratio
        logr <- tracep/sigma_k2 - P*logdetp - traces/sigma_k2 + P*logdets
        
        # accept or reject
        if (log(runif(1)) <= logr) {
            newU_ks[, , k] <- Up
        } else {
            newU_ks[, , k] <- Us
            accCount[k] <- FALSE
        }
    } # end of sampling for each k
    
    return(list(U_ks = newU_ks, accCount_U_k = accCount))
}

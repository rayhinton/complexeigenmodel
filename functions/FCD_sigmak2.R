# sigmak2 FCD samplers


##### function

sigmak2_gibbs_densCovar <- function(data_list, param_list) {
    result_sigmas <- vector(mode = "numeric", K)
    
    # PxP identity matrix
    IP <- diag(param_list$P)
    
    for (k in 1:K) {
        Pk <- data_list[[k]]
        
        # choose the most recent value of Uk
        Uk <- param_list$U_ks[, , k]
        # choose the most recent value of Lambda_k
        Lambdak <- diag(param_list$Lambda_ks[, k])
        # temporary matrix to ease calculations
        Mk <- Uk %*% Lambdak %*% t(Conj(Uk)) + IP
        
        # calculate parameters for the Inverse Gamma distribution
        # ak <- P * (nk[k] - 1)
        # ak <- P * (param_list$n_k[k])
        ak <- .ak_sigmak2_densCovar(P, param_list$n_k[k])
        # bk <- (solve(Mk) %*% Pk) |> diag() |> sum() |> Re()
        bk <- .bk_sigmak2_densCovar(Mk, Pk)
        
        # FCD is Inverse Gamma: 1 / Gamma, in terms of the rate parameter
        result_sigmas[k] <- 1/rgamma(1, ak, rate = bk)    
    }
    
    result_sigmas
}

.ak_sigmak2_densCovar <- function (P, nk) {
    return(P * nk)
}

.bk_sigmak2_densCovar <- function(Mk, Pk) {
    return(
        (solve(Mk) %*% Pk) |> diag() |> sum() |> Re()
    )
}
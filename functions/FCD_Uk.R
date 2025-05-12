Uk_gibbs_densCovar <- function(data_list, param_list) {
    VAVH <- param_list$Vs %*% diag(param_list$As) %*% t(Conj(param_list$Vs))
    
    U_ks <- param_list$U_ks
    
    for (k in 1:param_list$K) {
        # Ukc <- Uks[, , k, s-1]
        Ukc <- param_list$U_ks[, , k]
        
        for (j in sample(param_list$d)) {
            # bj <- B[j, j]
            bj <- param_list$Bs[j]
            # lambdajk <- Lambdaks[j, j, k]
            lambdajk <- param_list$Lambda_ks[j, k]
            
            # TODO VAV^H could be calculated beforehand - the value does not change
            # parj <- bj * VAVH + lambdajk/(1 + lambdajk) * Yks[, , k]
            parj <- bj * VAVH + lambdajk/(1 + lambdajk) * data_list[[k]]
            
            Ukc[, j] <- rcvb_LN(Ukc, j, parj)
        }
        
        U_ks[, , k] <- Ukc
    }
    
    return(U_ks)
}

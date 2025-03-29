# function to sample all j values, for each group k
Lambdak_gibbs_densCovar <- function(data_list, param_list) {
    d <- param_list$d
    K <- param_list$K
    
    # an array to store the sampled values
    result_Lambdas <- array(NA, c(d, K))
    
    # for each group, k
    for (k in 1:K) {
        # extract parameters and data from input lists
        Y <- data_list[[k]]
        n <- param_list$n_k[k]
        U <- param_list$U_ks[, , k]
        sigma2 <- param_list$sigma_k2s[k]
        
        # for each Lambdak value, j
        for (j in 1:d) {
            # xi_jk has a truncated Gamma(nk, tjk) distribution (shape, rate), 
            # and tjk is a temporary parameter defined as follows.
            # (note we take the Real part, since the quadratic form should be real, 
            # but may have a small complex part due to numerical issues.)
            tjk <- Re( t(Conj(U[, j])) %*% Y %*% U[, j] ) / sigma2
            
            ##### sample from a truncated Gamma distribution
            # lower and upper bounds for truncating the Gamma distribution
            lb <- 0
            ub <- 1
            
            # upper and lower percentile bounds
            lp <- pgamma(lb, shape = n, rate = tjk)
            up <- pgamma(ub, shape = n, rate = tjk)
            
            # generate the random sample
            u <- runif(1, lp, up)
            # invert the Uniform r.v. to a truncated Gamma r.v.
            # TODO if the data CW df is n-1, then shape should be n
            # if the data CW df is n, then shape should be n+1
            xi_jk <- qgamma(u, shape = n+1, rate = tjk)
            # convert xi to Lambda
            result_Lambdas[j, k] <- 1/xi_jk - 1
        }
    }
    
    return(result_Lambdas)
}
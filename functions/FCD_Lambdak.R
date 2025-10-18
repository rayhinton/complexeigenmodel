# function to sample all j values, for each group k
Lambdak_gibbs_densCovar <- function(data_list, param_list, doOrdered = FALSE) {
    d <- param_list$d
    K <- param_list$K
    
    # an array to store the sampled values
    #result_Lambdas <- array(NA, c(d, K))
    result_Lambdas <- param_list$Lambda_ks
    
    # for each group, k
    for (k in 1:K) {
        # extract parameters and data from input lists
        Y <- data_list[[k]]
        n <- param_list$n_k[k]
        U <- param_list$U_ks[, , k]
        sigma2 <- param_list$sigma_k2s[k]
        
        # for each Lambdak value, j
        for (j in sample(d)) {
            # xi_jk has a truncated Gamma(nk, tjk) distribution (shape, rate), 
            # and tjk is a temporary parameter defined as follows.
            # (note we take the Real part, since the quadratic form should be real, 
            # but may have a small complex part due to numerical issues.)
            tjk <- Re( t(Conj(U[, j])) %*% Y %*% U[, j] ) / sigma2
            
            ##### sample from a truncated Gamma distribution
            # lower and upper bounds for truncating the Gamma distribution
            if (!doOrdered) {
                lb <- 0
                ub <- 1
            } else {
                if (j == 1) lb <- 0 else lb <- 1 / (result_Lambdas[j-1, k] + 1)
                if (j == d) ub <- 1 else ub <- 1 / (result_Lambdas[j+1, k] + 1)
            }
            
            # upper and lower percentile bounds
            # TODO change shape to n+1 or n, depending on TODO note below
            #lp <- pgamma(lb, shape = n+1, rate = tjk)
            #up <- pgamma(ub, shape = n+1, rate = tjk)
            
            #print(c(lb, ub, lp, up, n+1, tjk))
            
            #print(c(lp, up))
            # generate the random sample
            #u <- runif(1, lp, up)
            # invert the Uniform r.v. to a truncated Gamma r.v.
            # TODO if the data CW df is n-1, then shape should be n
            # if the data CW df is n, then shape should be n+1
            # print(c(u, n+1, tjk))
            #xi_jk <- qgamma(u, shape = n+1, rate = tjk)
            
            #print(c(lb, ub, n+1, tjk))
            
            #mlog_lp <- -pgamma(ub, shape = n+1, rate = tjk, log.p = TRUE)
            #mlog_up <- -pgamma(lb, shape = n+1, rate = tjk, log.p = TRUE)
            
            #logu <- -truncdist::rtrunc(1, "exp", mlog_lp, mlog_up)
            #xi_jk <- qgamma(logu, shape = n+1, rate = tjk, log.p = TRUE)
            
            distr <- Runuran::udgamma(shape = n+1, scale = 1/tjk, lb, ub)
            #gen <- Runuran::pinvd.new(distr)
            gen <- Runuran::arsd.new(distr)
            xi_jk <- Runuran::ur(gen, 1)
            
            #xi_jk <- Runuran::urgamma(1, n+1, 1/tjk, lb, ub)
            
            # convert xi to Lambda
            result_Lambdas[j, k] <- 1/xi_jk - 1
        }
    }
    
    return(result_Lambdas)
}

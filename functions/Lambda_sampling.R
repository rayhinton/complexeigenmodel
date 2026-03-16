# sample_lambda_first

# if w_ss is not specified, then use the tuned values 
sample_Lambda_first <- function(k, data_list_w, thisUkl, thisSigmak2, 
                                thisLambda, d, LL, w_ss = NULL, m_ss = Inf) {
    ### Lambda: 1st frequency
    return_Lambda <- thisLambda[, 1]
    LSkw <- data_list_w[[1]][[k]]
    Ukl <- thisUkl[, , 1]
    
    # for each Lambdak value, j
    for (j in sample(d)) {
        # xi_jk has a truncated Gamma(nk, tjk) distribution (shape, rate), 
        # and tjk is a temporary parameter defined as follows.
        # (take the Real part, since the quadratic form should be real, 
        # but may have a small complex part due to numerical issues.)
        tjk <- Re( t(Conj(Ukl[, j])) %*% LSkw %*% Ukl[, j] ) / 
            thisSigmak2
        
        ##### sample from a truncated Gamma distribution
        # lower and upper bounds for truncating the Gamma distribution
        if (j == 1) lb <- 0 else lb <- 1 / (return_Lambda[j-1] + 1)
        if (j == d) ub <- 1 else ub <- 1 / (return_Lambda[j+1] + 1)
        
        if (is.null(w_ss)) {
            w_ss <- w_ss_tune[j, k, 1]
        }
        # by slice sampling, instead
        xi_jk <- uni.slice(1/(1 + return_Lambda[j]),
                           dgamma,
                           w = w_ss, m = m_ss, lower = lb, upper = ub,
                           shape = LL+1, scale = 1/tjk,
                           log = TRUE)
        
        # convert xi to Lambda
        # thisLambda[j, 1] <- 1/xi_jk - 1
        return_Lambda[j] <- 1/xi_jk - 1
    } # end of l = 1 Lambda sampling

    return(return_Lambda)
}

sample_Lambda_interior <- function(k, data_list_w, thisUkl, thisSigmak2, 
                                    thisLambda, thisTaujkl2, Lambda_prior,
                                    d, LL, num_freqs, unnorm_logPDF,
                                    w_ss = NULL, m_ss = Inf) {
    ### Lambda: frequencies in 2 to num_freqs-1
    if (Lambda_prior %in% c("1RW", "2RWPN")) {
        return_Lambda <- thisLambda
        
        for ( l in sample(2:(num_freqs-1)) ) {
            Ukw <- thisUkl[, , l]
            Dkw1 <- data_list_w[[l]][[k]]
            
            for (j in sample(d)) {
                ajk <- t(Conj(Ukw[, j])) %*% Dkw1 %*% Ukw[, j] / 
                    thisSigmak2
                # should be strictly real
                ajk <- Re(ajk)
                
                if (Lambda_prior == "2RWPN") {
                    # 2nd order RW, previous and next neighbors
                    lp <- return_Lambda[j, l-1]
                    ln <- return_Lambda[j, l+1]
                    mu <- (lp + ln)/2
                }
                
                else if (Lambda_prior == "1RW") {
                    # 1st order RW
                    mu <- return_Lambda[j, l-1]
                }
                
                # for single smoothing parameter
                varpar_jkl <- thisTaujkl2[j, l-1]
                
                if (is.null(w_ss)) {
                    w_ss <- w_ss_tune[j, k, l]
                }
                # by slice sampling, instead
                newdraw <- uni.slice(return_Lambda[j, l],
                                     unnorm_logPDF,
                                     w = w_ss, m = m_ss, lower = 0, 
                                     upper = Inf, 
                                     tau2 = varpar_jkl,
                                     ajk = ajk, mu = mu, N = LL, 
                                     logscale = TRUE)    
                # result_Lambdas[j, k, l] <- newdraw
                return_Lambda[j, l] <- newdraw
                
            } # end of sampling j in 1:d
        } # end of sampling over frequencies
    } # end of the conditional to check the type of prior for Lambda
    
    return(return_Lambda)
}

sample_Lambda_last <- function(k, data_list_w, thisUkl, thisSigmak2, thisLambda,
                          thisZetajk2, d, LL, num_freqs, unnorm_logPDF,
                          w_ss = NULL, m_ss = Inf) {
    ### Lambda: last frequency
    return_Lambda <- thisLambda[, num_freqs]
    Uw1 <- thisUkl[, , num_freqs]
    Dkw1 <- data_list_w[[num_freqs]][[k]]
    
    for (j in sample(d)) {
        ajk <- t(Conj(Uw1[, j])) %*% Dkw1 %*% Uw1[, j] / 
            thisSigmak2
        # should be strictly real - this is done in the existing Lambda
        # FCD sampler.
        ajk <- Re(ajk)
        
        # 1st order random walk
        mu <- thisLambda[j, num_freqs - 1]
        
        # by slice sampling, instead
        if (is.null(w_ss)) {
            w_ss <- w_ss_tune[j, k, num_freqs]
        }
        newdraw <- uni.slice(return_Lambda[j], 
                             unnorm_logPDF,
                             w = w_ss, m = m_ss, lower = 0, upper = Inf, 
                             tau2 = thisZetajk2[j],
                             ajk = ajk, mu = mu, N = LL, 
                             logscale = TRUE)    
        return_Lambda[j] <- newdraw
        
    } # end of sampling j in 1:d
    
    return(return_Lambda)
}

sample_taujkl2 <- function(thisLambda, tau2_a, tau2_b, Lambda_prior, d, num_freqs) {
    if (Lambda_prior == "2RWPN") {
        return_taujkl2 <- matrix(0, d, num_freqs - 2)
    } else {
        stop(paste0("No sampler implemented for taujkl when Lambda_prior is ", Lambda_prior))
    }
    
    for (j in 1:d) {
        if (Lambda_prior == "2RWPN") {
            # 2nd order RW, next and previous neighbors
            ### BEGIN STANDARD SAMPLING
            sumalljk <- (
                (thisLambda[j, 2:(num_freqs-1)] -
                     .5*thisLambda[j, 1:(num_freqs-2)] -
                     .5*thisLambda[j, 3:(num_freqs)])^2)
            
            return_taujkl2[j, ] <-
                1/rgamma(length(sumalljk), 
                         tau2_a + .5,
                         rate = tau2_b + .5*sumalljk)
            ### END STANDARD SAMPLING
        } else {
            stop(paste0("No sampler implemented for taujkl when Lambda_prior is ", Lambda_prior))
        }
    }
    
    return(return_taujkl2)
}

sample_zetajk2 <- function(thisLambda, tau2_a, tau2_b, d, num_freqs) {
    return_zetajk2 <- rep(0, d)
    
    for (j in 1:d) {
        sum_zeta <- (thisLambda[j, num_freqs] - 
                 thisLambda[j, num_freqs-1])^2
        return_zetajk2[j] <- 1/rgamma(1, tau2_a + .5,
                                       rate = tau2_b + .5*sum_zeta)
    }

    return(return_zetajk2)
}

# B-spline prior-based sampling -------------------------------------------

setup_bspline <- function(num_freqs, n_knots = 25, degree = 3, 
                          ridge_eps = 1e-6) {
    # Frequencies: omega_l for l = 1, ..., num_freqs
    # We don't need the actual omega values, just indices 1:num_freqs 
    # as the "covariate" for the spline basis.
    freq_indices <- 1:num_freqs
    
    # Place knots: equally spaced, with boundary knots at the endpoints.
    # splineDesign needs (n_knots + degree + 1) knots total, 
    # including (degree + 1) repeated boundary knots on each side.
    n_inner_knots <- n_knots
    inner_knots <- seq(from = 1, to = num_freqs, length.out = n_inner_knots + 2)
    # drop the endpoints — they become boundary knots
    inner_knots <- inner_knots[-c(1, n_inner_knots + 2)]
    
    boundary_knots <- c(1, num_freqs)
    
    # Full knot vector with repeated boundary knots
    knots <- c(rep(boundary_knots[1], degree + 1),
               inner_knots,
               rep(boundary_knots[2], degree + 1))
    
    # Number of basis functions: n_inner_knots + degree + 1
    M <- n_inner_knots + degree + 1
    
    # Design matrix: N_omega x M
    B <- splineDesign(knots, x = freq_indices, ord = degree + 1, outer.ok = FALSE)
    
    # Second-difference matrix: (M-2) x M
    D2 <- diff(diag(M), differences = 2)
    
    # Penalty matrix: M x M, rank M-2
    P <- t(D2) %*% D2
    
    # Optionally add small ridge for numerical stability
    P_ridge <- P + ridge_eps * diag(M)
    
    # Identify the pinning index: which basis function has the largest 
    # value at the first frequency (l = 1)?
    b1 <- B[1, ]  # basis vector at first frequency
    pin_idx <- which.max(abs(b1))
    free_idx <- setdiff(1:M, pin_idx)
    
    return(list(
        B = B,               # N_omega x M design matrix
        P = P,               # M x M penalty matrix (rank M-2)
        P_ridge = P_ridge,   # M x M penalty matrix with ridge
        D2 = D2,             # (M-2) x M second-difference matrix
        M = M,               # number of basis functions
        knots = knots,       # full knot vector
        degree = degree,
        num_freqs = num_freqs,
        pin_idx = pin_idx,   # index of pinned coefficient
        free_idx = free_idx, # indices of free coefficients
        b1 = b1,             # basis vector at first frequency
        ridge_eps = ridge_eps
    ))
}

reconstruct_beta <- function(beta_free, lambda_first_j, bspline_setup) {
    # Given free coefficients and the pinned lambda value at freq 1,
    # reconstruct the full beta vector.
    # 
    # Constraint: b1^T beta = log(lambda_{j,k,1})
    # => beta[pin] = (log(lambda_first_j) - sum(b1[free] * beta_free)) / b1[pin]
    
    b1 <- bspline_setup$b1
    pin <- bspline_setup$pin_idx
    free <- bspline_setup$free_idx
    
    beta_full <- numeric(bspline_setup$M)
    beta_full[free] <- beta_free
    beta_full[pin] <- (log(lambda_first_j) - sum(b1[free] * beta_free)) / b1[pin]
    
    return(beta_full)
}

beta_to_lambda <- function(beta_full, bspline_setup) {
    # Convert full beta vector to lambda values at all frequencies.
    # lambda_l = exp(B[l, ] %*% beta)
    as.numeric(exp(bspline_setup$B %*% beta_full))
}

elliptical_slice_sample <- function(current, log_lik_fn, prior_chol, ..., 
                                     max_iter = 200) {
    # Elliptical slice sampling (Murray, Adams & MacKay 2010).
    # 
    # current:    current state vector (length M-1)
    # log_lik_fn: function(x, ...) returning log-likelihood at x
    # prior_chol: Cholesky factor (upper triangular, from chol()) of the 
    #             prior covariance, so that prior_chol^T %*% prior_chol = Sigma_prior
    #             Draw from prior via: prior_chol^T %*% z, z ~ N(0, I)
    # ...:        additional arguments passed to log_lik_fn
    #
    # Returns: new state vector
    
    n <- length(current)
    
    # Draw from the prior: nu ~ N(0, Sigma_prior)
    nu <- as.numeric(t(prior_chol) %*% rnorm(n))
    
    # Log-likelihood threshold
    log_y <- log_lik_fn(current, ...) + log(runif(1))
    
    # Draw initial angle and define bracket
    theta <- runif(1, 0, 2 * pi)
    theta_min <- theta - 2 * pi
    theta_max <- theta
    
    for (iter in 1:max_iter) {
        # Propose on the ellipse
        proposal <- current * cos(theta) + nu * sin(theta)
        
        if (log_lik_fn(proposal, ...) > log_y) {
            return(proposal)
        }
        
        # Shrink the bracket
        if (theta < 0) {
            theta_min <- theta
        } else {
            theta_max <- theta
        }
        theta <- runif(1, theta_min, theta_max)
    }
    
    # If we exhaust iterations, return current (should be rare)
    warning("ESS reached max_iter without accepting; returning current state")
    return(current)
}

setup_conditional_prior <- function(bspline_setup, use_ridge = TRUE) {
    # One-time precomputation of quantities that depend only on the 
    # B-spline basis and penalty, NOT on tau2.
    # Call once, then use scale_conditional_prior() each iteration.
    
    pin <- bspline_setup$pin_idx
    free <- bspline_setup$free_idx
    b1 <- bspline_setup$b1
    M <- bspline_setup$M
    
    if (use_ridge) {
        P <- bspline_setup$P_ridge
    } else {
        P <- bspline_setup$P
    }
    
    # Substitution matrix: full beta = A %*% beta_free + v * log(lambda_first)
    b1_p <- b1[pin]
    b1_f <- b1[free]
    
    A <- matrix(0, M, M - 1)
    A[free, ] <- diag(M - 1)
    A[pin, ] <- -b1_f / b1_p
    
    v <- numeric(M)
    v[pin] <- 1 / b1_p
    
    AtPA <- t(A) %*% P %*% A        # (M-1) x (M-1), fixed
    Atpv <- as.numeric(t(A) %*% P %*% v)  # (M-1) vector, fixed
    vtPv <- as.numeric(t(v) %*% P %*% v)  # scalar, fixed
    
    # Precompute: solve(AtPA) and chol(solve(AtPA)) — the expensive parts
    AtPA_inv <- solve(AtPA)          # (M-1) x (M-1)
    AtPA_inv_chol <- chol(AtPA_inv)  # upper triangular: t(R) %*% R = AtPA_inv
    
    # mean_coeff (does not depend on tau2):
    # mean_free = -AtPA_inv %*% Atpv * log(lambda_first)
    mean_coeff <- -as.numeric(AtPA_inv %*% Atpv)
    
    return(list(
        A = A, v = v,
        AtPA = AtPA, Atpv = Atpv, vtPv = vtPv,
        AtPA_inv = AtPA_inv,
        AtPA_inv_chol = AtPA_inv_chol,
        mean_coeff = mean_coeff
    ))
}

scale_conditional_prior <- function(tau2_jk, cond_prior_setup) {
    # Per-iteration scaling: only O(M^2) multiplications.
    # cov_free = tau2 * AtPA_inv
    # cov_chol = sqrt(tau2) * AtPA_inv_chol
    sqrt_tau2 <- sqrt(tau2_jk)
    
    list(
        cov_chol = sqrt_tau2 * cond_prior_setup$AtPA_inv_chol,
        mean_coeff = cond_prior_setup$mean_coeff
    )
}

sample_Lambda_bspline <- function(k, data_list_w, thisUkl, thisSigmak2,
                                  thisLambda, thisBetaFree, new_tau2,
                                  bspline_setup, cond_prior_setup,
                                  d, LL, num_freqs) {
    # Sample Lambda curves at frequencies 2:num_freqs using B-spline + ESS,
    # for one observation k.
    #
    # thisLambda:      d x num_freqs matrix (column 1 is the pinned first frequency)
    # thisBetaFree:    list of length d, each element is a vector of length M-1
    # new_tau2:        vector of length d, current smoothing variances
    # cond_prior_setup: output of setup_conditional_prior() (precomputed once)
    #
    # Returns: list with updated Lambda, BetaFree, and BetaFull
    
    B <- bspline_setup$B
    return_Lambda <- thisLambda
    return_BetaFree <- thisBetaFree
    return_BetaFull <- vector("list", d)
    
    # Precompute a_{j,l} for all j and l (once per call)
    precomp_ajl <- matrix(NA, d, num_freqs)
    for (l in 1:num_freqs) {
        Ukw <- thisUkl[, , l]
        Dkw1 <- data_list_w[[l]][[k]]
        for (j in 1:d) {
            precomp_ajl[j, l] <- Re(t(Conj(Ukw[, j])) %*% Dkw1 %*% Ukw[, j]) / 
                thisSigmak2
        }
    }
    
    # Scale conditional priors by current tau2
    cond_priors <- vector("list", d)
    for (j in 1:d) {
        cond_priors[[j]] <- scale_conditional_prior(new_tau2[j], cond_prior_setup)
    }
    
    for (j in sample(d)) {
        lambda_first_j <- return_Lambda[j, 1]
        prior_mean_j <- cond_priors[[j]]$mean_coeff * log(lambda_first_j)
        ajl_j <- precomp_ajl[j, 2:num_freqs]
        
        # Log-likelihood: vectorized over frequencies 2:num_freqs
        log_lik_j <- function(x_centered) {
            bf <- x_centered + prior_mean_j
            beta_full <- reconstruct_beta(bf, lambda_first_j, bspline_setup)
            lambda_vals <- exp(as.numeric(B %*% beta_full))
            lv <- lambda_vals[2:num_freqs]
            sum(-LL * log(1 + lv) - ajl_j / (1 + lv))
        }
        
        current_centered <- return_BetaFree[[j]] - prior_mean_j
        new_centered <- elliptical_slice_sample(
            current = current_centered,
            log_lik_fn = log_lik_j,
            prior_chol = cond_priors[[j]]$cov_chol
        )
        
        new_beta_free <- new_centered + prior_mean_j
        beta_full <- reconstruct_beta(new_beta_free, lambda_first_j, bspline_setup)
        
        return_BetaFree[[j]] <- new_beta_free
        return_BetaFull[[j]] <- beta_full
        return_Lambda[j, ] <- exp(as.numeric(B %*% beta_full))
        return_Lambda[j, 1] <- lambda_first_j  # preserve exact pinned value
    }
    
    # Fill in BetaFull for any j not yet updated (if d > 1 and sample order matters)
    for (j in 1:d) {
        if (is.null(return_BetaFull[[j]])) {
            return_BetaFull[[j]] <- reconstruct_beta(return_BetaFree[[j]],
                                                      return_Lambda[j, 1],
                                                      bspline_setup)
        }
    }
    
    return(list(Lambda = return_Lambda, BetaFree = return_BetaFree, 
                BetaFull = return_BetaFull))
}

sample_tau2_bspline <- function(thisBetaFull, bspline_setup, tau2_a, tau2_b, d,
                                use_ridge = TRUE) {
    # Conjugate IG update for tau^2_{j,k} under the P-spline prior.
    #
    # thisBetaFull: list of length d, each element is full beta vector (length M)
    # tau2_a, tau2_b: IG(a, b) prior hyperparameters
    #
    # FCD: tau2 | beta ~ IG(a + rank(P)/2, b + beta^T P beta / 2)
    #
    # Returns: vector of length d with updated tau2 values
    
    if (use_ridge) {
        P <- bspline_setup$P_ridge
        P_rank <- bspline_setup$M  # full rank with ridge
    } else {
        P <- bspline_setup$P
        P_rank <- bspline_setup$M - 2
    }
    
    return_tau2 <- numeric(d)
    for (j in 1:d) {
        beta_j <- thisBetaFull[[j]]
        quad <- as.numeric(t(beta_j) %*% P %*% beta_j)
        return_tau2[j] <- 1 / rgamma(1, 
                                       shape = tau2_a + P_rank / 2,
                                       rate = tau2_b + quad / 2)
    }
    
    return(return_tau2)
}

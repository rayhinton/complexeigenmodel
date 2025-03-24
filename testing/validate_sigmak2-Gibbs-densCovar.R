# validate the sigmak2 Gibbs sampler for dense covariance matrix model

library(cmvnorm)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_sigmak2.R")

# dense covariance model

# depends on:
# - P, size of the covariance matrix (eventually, SDM)
# - n_k, degrees of freedom in CW
#   * number of rows in data matrix k (not sure about the equivalence for time series - maybe, the length of the time series?)
# - P_k, PxP, data matrices (either the estimated covariance or SDM matrix)
# - U_k, Pxd, eigenvector matrices
# - Lambda_k, eigenvalues

# inputs
# - data - list of observed data matrices (estimated covariances or SDMs)
# - params - list of currently sampled parameter values

# return
# - a vector of sampled sigmak2 values

demo_FCD <- function(Lambda_ks, P = 12, d = 4, K = 3, nk_scale = 1000, 
                     gibbs_its = 5000) {
    
    n_k <- nk_scale + 5*(1:K)
    # number of Gibbs iterations, for testing
    gibbs_its <- 5000
    
    data_list <- list()
    
    U_ks <- array(NA, c(P, d, K))

    sigma_k20 <- seq(5, 20, length.out = K)
    
    for (k in 1:K) {
        U_ks[, , k] <- runitary(P, d)
        
        # temp matrix
        M_k <- diag(P) + 
            U_ks[, , k] %*% diag(Lambda_ks[, k]) %*% t(Conj(U_ks[, , k]))
        M_k <- round(M_k, 9)
        print(paste0("Summary of eigenvalues for M_",k," matrix"))
        print(summary(eigen(M_k)$values))
        
        data_list[[k]] <- rcwis(n_k[k], sigma_k20[k] * M_k)
    }
    
    # stand-in for sigma_k2, the parameters to be sampled
    sigma_k2s <- vector("numeric", K)
    
    ##### Parameters list
    param_list <- list(P = P,
                       d = d,
                       K = K,
                       n_k = n_k,
                       U_ks = U_ks,
                       Lambda_ks = Lambda_ks,
                       sigma_k2s = sigma_k2s)
    
    # test function
    sigmak2_gibbs_densCovar(data_list, param_list)
    
    # test the function within a sampling loop
    # initialize an array to track the sigma_k2 values
    sigma_k2_S <- array(NA, c(K, gibbs_its))
    
    for (s in 1:gibbs_its) {
        sigma_k2s <- sigmak2_gibbs_densCovar(data_list, param_list)
        param_list$sigma_k2s <- sigma_k2s
        sigma_k2_S[, s] <- sigma_k2s
    }
    
    # in this case, we know that the distribution should be certain Inv Gammas
    
    # Exact posterior mean
    print("Exact posterior means")
    for (k in 1:K) {
        # temp matrix
        M_k <- diag(P) + 
            U_ks[, , k] %*% diag(Lambda_ks[, k]) %*% t(Conj(U_ks[, , k]))
        print(Re(sum(diag(solve(M_k) %*% data_list[[k]]))) / (P * n_k[k] - 1))
    }
    
    post_samples <- round(gibbs_its/2):gibbs_its
    
    # Sample posterior mean
    print("Sample posterior means")
    print(rowMeans(sigma_k2_S[, post_samples]))
    # Known sigmak2 values
    print("Known sigma_k2 values")
    print(sigma_k20)
    
    # plot observed and true densities
    for (k in 1:K) {
        xmin <- min(sigma_k2_S[k, post_samples]) * .8
        xmax <- max(sigma_k2_S[k, post_samples]) * 1.2
        plot(density(sigma_k2_S[k, post_samples]), type = "l", 
             main = paste0("Observed and true posterior densities, k = ",k))
        
        M_k <- diag(P) + 
            U_ks[, , k] %*% diag(Lambda_ks[, k]) %*% t(Conj(U_ks[, , k]))
        
        ak <- P*n_k[k]
        bk <- Re(sum(diag(solve(M_k) %*% data_list[[k]])))
        # TODO note the transformation required to plot inverse gamma, i.e. /x^2,
        curve(dgamma(1/x, shape = ak, rate = bk) / x^2,
          from = xmin, to = xmax, n = 501, add = TRUE,
          col = "red", ylab = "density")
    }
}

P = 12
d = 4
K = 3
nk_scale = 1000
gibbs_its = 5000

# test cases:

# Lambda_ks[, k] <- rgamma(d, 1, 1) |> sort(decreasing = TRUE)
# Lambda_ks[, k] <- seq(5, 10^(k/2), length.out = d) |> sort(decreasing = TRUE)

# expected results:
# - for all, the sample densities match the known FCD densities
# - with the random rgammas - the posterior expectations match the true sigmas

# expected issues:
# - with the sequences - the posterior expectations do not match the true sigmas

set.seed(13032025)
Lambda_ks <- array(NA, c(d, K))
# different Lambda values (eigenvalues) lead to different estimator quality
for(k in 1:K) {
    Lambda_ks[, k] <- rgamma(d, 1, 1) |> sort(decreasing = TRUE)
    # Lambda_ks[, k] <- seq(5, 10^(k/2), length.out = d) |> sort(decreasing = TRUE)
    # Lambda_ks[, k] <- seq(5, 10^(k), length.out = d) |> sort(decreasing = TRUE)
}

demo_FCD(Lambda_ks)


# Inverse gamma densities, CDFs, and quantiles ----------------------------

dinvgamma <- function(x, shape, scale) {
    return(dgamma(1/x, shape = shape, rate = scale) / x^2)
}

curve(dinvgamma(x, 5, 25),
      from = 0, to = 20, n = 501, ylab = "density")

q975 <- 1/qgamma(.025, shape = 5, rate = 25)
q025 <- 1/qgamma(.975, shape = 5, rate = 25)

abline(v = q975)
abline(v = q025)

integrate(dinvgamma, 0, q975, shape = 5, scale = 25)
integrate(dinvgamma, 0, q025, shape = 5, scale = 25)

# Choose different M_k values ---------------------------------------------

demo_random_post <- function(Mk, sigma_02, n, Mk_df, ex_its, hpd_round = 11) {

    P <- ncol(Mk)
    print("Summary of eigenvalues for M_k matrix: ")
    print(summary(eigen(Mk)$values))
    
    G0 <- round(sigma_02 * Mk, hpd_round)
    # isHermitian(G0)
    
    post_means <- rep(NA, ex_its)
    post_CIs <- matrix(NA, ex_its, 2)
    
    for (i in 1:ex_its) {
      Pk <- rcwis(n, G0)
      
      ak <- .ak_sigmak2_densCovar(P, n)
      bk <- .bk_sigmak2_densCovar(Mk, Pk)
      
      post_means[i] <- bk / (ak-1)
      post_CIs[i, 1] <- 1/qgamma(.975, shape = ak, rate = bk)
      post_CIs[i, 2] <- 1/qgamma(.025, shape = ak, rate = bk)
    }
    
    print(
      paste0("True sigma2 is within 95% cred. int. ", 
             round(100*mean((post_CIs[, 1] <= sigma_02) & (sigma_02 <= post_CIs[, 2])), 2),
             "% of the time"))
    
    print(paste0("Median .025 posterior quantiles: ",
             median(post_CIs[, 1]) |> round(3)))
    print(paste0("Median .975 posterior quantiles: ",
             median(post_CIs[, 2]) |> round(3)))
    print(paste0("Mean of posterior means: ",
             mean(post_means) |> round(3)))
    print(paste0("Median of posterior means: ",
             median(post_means) |> round(3)))
    
    plot(density(post_means),
         main = "Dens. of posterior means")
    abline(v = sigma_02, col = "red")
}

n <- 1000
sigma_02 <- 5
ex_its <- 1000
Mk_df <- 1000
hpd_round <- 11

set.seed(14032025)
Mk <- rcwis(Mk_df, diag(P))

Mk <- eigen(Mk)$vectors %*% 
  # Mk has eigenvalues of similar scale
  diag(seq(30, 20, length.out = P)) %*%
  # Mk has eigenvalues of different scales
  # diag(seq(1000, 1, length.out = P)) %*%
  t(Conj(eigen(Mk)$vectors))

demo_random_post(Mk, sigma_02, n, Mk_df, ex_its, hpd_round = hpd_round)

set.seed(14032025)
Mk <- rcwis(Mk_df, diag(P))

Mk <- eigen(Mk)$vectors %*% 
  # Mk has eigenvalues of similar scale
  # diag(seq(30, 20, length.out = P)) %*%
  # Mk has eigenvalues of different scales
  diag(seq(1000, 1, length.out = P)) %*%
  t(Conj(eigen(Mk)$vectors))

demo_random_post(Mk, sigma_02, n, Mk_df, ex_its, hpd_round = hpd_round)

# compare mean traces from rcwis ------------------------------------------

# sanity check the rcwis function
# show that it generates random matrices with traces that seem to match the expectation of its trace

trace_its <- 100
traces <- rep(NA, trace_its)

# Mk <- diag(P)

set.seed(14032025)
Mk <- rcwis(1000, diag(P))
Mk <- eigen(Mk)$vectors %*% 
  # Mk has eigenvalues of similar scale
  # diag(seq(30, 20, length.out = P)) %*%
  # Mk has eigenvalues of different scales
  diag(seq(1000, 1, length.out = P)) %*%
  t(Conj(eigen(Mk)$vectors))
Mk <- round(Mk, 11)

set.seed(14032025)
for (i in 1:trace_its) {
  traces[i] <- sum(Re(diag(rcwis(n, Mk))))
}

summary(traces)
n * sum(Re(diag(Mk)))

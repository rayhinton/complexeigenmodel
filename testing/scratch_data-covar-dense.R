# generate sample data based on the dense covariance data model

# based on inputs:

# - K, number of groups
# - nk, number rows per observation matrix
# - P, dimension of the observed data vectors
# - d, reduced dimension of interest
# - alphaBetaRange - list of min and max of alpha and beta; (0, 1) or (1, 2)?
# - bing_its, number of iterations for sampling from matrix Bingham distribution
# - simdataseed, a seed provided by the calling script
# note: U_k matrices will be Pxd

set.seed(simdataseed)

# generate true underlying covariance matrices ----------------------------
# Sigma_k = sigma_k^2 * (U_k * Lambda_k * U_k^H + I_P)
# need to generate sigma_k, U_k, Lambda_k

##### generate sigma_k
sigma_k2_0 <- 1/rgamma(K, 1, 1)

##### V matrix, PxP
V_0 <- runitary(P, P)

##### A and B matrices
# need to sample in terms of alpha, beta vectors
# and w scalar
minab <- min(alphaBetaRange)
maxab <- max(alphaBetaRange)
alpha_0 <- c(maxab, 
             runif(P-2, min = minab, max = maxab) |> sort(decreasing = TRUE), 
             minab)
beta_0 <- c(maxab, 
            runif(d-2, min = minab, max = maxab) |> sort(decreasing = TRUE), 
            minab)

# prior for w scalar 
# - I believe Hoff is recommending these parameters in terms of the shape and SCALE for the Gamma distribution. Otherwise, a large tau parameter does not seem to lead to a diffuse prior. 
# - Also, comparing some of his FCDs that are Inverse Gamma distributions, he seems to have the parameters directly multiplying the inverse of the r.v. value
eta_0 <- 2
tau2_0 <- 100

w_0 <- rgamma(1, eta_0/2, scale = tau2_0)

A_0 <- diag(sqrt(w_0) * alpha_0)
B_0 <- diag(sqrt(w_0) * beta_0)

##### generate U_k
# generating U_k, Lambda_k - these are just the eigenvectors, eigenvalues of some complex symmetric matrix
# these are the eigenvectors of part of the covariance matrix
U_k_0 <- array(NA, c(P, d, K))
# generate the U_k matrices from the prior, complex gen. Bingham (A, B, V)
A_0_gb <- V_0 %*% A_0 %*% t(Conj(V_0))
for (k in 1:K) {
    U_k_0[, , k] <- runitary(P, d)
    for (i in 1:bing_its) {
        U_k_0[, , k] <- rcmb(U_k_0[, , k], A_0_gb, B_0)
    }
}

###### generate Lambda_k matrices
Lambda_k_0 <- array(NA, c(d, d, K))
for (k in 1:K) {
    omega_k <- runif(d)
    lambda_k <- omega_k / (1 - omega_k)
    Lambda_k_0[, , k] <- diag(lambda_k)
}

# Sigma_k covariance matrices, in terms of the previous parameters
Sigma_k_0 <- array(NA, c(P, P, K))
for (k in 1:K) {
    Sigma_k_0[, , k] <- sigma_k2_0[k] * 
        (U_k_0[, , k] %*% Lambda_k_0[, , k] %*% t(Conj(U_k_0[, , k])) + 
             diag(P))
}

# generate data, Y_k and P_k ----------------------------------------------

# generate true mean vectors, mu_k
# mu_k_0 <- rcmvnorm(K, rep(0, P))
mu_k_0 <- matrix(0, K, P)

###### Y_k, List of observed data matrices
Y_k <- list()

for (k in 1:K) {
    Y_k[[k]] <- rcmvnorm(nk[k], mu_k_0[k, ], Sigma_k_0[, , k])
}

###### P_k, sums of squares matrices
P_k <- array(NA, c(P, P, K))

for (k in 1:K) {
    # P_k[, , k] <- t(Conj(Y_k[[k]])) %*% 
    #     (diag(nk[k]) - matrix(1/nk[k], nk[k], nk[k])) %*% Y_k[[k]]
    # P_k [, , k] <- t(Conj(Y_k[[k]])) %*% Y_k[[k]]
    P_k[, , k] <- rcwis(nk[k], Sigma_k_0[, , k])
}
# scratch and testing for FCDs

library(cmvnorm)

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")

# generate testing data ---------------------------------------------------

# number of groups
K <- 3

# number rows per observation matrix
nk <- c(20, 25, 30)
stopifnot(length(nk) == K)

# U_k matrices will be Pxd
# dimension of the observed vectors
P <- 6
# reduced dimension of interest
d <- 2

# generate true underlying covariance matrices ----------------------------
# Sigma_k = sigma_k^2 * (U_k * Lambda_k * U_k^H + I_P)
# need to generate sigma_k, U_k, Lambda_k

##### generate sigma_k
set.seed(15022025)
sigma_k2_0 <- 1/rgamma(K, 1, 1)

##### V matrix, PxP
V_0 <- runitary(P, P)

##### A and B matrices
# need to sample in terms of alpha, beta vectors
# and w scalar
alpha_0 <- c(1, runif(P-2) |> sort(decreasing = TRUE), 0)
beta_0 <- c(1, runif(d-2) |> sort(decreasing = TRUE), 0)

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
    for (i in 1:1000) {
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
mu_k_0 <- rcmvnorm(K, rep(0, P))

###### Y_k, List of observed data matrices
Y_k <- list()

for (k in 1:K) {
    Y_k[[k]] <- rcmvnorm(nk[k], mu_k_0[k, ], Sigma_k_0[, , k])
}

###### P_k, sums of squares matrices
P_k <- array(NA, c(P, P, K))

for (k in 1:K) {
    P_k[, , k] <- t(Conj(Y_k[[k]])) %*% 
        (diag(nk[k]) - matrix(1/nk[k], nk[k], nk[k])) %*% Y_k[[k]]
}


# sampling initialization -------------------------------------------------
# number of MCMC iteratons
S <- 10

# initialize arrays to store values
V_s <- array(NA, c(P, P, S))
V_s[, , 1] <- runitary(P, P)

##### A and B matrices
# need to sample in terms of alpha, beta vectors
# and w scalar
alpha_s <- matrix(NA, P, S)
beta_s <- matrix(NA, d, S)
w_s <- vector("numeric", S)

alpha_s[, 1] <- c(1, runif(P-2) |> sort(decreasing = TRUE), 0)
beta_s[, 1] <- c(1, runif(d-2) |> sort(decreasing = TRUE), 0)

# prior for w scalar 
# - I believe Hoff is recommending these parameters in terms of the shape and SCALE for the Gamma distribution. Otherwise, a large tau parameter does not seem to lead to a diffuse prior. 
# - Also, comparing some of his FCDs that are Inverse Gamma distributions, he seems to have the parameters directly multiplying the inverse of the r.v. value
eta_0 <- 2
tau2_0 <- 100

w_s[1] <- rgamma(1, eta_0/2, scale = tau2_0)

##### U_k matrices
# initialize an array to store the U_k over iterations
U_k_s <- array(NA, c(P, d, K, S))

# calculate temporary A and B matrices, and transformation of A for Bingham
A_s <- diag(sqrt(w_s[1]) * alpha_s[, 1])
A_gb <- V_s[, , 1] %*% A_s %*% t(Conj(V_s[, , 1]))
B_s <- diag(sqrt(w_s[1]) * beta_s[, 1])
for (k in 1:K) {
    U_k_s[, , k, 1] <- rcmb(runitary(P, d), A_gb, B_s)
}

##### Lambda_k matrices
# initialize an array to store Lambda_k diagonals
Lambda_k_s <- array(NA, c(d, K, S))

for (k in 1:K) {
    omega_k <- runif(d)
    Lambda_k_s[, k, 1] <- omega_k / (1 - omega_k)
}

##### sigma_k2 covariance matrix scales
sigma_k2_s <- matrix(NA, K, S)

sigma_k2_s[, 1] <- 1/rgamma(K, 1, 1)

# temporary sampling fill-in ----------------------------------------------

# TODO remove s, the explicit MCMC iteration value
s <- 2

# TODO remove these parameter fill-in values
# fill in the U_k and Lambda_k values, for now
U_k_s[, , 1:K, s] <- U_k_s[, , 1:K, 1]
Lambda_k_s[, 1:K, s] <- Lambda_k_s[, 1:K, 1]
w_s[s] <- w_s[1]
beta_s[, s] <- beta_s[, 1]
alpha_s[, s] <- alpha_s[, 1]
V_s[, , s] <- V_s[, , 1]

# sigma_k2 sampling -------------------------------------------------------

# sigma_k^2: an Inverse Gamma distribution

# depends on parameters U_k, Lambda_k, and observed P_k
# also depends on the complex Wishart degrees of freedom
# and the dimension the observation vectors, i.e. Px1

IP <- diag(P)

for (k in 1:K) {
    # set up temp matrices based on previously sampled values
    # TODO what order am I sampling, within an iteration? possibly change s to s-1
    Uk <- U_k_s[, , k, s]
    Lambdak <- diag(Lambda_k_s[, k, s])
    
    # calculate parameters for the Inverse Gamma distribution
    ak <- P * (nk[k] - 1)
    bk <- Re(sum(diag(
        solve(Uk %*% Lambdak %*% t(Conj(Uk)) + IP) %*% P_k[, , k]
    )))
    
    # FCD is Inverse Gamma: 1 / Gamma, in terms of the rate parameter
    sigma_k2_s[k, s] <- 1/rgamma(1, ak, rate = bk)
}

# V matrix sampling -------------------------------------------------------
# FCD is a complex matrix Bingham distribution
# depending on U_k, B, and A

# set up temporary vectors and matrices for a, A, b, B, for convenience
a_s <- sqrt(w_s[s]) * alpha_s[, s]
A_s <- diag(a_s)
b_s <- sqrt(w_s[s]) * beta_s[, s]
B_s <- diag(b_s)

# TODO could possibly do faster with sort of `apply`
sumMat <- matrix(0 + 0i, P, P)
for (k in 1:K) {
    # TODO since B_s is a diagonal matrix, could I do this multiplication more
    # quickly using recycling?
    sumMat <- sumMat + U_k_s[, , k, s] %*% B_s %*% t(Conj(U_k_s[, , k, s]))
}

# TODO need to make sure rcmb can sample **square matrices**
V_s[, , s] <- rcmb(V_s[, , s-1], sumMat, B_s)

# U_k, eigenvector matrices -----------------------------------------------

# there will be two loops:
# - over all the K U_k matrices
# - over the columns in a U_k matrix

# set k to 1, to just start with sampling one matrix
for (k in 1:K) {
    # make a temporary matrix to sample; starting with the previous sample s-1
    Uk <- U_k_s[, , k, s-1]
    # number of columns in U_k is d
    for (j in sample(1:d)) {
        # b_j is the jth diagonal entry of matrix B
        b_j <- b_s[j]
        omega_j <- Lambda_k_s[j, k, s] / (Lambda_k_s[j, k, s] + 1)
        tempmat <- V_s[, , s] %*% A_s %*% t(Conj(V_s[, , s]))

        Uk[, j] <- rcvb_LN(Uk, j, b_j * tempmat + omega_j * P_k[, , k])
    }
    # put the newly sampled matrix into the larger array at index s
    U_k_s[, , k, s] <- Uk
}
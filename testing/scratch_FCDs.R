# scratch and testing for FCDs

library(cmvnorm)

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")

# generate testing data ---------------------------------------------------

# we must generate "true" parameters, and then data that has a distribution
# based on these parameters

# number of groups
K <- 3

# number rows per observation matrix
nk <- c(20, 25, 30)
stopifnot(length(nk) == K)

# U_k matrices will be Pxd
# dimension of the observed vectors
P <- 8
# reduced dimension of interest
d <- 4

# min and max of alpha, beta vectors
alphaBetaRange <- c(1, 2)

# number of iterations for sampling from matrix Bingham distribution
bing_its <- 100

# generate parameters and data (Y_k, P_k) using this script
simdataseed <- 15022025
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_data-covar-dense.R")

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

alpha_s[, 1] <- c(max(alphaBetaRange), 
                  runif(P-2, min(alphaBetaRange), max(alphaBetaRange)) |> 
                      sort(decreasing = TRUE), 
                  min(alphaBetaRange))
beta_s[, 1] <- c(max(alphaBetaRange), 
                 runif(d-2, min(alphaBetaRange), max(alphaBetaRange)) |> 
                     sort(decreasing = TRUE), 
                 min(alphaBetaRange))

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
V_s[, , s] <- rcBingUP_gibbs(V_s[, , s-1], sumMat, A_s)

# U_k, eigenvector matrices -----------------------------------------------

# there will be two loops:
# - over all the K U_k matrices
# - over the columns in a U_k matrix

# this is a matrix quantity which is the same for all groups and columns
tempmat <- V_s[, , s] %*% A_s %*% t(Conj(V_s[, , s]))
# set k to 1, to just start with sampling one matrix
for (k in 1:K) {
    # make a temporary matrix to sample; starting with the previous sample s-1
    Uk <- U_k_s[, , k, s-1]
    # number of columns in U_k is d
    for (j in sample(1:d)) {
        # extract and calculate some temporary values
        # b_j is the jth diagonal entry of matrix B
        b_j <- b_s[j]
        # TODO possibly change to the (s-1) iteration of Lambda_k_s. If I sample
        # Lambda before this, then it should be s. If I sample Lambda after,
        # then the most recent iteration is s-1.
        omega_j <- Lambda_k_s[j, k, s] / (Lambda_k_s[j, k, s] + 1)

        # sample column j as a Bingham vector, orthogonal to other columns in Uk
        Uk[, j] <- rcvb_LN(Uk, j, b_j * tempmat + omega_j * P_k[, , k])
    }
    # put the newly sampled matrix into the larger array at index s
    U_k_s[, , k, s] <- Uk
}

# Lambda_k, eigenvector matrices ------------------------------------------

# If each omega_jk has a Uniform prior distribution, then the FCD of 1 - omega
# is a Gamma distribution, truncated to (0, 1).

for (k in 1:K) {
    for (j in 1:d) {
        # xi_jk Gamma(nk, tjk), shape, rate, 
        # and tjk is a temporary parameter defined as:
        # (note we take the Real part, since the quadratic form should be real, 
        # but may have a small complex part due to numerical issues.)
        # TODO change s, possibly: If I sample U_k before, then s. Otherwise, s-1. 
        tjk <- Re(
            t(Conj(U_k_s[, j, k, s])) %*% P_k[, , k] %*% U_k_s[, j, k, s]
            ) / sigma_k2_s[k, s]
        
        ##### sample from a truncated Gamma distribution
        # lower and upper bounds for truncating the Gamma distribution
        lb <- 0
        ub <- 1
        
        # upper and lower percentile bounds
        lp <- pgamma(lb, shape = nk[k], rate = tjk)
        up <- pgamma(ub, shape = nk[k], rate = tjk)
        
        # generate the random sample
        u <- runif(1, lp, up)
        # invert the Uniform r.v. to a truncated Gamma r.v.
        xi_jk <- qgamma(u, shape = nk[k], rate = tjk)
        # convert xi to Lambda
        Lambda_k_s[j, k, s] <- 1/xi_jk - 1
    }
}

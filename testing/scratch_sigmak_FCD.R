# validate sigma_k FCD

library(cmvnorm)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")

# generate data -----------------------------------------------------------

# - K, number of groups
# - nk, number rows per observation matrix
# - P, dimension of the observed data vectors
# - d, reduced dimension of interest
# - bing_its, number of iterations for sampling from matrix Bingham distribution
# - simdataseed, a seed provided by the calling script
# note: U_k matrices will be Pxd

P <- 8
d <- 4
K <- 3
nk <- 1000 + (1:K)*5
alphaBetaRange <- c(1, 2)
bing_its <- 100
simdataseed <- 8
# simdataseed <- 10032025
set.seed(simdataseed)

# source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_data-covar-dense.R")

sigma_k2_0 <- 1/rgamma(K, 1, 1)

# set up for the sampler --------------------------------------------------

M_k <- array(NA, c(P, P, K))
P_k <- array(NA, c(P, P, K))
for (k in 1:K) {
    M_k[, , k] <- rcwis(P+2, diag(P))
    P_k[, , k] <- rcwis(nk[k], sigma_k2_0[k]*M_k[, , k])
}

# number of sampling iterations
S <- 5000

##### some initial values of other parameters

##### array for sigma_k2 values
sigma_k2_s <- matrix(NA, K, S)

sigma_k2_s[, 1] <- 1/rgamma(K, 1, 1)

# sigmak sampler ----------------------------------------------------------

# sigma_k^2: an Inverse Gamma distribution

# depends on parameters U_k, Lambda_k, and observed P_k
# also depends on the complex Wishart degrees of freedom
# and the dimension the observation vectors, i.e. Px1

IP <- diag(P)

for (s in 1:S) {
    for (k in 1:K) {
        # set up temp matrices based on previously sampled values
        # TODO what order am I sampling, within an iteration? possibly change s to s-1
        
        # choose the most recent value of Uk
        # Uk <- U_k_s[, , k, s]
        # TODO testing case, just use the known value U_k_0
        # Uk <- U_k_0[, , k]
        
        # choose the most recent value of Lambda_k
        # Lambdak <- diag(Lambda_k_s[, k, s])
        # TODO testing case, just use the known value Lambda_k_0
        # Lambdak <- Lambda_k_0[, , k]
    
        # Mk <- Uk %*% Lambdak %*% t(Conj(Uk)) + IP
        Mk <- M_k[, , k]
        
        # calculate parameters for the Inverse Gamma distribution
        # ak <- P * (nk[k] - 1)
        ak <- P * (nk[k])
        bk <- Re(sum(diag(
            solve(Mk) %*% P_k[, , k]
        )))
        
        # FCD is Inverse Gamma: 1 / Gamma, in terms of the rate parameter
        sigma_k2_s[k, s] <- 1/rgamma(1, ak, rate = bk)
    }
}

rowMeans(sigma_k2_s)
apply(sigma_k2_s, 1, median)
sigma_k2_0

quantile(sigma_k2_s[1, ], c(.025, .975))
quantile(sigma_k2_s[2, ], c(.025, .975))
quantile(sigma_k2_s[3, ], c(.025, .975))

# note that the mean of such a Inv-Gamma is bk/(ak-1)
bk/sigma_k2_0[k] + 1
# ak is essentially fixed, for a certain data vector size and sample size
# so, what I should really be seeing is, how far off is bk, from the parameter value which would result in a mean which exactly matches the true sigmak?
sigma_k2_0[k] * (ak - 1)
bk/(ak - 1)

eigen(Mk)$values |> log10()

plot(sigma_k2_s[1, round(S/2):S])

curve(dgamma(1/x, ak, bk)/x^2, .8, 1.2)


# sampling distribution of this estimator ---------------------------------

# source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

# rcomplex_wishart(6, 4)

P1 <- 12
n1 <- 1000
set.seed(10032025)
M1 <- rcwis(n1, diag(P1))
M1 <- eigen(M1)$vectors %*% diag(seq(30, 0.5, length.out = P1)) %*%
    t(Conj(eigen(M1)$vectors)); diag(M1) <- Re(diag(M1))
# M1 <- rcomplex_wishart(P1+2, diag(P1))
# M1 <- diag(P1)
invM1 <- solve(M1)
# invM1 <- M1

# True sigma2 value
sigma2_1 <- 5
# True CW parameter
G1 <- sigma2_1 * M1
isSymmetric(G1)
isHermitian(G1)
G1 <- round(G1, digits = 9)

ests <- 1000
sigma2_ests <- rep(NA, ests)

for (i in 1:ests) {
    Pi <- rcwis(n1, G1)
    # Pi <- rcomplex_wishart(n1, G1)
    sigma2_ests[i] <- Re(sum(diag(invM1 %*% Pi))) / (P1*n1 - 1)
}

mean(sigma2_ests)
((P1*n1) / (P1*n1 - 1) * sigma2_1)
plot(density(sigma2_ests))

# if the diagonal elements are large, then the estimator seems to be good
M1diags <- diag(M1)
M1offdiags <- M1[!diag(TRUE, P1)]
sum(abs(M1offdiags)) / sum(abs(M1diags))

# does this function even give matrices with the right trace?
check_ests <- rep(NA, ests)
for (i in 1:ests) {
    # Xi <- rcwis(n1, M1)
    Xi <- rcomplex_wishart(n1, M1)
    check_ests[i] <- Re(sum(diag(Xi)))
}

n1 * P1
mean(check_ests)
sum(diag(M1)) |> Re()
plot(density(check_ests))


# different credible intervals --------------------------------------------

P1 <- 12
n1 <- 1000
set.seed(10032025)
M1 <- rcwis(n1, diag(P1))
M1 <- eigen(M1)$vectors %*% diag(seq(1, 0.5, length.out = P1)) %*% t(Conj(eigen(M1)$vectors)); diag(M1) <- Re(diag(M1))
summary(eigen(M1)$values)

1/eigen(M1)$values

sigma2_1 <- 5

G1 <- sigma2_1 * M1
G1 <- round(G1, 9)
invM1 <- solve(M1)
summary(eigen(invM1)$values)

set.seed(12032025)
Pk1 <- rcwis(n1, G1)
summary(eigen(Pk1)$values)

summary(abs(eigen(invM1 %*% Pk1)$values))

(a1 <- P1*n1)
Im(sum(diag(invM1 %*% Pk1)))
(b1 <- Re(sum(diag(invM1 %*% Pk1))))

curve(dgamma(1/x, shape = a1, rate = b1), from = 1, to = 20, n = 1001, add = FALSE)
abline(v = c(1/qgamma(.025, shape = a1, rate = b1),
        1/qgamma(.975, shape = a1, rate = b1))
)


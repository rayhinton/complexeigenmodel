# Complex Wishart using Cholesky decomposition

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

# sample the Cholesky decomposition components

# Goal: A ~ CW(n, Im), 
# - where A is mxm 
# - and d.f. = n

# functions ---------------------------------------------------------------

rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

# esp

# experiment with square roots and hermitian ------------------------------

Sigma <- rcomplex_wishart(100, diag(8))

microbenchmark::microbenchmark(
    rcomplex_wishart(1000, Sigma, useEigenR = FALSE, byCholesky = FALSE),
    rcomplex_wishart(1000, Sigma, useEigenR = FALSE, byCholesky = TRUE),
    rcomplex_wishart(1000, Sigma, useEigenR = TRUE, byCholesky = FALSE),
    rcomplex_wishart(1000, Sigma, useEigenR = TRUE, byCholesky = TRUE),
    times = 100
)

# ii, the diagonal components ---------------------------------------------
# assume tij are the components of a upper triangular matrix Tt

alt_rcomplex_wishart <- function(df, Sigma, useEigenR = FALSE) {
    
    # the square root can be faster to compute for small matrices via Eigen_sqrt
    if (useEigenR) {
        C <- EigenR::Eigen_sqrt(Sigma)
    } else {
        Sigmaevd <- eigen(Sigma)
        C <- Sigmaevd$vectors %*% diag(sqrt(Sigmaevd$values)) %*%
            t(Conj(Sigmaevd$vectors))
    }

    m <- nrow(Sigma)
    Tt <- matrix(0 + 0i, m, m)
    diag(Tt) <- sqrt(rgamma(m, df - m + 1:m))
    Tt[upper.tri(Tt)] <- rscnorm(choose(m, 2))
    WI <- Tt %*% t(Conj(Tt))
    
    # TODO get rid of t Conj, since C is Hermitian?
    # for small matrices, it can be faster
    A <- C %*% WI %*% C
    
    diag(A) <- Re(diag(A))
    
    return(A)
}

# Tt <- matrix(0 + 0i, m, m)

# tii ~ Gamma(n - m + i, 1)
# diag(Tt) <- sqrt(rgamma(m, n - m + 1:m))
# 
# Tt[upper.tri(Tt)] <- rscnorm(choose(m, 2))
# 
# A <- Tt %*% t(Conj(Tt))

m <- 8
n <- 100

A <- alt_rcomplex_wishart(n, diag(m))
A
# is A full rank? all non-zero eigenvalues
eigen(A)$values
# is A Hermitian?
isSymmetric(A)
# are the diagonals real-valued?
diag(A)

microbenchmark::microbenchmark(
    EigenR::Eigen_sqrt(A),
    {
        Aevd <- eigen(A)
        C <- Aevd$vectors %*% diag(sqrt(Aevd$values)) %*% 
            t(Conj(Aevd$vectors))
    },
    times = 1000
)

As <- array(NA, c(m, m, 1e4))

for (s in 1:1e4) {
    As[, , s] <- alt_rcomplex_wishart(n, diag(m))
}

avgA <- apply(As, c(1, 2), mean)
View(avgA)

# Frobenius distance between true mean and sample mean.
# I believe that the Frobenius distance is unbounded in the space of HPD matrices.
norm(diag(n, nrow = m) - avgA, "F")

microbenchmark::microbenchmark(
    alt_rcomplex_wishart(n, diag(m)),
    rcomplex_wishart(n, diag(m)),
    times = 1000
)

# CW matrices with arbitrary parameter ------------------------------------

W <- alt_rcomplex_wishart(n, A)

View(W)
eigen(W)$values
isSymmetric(W)
diag(W)

Ws <- array(NA, c(m, m, 1e4))

for (s in 1:1e4) {
    Ws[, , s] <- alt_rcomplex_wishart(n, A, useEigenR = TRUE)
}

avgW <- apply(Ws, c(1, 2), mean)

View(avgW)

A / avgW

cbind(n*A[, 2], avgW[, 2])

# distance between true mean and sample mean
norm(n*A - avgW, "F")
# note that the Frobenius distance upper bound depends on the eigenvalues, which can be large
cbind(n*eigen(A)$values, eigen(avgW)$values)


microbenchmark::microbenchmark(
    alt_rcomplex_wishart(n, A),
    rcomplex_wishart(n, A),
    times = 1000
)



# test Hermitian solution -------------------------------------------------

W <- rcomplex_wishart(5, diag(4), byCholesky = TRUE)
isSymmetric(W)

Wavg <- (W + t(Conj(W))) / 2
norm(W - Wavg, "F")
max(abs(W - Wavg))
cbind(W[, 1], Wavg[, 1])

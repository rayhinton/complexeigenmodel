# test rcgb, sampling from the Complex matrix generalized Bingham distribution

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")

# parameters
# random matrix U is NxP
# N, dimension of observed vectors (number of parameters) 
N <- 5
# P, dimension reduction, P < N
P <- 3
# n, number of observations, related to degrees of freedom for Wishart matrix
n <- 10

# test the real vector sampler first --------------------------------------

set.seed(12022025)
# Ar, a N x N positive semi-definite symmetric matrix
Ar <- rWishart(1, n, diag(N))[, , 1]

# sample from real valued vector Bingham Langevin distribution.
# According to Algorithm 4 in Appendix 1 of Abdallah

# expected output: 
# - random vector that is Nx1
# - that has norm 1
ur <- rvB(Ar)
dim(ur)
crossprod(ur) == 1

# length should be N
testthat::expect_length(ur, N)
# should be numeric, or double
testthat::expect_type(ur, "double")
# should not be complex
testthat::expect_false(is.complex(ur))
# output should have unit norm
testthat::expect_equal(crossprod(ur)[1,1], 1)

##########
# test the complex vector sampler -----------------------------------------

# complex parameters
# c, a N x 1 complex vector
# A, a N x N positive semi-definite Hermitian matrix


set.seed(3012025)
Ac <- rcwis(N, diag(N))

# A <- rWishart(1, N, diag(N))[, , 1]
# c <- t(rmvnorm(1, mean = rep(0, N)))

uc <- rcvb(Ac)

# length should be N
testthat::expect_length(uc, N)
# should be numeric, or double
testthat::expect_type(uc, "complex")
# should not be complex
testthat::expect_false(is.double(uc))
# output should have unit norm
# - real part should be 1
# - complex part should be 0
testthat::expect_equal(Re(t(Conj(uc)) %*% uc)[1, 1], 1)
testthat::expect_equal(Im(t(Conj(uc)) %*% uc)[1, 1], 0)

##########
# matrix complex generalized Bingham Langevin distribution ----------------

# parameters
# C, arbitrary NxP matrix
# C <- matrix(rnorm(N * P), N, P)
# set of matrices {A_p}, i.e. potentially one matrix for each column
# AP <- array(NA, dim = c(N, N, P))
# for (p in 1:P) {
    # AP[, , p] <- rcwis(N, diag(N))
# }

# parameters
# B, a PxP diagonal matrix
# A, a NxN positive semi-definite Hermitian matrix
set.seed(12022025)
Bc <- rgamma(P, 1, 0.5) |> sort(decreasing = TRUE)


# kappa, scalar

# generate a random initialization matrix, U_init
Uc <- rustiefel(N, P)
# p <- 1

p_test <- 1
up_LN <- rcvb_LN(Uc, p_test, Ac, Bc)

# test: the sampled u_p should be orthogonal to the other columns in U
# meaning: the dot product of the columns of U[, -p] with u_p should be a 0 vector.
# since the vectors are complex, need to decompose these into real and complex parts. 
# - real part should be 0
# - complex part should be 0
testthat::expect_equal(Re(t(Conj(Uc[, -p_test])) %*% up_LN), 
                       matrix(0, nrow = P-1, ncol = 1))
testthat::expect_equal(Im(t(Conj(Uc[, -p_test])) %*% up_LN), 
                       matrix(0, nrow = P-1, ncol = 1))

# test: new vector should have unit norm
testthat::expect_equal((t(Conj(up_LN)) %*% up_LN)[1, 1], 
                       (1 + 0i))

# test: new vector should have correct dimension
# - number of rows should equal number of rows of the original matrix
# - should be a column vector, i.e. a matrix with 1 column
testthat::expect_equal(dim(up_LN), 
                       c(nrow(Uc), 1))

# test full matrix sampling -----------------------------------------------

# store the old dimension, to compare to dimensions of newly sampled matrix
old_dim <- dim(Uc)

# draw a sample from complex matrix Bingham distribution
Uc <- rcmb(Uc, Ac, Bc)

# expect an error when using a square U matrix
testthat::expect_error(rcmb(diag(N), diag(N), diag(P)),
                       regexp = "nrow(U) must be greater than ncol(U)",
                       fixed = TRUE)

# the dimensions of the newly sampled matrix should be the same as the old matrix
testthat::expect_equal(old_dim, dim(Uc))

# columns of the new matrix should be orthonormal
# - real part of matrix should be I_p identity matrix
# - imaginary part of matrix should be zeros
testthat::expect_equal(Re(t(Conj(Uc)) %*% Uc), diag(P))
testthat::expect_equal(Im(t(Conj(Uc)) %*% Uc), matrix(0, P, P))

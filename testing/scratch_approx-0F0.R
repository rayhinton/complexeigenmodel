# approximating 0F0(A, B) by Nasuda

# problems with rcmvnorm: it is generating an error saying that the input Sigma_k matrices are not Hermitian positive definite
# the function isHermitian says that some of them are not Hermitian
diagComplex <- diag(1 + 1i, nrow = 4)
isHerm

# based on inputs:

# - K, number of groups
# - nk, number rows per observation matrix
# - P, dimension of the observed data vectors
# - d, reduced dimension of interest
# - bing_its, number of iterations for sampling from matrix Bingham distribution
# - simdataseed, a seed provided by the calling script
# note: U_k matrices will be Pxd

# calculating sigma(j) ----------------------------------------------------

# A, PxP, Hermitian
# B dxd, d < P, diagonal (real)
# B_1, PxP, block diagonal with B in upper left and 0s elsewhere

# assumption: values in b are distinct, thus have multiplicity 1 each
# thus we have multiplicities of b_1 that are:
# b_1 = (1, 1, ..., 1, P-d), length d+1

P <- 10
d <- 4
K <- 3
nk <- c(35, 40, 45)
alphaBetaRange <- c(1, 2)
bing_its <- 100
# simdataseed <- 6
simdataseed <- 6032025

# grid size to use for the discrete density samplers
gs <- 201

ms_1 <- c(rep(1, d), P-d) 
cumsum(ms_1)

# we need to calculate sigma(P), i.e. not just up to d or d+1, but for all values
sigmas <- c(1:d, rep(d+1, P-d))
length(sigmas)

# all i < j
firstpairs <- combn(1:P, 2) |> t()
# such that i < d + 1, equivalent to i <= d
Ppairs <- firstpairs[firstpairs[, 1] <= d, ]

# generate some random stand-in parameter values --------------------------

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_data-covar-dense.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")

# V matrix parameter
# V_0

# U_k matrix parameters
# U_k_0

# w scalar
# w_0

# "previous" alpha, beta samples
alpha_0
beta_0

# calculate matrix M, sum of function of V and U_k matrices
M <- matrix(0, nrow = P, ncol = d)
for (k in 1:K) {
    M <- M + Conj( t(Conj(V_0))%*%U_k_0[, , k] ) * t(Conj(V_0))%*%U_k_0[, , k]
}
# M should be exactly real, but it comes out complex with 0 imaginary values
M <- Re(M)

# TODO calculate vector (alpha' * M), the vector alpha times M


# sample alpha values -----------------------------------------------------

# this vector stays constant, when sampling alpha=
Mbeta <- M %*% beta_0
# start a new alpha vector, to collect the newly sampled values
alpha_s <- alpha_0

# TODO sample a_ind randomly, so we do it in a different order each time

# for the a index: select all the rows where the 1st or 2nd column is a_ind
a_ind <- 3
# these are the combined columns
# Ppairs[Ppairs[, 1] == a_ind | Ppairs[, 2] == a_ind, ]
# these have the indices I actually need
others <- Ppairs[Ppairs[, 1] == a_ind, 2]
some <- Ppairs[Ppairs[, 2] == a_ind, 1]

# for each value in grid xs = seq(0, 1, .1) (or some other sequence), calculate
# the product of (a[some] - x) and (x - a[others]). the sampled alpha value
# should be between the surrounding alpha values. since the alphas are sorted
# decreasing, that means it should be between the alpha values with a_ind+1 and
# a_ind-1, minimum and maximum respectively. However, since the prior is that
# the alpha values are order statistics from a Uniform(1, 2) distribution, and
# thus they are equal with probability 0, and further, since the approximation
# we are using requires that the alpha values are not equal, we should not give
# positive probability to the values possibly being equal in the posterior.
# Therefore, we should not use a sequence starting and ending exactly with the
# respective min and max alpha values. Instead, we will begin and end the
# sequences on a small offset from the exact surrounding alpha values.

# xs <- seq(0, 1, .1)
# create the min and max values for the sequence. Give them an offset, so that we do not possibly sample equal alpha values.
min_al <- alpha_s[a_ind + 1]
max_al <- alpha_s[a_ind - 1]
(max_al - min_al)/(gs-1)/2
gint <- (max_al - min_al)/(gs-1)
min_al <- min_al + gint/2
max_al <- max_al - gint/2

xs <- seq(min_al, max_al, length.out = gs)

# TODO consider making these logarithms
# TODO then doing Gumbel trick sampling? or whatever
aprod1 <- sapply(xs, function(x) {
    sum(log((alpha_s[some] - x))) + sum(log(x - alpha_s[others]))
    })

# aprod2
aprod2 <- w_0 * xs * Mbeta[a_ind]
# calculate the proportional density, with a possible additional factor
proplogdens <- K * aprod1 + aprod2
if (2 <= a_ind & a_ind <= d) {
    # optional: aprod3 (only needed if 2 <= a_ind <= d)
    aprod3 <- -K * w_0 * alpha_s[a_ind] * beta_0[a_ind]
    proplogdens <- proplogdens + aprod3
}

summary(proplogdens)

# proplogdens <- proplogdens - max(proplogdens)
# summary(proplogdens)

# draw one value from xs (i.e. an alpha value), using the logweights
sample_gumbel(xs, 1, proplogdens)
c(min_al, max_al)

# for the b index: since we are assuming b values are distinct, and we only need
# b values from 2 to d-1, we in fact have that sigma(k) = k for 2 <= k <= d-1.

# select all the rows where the 1st or 2nd column is b_ind
b_ind <- 2
Ppairs[Ppairs[, 1] == b_ind | Ppairs[, 2] == b_ind, ]


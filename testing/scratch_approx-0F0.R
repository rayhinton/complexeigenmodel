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

P <- 8
d <- 4
K <- 3
nk <- 30 + (1:K)*5
alphaBetaRange <- c(1, 2)
bing_its <- 100
# simdataseed <- 6
simdataseed <- 8032025

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
# TODO the choice of pairs should depend on the smallest value in A, B matrices
# - if alpha, beta range from 2 to 1: use i <= d
# - if alpha, beta range from 1 to 0: use i <= d - 1
Ppairs <- firstpairs[firstpairs[, 1] <= d, ]

nrow(firstpairs)
nrow(Ppairs)
choose(P, 2) - choose(P-d, 2)

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
# alpha_s <- alpha_0

S_its <- 1000

alpha_S <- matrix(NA, P, S_its)
alpha_S[c(1, P), ] <- rev(alphaBetaRange)
alpha_S[, 1:5]

# initialize the first sample
set.seed(8032025)
# alpha_S[2:(P-1), 1] <- sort(runif(P-2, 1, 2), decreasing = TRUE)
alpha_S[, 1] <- seq(alphaBetaRange[2], alphaBetaRange[1], length.out = P)
alpha_S[, 1:5]
for (s in 2:S_its) {
    
    alpha_s <- alpha_S[, s-1]
    
    # TODO sample a_ind randomly, so we do it in a different order each time
    for (a_ind in sample(2:(P-1))) {
        
        # for the a index: select all the rows where the 1st or 2nd column is a_ind
        # a_ind <- 3
        # these are the combined columns
        # Ppairs[Ppairs[, 1] == a_ind | Ppairs[, 2] == a_ind, ]
        # these have the indices I actually need
        others <- Ppairs[Ppairs[, 1] == a_ind, 2]
        some <- Ppairs[Ppairs[, 2] == a_ind, 1]
        
        # for each value in grid xs = seq(0, 1, .1) (or some other sequence),
        # calculate the product of (a[some] - x) and (x - a[others]). the
        # sampled alpha value should be between the surrounding alpha values.
        # since the alphas are sorted decreasing, that means it should be
        # between the alpha values with a_ind+1 and a_ind-1, minimum and maximum
        # respectively. However, since the prior is that the alpha values are
        # order statistics from a Uniform(1, 2) distribution, and thus they are
        # equal with probability 0, and further, since the approximation we are
        # using requires that the alpha values are not equal, we should not give
        # positive probability to the values possibly being equal in the
        # posterior. Therefore, we should not use a sequence starting and ending
        # exactly with the respective min and max alpha values. Instead, we will
        # begin and end the sequences on a small offset from the exact
        # surrounding alpha values.
        
        # xs <- seq(0, 1, .1)
        # create the min and max values for the sequence. Give them an offset, so that we do not possibly sample equal alpha values.
        min_al <- alpha_s[a_ind + 1]
        max_al <- alpha_s[a_ind - 1]
        gint <- (max_al - min_al)/(gs-1)
        min_al <- min_al + gint/2
        max_al <- max_al - gint/2
        
        xs <- seq(min_al, max_al, length.out = gs)
        
        # TODO check that this sum of logs is actually calculating the right log density
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
        
        # draw one value from xs (i.e. an alpha value), using the logweights
        alpha_s[a_ind] <- sample_gumbel(xs, 1, proplogdens)
    }
    
    alpha_S[, s] <- alpha_s
}

rowMeans(alpha_S[, floor(S_its/2):S_its])
alpha_0

# for the b index: since we are assuming b values are distinct, and we only need
# b values from 2 to d-1, we in fact have that sigma(k) = k for 2 <= k <= d-1.

# compare - previous sample means
# [1] 2.000000 1.990269 1.965859 1.928389
# [5] 1.814975 1.770744 1.554632 1.306260
# [9] 1.053582 1.000000

# recent sample means
# [1] 2.000000 1.990328 1.965964 1.928305
# [5] 1.813397 1.768915 1.562277 1.335526
# [9] 1.129983 1.000000


# B sampling --------------------------------------------------------------

alphaM <- t(alpha_s) %*% M

S_its <- 1000

beta_S <- matrix(NA, d, S_its)
beta_S[c(1, d), ] <- c(2, 1) 
beta_S[, 1:5]

# initialize the first sample
set.seed(8032025)
# alpha_S[2:(P-1), 1] <- sort(runif(P-2, 1, 2), decreasing = TRUE)
beta_S[, 1] <- seq(2, 1, length.out = d)
beta_S[, 1:5]

b_ind <- 2
for (s in 2:S_its) {
    
    beta_s <- beta_S[, s - 1]
    
    for (b_ind in sample(2:(d-1))) {
        
        # select all the rows where the 1st or 2nd column is b_ind
        # indexes that are after b_ind (i.e. 2nd)
        b2s <- Ppairs[Ppairs[, 1] == b_ind & Ppairs[, 2] <= d, 2]
        # indexes that are before b_ind (i.e. 1st)
        b1s <- Ppairs[Ppairs[, 2] == b_ind, 1]
        
        # create a grid of beta values, depending on the next highest and lowest
        # beta_s values. first, determine the upper and lower limits of the
        # grid.
        min_bl <- beta_s[b_ind + 1]
        max_bl <- beta_s[b_ind - 1]
        gint <- (max_bl - min_bl)/(gs-1)
        min_bl <- min_bl + gint/2
        max_bl <- max_bl - gint/2
        # make the grid
        xs <- seq(min_bl, max_bl, length.out = gs)
        
        # calculate terms in the log density
        bprod1 <- sapply(xs, function(x) {
            sum(log((beta_s[b1s] - x))) + sum(log(x - beta_s[b2s]))
            })
        bprod2 <- (P-d)*log(beta_s[b_ind])
        bprod3 <- xs * (-K * w_0 * alpha_s[b_ind] + alphaM[b_ind])
        
        # calculate the log density
        proplogdens_b <- K*bprod1 + K*bprod1 + bprod3
        
        # draw one value from xs (i.e. beta value), using the logweights
        beta_s[b_ind] <- sample_gumbel(xs, 1, proplogdens_b)
    }
    
    beta_S[, s] <- beta_s
}

rowMeans(beta_S[, floor(S_its/2):S_its])
beta_0

# sample w ----------------------------------------------------------------

fstar <- choose(P, 2) - choose(P-d, 2)
betatilde <- c(beta_s, rep(0, P-d))

w_shape <- K*fstar + eta_0/2
w_rate <- K * (alpha_s %*% betatilde) - t(alpha_s) %*% M %*% beta_s + tau2_0

w_s <- rgamma(1, shape = w_shape, rate = w_rate)
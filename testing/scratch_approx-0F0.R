# approximating 0F0(A, B) by Nasuda

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
nk <- 1000 + (1:K)*5
alphaBetaRange <- c(1, 2)
bing_its <- 5000
# simdataseed <- 6
# simdataseed <- 8032025


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
# source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_data-covar-dense.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")

# set.seed(8052025) 
set.seed(3145)

# V matrix parameter
# V_0
V_0 <- runitary(P, P)

# w scalar
# w_0
w_0 <- rgamma(1, 1, 1)

# "previous" alpha, beta samples
alpha_0 <- c(alphaBetaRange[2],
             runif(P-2, alphaBetaRange[1], alphaBetaRange[2]) |> sort(decreasing = TRUE),
             alphaBetaRange[1])
beta_0 <- c(alphaBetaRange[2],
            runif(d-2, alphaBetaRange[1], alphaBetaRange[2]) |> sort(decreasing = TRUE),
            alphaBetaRange[1])
alpha_0
beta_0

A_0 <- diag(sqrt(w_0) * alpha_0)
B_0 <- diag(sqrt(w_0) * beta_0)

# U_k matrix parameters
# U_k_0
U_k_0 <- array(NA, c(P, d, K))
Ukinit <- runitary(P, d)
for (i in 1:bing_its) {
    if (i %% 500 == 0) {
        print(paste0("i = ", i))
    }
    Ukinit <- rcmb(Ukinit, A_0, B_0)
}
for (s in 1:(K*100)) {
    Ukinit <- rcmb(Ukinit, A_0, B_0)
    if (s %% 100 == 0) {
        print(paste0("k = ", s/100))
        U_k_0[, , s/100] <- Ukinit
    }
}

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

S_its <- 10000

alpha_S <- matrix(NA, P, S_its)
alpha_S[c(1, P), ] <- rev(alphaBetaRange)
alpha_S[, 1:5]

# initialize the first sample
# set.seed(8032025)
set.seed(3145)
# alpha_S[2:(P-1), 1] <- sort(runif(P-2, 1, 2), decreasing = TRUE)
alpha_S[, 1] <- seq(alphaBetaRange[2], alphaBetaRange[1], length.out = P)
# TODO unrealistic starting values, but testing
# alpha_S[, 1] <- alpha_0
alpha_S[, 1:5]
for (s in 2:S_its) {
    if (s %% 500 == 0) {print(paste0("s = ", s))}
    
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

quantile(alpha_S[7, floor(S_its/2):S_its], probs = c(.025, .975))

apply(alpha_S[2:(P-1), floor(S_its/2):S_its], 1, 
      quantile, probs = c(0.025, 0.975)) |> 
    t() |>
    cbind("true" = alpha_0[2:(P-1)])

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
        proplogdens_b <- K*bprod1 + K*bprod2 + bprod3
        
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

# test approximations of 0F0 with known values ----------------------------

Betaf <- 2
ms <- c(6, 2)
m <- sum(ms)
C <- m*Betaf/2

Bmvgamma <- function(C, m, Betaf) {
    term1 <- pi^(.25*m*(m-1)*Betaf)
    term2 <- prod(gamma(C - ((1:m) - 1)/2*Betaf))
    return(term1*term2)
}

# Bmvgamma(C, m, Betaf)

omegamBeta <- function(m, Betaf) {
    term1 <- m*log(2)
    term2 <- (m^2 * Betaf / 2) * log(pi)
    term3 <- log(Bmvgamma(m*Betaf/2, m, Betaf))
    
    return(exp(term1 + term2 - term3))
}

# omegamBeta(m, Betaf)

Omegams <- function(ms, Betaf) {
    term2 <- 0
    for(j in 1:length(ms)) {
        term2 <- term2 + log(omegamBeta(ms[j], Betaf))
    }
    
    logresult <- term2 - log(omegamBeta(sum(ms), Betaf))
    
    return(exp(logresult))
}

###
# pi
###
# simplifying case: when B is I_4, then the multiplicities are 4, 4 and thus
# s is 4*4 = 16
ss <- prod(ms)

pi^(Betaf/2 * ss)

###
# Omega
###

Omegams(ms, 2)

### 
# J(A, B)
###

as <- seq(2, 1, length.out = m)

ijs <- cbind(rep(1:(ms[1]), each = ms[2]), 
             rep((ms[1]+1):m, times = ms[2]))
lJAB <- sum(log(as[ijs[, 1]] - as[ijs[, 2]]))

###
# etr
###
exp(sum(as[1:ms[1]]))

###
# approximation
###

(log(pi)*(Betaf/2 * ss) +
    log(Omegams(ms, 2)) +
    lJAB*(-Betaf/2) +
    sum(as[1:ms[1]])) |> exp()

###
# True
###
exp(sum(as))


# Monte Carlo approximation of 0F0 ----------------------------------------

# simulate uniform matrix from U(P)

rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}
rnorm(1, 0, 1/sqrt(2)) + 1i * rnorm(1, 0, 1/sqrt(2))

rcstiefel <- function(P, d) {
    X1 <- matrix(rscnorm(P*d), ncol = d)
    return(qr(X1) |> qr.Q())
}

P <- 4

M <- 1e4
intgd <- rep(NA, M)

set.seed(10052025)
as <- c(2, 
       runif(P-2, 1, 2) |> sort(decreasing = TRUE),
       1)
bs <- c(2, 
       runif(P-2, 1, 2) |> sort(decreasing = TRUE),
       1)

A <- diag(as)
B <- diag(bs)

for (i in 1:M) {
    U <- rcstiefel(P, P)
    # intgd[i] <- (A %*% U %*% B %*% t(Conj(U))) |> 
    #     diag() |> sum() |> Re() |> exp()
    intgd[i] <- exp(Re(t(as) %*% (Conj(U) * U) %*% bs))
}

mean(intgd)
summary(intgd)

# approximation

ss <- P*(P-1)/2

ijs <- t(combn(1:P, 2))

JAB <- prod((as[ijs[, 1]] - as[ijs[, 2]]) * (bs[ijs[, 1]] - bs[ijs[, 2]]))

Oms <- Omegams(rep(1, P), 2)

exp(sum(as * bs))

pi^ss *
    Oms *
    (JAB^-1) *
    exp(sum(as * bs))


# approximation by truncated sum ------------------------------------------

Is <- rep(1, P)

# this calculates kappa, the partitions for a particular k
kappa <- partitions::parts(5)
kappa1 <- kappa[, 1]

truncK <- 5
hgfsum <- 0

for (k in 1:truncK) {
    logk <- lfactorial(k)
    partsk <- partitions::parts(k)
    
    for(j in 1:ncol(partsk)) {
        kappa <- partsk[, j]
        logCA <- log(jack::Schur(as, kappa))
        logCB <- log(jack::Schur(bs, kappa))
        logCI <- log(jack::Schur(Is, kappa))
        sumadd <- exp(logCA + logCB - logk - logCI)
        
        hgfsum <- hgfsum + ifelse(is.nan(sumadd), 0, sumadd)
        # hgfsum <- hgfsum + exp(logCA)*exp(logCB)/ exp(logk) / exp(logCI)
    }
}

hgfsum


# other packages ----------------------------------------------------------


# Bmvgamma <- function(C, m, Betaf) {
Bmvgamma(10, 8, 1)
HypergeoMat::mvgamma(10+0i, 8)

ijs <- t(combn(1:P, 2))

set.seed(10052025)
as <- c(2, 
        runif(P-2, 1, 2) |> sort(decreasing = TRUE),
        1)
bs <- c(2, 
        runif(P-2, 1, 2) |> sort(decreasing = TRUE),
        1)

cij <- (as[ijs[, 1]] - as[ijs[, 2]]) * (bs[ijs[, 1]] - bs[ijs[, 2]])

2^P *
    exp(sum(as*bs)) *
    prod(sqrt(pi/cij))


# truncated sum

Is <- rep(1, P)

# this calculates kappa, the partitions for a particular k
kappa <- partitions::parts(5)
kappa1 <- kappa[, 1]

truncK <- 5
hgfsum <- 0

for (k in 1:truncK) {
    logk <- lfactorial(k)
    partsk <- partitions::parts(k)
    
    for(j in 1:ncol(partsk)) {
        kappa <- partsk[, j]
        logCA <- log(jack::Zonal(as, kappa))
        logCB <- log(jack::Zonal(bs, kappa))
        logCI <- log(jack::Zonal(Is, kappa))
        sumadd <- exp(logCA + logCB - logk - logCI)
        
        hgfsum <- hgfsum + ifelse(is.nan(sumadd), 0, sumadd)
        # hgfsum <- hgfsum + exp(logCA)*exp(logCB)/ exp(logk) / exp(logCI)
    }
}

hgfsum


# HGF of a single matrix argument -----------------------------------------

# 0F0(A) = etr(A)

P <- 4

library(foreach)
library(doParallel)

parallel::detectCores()

cluster <- makeCluster(12)
registerDoParallel(cluster)

# set.seed(10052025)
set.seed(12052025)
as <- c(2, 
        runif(P-2, 1, 2) |> sort(decreasing = TRUE),
        1)

truncK <- 15
hgfsum <- 0

for (k in 1:truncK) {
    logkfac <- lfactorial(k)
    partsk <- partitions::parts(k)
    
    print(paste0("k = ", k, "; num. parts = ", ncol(partsk)))
    # partsums <- rep(0, ncol(partsk))
    partsums <- vector("list", ncol(partsk))
    # for(j in 1:ncol(partsk)) {
    partsums <- foreach(j = 1:ncol(partsk)) %dopar% {
        kappa <- partsk[, j]
        # Zonal, real
        # logCA <- log(jack::Zonal(as, kappa))
        # Schur, complex
        logCA <- log(jack::Schur(as, kappa))
        sumadd <- exp(logCA - logkfac)
        
        partsums[j] <- ifelse(is.nan(sumadd), 0, sumadd)
        # hgfsum <- hgfsum + exp(logCA)*exp(logCB)/ exp(logk) / exp(logCI)
    }
    hgfsum <- hgfsum + Reduce("+", partsums)
}

hgfsum
exp(sum(as))
hgfsum/exp(sum(as))

# alpha 1 is complex, Schur
HypergeoMat::hypergeomPFQ(20, NULL, NULL, as, alpha = 1)
# alpha 2 is real, Zonal
HypergeoMat::hypergeomPFQ(20, NULL, NULL, as, alpha = 2)

stopCluster(cl = cluster)

partsk <- partitions::parts(3)
kappa <- partsk[, 1]
jack::Schur(as, kappa)
jack::Zonal(as, kappa)
jack::Jack(as, kappa, 1/4)
jack::Jack(c(1, 3/2, -2/3), lambda = c(3, 1), alpha = 1/4, which = "Z")

jack::Zonal(as, 0)


# testing the normalization of the Jack polynomials -----------------------

k <- 4
partsk <- partitions::parts(k)
as <- 1:4

sumkappa <- 0
for (j in 1:ncol(partsk)) {
    kappa <- partsk[, j]
    # Zonal, real
    # sumkappa <- sumkappa + jack::Zonal(as, kappa)
    # Schur, complex
    # sumkappa <- sumkappa + jack::Schur(as, kappa)
    sumkappa <- sumkappa + (jack::JackPol(P, kappa, 1, which = "C") |> qspray::evalQspray(tas))
    
}

sumkappa
sum(as)^k

# testing the proportions of Jack polynomials
# a test vector, that is compatible with qspray
tas <- 1:4

kappa <- partsk[, 5]

jack::Zonal(tas, kappa)
jack::Jack(tas, kappa, alpha = 2)

# ratio
jack::Zonal(tas, kappa) / jack::Jack(tas, kappa, alpha = 2)

# which can be J, P, Q, C: seems like J is default
jack::JackPol(4, kappa, 2, which = "C") |> print()
# alpha 2 and which = C seems to match with Zonal; which matches with documentation
jack::JackPol(4, kappa, 2, which = "C") |> qspray::evalQspray(tas)

# real, 2 matrix arguments ------------------------------------------------

cluster <- makeCluster(12)
registerDoParallel(cluster)

Is <- rep(1, P)

# set.seed(10052025)
set.seed(10052025)
as <- c(2, 
        runif(P-2, 1, 2) |> sort(decreasing = TRUE),
        1)
bs <- c(2, 
        runif(P-2, 1, 2) |> sort(decreasing = TRUE),
        1)

truncK <- 18
hgfsum <- 0

for (k in 1:truncK) {
    logkfac <- lfactorial(k)
    partsk <- partitions::parts(k)
    
    print(paste0("k = ", k, "; num. parts = ", ncol(partsk)))
    # partsums <- rep(0, ncol(partsk))
    partsums <- vector("list", ncol(partsk))
    # for(j in 1:ncol(partsk)) {
    partsums <- foreach(j = 1:ncol(partsk)) %dopar% {
        kappa <- partsk[, j]
        
        # Zonal, real
        logCA <- log(jack::Zonal(as, kappa))
        logCB <- log(jack::Zonal(bs, kappa))
        logCI <- log(jack::Zonal(Is, kappa))
        
        # Schur, complex
        # logCA <- log(jack::Schur(as, kappa))
        
        sumadd <- exp(logCA + logCB - logkfac - logCI)
        
        partsums[j] <- ifelse(is.nan(sumadd), 0, sumadd)
    }
    hgfsum <- hgfsum + Reduce("+", partsums)
}

hgfsum

stopCluster(cl = cluster)


# real Monte Carlo integral approximation ---------------------------------

M <- 1e5
intgd <- rep(NA, M)

for (i in 1:M) {
    # U <- rcstiefel(P, P)
    U <- rstiefel::rustiefel(P, P)
    # intgd[i] <- (A %*% U %*% B %*% t(Conj(U))) |> 
    #     diag() |> sum() |> Re() |> exp()
    intgd[i] <- exp(Re(t(as) %*% (Conj(U) * U) %*% bs))
}

mean(intgd)
summary(intgd)

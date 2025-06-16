# scratch CMACG prior for Uk

library(EigenR)
library(MASS)

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

# functions ---------------------------------------------------------------

mvGamma <- function(c, m, betaf, logscale = FALSE) {
    term1 <- .25*m*(m-1)*betaf*log(pi)
    term2 <- sum(lgamma(c - ((1:m)-1)/2 * betaf))
    
    if (logscale) {
        return(term1 + term2)
    } else {
        return(exp(term1 + term2))
    }
}

rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

runif_stiefel <- function(P, d, betaf) {
    # real
    if (betaf == 1) {
        X1 <- matrix(rnorm(P*d), ncol = d)
        # complex
    } else if (betaf == 2) {
        X1 <- matrix(rscnorm(P*d), ncol = d)        
    }
    
    # take QR decomposition - not guaranteed to be unique due to numerical methods
    X1qr <- qr(X1)
    QX1 <- qr.Q(X1qr)
    # extract sign of the diagonal elements of R
    D <- sign(Re(diag(qr.R(X1qr))))
    # transform by Q = QS, where S = diag(D). this transformation guarantees Q is unique.
    
    # Qfinal <- t(t(QX1) * D) # faster than matrix multiplication, for larger matrices
    Qfinal <- QX1 %*% diag(D)
    
    return(Qfinal)
}

rCMACG <- function(nrow, ncol, Sigma) {
    # simulate X ~ CMN(0, I_P, I_d),
    # calculate Z by transforming X with Cholesky decomposition of Sigma0
    # let H_Z be orthogonal portion of polar decomposition of Z
    X <- matrix(rscnorm(nrow*ncol), nrow, ncol)
    
    # check if this is the upper or lower triang version
    # want the matrix such that Sigma = cholSigma %*% t(Conj(cholSigma))
    # Eigen_chol returns the upper triangular portion
    cholSigma <- EigenR::Eigen_chol(Sigma)
    
    Z <- t(Conj(cholSigma)) %*% X 
    sqrtZHZ <- EigenR::Eigen_sqrt( t(Conj(Z)) %*% Z )
    return(Z %*% solve(sqrtZHZ))
}

dCMACG <- function(X, Sigma, logscale = TRUE) {
    m <- nrow(X)
    r <- ncol(X)
    
    invSigma <- solve(Sigma)
    
    term1 <- -r * log( EigenR::Eigen_det(Sigma) )
    term2 <- -m * log( EigenR::Eigen_det( t(Conj(X)) %*% invSigma %*% X ) )
    
    log_density <- Re(term1 + term2)
    
    if (logscale) {
        return(log_density)
    } else {
        return(exp(log_density))
    }
}

firstRowPos <- function(X) {
    rowsigns <- sign( Re(X[1, ]) )
    colflip <- ifelse(rowsigns == 0, 1, rowsigns)
    return( t(t(X) * colflip) )
}

std_cmpx_evec <- function(X) {
    Xc <- matrix(NA, nrow(X), ncol(X))
    # input: matrix X, Pxd
    for (j in 1:ncol(X)) {
        k <- which.max(Mod(X[, j]))
        phik <- Arg(X[k, j])
        ck <- exp(-1i * phik)
        Xc[, j] <- X[, j] * ck
    }
    
    return(Xc)
}

frob_dist <- function(A, B, returnDists = FALSE) {
    diffAB <- A - B
    sqDists <- diag( t(Conj(diffAB)) %*% diffAB )
    Fdist <- Re(sum(sqDists))
    if (returnDists) {
        return(list(Fdist = sqrt(Fdist),
                    sqDists = sqDists))
    } else {
        return(sqrt(Fdist))
    }
}

# parameter generation ----------------------------------------------------

P <- 8
d <- 4
nk <- 1000
S_its <- 1e5
tau_U <- .001 # SD of the entries of the Cayley transformation

doCayleyZeros <- FALSE
CayleyZeroProb <- 0.15

tracenorm <- TRUE
customEvals <- FALSE

set.seed(27052025)
# Sigma0 ~ diffuse Complex Wishart, PxP
Sigma0 <- rcomplex_wishart(P+1, P, diag(P))
if (tracenorm) {
    Sigma0 <- P * Sigma0 / Re(sum(diag(Sigma0)))
}
if (customEvals) {
    SigmaV <- eigen(Sigma0)$vectors
    sigmaevals <- c(rep(1.5, d),
                    rep(.5, d))
    Sigma0 <- SigmaV %*% diag(sigmaevals) %*% t(Conj(SigmaV))
}

invSigma0 <- solve(Sigma0)

eigen(Sigma0)$values
eigen(Sigma0)$values |> cumsum() / P

###
# Uk0 ~ CMACG(Sigma0)
###
Uk0 <- rCMACG(P, d, Sigma0)

t(Conj(Uk0)) %*% Uk0

# Lambdak0: diagonal matrix of decreasing positive entries
Lambdak0 <- diag(d:1)
Omegak0 <- solve( solve(Lambdak0) + diag(d) )

# sigmak02: random positive scalar
sigmak02 <- rgamma(1, 1, 1)

# Yk ~ CW(nk, sigmak02 * (Uk0 %*% Lambdak0 %*% t(Conj(Uk0))) )
Gamma0 <- sigmak02 * (Uk0 %*% Lambdak0 %*% t(Conj(Uk0)) + diag(P))
Yk <- rcomplex_wishart(nk, P, Gamma0)

# understand U0 and UE, true eigenvectors ---------------------------------

set.seed(6062025)
U0 <- runif_stiefel(P, d, 2)

Cov0 <- U0 %*% Lambdak0 %*% t(Conj(U0))
UE <- eigen(Cov0)$vectors[, 1:d]

# not the same
cbind(U0[, 1], UE[, 1])

# appear equal
(U0 %*% t(Conj(U0)))[, 1] 
(UE %*% t(Conj(UE)))[, 1]

# equal, within numerical tolerance
all.equal(
    (U0 %*% t(Conj(U0))),
    (UE %*% t(Conj(UE)))
    )

Rmat <- t(Conj(UE)) %*% U0
Rmatz <- zapsmall(Re(Rmat)) + 1i * zapsmall(Im(Rmat))

t(Conj(Rmatz)) %*% Rmatz
Rmatz %*% t(Conj(Rmatz))

t(Conj(Rmat)) %*% Rmat
Rmat %*% t(Conj(Rmat))

all.equal(
    (U0 %*% Lambdak0 %*% t(Conj(U0))),
    (UE %*% Lambdak0 %*% t(Conj(UE)))
)

U0[, 1]

# find a complex scalar so that the first real entry of the row is 0.5
a <- Re(U0[1, 1])
b <- Im(U0[1, 1])

C <- (a - sqrt(4*a^2 * b^2 + 4*b^4 - b^2 + 0i)) / (2*a^2 + 2*b^2)
D <- sqrt(1 - C^2)


# understand scaling of columns -------------------------------------------

X <- runif_stiefel(P, d, 2)

Xc <- std_cmpx_evec(X)

# a test that the resulting matrices are "equivalent"
# 1. They create the same projection matrix
(X %*% t(Conj(X)))[, 1]
(Xc %*% t(Conj(Xc)))[, 1]

# 2. The right-rotation matrix that transforms one to the other is a diagonal matrix
Rc <- t(Conj(Xc)) %*% X

zapsmall(Re(Rc))
zapsmall(Im(Rc))

# 3. The output for a standardized matrix is the same thing
Xc2 <- std_cmpx_evec(Xc)

zapsmall(Xc) == zapsmall(Xc2)
all.equal(Xc, Xc2)

# 4. Is the output the same for a -1 flipped column?
Xc3 <- std_cmpx_evec(cbind(-X[, 1], X[, 2:d]))

cbind(Xc3[, 1], Xc2[, 1])
all.equal(Xc3, Xc)

# 5. The columns are still orthogonal and unit norm
(t(Conj(Xc)) %*% Xc)
(t(Conj(Xc2)) %*% Xc2)
(t(Conj(Xc3)) %*% Xc3)

# MH proposal: Uniform on Stiefel(d, P) -----------------------------------

Uk_S <- array(NA, c(P, d, S_its))
Uk_S[, , 1] <- runif_stiefel(P, d, betaf = 2)
Us <- Uk_S[, , 1]

accCount <- 0
for (s in 2:S_its) {
    Up <- runif_stiefel(P, d, betaf = 2)
    
    tracep <- Re(sum(diag( Up %*% Omegak0 %*% t(Conj(Up)) %*% Yk )))
    # there is a separate determinant function that has a logarithm option
    logdetp <- log(Re(EigenR::Eigen_det( t(Conj(Up)) %*% invSigma0 %*% Up )))
    
    traces <- Re(sum(diag( Us %*% Omegak0 %*% t(Conj(Us)) %*% Yk )))
    # there is a separate determinant function that has a logarithm option
    logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*% invSigma0 %*% Us )))
    
    (logr <- tracep/sigmak02 - P*logdetp - traces/sigmak02 + P*logdets)
    
    if (log(runif(1)) <= logr) {
        Uk_S[, , s] <- Up
        Us <- Up
        accCount <- accCount + 1
    } else {
        Uk_S[, , s] <- Us
    }
}
accCount/(S_its-1)


# MH proposal: CMACG based on the previous sample -------------------------

# essentially a tuning parameter
# Dp <- diag(rep(c(100, 1), times = c(d, P-d)))
# Dp <- diag(rep(c(2, 1), times = c(P-d, d)))
Dp <- diag(P:1)
# tau_U <- 0.0000001

# Sampler based on CMACG proposal
Uk_S <- array(NA, c(P, d, S_its))
Uk_S[, , 1] <- runif_stiefel(P, d, betaf = 2)
# Uk_S[, , 1] <- Uk0
Us <- Uk_S[, , 1]

accCount <- 0
s <- 2
for (s in 2:S_its) {
    # calculate parameters for proposal
    Us_perp <- MASS::Null(Us)
    Ss <- cbind(Us, Us_perp)
    Sigmas <- Ss %*% Dp %*% t(Conj(Ss))
    # Sigmas <- Us %*% t(Conj(Us)) + tau_U * Us_perp %*% t(Conj(Us_perp))
    
    # propose Up
    Up <- rCMACG(P, d, Sigmas)
    
    # calculate terms based on priors
    tracep <- Re(sum(diag( Up %*% Omegak0 %*% t(Conj(Up)) %*% Yk )))
    logdetp <- log(Re(EigenR::Eigen_det( t(Conj(Up)) %*% invSigma0 %*% Up )))
    traces <- Re(sum(diag( Us %*% Omegak0 %*% t(Conj(Us)) %*% Yk )))
    logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*% invSigma0 %*% Us )))

    # calculate terms based on proposal    
    Up_perp <- MASS::Null(Up)
    Sp <- cbind(Up_perp, Up)
    Sigmap <- Sp %*% Dp %*% t(Conj(Sp))
    # Sigmap <- Up %*% t(Conj(Up)) + tau_U * Up_perp %*% t(Conj(Up_perp))
    
    (logr <- tracep/sigmak02 - P*logdetp - traces/sigmak02 + P*logdets +
            dCMACG(Us, Sigmap) - dCMACG(Up, Sigmas))
    
    if (log(runif(1)) <= logr) {
        Uk_S[, , s] <- Up
        Us <- Up
        accCount <- accCount + 1
    } else {
        Uk_S[, , s] <- Us
    }
} # end of s sampling loop
    
accCount/(S_its-1)


# Cayley transform --------------------------------------------------------

m <- 4

S <- matrix(0, m, m)
# S[upper.tri(S)] <- c(1, rep(0, choose(m, 2)-1))
# S[upper.tri(S)] <- c(1.5, -0.5, rep(0, choose(m, 2)-2))
# S[upper.tri(S)] <- rep(-0.1, choose(m, 2))
# S[upper.tri(S)] <- runif(choose(m, 2), -1, 1)
S[upper.tri(S)] <- runif(choose(m, 2), -1, 1) + runif(choose(m, 2), -1, 1) * 1i

for(j in 1:(m-1)) {
    for (i in (j+1):m) {
        S[i, j] <- -Conj(S[j, i])
    }
}
S
S %*% S
eigen(S)
EigenR::Eigen_det(S)
eigen(diag(m) + S)
EigenR::Eigen_det(diag(m) + S)
EigenR::Eigen_det(diag(m) - S)

U <- (diag(m) - S) %*% solve(diag(m) + S)
U
EigenR::Eigen_det(U)
# crossprod(U)
norm(diag(m) - U, type = "F")

Uother <- solve(diag(m) - S) %*% (diag(m) + S)
Uother

# U <- runif_stiefel(m, m, 2)
SU <- (diag(m) - U) %*% solve(diag(m) + U)
SU[upper.tri(SU)]


SUt <- (diag(m) - t(U)) %*% solve(diag(m) + t(U))
SUt[upper.tri(SUt)]
S[upper.tri(S)]

dnorm(S[upper.tri(S)])
dnorm(SUt[upper.tri(SUt)])


Uk <- runif_stiefel(m, 2, betaf = 2)

frob_dist(Uk, U %*% Uk)
frob_dist(diag(1, nrow = 4, ncol = 2), U %*% diag(1, nrow = m, ncol = 2))

# rotation matrix solutions -----------------------------------------------

set.seed(28052025)

# - solve for R in mat2 = Rs %*% mat1
# - t(mat2) = t(R %*% mat1) = t(mat1) %*% t(R)
# - solve(t(mat1), t(mat2)), then take transpose

mat1 <- runif_stiefel(4, 2, 1)
mat2 <- runif_stiefel(4, 2, 1)

R <- t( solve(t(mat1), t(mat2)) )

# Frobenius inner products ------------------------------------------------

U1 <- runif_stiefel(4, 2, 2)
U2 <- runif_stiefel(4, 2, 2)

sum(diag( t(Conj(U1 - U2)) %*% (U1 - U2) ))

sum(diag( t(Conj(U1)) %*% U2 ))

2*2 - 2*Re(sum(diag( t(Conj(U1)) %*% U2 )))


# Understand invariance of CMACG and/or Wishart ---------------------------

# calculate CMACG density for U, UQ, and RU
# where U is P by d
# Q is unitary d by d
# R is unitary P by P

P <- 4
d <- 2

set.seed(1062025)
U <- runif_stiefel(P, d, betaf = 2)
Q <- runif_stiefel(d, d, betaf = 2)
R <- runif_stiefel(P, P, betaf = 2)

UQ <- U %*% Q
RU <- R %*% U

Sigma <- rcomplex_wishart(P+1, P, diag(P))

Y <- rcomplex_wishart(P+1, P, diag(P))
LA <- diag(d:1)

# dCMACG <- function(X, Sigma, logscale = TRUE) {
dCMACG(U, Sigma)
dCMACG(UQ, Sigma) # equal to density for U, as expected (right unitary transformation invariant)
dCMACG(RU, Sigma) # not equal

# also, look at the projection matrices

U %*% t(Conj(U))
UQ %*% t(Conj(UQ))
RU %*% t(Conj(RU))

# what about sign flipping?
Sf <- diag(c(-1, 1))

dCMACG(U %*% Sf, Sigma) # equal to the density for U; antipodally symmetric?
dCMACG(UQ %*% Sf, Sigma) # equal to density for U, as expected
dCMACG(RU %*% Sf, Sigma) # not equal

# Complex Wishart densities
# but I should actually include the normalizing constants here! Since if the parameter changes, then the normalizing constant changes
dcomplex_wishart <- function(X, covPar, df, logscale = FALSE) {
    p <- nrow(X)
    invcovPar <- solve(covPar)
    
    logdens <- -mvGamma(df, p, 2, logscale = TRUE) +    
        -df * log( Re(EigenR::Eigen_det(covPar)) ) +
        (df-p) * log( Re(EigenR::Eigen_det(X)) ) + 
        -Re(sum(diag(invcovPar %*% X)))
        
    if (logscale) {
        return(logdens)
    } else {
        return(exp(logdens))
    }
}

Uf <- U %*% Sf
sign( Re(Uf[1, ]) )

Sc <- diag(complex(1, modulus = 1, argument = rnorm(2)))
Uc <- U %*% Sc

dcomplex_wishart(Y, U %*% LA %*% t(Conj(U)) + diag(P), P+1)
dcomplex_wishart(Y, UQ %*% LA %*% t(Conj(UQ)) + diag(P), P+1)
dcomplex_wishart(Y, RU %*% LA %*% t(Conj(RU)) + diag(P), P+1)

dcomplex_wishart(Y, Uf %*% LA %*% t(Conj(Uf)) + diag(P), P+1)
dcomplex_wishart(Y, Uc %*% LA %*% t(Conj(Uc)) + diag(P), P+1)

rowsigns <- sign( Re(Uf[1, ]) )

colflip <- ifelse(rowsigns == 0, 1, rowsigns)

t( t(Uf) * colflip )

firstRowPos(Uf) == U
Uf
U

# MH proposal: small rotation with Cayley transform and CMACG prior -------

###
# generate Up and calculate J(Up | Us)
###

# generate the upper triangular portions of S, as independent complex normal variables
# need to be able to modify their variance or standard deviation

# S upper-tris = rscnorm( choose(P, 2) ) * U_prop_sd
# apply the steps above to transform Sp to Rp, a rotation matrix

# the proposal is now Up <- Rp %*% Us

# J(Up | Us): calculate the complex normal density values of the upper-tri part of Sp

# complex normal densities


###
# calculate J(Us | Up)
###

# then, need to determine the density of the rotation matrix that would have generated Us, from Up
# - solve for Rs in Us = Rs %*% Up
# - t(Us) = t(Rs %*% Up) = t(Up) %*% t(Rs)
# - solve(t(Up), t(Us)), then take transpose

# transform Rs to Us
# extract upper triangular elements
# J(Us | Up): compute complex normal densities

# Sampler based on CMACG proposal
Uk_S <- array(NA, c(P, d, S_its))
Uk_S[, , 1] <- runif_stiefel(P, d, betaf = 2)
# Uk_S[, , 1] <- Uk0
Us <- Uk_S[, , 1]

accCount <- rep(TRUE, S_its)
s <- 2
for (s in 2:S_its) {
    if (s %% (S_its/100) == 0) print(paste0("s = ", s))
    
    # generate proposal rotation, Rp
    Sp <- matrix(0, P, P)
    Sp[upper.tri(Sp)] <- tau_U*rscnorm(choose(P, 2))
    
    if (doCayleyZeros) {
        Sp[upper.tri(Sp)] <- 
            Sp[upper.tri(Sp)] * rbinom(choose(P, 2), 1, CayleyZeroProb)
    }
    
    for(j in 1:(P-1)) {
        for (i in (j+1):P) {
            Sp[i, j] <- -Conj(Sp[j, i])
        }
    }
    Rp <- (diag(P) - Sp) %*% solve(diag(P) + Sp)
    Up <- Rp %*% Us
    
    # calculate terms based on priors
    tracep <- Re(sum(diag( Up %*% Omegak0 %*% t(Conj(Up)) %*% Yk )))
    logdetp <- log(Re(EigenR::Eigen_det( t(Conj(Up)) %*% invSigma0 %*% Up )))
    traces <- Re(sum(diag( Us %*% Omegak0 %*% t(Conj(Us)) %*% Yk )))
    logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*% invSigma0 %*% Us )))
    
    # calculate acceptance ratio
    (logr <- tracep/sigmak02 - P*logdetp - traces/sigmak02 + P*logdets)
    
    if (log(runif(1)) <= logr) {
        Uk_S[, , s] <- Up
        Us <- Up
        # accCount <- accCount + 1
        # accCount[s] <- TRUE
    } else {
        Uk_S[, , s] <- Us
        accCount[s] <- FALSE
    }
} # end of s sampling loop

# accCount/(S_its-1)
mean(accCount)
plot(accCount)
plot(cumsum(accCount))

# show the sampled matrices have orthonormal columns
t(Conj(Uk_S[, , s])) %*% Uk_S[, , s]

itsKeep <- (S_its/2 + 1):S_its
# itsKeep <- 5000:S_its

projSum <- matrix(0, P, P)
covSum <- matrix(0, P, P)
avgUkSum <- matrix(0, P, d)
for (s in itsKeep) {
    N <- length(itsKeep)
    projSum <- projSum + 
        (1/N) * Uk_S[, , s] %*% t(Conj(Uk_S[, , s]))
    
    covSum <- covSum +
        (1/N) * Uk_S[, , s] %*% Lambdak0 %*% t(Conj(Uk_S[, , s]))
    
    avgUkSum <- avgUkSum + 
        (1/N) * std_cmpx_evec(Uk_S[, , s])
}
avgUkSum
t(Conj(avgUkSum)) %*% avgUkSum

avgUkSum_Q <- qr(avgUkSum) |> qr.Q() |> std_cmpx_evec()
cbind(avgUkSum[, 1], avgUkSum_Q[, 1])

# eigenvectors of average projection and covariance matrices
evec_projSum <- eigen(projSum)$vectors[, 1:d]
evec_covSum <- eigen(covSum)$vectors[, 1:d]

Uk0_true <- eigen(Uk0 %*% Lambdak0 %*% t(Conj(Uk0)))$vectors[, 1:d]
Uk0_std <- std_cmpx_evec(Uk0)

frob_dist(avgUkSum_Q, Uk0_std, returnDists = TRUE)

# actual average of Uks?
avgUk <- apply(Uk_S[, , itsKeep], c(1, 2), mean)

std_Uk_S <- array(NA, c(P, d, length(itsKeep)))

# it is not even truly semi-unitary
t(Conj(avgUk)) %*% avgUk

avgUkQ <- qr(avgUk) |> qr.Q()
t(Conj(avgUkQ)) %*% avgUkQ

# the "true" eigenvectors are different
cbind(Uk0[, 1], Uk0_true[, 1])

# the "true" projection matrices are the same!
cbind((Uk0 %*% t(Conj(Uk0)))[, 1],
      (Uk0_true %*% t(Conj(Uk0_true)))[, 1])

coli <- 3
# look at "average" eigenvectors and true eigenvectors
cbind(Uk0[, coli], evec_covSum[, coli], evec_projSum[, coli])
cbind(Uk0_true[, coli], evec_covSum[, coli], evec_projSum[, coli])
cbind(Uk0_true[, coli], evec_covSum[, coli], evec_projSum[, coli], avgUk[, coli])

### compare Frobenius distance between Projection eigenvectors
# standardized by first row
frob_dist(std_cmpx_evec(evec_projSum), std_cmpx_evec(Uk0))
frob_dist(std_cmpx_evec(evec_projSum), std_cmpx_evec(Uk0_true), returnDists = TRUE)
# raw eigenvectors
frob_dist(evec_projSum, Uk0)

### compare Frobenius distance between covariance eigenvectors
# standardized by first row
frob_dist(std_cmpx_evec(evec_covSum), std_cmpx_evec(Uk0), returnDists = TRUE)
frob_dist(std_cmpx_evec(evec_covSum), std_cmpx_evec(Uk0_true))

cbind(
    std_cmpx_evec(evec_covSum)[, 2], 
    std_cmpx_evec(Uk0)[, 2]
)

# raw eigenvectors
frob_dist(evec_covSum, Uk0)

### compare Frobenius distance between covariance matrices
# same for either subspace representation of Uk0
frob_dist(covSum, Uk0 %*% Lambdak0 %*% t(Conj(Uk0)))
frob_dist(covSum, Uk0_true %*% Lambdak0 %*% t(Conj(Uk0_true)))

### compare Frobenius distance between projection matrices
# same, since Uk0 and Uk0_true form the same subspace
frob_dist(projSum, Uk0 %*% t(Conj(Uk0)))
frob_dist(projSum, Uk0_true %*% t(Conj(Uk0_true)))
# the "average" of projection matrices is not technically a proper projection matrix.
# what if I take the "subspace" of the average projection matrix?
frob_dist(evec_projSum %*% t(Conj(evec_projSum)), Uk0 %*% t(Conj(Uk0)))
frob_dist(evec_projSum %*% t(Conj(evec_projSum)), Uk0_true %*% t(Conj(Uk0_true)))

# how far is the "average subspace" from the "average projection matrix"
frob_dist(evec_projSum %*% t(Conj(evec_projSum)), projSum)

# the "average subspace" seems to be rank d
eigen(evec_projSum %*% t(Conj(evec_projSum)))$values

# look at columns of the projection matrices
cbind(
    (evec_projSum %*% t(Conj(evec_projSum)))[, 1],
    (Uk0 %*% t(Conj(Uk0)))[, 1]
)

# Procrustes average of eigenvectors --------------------------------------

### aligned to true Uk0 - but is that step legitimate?
Eprimek <- array(NA, c(P, d, length(itsKeep)))

Ehat_c <- Uk_S[, , itsKeep[1]]
s <- 2

for (s in 1:length(itsKeep)) {
    Mk <- t(Conj(Uk_S[, , itsKeep[s]])) %*% Ehat_c
    svdMk <- svd(Mk)
    Rk <- svdMk$u %*% t(Conj(svdMk$v))
    
    Eprimek[, , s] <- Uk_S[, , itsKeep[s]] %*% Rk
}

Araw <- apply(Eprimek, c(1, 2), mean)
svdAraw <- svd(Araw)
Ehat_n <- svdAraw$u %*% t(Conj(svdAraw$v))

frob_dist(Ehat_c, Ehat_n)

Ehat_c <- Ehat_n

cbind(std_cmpx_evec(Ehat_c)[, 1], Uk0_std[, 1])

frob_dist((Ehat_c), Uk0_std)
frob_dist(std_cmpx_evec(Ehat_c), Uk0_std)

Maligned <- t(Conj(Ehat_c)) %*% Uk0_std
svdMaligned <- svd(Maligned)
Raligned <- svdMaligned$u %*% t(Conj(svdMaligned$v))
Ealigned <- Ehat_c %*% Raligned

frob_dist(Ealigned, Uk0_std)
frob_dist(std_cmpx_evec(Ealigned), Uk0_std)

### natural mean, from Chikuse 2012, p. 237

Uk_S_std <- array(NA, c(P, d, S_its))
for (s in 1:S_its) {
    Uk_S_std[, , s] <- std_cmpx_evec(Uk_S[, , s])
}

S <- apply(Uk_S[, , itsKeep], c(1, 2), mean)
# S <- apply(Uk_S_std[, , itsKeep], c(1, 2), mean)
HS <- S %*% EigenR::Eigen_sqrt(solve(t(Conj(S)) %*% S))
QS <- qr(S) |> qr.Q()

frob_dist(HS, Uk0, returnDists = TRUE)
frob_dist(std_cmpx_evec(HS), Uk0_std, returnDists = TRUE)

cbind(HS[, 2], Uk0[, 2])
cbind(std_cmpx_evec(HS)[, 2], Uk0_std[, 2])

# compare distances to the "mean"
dists_to_mean <- rep(NA, S_its)
for (s in 1:S_its) {
    dists_to_mean[s] <- frob_dist(HS, Uk_S[, , s])
    # dists_to_mean[s] <- frob_dist(diag(1, P, d), std_cmpx_evec(Uk_S[, , s]))
}
quantile(dists_to_mean[itsKeep], c(.025, .25, .5, .75, .975))

plot(dists_to_mean[seq(1, S_its, by = 50)], type = "l")

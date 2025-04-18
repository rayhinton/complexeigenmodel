# issue: 2x2 complex Bingham sampler sometimes has high rejection rates

library(rstiefel)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

its <- 1000

C <- diag(c(2, 1))

Aevals <- c(80, 20)

(Av <- 1/sqrt(2) * matrix(c(1, 1, -1, 1), ncol = 2))
(A <- Av %*% diag(Aevals) %*% t(Av))
eigen(A)

B <- diag(c(4, 1))

# observe that the number of rejections is high
(X <- my.rCbing.Op(A, B, istatus = 100))

X$X %*% t(Conj(X$X))
X$X %*% C %*% t(Conj(X$X))

Ucs <- array(NA, dim = c(2, 2, its))
Uc_fs <- array(NA, dim = c(2, 2, its))
Covs <- array(NA, dim = c(2, 2, its))
nrejs <- rep(NA, its)

# reference vector
ur <- matrix(c(1 +0i, 0+0i), ncol = 1)
set.seed(3042025)
for (i in 1:its) {
    # generate a random matrix
    Usamp <- my.rCbing.Op(A, B)
    U <- Usamp$X
    nrejs[i] <- Usamp$nrej
    # store the raw sample
    Ucs[, , i] <- U
    
    # flip the column signs, relative to reference vector
    if (Re(t(U[, 1]) %*% ur) < 0) {
        U[, 1] <- -U[, 1]
    }
    if (Re(t(U[, 2]) %*% ur) < 0) {
        U[, 2] <- -U[, 2]
    }
    
    # store the modified sample
    Uc_fs[, , i] <- U
    
    # store a "covariance" matrix from the sampled value
    Covs[, , i] <- U %*% C %*% t(Conj(U))
}

(Uc_fs_mean <- apply(Uc_fs, c(1, 2), mean))
Covs_mean <- apply(Covs, c(1, 2), mean)
eigen(Covs_mean)
Av

t(Conj(Uc_fs_mean)) %*% Uc_fs_mean

summary(nrejs)
quantile(nrejs, c(.95, .99))

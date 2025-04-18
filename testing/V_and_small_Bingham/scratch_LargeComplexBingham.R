# sampling for V matrices

# Assuming that V has even dimensions

# library(cmvnorm)
library(rstiefel)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")

# an example of orthogonal 2D complex vectors
x <- matrix(c(1, 1i), ncol = 1)
y <- matrix(c(-1, 1i), ncol = 1)

t(Conj(x)) %*% y
t(Conj(y)) %*% x

# begin sampling

P <- 8

set.seed(9032025)
Av <- runitary(P, P)
Aevals <- (P:1)
# A <- rcomplex_wishart(P, P, diag(P))
A <- Av %*% diag(Aevals) %*% t(Conj(Av))
# B <- runif(P, 0, 10) |> sort(decreasing = TRUE)
B <- (P:1)

eigen(A)$values
eigen(A)
B

# generate a random initial matrix, X, which should have orthonormal columns
X <- (rcomplex_wishart(P, P, diag(P)) |> eigen())$vector
# the matrix has orthonormal columns
t(Conj(X)) %*% X

# sample matrix using the function
newX <- rcBingUP_gibbs(X, A, B)
    
# the matrix is orthonormal
zapsmall(t(Conj(newX)) %*% newX)
zapsmall(newX %*% t(Conj(newX)))

# compare sampled columns to a reference vector
ur <- matrix(c(1, rep(0, P-1)), nrow = 1)
Re(ur %*% newX)
sign(Re(ur %*% newX))

# sample an X many times
gibbsIts <- 2500
Xs <- array(NA, c(P, P, gibbsIts))
flipped_Xs <- array(NA, c(P, P, gibbsIts))
# array to hold some projection matrices
projs <- array(NA, c(P, P, gibbsIts))
Covs <- array(NA, c(P, P, gibbsIts))

set.seed(25032025)
# generate a random initial matrix, X, which should have orthonormal columns
Xs[, , 1] <- (rcomplex_wishart(P, P, diag(P)) |> eigen())$vector
for (i in 2:gibbsIts) {
    if (i %% 100 == 0) print(i)
    
    X <- rcBingUP_gibbs(Xs[, , i-1], A, B, istatus = 0, 
                        Imtol = 10^ceiling(log10(.Machine$double.eps)))
    Xs[, , i] <- X
    
    projs[, , i] <- X[, 1:(P-1)] %*% t(Conj(X[, 1:(P-1)]))
    # flipped_Xs
    colmod <- sign(Re(ur %*% X))
    flipped_Xs[, , i] <- t(t(X) * c(colmod))
    
    Covs[, , i] <- X %*% diag(P:1) %*% t(Conj(X))
    
    # colmod <- sign(t(ur) %*% U)
    # U <- t(t(U) * c(colmod))
    # 
    # U_flipped[, i] <- colmod == -1
    # 
    # # store the modified sample
    # U_fs[, , i] <- U
}

(eigen(A)$vectors[, 1:(P-1)] %*% t(Conj(eigen(A)$vectors[, 1:(P-1)])))[, 3]
apply(projs[, , (gibbsIts/2):gibbsIts], c(1, 2), mean)[, 3]

(meanflipped <- apply(flipped_Xs[, , (gibbsIts/2):gibbsIts], c(1, 2), mean))
eigen(A)

# I think this is the best to compare - "average" eigenvectors
avgCovs <- eigen(apply(Covs[, , (gibbsIts/2):gibbsIts], c(1, 2), mean))$vectors
Avecs <- eigen(A)$vectors

coli <- 4
cbind(avgCovs[, coli], Avecs[, coli])

Re(diag(t(Conj(meanflipped)) %*% meanflipped))

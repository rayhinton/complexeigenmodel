# scratch - matrix distances or statistics

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")

frob_norm <- function(X) {
    outdiag <- Re(diag(crossprod(Conj(X), X)))
    outlist <- list(fnorm = sum(outdiag), diags = outdiag)
    
    return(outlist)
}

frame_distance <- function(A, B) {
    # return(norm(A - B, "F"))
    return(frob_norm(A - B))
}

procrustes_distance <- function(A, B) {
    # Get orthogonal bases
    Q_A <- qr.Q(qr(A))
    Q_B <- qr.Q(qr(B))
    
    # Find optimal rotation between the frames
    svd_result <- svd(t(Conj(Q_B)) %*% Q_A)
    R_opt <- svd_result$v %*% t(Conj(svd_result$u))
    
    # Compute distance after optimal alignment
    # return(norm(Q_A - Q_B %*% R_opt, "F"))
    return(frob_norm(Q_A - Q_B %*% R_opt))
}

# Description: calculate Grassmann Distance between subspaces

# Inputs:
# A, p x k matrix, semi-unitary
# B, p x k matrix, semi-unitary
# r, the number of columns of A and B to use for calculating subspaces
# s_tol, tolerance to use for checking if singular values are greater than 1

# Output:
# sub_dist, the Grassmann distance between the spaces spanned by A[, 1:d] and B[, 1:d]

grass_dist <- function (A, B, r = ncol(A), s_tol = 2*.Machine$double.eps) {
    stopifnot("A and B must have the same dimensions" = all(dim(A) == dim(B)))
    
    Sigma <- t(Conj(A)) %*% B
    svd_Sigma <- svd(Sigma)
    
    # check for numerical issues with SVDs
    max_diff <- max(svd_Sigma$d[1:r] - 1)
    if (any(svd_Sigma$d[1:r] > (1 + s_tol))) {
        # warning()
        stop(paste0("Singular values exceed 1 by more than s_tol = ", s_tol, 
                    ". Maximum excess: ", max_diff))
    } else {
        sigmas <- pmin(svd_Sigma$d[1:r], 1)
    }
    
    thetas <- acos(sigmas)
    sq_thetas <- thetas^2
    sub_dist <- sqrt(sum(sq_thetas))
    
    outlist <- list(sq_thetas = sq_thetas, sub_dist = sub_dist)
    
    return(outlist)
}

# test the distance between two 2D subspaces

# this is A from ex. site below
twoD_1 <- matrix(c(1,0,0,0,1,0), ncol = 2, byrow = FALSE)
twoD_2 <- 1/sqrt(2) * matrix(c(1,1,0,1,-1,0), ncol = 2, byrow = FALSE)
twoD_3 <- 1/sqrt(2) * matrix(c(0,1,1,0,-1,1), ncol = 2, byrow = FALSE)

# both matrices have orthonormal columns
crossprod(twoD_1)
crossprod(twoD_2)
crossprod(twoD_3)

# two of the distances are the same, one is different
grass_dist(twoD_1, twoD_2)
grass_dist(twoD_1, twoD_3)
grass_dist(twoD_2, twoD_3)

# are the other distances the same, if one of the column signs are flipped?
t(t(twoD_1) * c(-1, 1))
grass_dist(t(t(twoD_1) * c(-1, 1)), twoD_3)

# the Grassmannian distance between 1 and 2 was essentially 0
# but it is nonzero for the Procrustes Distance
procrustes_distance(twoD_1, twoD_2)
procrustes_distance(t(t(twoD_1) * c(-1, 1)), twoD_2)

procrustes_distance(twoD_1, twoD_3)
procrustes_distance(t(t(twoD_1) * c(-1, 1)), twoD_3)

# Frobenius Norm - not the same when a sign of one column is flipped
norm(twoD_1 - twoD_2, "F")
norm(t(t(twoD_1) * c(-1, 1)) - twoD_2, "F")

frob_norm(twoD_1 - twoD_2)
frob_norm(t(t(twoD_1) * c(-1, 1)) - twoD_2)

# this is B from ex. site below
twoD_4 <- matrix(c(1/sqrt(2), 0, 1/sqrt(2), 0, 1, 0), ncol = 2, byrow = FALSE)
crossprod(twoD_4)

grass_dist(twoD_1, twoD_4)

twoD_5 <- matrix(c(0,0,1,0,1,0), ncol = 2, byrow = FALSE)
crossprod(twoD_5)

grass_dist(twoD_4, twoD_5)
grass_dist(twoD_1, twoD_5)

grass_dist(twoD_4, twoD_4)

# examples from https://jyopari.github.io/posts/grassman

# this is C from example previously
twoD_6 <- matrix(c(0, 1, 0, 0, 0, 1), ncol = 2, byrow = FALSE)

grass_dist(twoD_1, twoD_6)

# generate parameters to the Complex Matrix Bingham distribution

sum(diag(crossprod(twoD_1 - twoD_5)))
sum(diag(crossprod(twoD_1 - twoD_5)))

#####
# 1d projection matrices
#####

# is the average of projection matrices, a projection matrix?

# show two different semi-orthogonal matrices
twoD_1
twoD_2

# compute the projection matrices of their first columns
A1 <- tcrossprod(twoD_1[, 1])
A2 <- tcrossprod(twoD_2[, 1])

# take the average of the two different projection matrices
avgA <- (A1 + A2)/2
    
    # is it symmetric? Yes
isSymmetric(avgA)
    # is it idempotent? No

avgA %*% avgA

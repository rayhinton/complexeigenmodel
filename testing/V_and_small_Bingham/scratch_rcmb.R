# test rcmb function

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")


# distance functions ------------------------------------------------------

frob_norm <- function(X) {
    return(Re(sum(diag(crossprod(Conj(X), X)))))
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
    sub_dist <- sqrt(sum(thetas^2))
    
    return(sub_dist)
}

# set up sampler testing --------------------------------------------------

# rcmb(U, A, B)

# - rcmb is made for steps in a Gibbs sampler. 
# - you must supply U, the previous sampled value, and then columns are sampled from FCD

P <- 8
d <- 4
its <- 2500

set.seed(9042025)
# generate an 8x8 parameter G
G <- rcomplex_wishart(P+1, P, diag(P))
Gother <- rcomplex_wishart(P+1, P, diag(P:1))

# distance between G and Gother
grass_dist(eigen(G)$vectors[, 1:d], eigen(Gother)$vectors[, 1:d])
procrustes_distance(eigen(G)$vectors[, 1:d], eigen(Gother)$vectors[, 1:d])
frame_distance(eigen(G)$vectors[, 1:d], eigen(Gother)$vectors[, 1:d])

eigen(G)$values

eigen(G)$vectors[, 1]
eigen(Gother)$vectors[, 1]

S <- solve(diag(ceiling(max(eigen(G)$values)), nrow = P) - G)
eigen(S)$values

W <- rcomplex_wishart(P+1, P, S)
U <- eigen(W)$vectors
L <- diag(eigen(W)$values)

U1 <- U[, 1:d]
L1 <- L[1:d, 1:d]

# sampling ----------------------------------------------------------------

# Generate multiple samples from rcmb(previousU, G, L1)
Us <- array(NA, c(P, d, its))
Covs <- array(NA, c(P, P, its))

# first test - if the sampler is initialized with a sample from the distribution, does it converge well?
# Us[, , 1] <- U1
# Us[, , 1] <- runitary(P, d)
Us[, , 1] <- diag(1, nrow = P, ncol = d)
Covs[, , 1] <- Us[, , 1] %*% diag(4:1) %*% t(Conj(Us[, , 1]))
# TODO later test - if the sample is initialized with a random semi-unitary matrix, does it converge well?

set.seed(8042025)
for (s in 2:its) {
    if (s %% 500 == 0) print(paste0("s = ", s))
    # Us[, , s] <- rcmb(Us[, , s-1], G, L1)
    Us[, , s] <- rcmb(Us[, , s-1], Gother, L1)
    Covs[, , s] <- Us[, , s] %*% diag(4:1) %*% t(Conj(Us[, , s]))
    
}

# compare a summary of samples to certain known parameter -----------------

avgCovs <- apply(Covs[, , (its/2) : its], c(1, 2), mean)

Vhat <- eigen(avgCovs)$vectors[, 1:d]
Vtilde <- eigen(G)$vectors[, 1:d]
Vtilde_other <- eigen(Gother)$vectors[, 1:d]

# distances from observed "average" to incorrect or correct parameters
dists <- data.frame(
    eigenG = c(grass_dist(Vhat, Vtilde), procrustes_distance(Vhat, Vtilde),
               frame_distance(Vhat, Vtilde)),
    eigenGother = c(grass_dist(Vhat, Vtilde_other), 
                    procrustes_distance(Vhat, Vtilde_other),
                    frame_distance(Vhat, Vtilde_other)))
rownames(dists) <- c("grass", "procrustes", "frobenius")
dists

# compare individual columns
coli <- 1
cbind(
    Vhat[, coli],
    Vtilde[, coli],
    Vtilde_other[, coli]
)

# TODO compare the projection matrices too, 
# either full projection matrices, or those created from a single column at a time

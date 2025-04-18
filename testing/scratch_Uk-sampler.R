# test Uk FCD sampler

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")

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

# the columns, conditional on others, have certain vector Bingham distribution

# Uk | A, B, V, Lambdak, Pk
# where Pk | Uk, Lambdak, sigmak2 ~ Complex Wishart

# - need true Uk0 matrices, not necessarily from any particular distribution
# - do need Pk generated from correct CW distribution

P <- 8
d <- 4
its <- 10000

nk <- 1000

set.seed(12042025)

A <- rgamma(P, 1, 1) |> sort(decreasing = TRUE) |> diag()
B <- rgamma(d, 1, 1) |> sort(decreasing = TRUE) |> diag()
V <- runitary(P, P)
VAVH <- V %*% A %*% t(Conj(V))

Lambdak <- rgamma(d, 1, 1) |> sort(decreasing = TRUE) |> diag()
sigmak2 <- rgamma(1, 1, 1)

Uk0 <- runitary(P, d)
Gammak0 <- sigmak2 * ( Uk0 %*% Lambdak %*% t(Conj(Uk0)) + diag(P) )

Uk0_eigenvectors <- eigen(Uk0 %*% Lambdak %*% t(Conj(Uk0)))$vectors[, 1:d]

cbind(Uk0[, 1], Uk0_eigenvectors[, 1])

Yk <- rcomplex_wishart(nk, P, Gammak0)

# initialize a sample for Uk
set.seed(13042025)

Uks <- array(NA, c(P, d, its))
avgCovs <- array(NA, c(P, P, its))
# Uks[, , 1] <- runitary(P, d)
# Uks[, , 1] <- diag(1+0i, nrow = P, ncol = d)
# try an initialization that is based on the data
Uks[, , 1] <- eigen((nk*sigmak2)^-1 * Yk - diag(P))$vectors[, 1:d]

for (s in 2:its) {
    if (s %% 500 == 0) print(paste0("s = ", s))
    Ukc <- Uks[, , s-1]
    
    for (j in sample(d)) {
        bj <- B[j, j]
        lambdajk <- Lambdak[j, j]
        
        # TODO VAV^H could be calculated beforehand - the value does not change
        parj <- bj * VAVH + lambdajk/(1 + lambdajk) * Yk
        
        Ukc[, j] <- rcvb_LN(Ukc, j, parj)
    }
    
    Uks[, , s] <- Ukc
    avgCovs[, , s] <- Ukc %*% Lambdak %*% t(Conj(Ukc))
}

# nonburn <- seq(its/2, its, by = 10)
nonburn <- seq(its/2, its, by = 1)

Ukhat <- eigen(apply(avgCovs[, , nonburn], c(1, 2), mean))$vectors[, 1:d]

coli <- 1
cbind(Ukhat[, coli],
      Uk0[, coli])
cbind(Ukhat[, coli],
      Uk0_eigenvectors[, coli])

# Grassmannian distance between the subspaces
grass_dist(Ukhat, Uk0)
# Frobenius norm of difference between the projection matrices
frob_norm(Ukhat %*% t(Ukhat) - Uk0 %*% t(Conj(Uk0)))
# Distance between the eigenvector matrices
frame_distance(Ukhat, Uk0)
# Distance between the eigenvectors in terms of Procustes statistic
procrustes_distance(Ukhat, Uk0)

# Grassmannian distance between the subspaces
grass_dist(Ukhat, Uk0_eigenvectors)
# Frobenius norm of difference between the projection matrices
frob_norm(Ukhat %*% t(Ukhat) - Uk0_eigenvectors %*% t(Conj(Uk0_eigenvectors)))
# Distance between the eigenvector matrices
frame_distance(Ukhat, Uk0_eigenvectors)
# Distance between the eigenvectors in terms of Procustes statistic
procrustes_distance(Ukhat, Uk0_eigenvectors)

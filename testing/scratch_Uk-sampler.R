# test Uk FCD sampler

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/matrix_distances.R")

Uk_gibbs_densCovar <- function(data_list, param_list) {
    VAVH <- param_list$Vs %*% param_list$As %*% t(Conj(param_list$Vs))
    
    U_ks <- param_list$U_ks
    
    for (k in 1:K) {
        # Ukc <- Uks[, , k, s-1]
        Ukc <- param_list$U_ks[, , k]
        
        for (j in sample(d)) {
            # bj <- B[j, j]
            bj <- param_list$Bs[j, j]
            # lambdajk <- Lambdaks[j, j, k]
            lambdajk <- param_list$Lambda_ks[j, j, k]
            
            # TODO VAV^H could be calculated beforehand - the value does not change
            # parj <- bj * VAVH + lambdajk/(1 + lambdajk) * Yks[, , k]
            parj <- bj * VAVH + lambdajk/(1 + lambdajk) * data_list[[k]]
            
            Ukc[, j] <- rcvb_LN(Ukc, j, parj)
        }
        
        U_ks[, , k] <- Ukc
    }
    
    return(U_ks)
} 

# frob_norm <- function(X) {
#     return(Re(sum(diag(crossprod(Conj(X), X)))))
# }
# 
# frame_distance <- function(A, B) {
#     # return(norm(A - B, "F"))
#     return(frob_norm(A - B))
# }
# 
# procrustes_distance <- function(A, B) {
#     # Get orthogonal bases
#     Q_A <- qr.Q(qr(A))
#     Q_B <- qr.Q(qr(B))
#     
#     # Find optimal rotation between the frames
#     svd_result <- svd(t(Conj(Q_B)) %*% Q_A)
#     R_opt <- svd_result$v %*% t(Conj(svd_result$u))
#     
#     # Compute distance after optimal alignment
#     # return(norm(Q_A - Q_B %*% R_opt, "F"))
#     return(frob_norm(Q_A - Q_B %*% R_opt))
# }
# 
# grass_dist <- function (A, B, r = ncol(A), s_tol = 2*.Machine$double.eps) {
#     stopifnot("A and B must have the same dimensions" = all(dim(A) == dim(B)))
#     
#     Sigma <- t(Conj(A)) %*% B
#     svd_Sigma <- svd(Sigma)
#     
#     # check for numerical issues with SVDs
#     max_diff <- max(svd_Sigma$d[1:r] - 1)
#     if (any(svd_Sigma$d[1:r] > (1 + s_tol))) {
#         # warning()
#         stop(paste0("Singular values exceed 1 by more than s_tol = ", s_tol, 
#                     ". Maximum excess: ", max_diff))
#     } else {
#         sigmas <- pmin(svd_Sigma$d[1:r], 1)
#     }
#     
#     thetas <- acos(sigmas)
#     sub_dist <- sqrt(sum(thetas^2))
#     
#     return(sub_dist)
# }

# the columns, conditional on others, have certain vector Bingham distribution

# Uk | A, B, V, Lambdak, Pk
# where Pk | Uk, Lambdak, sigmak2 ~ Complex Wishart

# - need true Uk0 matrices, not necessarily from any particular distribution
# - do need Pk generated from correct CW distribution

K <- 3
P <- 8
d <- 4
its <- 5000

nks <- rep(1000, K)

set.seed(12042025)

A <- rgamma(P, 1, 1) |> sort(decreasing = TRUE) |> diag()
B <- rgamma(d, 1, 1) |> sort(decreasing = TRUE) |> diag()
V <- runitary(P, P)
# VAVH <- V %*% A %*% t(Conj(V))

# Lambdak <- rgamma(d, 1, 1) |> sort(decreasing = TRUE) |> diag()
Lambdaks <- array(NA, c(d, d, K))
# sigmak2 <- rgamma(1, 1, 1)
sigmak2s <- rep(NA, K)

# Uk0 <- runitary(P, d)
Uk0s <- array(NA, c(P, d, K))
# Gammak0 <- sigmak2 * ( Uk0 %*% Lambdak %*% t(Conj(Uk0)) + diag(P) )

# Uk0_eigenvectors <- eigen(Uk0 %*% Lambdak %*% t(Conj(Uk0)))$vectors[, 1:d]
Uk0s_evecs <- array(NA, c(P, d, K))

# cbind(Uk0[, 1], Uk0_eigenvectors[, 1])

set.seed(22042025)
# Yk <- rcomplex_wishart(nk, P, Gammak0)
Yks <- array(NA, c(P, P, K))
data_list <- list()

Uks <- array(NA, c(P, d, K, its))

for (k in 1:K) {
    Lambdaks[, , k] <- rgamma(d, 1, 1) |> sort(decreasing = TRUE) |> diag()
    sigmak2s[k] <- rgamma(1, 1, 1)
    Uk0s[, , k] <- runitary(P, d)
    
    Gammak0 <- sigmak2s[k] * 
        ( Uk0s[, , k] %*% Lambdaks[, , k] %*% t(Conj(Uk0s[, , k])) + diag(P) )
    
    Yks[, , k] <- rcomplex_wishart(nks[k], P, Gammak0)
    data_list[[k]] <- Yks[, , k]
    
    # random initialization for each
    # Uks[, , k, 1] <- runitary(P, d)
    # ideal initialization at the exact true Uk0
    Uks[, , k, 1] <- Uk0s[, , k]
    
    Uk0s_evecs[, , k] <- eigen(Uk0s[, , k] %*% Lambdaks[, , k] %*% t(Conj(Uk0s[, , k])) )$vectors[, 1:d]
}

param_list <- list(
    P = P,
    d = d,
    K = K,
    n_k = nks,
    U_ks = Uks[, , , 1],
    Lambda_ks = Lambdaks,
    sigma_k2s = sigmak2s,
    Vs = V,
    As = A,
    Bs = B
)

# initialize a sample for Uk
set.seed(13042025)

# Uks <- array(NA, c(P, d, K, its))
avgCovs <- array(NA, c(P, P, K, its))
# Uks[, , 1] <- runitary(P, d)
# Uks[, , 1] <- diag(1+0i, nrow = P, ncol = d)
# try an initialization that is based on the data
# Uks[, , 1] <- eigen((nk*sigmak2)^-1 * Yk - diag(P))$vectors[, 1:d]

{
print(Sys.time())
for (s in 2:its) {
    if (s %% 500 == 0) print(paste0("s = ", s))
    
    # VAVH <- param_list$Vs %*% param_list$As %*% t(Conj(param_list$Vs))
    # 
    # for (k in 1:K) {
    #     # Ukc <- Uks[, , k, s-1]
    #     Ukc <- param_list$U_ks[, , k]
    #     
    #     for (j in sample(d)) {
    #         # bj <- B[j, j]
    #         bj <- param_list$Bs[j, j]
    #         # lambdajk <- Lambdaks[j, j, k]
    #         lambdajk <- param_list$Lambda_ks[j, j, k]
    #         
    #         # TODO VAV^H could be calculated beforehand - the value does not change
    #         # parj <- bj * VAVH + lambdajk/(1 + lambdajk) * Yks[, , k]
    #         parj <- bj * VAVH + lambdajk/(1 + lambdajk) * data_list[[k]]
    #         
    #         Ukc[, j] <- rcvb_LN(Ukc, j, parj)
    #     }
    #     
    #     Uks[, , k, s] <- Ukc
    #     param_list$U_ks[, , k] <- Ukc
    #     avgCovs[, , k, s] <- Ukc %*% Lambdaks[, , k] %*% t(Conj(Ukc))
    #     
    # }
    
    Uks[, , , s] <- Uk_gibbs_densCovar(data_list, param_list)
    param_list$U_ks <- Uks[, , , s]
    for (k in 1:K) {
        avgCovs[, , k, s] <- Uks[, , k, s] %*% Lambdaks[, , k] %*% t(Conj(Uks[, , k, s]))
    }
}
print(Sys.time())
}

# nonburn <- seq(its/2, its, by = 10)
# nonburn <- seq(its/2, its, by = 1)
nonburn <- seq(its/2, its, by = 1)

Ukhats <- array(NA, c(P, d, K))
for (k in 1:K) {
    Ukhats[, , k] <- eigen(apply(avgCovs[, , k, nonburn], c(1, 2), mean))$vectors[, 1:d]
}

# Ukhat <- eigen(apply(avgCovs[, , nonburn], c(1, 2), mean))$vectors[, 1:d]

k <- 3
coli <- 4
cbind(Ukhats[, coli, k],
      Uk0s[, coli, k])
cbind(Ukhats[, coli, k],
      Uk0s_evecs[, coli, k])

Ukhat <- Ukhats[, , k]
Uk0 <- Uk0s[, , k]
Uk0_evecs <- Uk0s_evecs[, , k]

# Grassmannian distance between the subspaces
grass_dist(Ukhat, Uk0)
# Frobenius norm of difference between the projection matrices
frob_norm(Ukhat %*% t(Ukhat) - Uk0 %*% t(Conj(Uk0)))
# Distance between the eigenvector matrices
frame_distance(Ukhat, Uk0)
# Distance between the eigenvectors in terms of Procustes statistic
procrustes_distance(Ukhat, Uk0)

# Grassmannian distance between the subspaces
grass_dist(Ukhat, Uk0_evecs)
# Frobenius norm of difference between the projection matrices
frob_norm(Ukhat %*% t(Ukhat) - Uk0_evecs %*% t(Conj(Uk0_evecs)))
# Distance between the eigenvector matrices
frame_distance(Ukhat, Uk0_evecs)
# Distance between the eigenvectors in terms of Procustes statistic
procrustes_distance(Ukhat, Uk0_evecs)

procrustes_distance(Uk0, Uk0_evecs)
grass_dist(Uk0, Uk0_evecs)

# -------------------------------------------------------------------------

# a "credible" interval?
# calculate distance of sampled matrices to the "average" matrix

# in order to do this right, I need a "distance" which properly 

dists_Uks_Ukhat <- rep(NA, length(nonburn))

for (i in 1:length(nonburn)) {
    # dists_Uks_Ukhat[i] <- procrustes_distance(Uks[, , nonburn[i]], Ukhat)$fnorm
    dists_Uks_Ukhat[i] <- grass_dist(Uks[, , nonburn[i]], Ukhat)$sub_dist
}

quantile(dists_Uks_Ukhat, probs = c(.025, .975))

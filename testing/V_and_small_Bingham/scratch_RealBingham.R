# comparing real Bingham samplers

library(rstiefel)
library(simdd)

P <- 8
its <- 2500

set.seed(17042025)

H <- diag(seq(10, 1, length.out = P))

G <- rWishart(1, P + 1, diag(P))[, , 1]
A <- diag(eigen(G)$values)
V <- eigen(G)$vectors

X <- matrix(rnorm(P*P), ncol = P) |> qr() |> qr.Q()

crossprod(X)
tcrossprod(X)

# Gibbs sampler
rbing.matrix.gibbs(G, H, X)
# rejection sampler
if (P <= 4) {
    rbing.Op(G, H)
}


# test out 1d projection matrices -----------------------------------------

tcrossprod(V[, 1])

# set up sampling ---------------------------------------------------------

Us_rst <- array(NA, c(P, P, its))
Covs_rst <- array(NA, c(P, P, its))

if (P <= 4) {
    Us_rj_rst <- array(NA, c(P, P, its))
    Covs_rj_rst <- array(NA, c(P, P, its))
}
# initialize the first entry as a random orthnormal matrix
Us_rst[, , 1] <- matrix(rnorm(P*P), ncol = P) |> qr() |> qr.Q()

set.seed(17042025) 
for (i in 2:its) {
    if (i %% 500 == 0) {print(i)}
    
    Us_rst[, , i] <- rbing.matrix.gibbs(G, H, Us_rst[, , i-1])
    Covs_rst[, , i] <- Us_rst[, , i] %*% A %*% t(Us_rst[, , i])
    
    if (P <= 4) {
        Us_rj_rst[, , i] <- rbing.Op(G, H)
        Covs_rj_rst[, , i] <- Us_rj_rst[, , i] %*% A %*% t(Us_rj_rst[, , i])
    }
}

# look at the last few samples to confirm they are random, not repeated
Us_rst[, , (its-3):its]
# Us_rj_rst[, , (its-3):its]

avgCovs_rst <- apply(Covs_rst[, , (its/2) : its], c(1, 2), mean)
avgUs_rst <- eigen(avgCovs_rst)$vectors

if (P <= 4) {
    avgCovs_rj_rst <- apply(Covs_rj_rst[, , (its/2) : its], c(1, 2), mean)
    avgUs_rj_rst <- eigen(avgCovs_rj_rst)$vectors
}

V

# compare two specific columns
coli <- 7
cbind(avgUs_rst[, coli], V[, coli])

# Frobenius norm of the difference of the "average" covariance and G parameter
sum(diag(crossprod(avgCovs_rst - G)))
# Frobenius norm of the difference of the "average" eigenvector and V parameter
# TODO this could be a bad comparison, since the column signs might be flipped
sum(diag(crossprod(avgUs_rst - V)))

# experiments about (fixed trace) Complex Wishart dist

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

logdet <- function(X) {
    return( log(EigenR::Eigen_det(X)) )
}

rFTCW <- function(Sigma, n, a) {
    W <- rcomplex_wishart(n, Sigma)
    Y <- a * W / Re(sum(diag(W)))
    return(Y)
}

# calculate mean
V <- rFTCW(diag(8), 9, 8)
Ws <- array(NA, c(8, 8, 10000))

for (s in 1:10000) {
    Ws[, , s] <- rFTCW(V, 9, 8)
}

avgW <- apply(Ws, c(1, 2), mean)

cbind(avgW[, 1], V[, 1])

# calculate mode



# determinants of observed Uk and Sigma -----------------------------------

P <- 64

Uex <- runif_stiefel(P, 10, 2)
Wex <- rFTCW(diag(P), P+1, P)

microbenchmark::microbenchmark(
    t(Conj(Uex)) %*% Wex %*% Uex,
    logdet(Wex),
    logdet(t(Conj(Uex)) %*% Wex %*% Uex),
    rFTCW(Wex, 2000, P),
    times = 1000)

microbenchmark::microbenchmark(
    EigenR::Eigen_det(t(Conj(Uex)) %*% Wex %*% Uex) * 
        EigenR::Eigen_det(t(Conj(Uex)) %*% Wex %*% Uex) *
        EigenR::Eigen_det(t(Conj(Uex)) %*% Wex %*% Uex) *
        EigenR::Eigen_det(t(Conj(Uex)) %*% Wex %*% Uex) * 
        EigenR::Eigen_det(t(Conj(Uex)) %*% Wex %*% Uex) *
        EigenR::Eigen_det(t(Conj(Uex)) %*% Wex %*% Uex),
    EigenR::Eigen_det(t(Conj(Uex)) %*% Wex %*% Uex %*%
                          t(Conj(Uex)) %*% Wex %*% Uex %*%
                          t(Conj(Uex)) %*% Wex %*% Uex %*%
                          t(Conj(Uex)) %*% Wex %*% Uex %*%
                          t(Conj(Uex)) %*% Wex %*% Uex %*%
                          t(Conj(Uex)) %*% Wex %*% Uex),
    times = 10000)




microbenchmark::microbenchmark(
    rFTCW(diag(64), 2000, 64),
    times = 1000
)


# if else test ------------------------------------------------------------


if (2 > 1) rrr <- 1 else rrr <- 0
print(rrr)

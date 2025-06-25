# sampler for Sigma parameter of CMACG prior

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

# functions ---------------------------------------------------------------

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

rFTCW <- function(Sigma, n, a) {
    W <- rcomplex_wishart(n, Sigma)
    Y <- a * W / Re(sum(diag(W)))
    return(Y)
}

logdet <- function(X) {
    return( log(EigenR::Eigen_det(X)) )
}

# setup -------------------------------------------------------------------

# need existing Uk values from CMACG(Sigma_0) distribution

P <- 64
d <- 4
K <- 100
S_its <- 20000
# tuning parameter for the proposal distribution
n_Sig <- P+10000
# for P = 8, and K = 10, P + 250 and P + 500 seemed to give acc. rates around 18% and 33%, respectively.
# for K = 100, increased to P+2000 (~.22 acc. rate)

burnin <- 0.5

if (burnin < 1) {
    itsKeep <- ceiling(S_its * burnin) : S_its
} else {
    itsKeep <- (burnin+1) : S_its
}

# experiments about timing ------------------------------------------------

# which takes longer - generating a high d.f. Wishart matrix, or taking many determinants?

# Uex <- runif_stiefel(8, 4, 2)
# Wex <- rFTCW(diag(8), 9, 8)
# 
# microbenchmark::microbenchmark(
#     t(Conj(Uex)) %*% Wex %*% Uex,
#     logdet(Wex),
#     logdet(t(Conj(Uex)) %*% Wex %*% Uex),
#     rFTCW(Wex, 2000, 8),
#     times = 1000)
# 
# microbenchmark::microbenchmark(
#     rFTCW(diag(64), 2000, 64),
#     times = 1000
# )

# the logdet of Uk matrix products is on the order of the time to generate a high-df FTCW, if you need say K = 100 of the Uk logdet

# thus, I think improving both is important

# meanwhile, collecting Sigma logdet terms will not help as much, compared to the number of sample logdets that are needed. But, it should still be done, for simplicity.

# generate data -----------------------------------------------------------

Sigma_0 <- rFTCW(diag(P), P+1, P)

Uks <- array(NA, c(P, d, K))
for (k in 1:K) {
    Uks[, , k] <- rCMACG(P, d, Sigma_0)
}

# initialize array to store samples
Sigma_S <- array(NA, c(P, P, S_its))
Sigma_S[, , 1] <- rFTCW(diag(P), P+1, P)
Sigmas <- Sigma_S[, , 1]
invSigmas <- solve(Sigmas)
accCount <- 0

# det function
#EigenR::Eigen_det
{
    print(Sys.time())
for (s in 2:S_its) {
    
    if(s %% (S_its/10) == 0) print(paste0("s = ", s))
    
    Sigmap <- rFTCW(Sigmas, n_Sig, P)
    invSigmap <- solve(Sigmap)
        
    # sumlog_p
    sumlog_p <- 0
    sumlog_s <- 0
    for (k in 1:K) {
        sumlog_p <- sumlog_p + 
            logdet( t(Conj(Uks[, , k])) %*% invSigmap %*% Uks[, , k] )
        sumlog_s <- sumlog_s + 
            logdet( t(Conj(Uks[, , k])) %*% invSigmas %*% Uks[, , k] )
    }
    
    logdens_num <- -d*K* logdet(Sigmap) - P * sumlog_p - 
        n_Sig * logdet(Sigmap) +
        (n_Sig - P) * logdet(Sigmas) - 
        P*n_Sig * log( Re(sum(diag( invSigmap %*% Sigmas))) ) 
        
    logdens_den <- -d*K* logdet(Sigmas) - P * sumlog_s - 
        n_Sig * logdet(Sigmas) +
        (n_Sig - P) * logdet(Sigmap) - 
        P*n_Sig * log( Re(sum(diag( invSigmas %*% Sigmap ))) ) 
    
    logr <- Re(logdens_num - logdens_den)
    
    if (log(runif(1)) <= logr) {
        Sigmas <- Sigmap
        invSigmas <- invSigmap
        accCount <- accCount + 1
    }
    
    Sigma_S[, , s] <- Sigmas
}
    print(Sys.time())
}

accCount/S_its

avgSigma <- apply(Sigma_S[, , itsKeep], c(1, 2), mean)

cbind(avgSigma[, 1], Sigma_0[, 1])

norm(avgSigma - Sigma_0, "F")

# trace plot of distances to mean -----------------------------------------

d_to_avgSigma <- rep(NA, S_its)
for (s in 1:S_its) {
    d_to_avgSigma[s] <- norm(avgSigma - Sigma_S[, , s], "F")
}

plot(d_to_avgSigma, type = "l", ylab = "Frob. dist.",
     main = "Frobenius distance of samples to mean Sigma")

quantile(d_to_avgSigma[itsKeep], c(.025, .5, .975))

# trace plot of diagonals -------------------------------------------------

diagi <- 1

main_expr_S <- expression(paste(Sigma[11], " trace plot"))

ylab_S <- bquote(Sigma[.(diagi) * .(diagi)])

plot(Re(Sigma_S[diagi, diagi, ]), type = "l",
     ylab = ylab_S, 
     main = main_expr_S)

quantile(Re(Sigma_S[diagi, diagi, itsKeep]), c(.025, .5, .975))
mean(Re(Sigma_S[diagi, diagi, itsKeep]))
Re(Sigma_0[diagi, diagi])

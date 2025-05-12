# testing all the 0F0 approximations


# functions ---------------------------------------------------------------

### For Nasuda approximation
Bmvgamma <- function(C, m, Betaf) {
    term1 <- pi^(.25*m*(m-1)*Betaf)
    term2 <- prod(gamma(C - ((1:m) - 1)/2*Betaf))
    return(term1*term2)
}

omegamBeta <- function(m, Betaf) {
    term1 <- m*log(2)
    term2 <- (m^2 * Betaf / 2) * log(pi)
    term3 <- log(Bmvgamma(m*Betaf/2, m, Betaf))
    
    return(exp(term1 + term2 - term3))
}

Omegams <- function(ms, Betaf) {
    term2 <- 0
    for(j in 1:length(ms)) {
        term2 <- term2 + log(omegamBeta(ms[j], Betaf))
    }
    
    logresult <- term2 - log(omegamBeta(sum(ms), Betaf))
    
    return(exp(logresult))
}

### functions to simulate from complex Stiefel manifold
rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

rcstiefel <- function(P, d) {
    X1 <- matrix(rscnorm(P*d), ncol = d)
    X1qr <- qr(X1)
    QX1 <- qr.Q(X1qr)
    dRX1 <- diag(qr.R(X1qr))
    
    diagD <- sign(Re(dRX1))
    
    Qfinal <- t(t(QX1) * diagD)
    return(Qfinal)
}

### properly normalized Schur polynomial
normSchur <- function(x, lambda) {
    resultQ <- (jack::JackPol(length(x), lambda, 1, which = "C") |> 
                qspray::evalQspray(x)) 
    return(gmp::asNumeric(resultQ))
}

# real, one matrix parameter ----------------------------------------------

P <- 4

# eigenvalues as character and numeric vectors
cas <- c("2", "5/3", "4/3", "1")
# cas <- c("2", "13/7", "12/7", "11/7", "10/7", "9/7", "8/7", "1")
as <- sapply(cas, function(x) eval(parse(text = x))) |> unname()

cbs <- c("2", "7/4", "5/4", "1")
# cbs <- c("2", "15/8", "14/8", "13/8", "12/8", "11/8", "10/8", "1")
bs <- sapply(cbs, function(x) eval(parse(text = x))) |> unname()

###
# truth: 0F0(A) = etr(A)
###
true_real_onepar <- exp(sum(as))

###
# truncated sum, fast calculation via package
###
# alpha = 2 is the real case
m <- 15 # how many terms to calculate
real_onepar_fasttrunc <- HypergeoMat::hypergeomPFQ(m, NULL, NULL, as, alpha = 2)

###
# truncated sum, calculated directly with Jack (zonal) polynomials
###

hgfsum <- 1
for (k in 1:m) {
    # used in each term
    logkfac <- lfactorial(k)
    # the partitions of k
    partsk <- partitions::parts(k)
    # will hold the sum terms
    partsums <- vector("list", ncol(partsk))
    # print status
    print(paste0("k = ", k, "; num. parts = ", ncol(partsk)))
    
    for (j in 1:ncol(partsk)) {
        kappa <- partsk[, j]
        
        # Zonal, real
        logCA <- log(jack::Zonal(as, kappa))
        
        sumadd <- exp(logCA - logkfac)
        
        partsums[j] <- ifelse(is.nan(sumadd), 0, sumadd)
    }
    
    hgfsum <- hgfsum + Reduce("+", partsums)
}

real_onepar_trunc <- hgfsum

rbind(true_real_onepar, real_onepar_fasttrunc, real_onepar_trunc)

# complex, one matrix parameter -------------------------------------------

###
# truth: 0F0(A) = etr(A)
###
true_comp_onepar <- exp(sum(as))

###
# truncated sum, fast calculation via package
###
# alpha = 2 is the real case
m <- 15 # how many terms to calculate
comp_onepar_fasttrunc <- HypergeoMat::hypergeomPFQ(m, NULL, NULL, as, alpha = 1)

###
# truncated sum, calculated directly with Jack (zonal) polynomials
###

hgfsum <- 1
for (k in 1:m) {
    # used in each term
    logkfac <- lfactorial(k)
    # the partitions of k
    partsk <- partitions::parts(k)
    # will hold the sum terms
    partsums <- vector("list", ncol(partsk))
    # print status
    print(paste0("k = ", k, "; num. parts = ", ncol(partsk)))
    
    for (j in 1:ncol(partsk)) {
        kappa <- partsk[, j]
        
        # Schur, complex
        # logCA <- log(jack::Schur(as, kappa))
        logCA <- normSchur(cas, kappa) |> log()
        
        sumadd <- exp(logCA - logkfac)
        
        partsums[j] <- ifelse(is.nan(sumadd), 0, sumadd)
    }
    
    hgfsum <- hgfsum + Reduce("+", partsums)
}

comp_onepar_trunc <- hgfsum

rbind(true_comp_onepar, comp_onepar_fasttrunc, comp_onepar_trunc)
# real, two matrix parameters ---------------------------------------------

m <- 10

library(foreach)
library(doParallel)

parallel::detectCores()

cluster <- makeCluster(12)
registerDoParallel(cluster)

Is <- rep(1, P)
cIs <- rep("1", P)

###
# Monte Carlo approximation
###

M <- 1e4
# intgd <- rep(NA, M)
intgd <- vector("list", M)

# for (i in 1:M) {
intgd <- foreach(i = 1:M) %dopar% {
    # U <- rcstiefel(P, P)
    U <- rstiefel::rustiefel(P, P)
    # intgd[i] <- (A %*% U %*% B %*% t(Conj(U))) |> 
    #     diag() |> sum() |> Re() |> exp()
    intgd[i] <- exp(Re(t(as) %*% (Conj(U) * U) %*% bs))
}

real_twopar_MC <- c(Reduce("+", intgd)/M)

# mean(intgd)
# summary(intgd)

###
# truncated sum
###

hgfsum <- 1
for (k in 1:15) {
    # used in each term
    logkfac <- lfactorial(k)
    # the partitions of k
    partsk <- partitions::parts(k)
    # will hold the sum terms
    partsums <- vector("list", ncol(partsk))
    # print status
    print(paste0("k = ", k, "; num. parts = ", ncol(partsk)))
    
    # for (j in 1:ncol(partsk)) {
    partsums <- foreach(j = 1:ncol(partsk)) %dopar% {
        kappa <- partsk[, j]
        
        # Zonal, real
        logCA <- log(jack::Zonal(as, kappa))
        logCB <- log(jack::Zonal(bs, kappa))
        logCI <- log(jack::Zonal(Is, kappa))
        
        sumadd <- exp(logCA + logCB - logCI - logkfac)
        
        partsums[j] <- ifelse(is.nan(sumadd), 0, sumadd)
    }
    
    hgfsum <- hgfsum + Reduce("+", partsums)
}

real_twopar_trunc <- hgfsum

stopCluster(cl = cluster)

###
# asymptotic approximation
###

betaf <- 1
ss <- P*(P - 1)/2

ijs <- t(combn(1:P, 2))
cij <- (as[ijs[, 1]] - as[ijs[, 2]]) * (bs[ijs[, 1]] - bs[ijs[, 2]])

# 2^P *
#     exp(sum(as*bs)) *
#     prod(sqrt(pi/cij))

# Nasuda approximation
# JAB
JAB <- prod(2/betaf*cij)

real_twopar_asymp <- (2/betaf*pi)^(betaf/2*ss) *
    Omegams(rep(1, P), betaf) *
    JAB^(-betaf/2) *
    exp(sum(as*bs))

rbind(real_twopar_MC, 
      real_twopar_trunc,
      real_twopar_asymp)

# complex, two matrix parameters ------------------------------------------

parallel::detectCores()

cluster <- makeCluster(12)
registerDoParallel(cluster)

###
# Monte Carlo approximation
###

M <- 1e4
# intgd <- rep(NA, M)
intgd <- vector("list", M)

# for (i in 1:M) {
intgd <- foreach(i = 1:M) %dopar% {
    U <- rcstiefel(P, P)
    intgd[i] <- exp(Re(t(as) %*% (Conj(U) * U) %*% bs))
}

comp_twopar_MC <- c(Reduce("+", intgd)/M)

###
# truncated sum
###

hgfsum <- 1
for (k in 1:m) {
    # used in each term
    logkfac <- lfactorial(k)
    # the partitions of k
    partsk <- partitions::parts(k)
    # will hold the sum terms
    partsums <- vector("list", ncol(partsk))
    # print status
    print(paste0("k = ", k, "; num. parts = ", ncol(partsk),
                 "; hgfsum = ", round(hgfsum, 2)))
    
    # for (j in 1:ncol(partsk)) {
    partsums <- foreach(j = 1:ncol(partsk)) %dopar% {
        kappa <- partsk[, j]

        # complex: properly Normalized Schur
        logCA <- log(normSchur(cas, kappa))
        logCB <- log(normSchur(cbs, kappa))
        logCI <- log(normSchur(cIs, kappa))
        
        sumadd <- exp(logCA + logCB - logCI - logkfac)
        
        partsums[j] <- ifelse(is.nan(sumadd), 0, sumadd)
    }
    
    hgfsum <- hgfsum + Reduce("+", partsums)
}

comp_twopar_trunc <- hgfsum

stopCluster(cl = cluster)

###
# asymptotic approximation
###

betaf <- 2
ss <- P*(P - 1)/2

ijs <- t(combn(1:P, 2))
cij <- (as[ijs[, 1]] - as[ijs[, 2]]) * (bs[ijs[, 1]] - bs[ijs[, 2]])

# 2^P *
#     exp(sum(as*bs)) *
#     prod(sqrt(pi/cij))

# Nasuda approximation
# JAB
JAB <- prod(2/betaf*cij)

comp_twopar_asymp <- (2/betaf*pi)^(betaf/2*ss) *
    Omegams(rep(1, P), betaf) *
    JAB^(-betaf/2) *
    exp(sum(as*bs))

rbind(comp_twopar_MC, comp_twopar_trunc, comp_twopar_asymp) |> 
    format(scientific = FALSE)

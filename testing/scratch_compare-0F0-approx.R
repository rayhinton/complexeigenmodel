# testing all the 0F0 approximations

logScaleStats <- function(log_values, confidence_level = 0.95) {
    # Helper functions
    logsumexp <- function(x) {
        max_x <- max(x)
        return(max_x + log(sum(exp(x - max_x))))
    }
    
    logDiffExp <- function(a, b) {
        return(a + log(1 - exp(b - a)))
    }
    
    logAbsDiffExp <- function(a, b) {
        return(max(a, b) + log(1 - exp(min(a, b) - max(a, b))))
    }
    
    n <- length(log_values)
    
    # Log of the mean
    log_mean <- logsumexp(log_values) - log(n)
    
    # Compute squared differences from mean in log space
    logSqDiffs <- sapply(log_values, function(x) {
        2 * logAbsDiffExp(x, log_mean)
    })
    
    # Log of sample standard deviation
    log_sample_var <- logsumexp(logSqDiffs) - log(n - 1)
    log_sample_sd <- 0.5 * log_sample_var
    
    # Log of standard error
    log_std_err <- 0.5 * (logsumexp(logSqDiffs) - log(n - 1) - log(n))
    
    # Calculate z-score for confidence interval
    alpha <- 1 - confidence_level
    z_score <- qnorm(1 - alpha/2)
    
    # 95% confidence interval bounds in log space
    # Upper bound: log(mean + z * std_err)
    log_upper <- logsumexp(c(log_mean, log(z_score) + log_std_err))
    
    # Lower bound: log(mean - z * std_err)
    # Only compute if mean > z * std_err (otherwise CI would include negative values)
    log_margin <- log(z_score) + log_std_err
    
    if (log_mean > log_margin) {
        log_lower <- logDiffExp(log_mean, log_margin)
    } else {
        log_lower <- -Inf  # CI extends to 0, so log(0) = -Inf
        warning("Lower confidence bound would be negative or zero. Setting to -Inf.")
    }
    
    # Return results
    return(list(
        log_mean = log_mean,
        log_sd = log_sample_sd,
        log_std_err = log_std_err,
        confidence_level = confidence_level,
        log_ci_lower = log_lower,
        log_ci_upper = log_upper,
        # For convenience, also return the original scale values
        mean = exp(log_mean),
        sd = exp(log_sample_sd),
        ci_lower = exp(log_lower),
        ci_upper = exp(log_upper)
    ))
}

# Example usage:
# Assuming you have a vector of log-scale values called 'log_data'
# results <- logScaleStats(log_data)
# print(results)

# FCD by MC integration functions ----------------------------------------

### functions to simulate from complex Stiefel manifold
# return n independent samples from the standard complex normal distribution
rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

# draw matrix distributed uniformly on V_d^P(F_betaf), where 
# - betaf = 1 --> Real
# - betaf = 2 --> Complex
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

# calculate LogSumExp
logsumexp <- function(x) {
    xstar <- max(x)
    lse <- xstar + log( sum(exp( x - xstar )) )
    return(lse)
}

# Approximate 0F0 by MC integration
hgf0F0MC <- function(x, y, betaf, MCits = 1e4, logscale = FALSE) {
    P <- length(x)
    d <- length(y)
    
    # compute the log integrand
    logintgd <- rep(NA, MCits)
    for (i in 1:MCits) {
        U <- runif_stiefel(P, d, betaf = betaf)
        logintgd[i] <- Re(crossprod(x, (Conj(U) * U) %*% y))
    }
    
    # return result on linear or log scales
    if (logscale == FALSE) {
        return(mean(exp(logintgd)))
    } else {
        # use LogSumExp trick for log scale approximation
        return(logsumexp(logintgd) - log(MCits))
    }
}

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
# return n independent samples from the standard complex normal distribution
# rscnorm <- function(n) {
#     return(rnorm(n, 0, 1/sqrt(2)) + 
#                1i * rnorm(n, 0, 1/sqrt(2)))
# }
# 
# # draw matrix distributed uniformly on V_d^P(C)
# rcstiefel <- function(P, d) {
#     # matrix of independent standard complex Normals
#     X1 <- matrix(rscnorm(P*d), ncol = d)
#     # take QR decomposition - not guaranteed to be unique due to numerical methods
#     X1qr <- qr(X1)
#     QX1 <- qr.Q(X1qr)
#     # extract sign of the diagonal elements of R
#     D <- sign(Re(diag(qr.R(X1qr))))
#     # transform by Q = QS, where S = diag(D), which ensures uniqueness of Q
#     Qfinal <- t(t(QX1) * D) # faster than matrix multiplication
#     return(Qfinal)
# }

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

cluster <- makeCluster(10)
registerDoParallel(cluster)

Is <- rep(1, P)
cIs <- rep("1", P)

###
# Monte Carlo approximation
###

M <- 1e5
# intgd <- rep(NA, M)
intgd <- vector("list", M)

ab_mult <- 20

### parallel over each integrand
# for (i in 1:M) {
intgd <- foreach(i = 1:M) %dopar% {
    # U <- rcstiefel(P, P)
    U <- rstiefel::rustiefel(P, P)
    # intgd[i] <- (A %*% U %*% B %*% t(Conj(U))) |> 
    #     diag() |> sum() |> Re() |> exp()
    intgd[i] <- Re(t(ab_mult*as) %*% (Conj(U) * U) %*% (ab_mult*bs))
}

intgd <- unlist(intgd)
# (real_twopar_MC <- mean(intgd))
real_twopar_MC <- logsumexp(intgd) - log(length(intgd))
summary(intgd)
mean(intgd)
sd(intgd)

a <- intgd[1]
b <- real_twopar_MC

a - b
log(a - b)
log(a) + log(1 - exp(log(b) - log(a)))

logDiffExp <- function(a, b) {
    return(a + log(1 - exp(b - a)))
}

logAbsDiffExp <- function(a, b) {
    # return( max(a, b) + log1p(1 - exp(min(a, b) - max(a, b))) )
    return( max(a, b) + log1p(-exp(min(a, b) - max(a, b))) )
}

logDiffExp(log(a), log(b))
logDiffExp(log(b), log(a))

logAbsDiffExp(log(a), log(b))
logAbsDiffExp(log(b), log(a))

logSqDiffs <- rep(NA, length(intgd))
for (i in 1:length(logSqDiffs)) {
    logSqDiffs[i] <- 2 * logAbsDiffExp(intgd[i], real_twopar_MC)
}

summary(logSqDiffs)
unique(logSqDiffs)

# sample variance, on the log scale
(logSampSD <- 0.5*(logsumexp(logSqDiffs) - log(length(logSqDiffs) - 1)))
# standard error, on the log scale
(logStdErr <- 0.5*(logsumexp(logSqDiffs) - log(length(logSqDiffs) - 1) - log(length(logSqDiffs))))

sqrt(var(intgd)) |> log()
sqrt(var(intgd) / length(logSqDiffs)) |> log()

# upper bound, linear scale
real_twopar_MC + 2*sqrt(var(intgd) / length(logSqDiffs))
# upper bound, log scale
logsumexp(c(real_twopar_MC, log(2) + logStdErr)) |> exp()
# lower bound, linear scale
real_twopar_MC - 2*sqrt(var(intgd) / length(logSqDiffs))
logDiffExp(log(real_twopar_MC), log(2) + logStdErr) |> exp()

# 95% CI on log scale
c(logsumexp(c((real_twopar_MC), log(2) + logStdErr)), 
  logDiffExp((real_twopar_MC), log(2) + logStdErr))

alphacrit <- 1 - 0.95
z_score <- qnorm(1 - alphacrit/2)

logScaleStats(intgd)

### parallel over batches
batch_count <- 10
batch_means <- vector("list", batch_count)

{
    print(Sys.time())
    batch_means <- foreach(i = 1:batch_count, .inorder = FALSE) %dopar% {
        hgf0F0MC(as, bs, betaf = 2, MCits = M/batch_count, logscale = TRUE)
    }
    print(Sys.time())
}
microbenchmark::microbenchmark(
    {
    batch_means <- foreach(i = 1:batch_count, .inorder = FALSE) %dopar% {
        hgf0F0MC(as, bs, betaf = 2, MCits = M/batch_count, logscale = TRUE)
    }
    
    (logsumexp(unlist(batch_means)) - log(batch_count)) |> exp()
    },
    times = 50
)
### serial
{
    print(Sys.time())
    hgf0F0MC(as, bs, betaf = 2, MCits = M, logscale = FALSE)    
    print(Sys.time())
}
microbenchmark::microbenchmark(
    {
    hgf0F0MC(as, bs, betaf = 2, MCits = M, logscale = FALSE)    
    },
    times = 50
)


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
rbind(comp_twopar_MC, comp_twopar_asymp) |> 
    format(scientific = FALSE)


# ratio, MC vs. Laplace ---------------------------------------------------

as1 <- c(2, runif(P-2, 1, 2) |> sort(decreasing = TRUE), 1)
bs1 <- c(2, runif(P-2, 1, 2) |> sort(decreasing = TRUE), 1)

w1 <- 2
w2 <- 3

cij1 <- cij * (w1^ss)
cij2 <- cij * (w2^ss)

JAB1 <- prod(2/betaf*cij1)
JAB2 <- prod(2/betaf*cij2)

w1_dens_asymp <- (2/betaf*pi)^(betaf/2*ss) *
    Omegams(rep(1, P), betaf) *
    JAB1^(-betaf/2) *
    exp(sum(w1*as*bs))

w2_dens_asymp <- (2/betaf*pi)^(betaf/2*ss) *
    Omegams(rep(1, P), betaf) *
    JAB2^(-betaf/2) *
    exp(sum(w2*as*bs))

w1_dens_asymp/w2_dens_asymp

###
# MC approx. for w1
###
M <- 1e4
# intgd <- rep(NA, M)
intgd <- vector("list", M)

# for (i in 1:M) {
intgd <- foreach(i = 1:M) %dopar% {
    U <- rcstiefel(P, P)
    intgd[i] <- exp(Re(w1 * t(as) %*% (Conj(U) * U) %*% bs))
}

w1_dens_MC <- c(Reduce("+", intgd)/M)

###
# MC approx. for w1
###
M <- 1e4
# intgd <- rep(NA, M)
intgd <- vector("list", M)

Sys.time()
for (i in 1:M) {
# intgd <- foreach(i = 1:M) %dopar% {
    U <- rcstiefel(P, P)
    intgd[i] <- exp(Re(w2 * t(as) %*% (Conj(U) * U) %*% bs))
}
Sys.time()

w2_dens_MC <- c(Reduce("+", intgd)/M)

log(w1_dens_MC/w2_dens_MC)
log(.000115)
log(.0001179)


# experiment with log-scale integrals -------------------------------------

M <- 1e4
# intgd <- rep(NA, M)
intgd <- vector("list", M)

for (i in 1:M) {
# intgd <- foreach(i = 1:M) %dopar% {
    U <- rcstiefel(P, P)
    intgd[i] <- exp(Re(t(as) %*% (Conj(U) * U) %*% bs))
}

c(Reduce("+", intgd)/M)

exp(matrixStats::logSumExp(log(unlist(intgd))) - log(M))

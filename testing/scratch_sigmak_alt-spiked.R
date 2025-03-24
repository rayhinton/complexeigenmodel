# sigmak^2 from alternative spiked covariance model

library(cmvnorm)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complex-computation-issue-checking.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

n <- 1000
P <- 12
d <- 4
# V matrix, PxP
set.seed(18032025)
U <- runitary(P, d)
lambdas <- seq(500, 3, length.out = d)
# lambdas <- seq(50, 30, length.out = d)

firstmat <- U %*% diag(lambdas) %*% t(Conj(U))
eigen(firstmat)$values

# check U is semi-unitary
t(Conj(U)) %*% U

sigma_02 <- 2

Sigma_0 <- U %*% diag(lambdas) %*% t(Conj(U)) + diag(sigma_02, nrow = P)

U_r <- complex_to_real(U)
L_r <- complex_to_real(diag(lambdas))
sigIp_r <- complex_to_real(sigma_02 * diag(P))
# Sigma_0_r <- U_r %*% L_r %*% t(U_r) + sigIp_r
Sigma_0_c <- real_to_complex(U_r %*% L_r %*% t(U_r) + sigIp_r)

eigen(Sigma_0)$values

# Y <- rcwis(n, Sigma_0_c)
Y <- rcomplex_wishart(n, P, Sigma = .5*Sigma_0_c)
Y_r <- complex_to_real(Y)

trY <- Re(sum(diag(Y)))

# compute log proportional density ----------------------------------------

# helper values based on given parameters and data
# eigenY <- eigen(Y)
# a <- eigenY$values
# V <- eigenY$vectors
# Z <- t(Conj(V)) %*% U
# aZZ <- t(a) %*% ( Conj(Z) * Z )

# over a range of x = sigma2 values
nxs <- 100001
xs <- seq(.00001, 1000, length.out = nxs)
logpropden <- rep(NA, nxs)

for (i in 1:nxs) {
    x <- xs[i]
    phi_sigma <- lambdas / (lambdas + x)
    Om <- diag(phi_sigma)
    Om_r <- complex_to_real(Om)
    # TODO change the first n*P to either n*P or (n*P + 1), depending on prior
    # logpropden[i] <- Re(-(n*P + 2)*log(x) - n * sum(log((lambdas + x)/x)) -
    #     (1/x) * (sum(diag(Y)) - aZZ %*% phi_sigma)) - 2/x
    # logpropden[i] <- Re(-(n*P)*log(x) - n * sum(log((lambdas + x)/x)) -
                            # (1/x) * (sum(diag(Y)) - aZZ %*% phi_sigma)) - 100*x
    
    prod2_r <- U_r %*% Om_r %*% t(U_r) %*% Y_r
    prod2_c <- real_to_complex(prod2_r)
    tr2 <- Re(sum(diag(prod2_c)))
    
    # tr2 <- sum(diag(U %*% Om %*% t(Conj(U)) %*% Y))
    
    logpropden[i] <- -(n*P + 2)*log(x) - (n)*sum(log((lambdas + x)/x)) -
                            (1/x)*(trY - tr2) - 2/x
}

summary(logpropden)
xs[which.max(logpropden)]

sample_sigma2 <- sample_gumbel(xs, 1000, logpropden)
summary(sample_sigma2)
quantile(sample_sigma2, c(0.025, 0.975))

plot(density(sample_sigma2))
hist(sample_sigma2)


# try integrating the proportional density --------------------------------

one_dens <- function(x) {
    phi_sigma <- lambdas / (lambdas + x)
    # TODO change the first n*P to either n*P or (n*P + 1), depending on prior
    return(
        Re(-(n*P)*log(x) - n * sum(log((lambdas + x)/x)) -
               (1/x) * (sum(diag(Y)) - aZZ %*% phi_sigma))
    )
}

# plot Inverse Gamma densities
curve(dgamma(1/x, shape = 1, rate = 5)/x^2, from = 0, to = 300)
curve(dexp(x, 100), from = 0, to = 300)

eval_dens <- function(x) {
    return(
        vapply(x, one_dens, 0)
    )
}

integrate(function(x) return(exp(eval_dens(x))),
          lower = 0.0001, upper = Inf)
integrate(function(x) return(5*dnorm(x)),
          lower = 0,
          upper = Inf)

# vapply(xs[1:10], eval_dens, 0)
# eval_dens(xs[1:10])
# logpropden[1:10]


# expected value of Complex Wishart traces? -------------------------------

its <- 1000

trace_Is <- rep(NA, its)
col1_Is <- matrix(NA, P, its)

set.seed(20032025)
for (i in 1:its) {
    A <- rcwis(n, diag(P))
    # A <- rcomplex_wishart(n, P, Sigma = .5*diag(P))
    trace_Is[i] <- sum(diag(A))
    col1_Is[, i] <- A[, 1]
}

# true mean:
n * P
# observed mean
mean(trace_Is)
summary(Re(trace_Is))
# true mean first col:
n * diag(P)[, 1]
rowMeans(col1_Is)


trace_Sigma0s <- rep(NA, its)
col1_Sigma0s <- matrix(NA, P, its)

set.seed(20032025)
for (i in 1:its) {
    A <- rcwis(n, Sigma_0)
    # A <- rcomplex_wishart(n, P, Sigma = .5*Sigma_0)
    trace_Sigma0s[i] <- sum(diag(A))
    col1_Sigma0s[, i] <- A[, 1]
}
# true mean:
n * sum(diag(Sigma_0))
mean(trace_Sigma0s)
summary(Re(trace_Sigma0s))

# true mean 1st col:
n * Sigma_0[, 1]
rowMeans(col1_Sigma0s)

diag(Sigma_0) <- Re(diag(Sigma_0))

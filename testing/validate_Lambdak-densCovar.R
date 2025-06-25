# validate Lambda_k sampler

library(truncdist)

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
# source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_sigmak2.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Lambdak.R")

# depends on n_k (sample size), U_k, P_k (data), and sigma_k2

# define n, sample size

# to test for a single value: 
# - need to define n, sample size,
# - need to generate
#   - U, Pxd matrix of orthonormal eigenvectors
#   - sigma2, a scale term
# - generate a true Lambda_0
# - generate a data observation Y from CW(n, sigma2 (U L U^H + I_P))

# initial parameters
K <- 3
n_k <- 100 + (1:K)*10
P <- 8
d <- 4
maxL <- 100
minL <- 2

# number of sampling iterations
S <- 1000
# which samples to keep for posterior analysis
gibbskeep <- round(S/2):S

# true Lambda values
Lambda_k_0 <- array(NA, c(d, K))

# generate parameters
U_ks <- array(NA, c(P, d, K))
Lambda_ks <- array(NA, c(d, K))
sigma_k2s <- array(NA, c(K))

data_list <- list()

set.seed(26032025)
for (k in 1:K) {
    # generate random true parameters
    U <- runitary(P, d)
    sigma2 <- 0.5 + (k-1)*5
    Lambda_0 <- sort(runif(d, minL, maxL), decreasing = TRUE)
    
    # data generating parameter
    Sigma_0 <- sigma2 * ( U %*% diag(Lambda_0) %*% t(Conj(U)) + diag(P) )
    
    # generate data value
    data_list[[k]] <- rcomplex_wishart(n_k[k], Sigma_0)
    
    # store values in arrays to be put in param_list
    U_ks[, , k] <- U
    sigma_k2s[k] <- sigma2
    Lambda_k_0[, k] <- Lambda_0
}

# create param_list
param_list <- list(
    P = P,
    d = d,
    K = K,
    n_k = n_k,
    U_ks = U_ks,
    # placeholder values, 1, for Lambda
    Lambda_ks = array(1, c(d, K)),
    sigma_k2s = sigma_k2s
)

# generate a sample, as an example
Lambdak_gibbs_densCovar(data_list, param_list)

# prepare to do sampling over different data observations
Ytests <- 200
postmeans <- array(NA, c(d, K, Ytests))
inpostCIs <- array(NA, c(d, K, Ytests))

set.seed(26032025)
for (yi in 1:Ytests) {
    # generate new data matrices
    for (k in 1:K) {
        Sigma_0 <- sigma_k2s[k] * ( U_ks[, , k] %*% diag(Lambda_k_0[, k]) %*% 
                                  t(Conj(U_ks[, , k])) + diag(P) )
        data_list[[k]] <- rcomplex_wishart(n_k[[k]], Sigma_0)
    }

    # holds all sampled Lambda values
    Lambda_s <- array(NA, c(d, K, S))
    
    # sample the posterior with these observations
    for (s in 1:S) {
        Lambda_s[, , s] <- Lambdak_gibbs_densCovar(data_list, param_list)
        param_list$Lambda_ks <- Lambda_s[, , s]
    }

    # calculate some posterior summaries for this Y value
    for (k in 1:K) {
        # calculate posterior means
        postmeans[, k, yi] <- apply(Lambda_s[, k, gibbskeep, drop = FALSE],
                                    c(1), mean)
        
        # get a 95% credible interval for each Lambda
        CI95 <- apply(Lambda_s[, k, gibbskeep, drop = FALSE], 
                       c(1), quantile, probs = c(0.025, 0.975))
        # do the CIs contain Lambda_0?
        inpostCIs[, k, yi] <- 
            (Lambda_k_0[, k] > CI95[1, ]) & (Lambda_k_0[, k] < CI95[2, ])         
    }
}

# look at summaries over all the sampled posteriors
k <- 3
apply(inpostCIs[, k, ], c(1), mean)
apply(postmeans[, k, ], c(1), mean)
Lambda_k_0

# plot observed densities -------------------------------------------------
# which group and term to plot
k <- 1
j <- 3

# plot observed density on Lambda scale
plot(density(Lambda_s[j, k, gibbskeep]))
abline(v = Lambda_0[j, j])

### plot on the xi scale, i.e. in [0, 1]
# calculate helper terms and density values
ajk <- Re(t(Conj(U[, j])) %*% Y %*% U[, j]) / sigma2
xs <- seq(0, 1, length.out = 1001)
truncgammas <- dtrunc(xs, "gamma", a = 0, b = 1, shape = n, rate = ajk)

# plot observed and true density
plot(density(1/(1+Lambda_s[j, k, gibbskeep])))
abline(v = 1/(1+Lambda_0[j, j]))
lines(xs, truncgammas, col = "red")

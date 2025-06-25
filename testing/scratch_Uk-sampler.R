# test Uk FCD sampler

library(foreach)

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

# the columns, conditional on others, have certain vector Bingham distribution

# Uk | A, B, V, Lambdak, Pk
# where Pk | Uk, Lambdak, sigmak2 ~ Complex Wishart

# - need true Uk0 matrices, not necessarily from any particular distribution
# - do need Pk generated from correct CW distribution

# setup -------------------------------------------------------------------

numchains <- 5

K <- 1
P <- 8
d <- 4
S_its <- 15000

burnin <- 5000
thinBy <- 5

nks <- rep(1000, K)

# initialize parameters and generate data ---------------------------------

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

set.seed(22042025)
# Yk <- rcomplex_wishart(nk, Gammak0)
Yks <- array(NA, c(P, P, K))
data_list <- list()

Uks <- array(NA, c(P, d, K, S_its))

for (k in 1:K) {
    Lambdaks[, , k] <- rgamma(d, 1, 1) |> sort(decreasing = TRUE) |> diag()
    sigmak2s[k] <- rgamma(1, 1, 1)
    Uk0s[, , k] <- runitary(P, d)
    
    Gammak0 <- sigmak2s[k] * 
        ( Uk0s[, , k] %*% Lambdaks[, , k] %*% t(Conj(Uk0s[, , k])) + diag(P) )
    
    Yks[, , k] <- rcomplex_wishart(nks[k], Gammak0)
    data_list[[k]] <- Yks[, , k]
    
    # random initialization for each
    # Uks[, , k, 1] <- runitary(P, d)
    # ideal initialization at the exact true Uk0
    # Uks[, , k, 1] <- Uk0s[, , k]
    Uks[, , k, 1] <- eigen(Yks[, , k])$vectors[, 1:d]

}

param_list <- list(
    P = P,
    d = d,
    K = K,
    n_k = nks,
    # U_ks = Uks[, , , 1],
    Lambda_ks = Lambdaks,
    sigma_k2s = sigmak2s,
    Vs = V,
    As = A,
    Bs = B
)

if (K == 1) {
    param_list$U_ks <- array(Uks[, , , 1], c(P, d, 1))
} else {
    param_list$U_ks <- Uks[, , , 1]
}

# make different initializations ------------------------------------------

# for testing purposes - assuming K = 1

# estimate Ukest from data, using first d eigenvectors
Ukest <- eigen(Yks[, , 1])$vectors[, 1:d]

# attempts to generate "overdispersed" starting points
initialUks <- list(
    NULL, 
    Uk_MH_Cayley(Ukest, .1, FALSE, NULL),
    Uk_MH_Cayley(Ukest, .2, FALSE, NULL),
    Ukest, MASS::Null(Ukest)[, 1:d]
)

# run sampler -------------------------------------------------------------

cluster <- parallel::makeCluster(numchains)
doParallel::registerDoParallel(cluster)

manysamples <- list()

set.seed(13042025)
{
    print(Sys.time())
    manysamples <- foreach(j = 1:numchains) %dopar% {
        
        Uks <- array(NA, c(P, d, K, S_its))
        
        if(is.null(initialUks[[j]])) {
            Uks[, , 1, 1] <- runif_stiefel(P, d, 2)
        } else {
            Uks[, , 1, 1] <- initialUks[[j]]
        }
        param_list$U_ks <- array(Uks[, , , 1], c(P, d, 1))
        
        for (s in 2:S_its) {
            if (s %% 500 == 0) print(paste0("s = ", s))
            
            Uks[, , , s] <- Uk_gibbs_densCovar(data_list, param_list)
            if (K == 1) {
                param_list$U_ks <- array(Uks[, , , s], c(P, d, 1))
            } else {
                param_list$U_ks <- Uks[, , , s]
            }
    
        }
        
        whichk <- 1
        
        manysamples[[j]] <- list(Uk_S = Uks[, , whichk, ], accCount = 1)
    }
    print(Sys.time())
}

parallel::stopCluster(cluster)

# use summarize function --------------------------------------------------

# make a new param_list that has required info for summary functions

# param_list <- list(
#     Sigmas = Sigma0,
#     invSigmas = invSigma0,
#     Lambda_ks = diag(Lambdak0),
#     Omega_ks = diag(Omegak0),
#     sigma_k2s = sigmak02,
#     P = P,
#     d = d,
#     n_k = n_k,
#     Uk0 = Uk0
# )

new_param_list <- list(Lambda_ks = diag(Lambdaks[, , 1]),
                       sigma_k2s = sigmak2s,
                       P = P,
                       d = d,
                       n_k = nks,
                       Uk0 = Uk0s[, , 1])

summarystuff <- summarize_Uk(manysamples, new_param_list, burnin = burnin)

summarystuff$avgs_df

dev.new()
pdf("Uk_CGB_sampler_trace.pdf", width = 9, height = 7.5)
dist_vals <- avg_tracePlots(manysamples, summarystuff, new_param_list, 
                            burnin = burnin, tracePlotEvery = thinBy)
dev.off()

bigoutmcmcs <- make_mcmclist(dist_vals, S_its, burnin, thinBy)

dist_vals$chain_quantiles
with(dist_vals,
     chain_quantiles[rownames(chain_quantiles) == "d_to_avgUk", ])

coda::gelman.diag(bigoutmcmcs)

coda::effectiveSize(bigoutmcmcs) |> t()
lapply(bigoutmcmcs, coda::effectiveSize)

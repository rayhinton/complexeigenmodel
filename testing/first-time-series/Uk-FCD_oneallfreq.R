# functions ---------------------------------------------------------------

source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")

# setup -------------------------------------------------------------------

P <- 4
d <- 2
K <- 2

tau_Uk <- rep(.005, K)

# initialize arrays -------------------------------------------------------

U_ks_all <- array(NA, c(P, d, K, gibbsIts))
U_ks_all[, , , 1] <- U_k0

Lambdak_w_s <- array(NA, c(d, K, num_freqs, gibbsIts))
Lambdak_w_s[, , , 1] <- Lambdak0_w

taujk2_s <- array(NA, c(d, K, gibbsIts))
taujk2_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)

# sampling ----------------------------------------------------------------

s <- 2

# before sampling
newU_ks <- array(NA, c(P, d, K))

# during sample, 1:K, 1:Nw

for (k in 1:K) {
    Us <- U_ks[, , k]
    # propose U by small rotation of Us
    Up <- Uk_MH_Cayley(Us, tau_Uk[k], doCayleyZeros, CayleyZeroProb)
    # let Y_kl be the scaled data matrix (i.e. SDM estimate at frequency index l)

    invSigmas <- solve(Sigmas)

    tracep <- 0
    traces <- 0
    
    # calculate traces and other terms, for the MH acceptance ratio
    for (l in 1:Nw) {
        Skl <- data_list_w[[l]][[k]]
        lambdakl <- result_Lambdas[, k, l]
        Omegakl <- diag(lambdakl / (1 + lambdakl))
        
        termp <- Re(sum(diag( Up %*% Omegakl %*% t(Conj(Up)) %*% Skl )))
        terms <- Re(sum(diag( Us %*% Omegakl %*% t(Conj(Us)) %*% Skl )))
        
        tracep <- tracep + termp
        traces <- traces + terms
    } # end of summing over frequencies

    logdetp <- log(Re(EigenR::Eigen_det( t(Conj(Up)) %*% invSigmas %*% Up )))
    logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*% invSigmas %*% Us )))

    # calculate acceptance ratio
    logr <- tracep/sigma_k2 - P*logdetp - traces/sigma_k2 + P*logdets

    # accept or reject
    if (log(runif(1)) <= logr) {
        newU_ks[, , k] <- Up
    } else {
        newU_ks[, , k] <- Us
        accCount[k] <- FALSE
    }

} # end of sampling over 1:K

U_ks_all[, , , s] <- newU_ks[, , k]

# functions ---------------------------------------------------------------

source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Uk_CMACG.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcomplex_wishart.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/matrix_distances.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/first-time-series/geoSS.R")

closest_semiorth <- function(U) {
 
    A <- Re(U)
    svd_result <- svd(A)
    V <- svd_result$u %*% t(svd_result$v)
    
    return(V)
}

sdf_ar2 <- function(omega, phis, sigma2 = 1) {
    phi1 <- phis[1]
    phi2 <- phis[2]
    
    denom <- 1 + phi1^2 + phi2^2 + 2*phi1*(phi2 - 1) * cos(2*pi*omega) - 
        2*phi2*cos(4*pi*omega)
    
    return(sigma2 / denom)
}

# unnormalized PDF
unnorm_logPDF <- function(x, tau2, ajk, mu, N, logscale = FALSE) {
    logf1 <- -N * log(1 + x)
    logf2 <- -ajk / (1+x) - (.5/tau2) * (x - mu)^2
    
    result <- ifelse(x <= 0, -Inf, logf1 + logf2)
    
    if (logscale) {
        return(result)
    } else {
        return(exp(result))
    }
}

logd_Uk <- function(Us, lw, invSigmas, sigma_k2, result_Lambdas, data_list_w, 
                    num_freqs, k) {
    
    # invSigmas <- solve(Sigmas)
    
    # calculate traces and other terms, for the MH acceptance ratio
    # for (l in 1:num_freqs) {
    Skl <- data_list_w[[lw]][[k]]
    lambdakl <- result_Lambdas[, k, lw]
    Omegakl <- diag(lambdakl / (1 + lambdakl))
    
    # trace <- Re(sum(diag( Us %*% Omegakl %*% t(Conj(Us)) %*% Skl )))
    traces <- Re(sum(Conj(Us) * (Skl %*% Us %*% Omegakl)))
    
    # } # end of summing over frequencies
    
    logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*% invSigmas %*% Us )))
    
    logdens <- traces/sigma_k2 - P*logdets
    
    return(logdens)
}

# setup -------------------------------------------------------------------

# dataseed <- 21092025
# dataseed <- 22092025
# dataseed <- 10102025
dataseed <- 17102025
parseed <- 314

P <- 4
d <- 2
K <- 2
Tt <- 1024 # length of time series
LL <- round(sqrt(Tt))

num_freqs <- Tt/2 - 1

# Geodesic slice sampling parameters
w <- 10
m <- 2

# tau_Uk <- rep(.1, K)
tau_Ukl <- array(0.1, c(K, num_freqs))
num_tau_check <- 20
show_tau_tune_summ <- TRUE
doCayleyZeros <- FALSE
CayleyZeroProb <- 0.5

omegaw <- seq(1/Tt, by = 1/Tt, length.out = num_freqs)

gibbsIts <- 3000
burnin <- 0.5
gibbsPrint <- 100

# hyperparameters to the prior for tau2
tau2_a <- 1
tau2_b <- 1

### adapting MH tuning parameters
# what s index is the maximum burnin iteration?
burninS <- floor(gibbsIts * burnin)
# how often should the adaptation happen?
tau_numin <- floor(burninS / (num_tau_check))
# at which s iterations should the adaptation happen?
tau_s_check <- seq(tau_numin, burninS, tau_numin)

# adaptive version
# tau_s_check <- seq(50, gibbsIts, 50)

# generate some true parameters -------------------------------------------

set.seed(parseed)
sigmak02 <- c(5, 10)
Sigma0 <- rFTCW(diag(P), P+1, P)

# TODO perhaps figure out a way to sample this from CMACG, but still get a real-valued 
# semi-orthogonal matrix, rather than complex-valued
# rCMACG <- function(nrow, ncol, Sigma) {
U_k0 <- array(NA, c(P, d, K))
for (k in 1:K) {
    # U <- qr.Q(qr(matrix(rnorm(P*d), ncol = d)))
    # U_k0[, , k] <- U
    
    U_k0[, , k] <- rCMACG(P, d, Sigma0) |> closest_semiorth()
}

# avg_Uk0 <- apply(U_k0, c(1, 2), mean)
# avg_Uk0 <- avg_Uk0 %*% 
#     solve( EigenR::Eigen_sqrt( t(Conj(avg_Uk0)) %*% avg_Uk0 ) )
# 
# avg_Uk0_perp <- (qr(avg_Uk0, complete = TRUE) |> 
#                      qr.Q(complete = TRUE))[, (d+1):P]
# 
# V_Sigma0 <- cbind(avg_Uk0, avg_Uk0_perp)
# Sigma0 <- V_Sigma0 %*% diag(P*(P:1) / sum(P:1)) %*% t(Conj(V_Sigma0))

# AR(2) parameters
# being put into an array columnwise
# so the first 4 numbers are two sets of AR(2) parameters for k = 1, and so on
# one pair of SDFs starts out separated by 1
# phis <- c(1, -.25, -.75, -.5,
#             .75, -.65, -1, -.05) |>
#     array(c(2, d, K))
# all that is different is the final AR(2) parameter; now the SDF is less extreme
phis <- c(1, -.25, -.75, -.5,
            .75, -.65, -1, -.5) |>
    array(c(2, d, K))

Lambdak0_w <- array(NA, c(d, K, num_freqs))

for (w in 1:num_freqs) {
    for (k in 1:K) {
        for (j in 1:d) {
            Lambdak0_w[j, k, w] <- sdf_ar2(omegaw[w], phis[, j, k])
        }
    }
}

# simulate time series ----------------------------------------------------
TS_data <- list()
SDMests <- list()

len_freq <- Tt/2

Utp <- waveslim::sine.taper(Tt, LL)

set.seed(dataseed)
for (k in 1:K) {
    Xm <- matrix(NA, Tt, d)
    
    for (j in 1:d) {
        Xm[, j] <- arima.sim(list(ar = phis[, j, k]), Tt)
    }
    X <- ts(Xm, start = c(1, 1), frequency = 1)
    
    Ymat <- sqrt(sigmak02[k]) * (U_k0[, , k] %*% t(X) + 
                                     matrix(rnorm(P*Tt), ncol = Tt))
    Yts <- ts(t(Ymat), start = c(1, 1), freq = 1)
    
    TS_data[[k]] <- Yts
    
    # estimate SDMs with multitaper
    Y_tp <- apply(Utp, MARGIN = 2, function(u) u*Yts, simplify = FALSE)
    F_tp_list <- lapply(Y_tp, FUN = function(Y) astsa::mvspec(Y,plot = FALSE))
    F_tp1 <- array(0, c(P, P, len_freq))
    
    for (ell in 1:len_freq) {
        for(j in 1:length(F_tp_list)){
            F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
        }
        F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
    }
    
    SDMests[[k]] <- F_tp1 * Tt
}

# Data (scaled SDM estimates) ---------------------------------------------

# generate the K observed P-dim. time series

trueSDMs <- array(NA, c(K, P, P, num_freqs))

# calculate true SDMs as data
data_list_w <- list()

# I need to do this in the opposite order, because the way Dr. Namdari's code
# works, I get all frequencies for only 1 observation (k) at a time.

for (w in 1:num_freqs){
    data_list <- list()
    
    for (k in 1:K) {
        # use the exact, known SDM
        # data_list[[k]] <-
        #     LL * sigmak02[k] * ( U_k0[, , k] %*% diag(Lambdak0_w[, k, w]) %*%
        #                             t(Conj(U_k0[, , k])) + diag(P) )
        
        trueSDMs[k, , , w] <- 
            sigmak02[k] * ( U_k0[, , k] %*% diag(Lambdak0_w[, k, w]) %*%
                                t(Conj(U_k0[, , k])) + diag(P) )
        
        # use the estimated SDM
        data_list[[k]] <- LL * SDMests[[k]][, , w]
    }
    
    data_list_w[[w]] <- data_list
}

# create parameter lists for sampling functions ---------------------------

param_list_w0 <- list(
    P = P, d = d, K = K, n_k = rep(LL, K),
    U_ks = U_k0, Lambda_ks = Lambdak0_w[, , 1],
    sigma_k2s = sigmak02, Sigmas = Sigma0)

# initialize arrays -------------------------------------------------------

# U_ks_all <- array(NA, c(P, d, K, gibbsIts))
# U_ks_all[, , , 1] <- U_k0
U_kls_all <- array(NA, c(P, d, K, num_freqs, gibbsIts))
U_kls_all[, , , , 1] <- U_k0

Lambdak_w_s <- array(NA, c(d, K, num_freqs, gibbsIts))
Lambdak_w_s[, , , 1] <- Lambdak0_w

taujk2_s <- array(NA, c(d, K, gibbsIts))
taujk2_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)

zetajk2_s <- array(NA, c(d, K, gibbsIts))
zetajk2_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)

accCount_s <- array(NA, c(K, num_freqs, gibbsIts))
accCount_s[, , 1] <- TRUE

# start cluster -----------------------------------------------------------

# library(foreach)
# library(doParallel)
# 
# n_cores <- 2
# cluster <- makeCluster(n_cores)
# registerDoParallel(cluster)

# sampling ----------------------------------------------------------------

{
starttime <- Sys.time()
for (s in 2:gibbsIts) {
    if (s %% gibbsPrint == 0) cat(paste0(s, ": ", Sys.time(), "\n"))
        
    # a temporary matrices with the previous samples
    U_kls <- U_kls_all[, , , , s-1]
    result_Lambdas <- Lambdak_w_s[, , , s-1]
    
    ##### Lambda sampling
    
    # section for frequency l = 1
    for (k in 1:K) {
        LSkw <- data_list_w[[1]][[k]]
        Ukl <- U_kls[, , k, 1]
        
        # for each Lambdak value, j
        for (j in sample(d)) {
            # xi_jk has a truncated Gamma(nk, tjk) distribution (shape, rate), 
            # and tjk is a temporary parameter defined as follows.
            # (take the Real part, since the quadratic form should be real, 
            # but may have a small complex part due to numerical issues.)
            tjk <- Re( t(Conj(Ukl[, j])) %*% LSkw %*% Ukl[, j] ) / sigmak02[k]
            
            ##### sample from a truncated Gamma distribution
            # lower and upper bounds for truncating the Gamma distribution
            if (j == 1) lb <- 0 else lb <- 1 / (result_Lambdas[j-1, k, 1] + 1)
            if (j == d) ub <- 1 else ub <- 1 / (result_Lambdas[j+1, k, 1] + 1)
            
            # distr <- Runuran::udgamma(shape = LL+1, scale = 1/tjk, lb, ub)
            # gen <- Runuran::arsd.new(distr)
            # xi_jk <- Runuran::ur(gen, 1)
            
            # by slice sampling, instead
            xi_jk <- uni.slice(1/(1 + result_Lambdas[j, k, 1]),
                               dgamma,
                               w = 1, m = Inf, lower = lb, upper = ub,
                               shape = LL+1, scale = 1/tjk,
                               log = TRUE)
                                    
            # convert xi to Lambda
            result_Lambdas[j, k, 1] <- 1/xi_jk - 1
        }
    }
    
    for ( w in sample(2:(num_freqs-1)) ) {
        for (k in 1:K) {
            Ukw <- U_kls[, , k, w]
            Dkw1 <- data_list_w[[w]][[k]]
            
            for (j in sample(d)) {
                ajk <- t(Conj(Ukw[, j])) %*% Dkw1 %*% Ukw[, j] / 
                    sigmak02[k]
                # should be strictly real
                ajk <- Re(ajk)
                
                # 2nd order RW, previous and next neighbors
                lp <- result_Lambdas[j, k, w-1]
                ln <- result_Lambdas[j, k, w+1]
                mu <- (lp + ln)/2
                
                # 2nd order RW, previous two
                # lm1 <- result_Lambdas[j, k, w-1]
                # lm2 <- result_Lambdas[j, k, w-2]
                # mu <- 2*lm1 - lm2
                
                # by slice sampling, instead
                newdraw <- uni.slice(result_Lambdas[j, k, w], unnorm_logPDF,
                                     w = 1, m = Inf, lower = 0, upper = Inf, 
                                     tau2 = taujk2_s[j, k, s-1], 
                                     ajk = ajk, mu = mu, N = LL, 
                                     logscale = TRUE)    
                result_Lambdas[j, k, w] <- newdraw
    
            } # end of sampling j in 1:d
        } # end of sampling k in 1:K
    } # end of sampling over frequencies
    
    ##### for the last frequency, l = num_freq
    # 1st order RW
    for (k in 1:K) {
        Uw1 <- U_kls[, , k, num_freqs]
        Dkw1 <- data_list_w[[num_freqs]][[k]]

        for (j in sample(d)) {
            ajk <- t(Conj(Uw1[, j])) %*% Dkw1 %*% Uw1[, j] / 
                sigmak02[k]
            # should be strictly real - this is done in the existing Lambda
            # FCD sampler.
            ajk <- Re(ajk)
            
            # 1st order random walk
            mu <- result_Lambdas[j, k, num_freqs - 1]
            
            # by slice sampling, instead
            newdraw <- uni.slice(result_Lambdas[j, k, num_freqs], unnorm_logPDF,
                                 w = 1, m = Inf, lower = 0, upper = Inf, 
                                 # TODO change to a variable, not fixed value
                                 tau2 = zetajk2_s[j, k, s-1], 
                                 ajk = ajk, mu = mu, N = LL, 
                                 logscale = TRUE)    
            result_Lambdas[j, k, num_freqs] <- newdraw
            
        } # end of sampling j in 1:d
    } # end of sampling k in 1:K
    
    # save the sampled Lambdas    
    Lambdak_w_s[, , , s] <- result_Lambdas
        
    # sample the separate taujk2 and zetajk2 parameters
    for (j in 1:d) {
        for (k in 1:K) {
            
            # 2nd order RW, next and previous neighbors
            sumalljk <- sum(
                (result_Lambdas[j, k, 2:(num_freqs-1)] -
                     .5*result_Lambdas[j, k, 1:(num_freqs-2)] -
                     .5*result_Lambdas[j, k, 3:(num_freqs)])^2)

            taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-2),
                                          rate = tau2_b + .5*sumalljk)

            # 2nd order RW, previous two
            # sumalljk <- sum(
            #     (result_Lambdas[j, k, 3:num_freqs] -
            #          result_Lambdas[j, k, 2:(num_freqs-1)] +
            #          .5*result_Lambdas[j, k, 1:(num_freqs-2)])^2)
            # 
            # taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-2),
            #                               rate = tau2_b + .5*sumalljk)
        
            # 1st order RW
            # sumalljk <- sum(
            #     (result_Lambdas[j, k, 2:(num_freqs-1)] -
            #          .5*result_Lambdas[j, k, 1:(num_freqs-2)])^2)
            # 
            # taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-1),
            #                               rate = tau2_b + .5*sumalljk)
            
            # zetajk2 sampling
            # based on 1st order random walk for the last Lambda
            sum_zeta <- (result_Lambdas[j, k, num_freqs] - 
                             result_Lambdas[j, k, num_freqs-1])^2
            zetajk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5,
                                           rate = tau2_b + .5*sum_zeta)
        }
    }
    
    ##### Uk sampling 
    
    # before sampling
    accCount <- array(TRUE, c(K, num_freqs))
    newU_kls <- array(NA, c(P, d, K, num_freqs))

    # during sample, 1:K, 1:num_freqs
    
    # result_list_Uk <- foreach(k = 1:K) %dopar% {
    for (k in 1:K) {
        sigma_k2 <- sigmak02[k]
        Sigmas <- Sigma0
        invSigmas <- solve(Sigmas)

        # newU_k <- array(NA, c(P, d, num_freqs))
        # accCount_k <- rep(TRUE, num_freqs)

        # calculate traces and other terms, for the MH acceptance ratio
        for (l in 1:num_freqs) {
            Us <- U_kls[, , k, l]
            
            ### MH sampling
            # propose U by small rotation of Us
            Up <- Uk_MH_Cayley(Us, tau_Ukl[k, l], doCayleyZeros, CayleyZeroProb)
            # let Skl be the scaled data matrix (i.e. SDM est. at freq. index l)
            Skl <- data_list_w[[l]][[k]]
            lambdakl <- result_Lambdas[, k, l]
            Omegakl <- diag(lambdakl / (1 + lambdakl))

            # tracep <- Re(sum(diag( Up %*% Omegakl %*% t(Conj(Up)) %*% Skl )))
            # traces <- Re(sum(diag( Us %*% Omegakl %*% t(Conj(Us)) %*% Skl )))
            tracep <- Re(sum(Conj(Up) * (Skl %*% Up %*% Omegakl)))
            traces <- Re(sum(Conj(Us) * (Skl %*% Us %*% Omegakl)))

            # remaining terms
            logdetp <- log(Re(EigenR::Eigen_det( t(Conj(Up)) %*%
                                                     invSigmas %*% Up )))
            logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*%
                                                     invSigmas %*% Us )))

            # calculate acceptance ratio
            logr <- tracep/sigma_k2 - P*logdetp - traces/sigma_k2 + P*logdets

            # accept or reject
            if (log(runif(1)) <= logr) {
                newU_kls[, , k, l] <- Up
                # newU_k[, , l] <- Up
            } else {
                newU_kls[, , k, l] <- Us
                accCount[k, l] <- FALSE
                # newU_k[, , l] <- Us
                # accCount_k <- FALSE
            }
            
            ### Geodesic slice sampling
            # newU_kls[, , k, l] <-
            #     geoSS(Us, logd_Uk, w = w, m = m,
            #           invSigmas = invSigmas, lw = l, sigma_k2 = sigma_k2,
            #           result_Lambdas = result_Lambdas,
            #           data_list_w = data_list_w, num_freqs = num_freqs, k = k)
            
        } # end of summing over frequencies
        
        # return this value for the dopar for loop
        # list(U = newU_k, acc = accCount_k)
    } # end of sampling over 1:K
    
    # save sampled values from the dopar list
    # for (k in 1:K) {
    #     U_kls_all[, , k, , s] <- result_list_Uk[[k]]$U
    #     accCount_s[k, , s] <- result_list_Uk[[k]]$acc
    # }
    
    U_kls_all[, , , , s] <- newU_kls
    accCount_s[, , s] <- accCount
    
    ### do adaptation of the Ukl MH tuning parameter
    if (s %in% tau_s_check) {
        cat(paste0("\n", s, ": Check for tau_Ukl adaptation\n"))
        
        curr_Ukl_acc_rate <-
            apply(accCount_s[, , (s - tau_numin + 1):s], c(1, 2), mean)

        # print(tau_Ukl[curr_Ukl_acc_rate == 0])
        print(which(curr_Ukl_acc_rate == 0))

        tau_Ukl[curr_Ukl_acc_rate >= .45] <-
            tau_Ukl[curr_Ukl_acc_rate >= .45] * 2
        tau_Ukl[curr_Ukl_acc_rate <= .15] <-
            tau_Ukl[curr_Ukl_acc_rate <= .15] / 4

        # print(tau_Ukl[curr_Ukl_acc_rate == 0])

        if (show_tau_tune_summ) {
            apply(accCount_s[, , (s - tau_numin + 1):s], c(1, 2), mean) |>
                t() |>
                summary() |>
                print()
        }
        
        # Roberts and Rosenthal, Adaptive Metropolis within Gibbs
        # modify the parameter by this amount
        # deltan <- min(0.01, (s%/%50)^(-.5))
        # 
        # # calculate the recent acceptance rates
        # curr_Ukl_acc_rate <-
        #     apply(accCount_s[, , (s - 50 + 1):s], c(1, 2), mean)
        # 
        # # add or subtract a factor to the log st. dev. tuning parameter
        # tau_Ukl[curr_Ukl_acc_rate >= .234] <-
        #     exp(log(tau_Ukl[curr_Ukl_acc_rate >= .234]) + deltan)
        # tau_Ukl[curr_Ukl_acc_rate < .234] <-
        #     exp(log(tau_Ukl[curr_Ukl_acc_rate < .234]) - deltan)
        # 
        # if (show_tau_tune_summ) {
        #     apply(accCount_s[, , (s - 50 + 1):s], c(1, 2), mean) |>
        #         t() |>
        #         summary() |>
        #         print()
        # }
    }
    
} # end of s, Gibbs sampling
endtime <- Sys.time()
}

# evaluate ----------------------------------------------------------------
stopCluster(cl = cluster)

print(endtime - starttime)

gibbsPostBurn <- round(gibbsIts * burnin):gibbsIts

apply(accCount_s[, , 2:gibbsIts], c(1, 2), mean) |> t() |> summary()
apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) |> t() |> summary()

apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) |>
    as.vector() |> 
    quantile(probs = c(0, 0.025, .25, .5, .75, .975, 1))

apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) |> t() |> 
    as.vector() |> plot()
abline(v = num_freqs, h = c(.15, .45), lty = 2)

# get matrix indices of the Ukls with certain acceptance rates
(apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) <= 0.15) |> 
    which(arr.ind = TRUE)

(apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) >= 0.5) |> 
    which(arr.ind = TRUE)

# what are the indices of the Ukls with the LOWEST acc. rates?
(as.vector(apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean)) |> 
    order())[1:10]

# what are the indices of the Ukls with the HIGHEST acc. rates?
(as.vector(apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean)) |> 
        order(decreasing = TRUE))[1:10]

which(apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) == 0)
tau_Ukl[which(apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) == 0)]

plot(Lambdak_w_s[1, 1, 3, 2:gibbsIts], type = "l")
plot(Lambdak_w_s[2, 1, 3, 2:gibbsIts], type = "l")

plot(Lambdak_w_s[1, 2, 3, 2:gibbsIts], type = "l")
plot(Lambdak_w_s[2, 2, 3, 2:gibbsIts], type = "l")

# evaluate Lambdas --------------------------------------------------------

upper_q <- apply(Lambdak_w_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.975)
lower_q <- apply(Lambdak_w_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.025)
# Lambda_means <- apply(Lambdak_w_s, c(1, 2, 3), mean)
Lambda_means <- apply(Lambdak_w_s[, , , gibbsPostBurn], c(1, 2, 3), median)

k <- 1

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, ylim = c(0, 17),
     main = paste0("k = ", k))
lines(lower_q[1, k, ], type = "l", lty = 2)
lines(Lambda_means[1, k, ])
lines(Lambdak0_w[1, k, ], lty = 3)

lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
lines(Lambda_means[2, k, ], , col = "red")
lines(Lambdak0_w[2, k, ], lty = 3, , col = "red")

k <- 2

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, 
     # xlim = c(400, 510),
     ylim = c(0, 11),
     main = paste0("k = ", k)
)
lines(lower_q[1, k, ], type = "l", lty = 2)
lines(Lambda_means[1, k, ])
lines(Lambdak0_w[1, k, ], lty = 3)

lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
lines(Lambda_means[2, k, ], , col = "red")
lines(Lambdak0_w[2, k, ], lty = 3, , col = "red")

# look at densities for one set of Lambda values --------------------------

k <- 1
l <- 273

for (j in 1:d) {
    if (j == 1) {
        plot(density(Lambdak_w_s[j, k, l, gibbsPostBurn]),
             main = paste0("Lambda: k = ", k, ", l = ", l),
             col = j)
    } else {
        lines(density(Lambdak_w_s[j, k, l, gibbsPostBurn]), col = j)
    }
    abline(v = c(lower_q[j, k, l], upper_q[j, k, l]), lty = 2, col = j)
    abline(v = Lambda_means[j, k, l], col = j)
    abline(v = Lambdak0_w[j, k, l], lty = 3, col = j)
}

legend(x = "topright", legend = c("density", "95% CI", "post. mean", "truth"),
       lty = c(1, 2, 1, 3), lwd = 2,
       y.intersp = 0.8)
# distances to true Ukl0 for all Ukl --------------------------------------

ds_to_true <- array(NA, c(K, num_freqs))

for (k in 1:K) {
    for (l in 1:num_freqs) {
        avgUkl <- apply(U_kls_all[, , k, l, gibbsPostBurn],
                        c(1, 2), mean)
        avgUkl <- avgUkl %*% 
            solve( EigenR::Eigen_sqrt( t(Conj(avgUkl)) %*% avgUkl ) )
        
        ds_to_true[k, l] <- fast_evec_Frob_stat(U_k0[, , k], avgUkl)
    }
}

summary(as.vector(ds_to_true))
quantile(as.vector(ds_to_true), probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
sum(ds_to_true >= sqrt(2))
which(ds_to_true >= sqrt(2), arr.ind = TRUE)

which(ds_to_true <= .3, arr.ind = TRUE)

which.max(ds_to_true)
which.min(ds_to_true)

# distance and trace plot for one Ukl -------------------------------------

k <- 1
l <- 273

apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean)[k, l]

avgUkl <- apply(U_kls_all[, , k, l, gibbsPostBurn],
               c(1, 2), mean)
avgUkl <- avgUkl %*% solve( EigenR::Eigen_sqrt( t(Conj(avgUkl)) %*% avgUkl ) )
print(fast_evec_Frob_stat(U_k0[, , k], avgUkl))

# avgUks[, , k] <- avgUk

d_to_avgUkl <- rep(NA, gibbsIts)
for (s in 1:gibbsIts) {
    d_to_avgUkl[s] <- fast_evec_Frob_stat(U_kls_all[, , k, l, s], avgUkl)
    # d_to_avgUkl[s] <- fast_evec_Frob_stat(U_kls_all[, , k, l, s], U_k0[, , k])
}

plot(d_to_avgUkl, type = "l",
     main = paste0("k = ", k, ", l = ", l))

evec_Frob_stat(U_k0[, , k], avgUkl, returnMats = TRUE)
evec_Frob_stat(avgUkl, U_k0[, , k], returnMats = TRUE)

# evaluate tau and zeta ---------------------------------------------------

summary(taujk2_s[1, 1, ])
summary(taujk2_s[1, 2, ])
summary(taujk2_s[2, 1, ])
summary(taujk2_s[2, 2, ])

summary(zetajk2_s[1, 1, ])
summary(zetajk2_s[2, 1, ])
summary(zetajk2_s[1, 2, ])
summary(zetajk2_s[2, 2, ])


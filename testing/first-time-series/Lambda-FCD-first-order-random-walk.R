source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Lambdak.R")

# functions ---------------------------------------------------------------

sdf_ar2 <- function(omega, phis, sigma2 = 1) {
    phi1 <- phis[1]
    phi2 <- phis[2]
    
    denom <- 1 + phi1^2 + phi2^2 + 2*phi1*(phi2 - 1) * cos(2*pi*omega) - 
        2*phi2*cos(4*pi*omega)
    
    return(sigma2 / denom)
}

# True parameters ---------------------------------------------------------

# dataseed <- 21092025
dataseed <- 22092025
parseed <- 314

P <- 4
d <- 2
K <- 2
Tt <- 1024 # length of time series
LL <- round(sqrt(Tt))

num_freqs <- 384

omegaw <- seq(1/Tt, by = 1/Tt, length.out = num_freqs)

gibbsIts <- 2000

tau2_Lambda <- 0.035

set.seed(parseed)
sigmak02 <- c(5, 10)
Sigma0 <- rFTCW(diag(P), P+1, P)

U_k0 <- array(NA, c(P, d, K))
for (k in 1:K) {
    U <- qr.Q(qr(matrix(rnorm(P*d), ncol = d)))
    U_k0[, , k] <- U
}

# AR(2) parameters
# one pair of SDFs starts out separated by 1
phis <- c(1, -.25, -.75, -.5,
            .75, -.65, -1, -.05) |>
    array(c(2, d, K))

# one pair of SDFs starts out separated by 0.1
# phis <- c(1, -.25, -.75, -.5,
#             .15, -.85, -1, -.05) |>
#     array(c(2, d, K))

Lambdak0_w <- array(NA, c(d, K, num_freqs))

for (w in 1:num_freqs) {
    for (k in 1:K) {
        for (j in 1:d) {
            Lambdak0_w[j, k, w] <- sdf_ar2(omegaw[w], phis[, j, k])
        }
    }
}

# plot true SDFs for the factors ------------------------------------------

# k = 2
curve(sdf_ar2(x, phis[, 1, 1]), 
      from = 0, to = 1/2, n = 501, ylab = "spec. density", xlab = "frequency")
curve(sdf_ar2(x, phis[, 2, 1]), 
      from = 0, to = 1/2, add = TRUE, n = 501, col = "red")

Lambdak0_w[1:2, , 1]

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

# Data (exact SDMs, for now) ----------------------------------------------

# generate the K observed P-dim. time series

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
    
#param_list_w1 <- list(
 #   P = P, d = d, K = K, n_k = rep(LL, K),
 #   U_ks = U_k0, Lambda_ks = Lambdak0_w1,
 #   sigma_k2s = sigmak02, Sigmas = Sigma0)

# initialize arrays -------------------------------------------------------

Lambdak_w_s <- array(NA, c(d, K, num_freqs, gibbsIts))

# sampling loop -----------------------------------------------------------

# sampling for lambda, omega_0: can use current FCD sampler, with ordering TRUE

for (s in 1:gibbsIts) {
    
    if (s %% 100 == 0) cat(paste0(s, ": ", Sys.time(), "\n"))
    
    param_list_w0$Lambda_ks <-
        Lambdak_gibbs_densCovar(data_list_w[[1]], param_list_w0, 
                                doOrdered = TRUE)
    Lambdak_w_s[, , 1, s] <- param_list_w0$Lambda_ks
    
    # sampling for lambda, omega_1: need new code, with a Gamma distribution
    # and in fact, we're sampling xijk, for omega_1, so be sure to do the transformation
    
    result_Lambdas <- Lambdak_w_s[, , , s]
    
    for (w in 2:num_freqs) {
        for (k in 1:K) {
            # TODO need to change this to be the actual updated values for a particular frequency.
            # for now, I am not sampling them, so they are constant and do not need proper indexing.
            # Uw1 <- param_list_w1$U_ks[, , k]
            Uw1 <- U_k0[, , k]
            Dkw1 <- data_list_w[[w]][[k]]
            
            for (j in 1:d) {
                xi_jks_w0 <- 1 / (1 + result_Lambdas[j, k, w-1])
                
                ajk <- t(Conj(Uw1[, j])) %*% Dkw1 %*% Uw1[, j] / 
                    sigmak02[k]
                # should be strictly real - this is done in the existing Lambda
                # FCD sampler.
                ajk <- Re(ajk)
                
                # Sample a new value from a Gamma distribution, constrained to (0, 1)
                # g_shape <- param_list_w1$n_k[k] + xi_jks_w0^2 / tau2_Lambda
                g_shape <- LL + xi_jks_w0^2 / tau2_Lambda
                g_scale <- (ajk + xi_jks_w0 / tau2_Lambda)^-1
                distr <- Runuran::udgamma(shape = g_shape, scale = g_scale, 0, 1)
                #gen <- Runuran::pinvd.new(distr)
                gen <- Runuran::arsd.new(distr)
                xi_jk <- Runuran::ur(gen, 1)
                
                result_Lambdas[j, k, w] <- 1/xi_jk - 1
            }
        }
        
        # eventually, update the param_list entry
        # or maybe, something slightly different 
        # - I should just return the d x K array of sampled Lambda values
        #param_list_w1$Lambda_ks <- result_Lambdas
        #Lambdak_w1_s[, , s] <- param_list_w1$Lambda_ks
    }    
    Lambdak_w_s[, , , s] <- result_Lambdas
}

# assess ------------------------------------------------------------------

w <- 2
gibbsPostBurn <- (gibbsIts/2):gibbsIts

# show the true values
# Lambdak0_w
Lambdak0_w[, , 1]

plot(Lambdak_w_s[1, 2, w, gibbsPostBurn], type = "l")
plot(Lambdak_w_s[2, 2, w, gibbsPostBurn], type = "l")

mean(Lambdak_w_s[1, 2, w, gibbsPostBurn])
quantile(Lambdak_w_s[1, 2, w, gibbsPostBurn], c(.025, .5, .975))

mean(Lambdak_w_s[2, 2, w, gibbsPostBurn])
quantile(Lambdak_w_s[2, 2, w, gibbsPostBurn], c(.025, .5, .975))

# plot density - not multi-modal (good), but overlapping (understandable)
plot(density(Lambdak_w_s[1, 2, w, gibbsPostBurn]))
lines(density(Lambdak_w_s[2, 2, w, gibbsPostBurn]), col = "red")
abline(v = Lambdak0_w[1:2, 2, w], col = 1:2)

# apply(Lambdak_w_s, c(1, 2, 3), mean)

upper_q <- apply(Lambdak_w_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.975)
lower_q <- apply(Lambdak_w_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.025)
Lambda_means <- apply(Lambdak_w_s[, , , gibbsPostBurn],
                      c(1, 2, 3), mean)

k <- 2

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, ylim = c(0, 13))
lines(lower_q[1, k, ], type = "l", lty = 2)
lines(Lambda_means[1, k, ])
lines(Lambdak0_w[1, k, ], lty = 3)

lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
lines(Lambda_means[2, k, ], , col = "red")
lines(Lambdak0_w[2, k, ], lty = 3, , col = "red")

mean(Lambdak_w_s[2, 2, 32, ])


# the other k -------------------------------------------------------------

k <- 1

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, ylim = c(0, 18))
lines(lower_q[1, k, ], type = "l", lty = 2)
lines(Lambda_means[1, k, ])
lines(Lambdak0_w[1, k, ], lty = 3)

lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
lines(Lambda_means[2, k, ], , col = "red")
lines(Lambdak0_w[2, k, ], lty = 3, , col = "red")

# compare densities for all frequencies of one entry ----------------------

plot(density(Lambdak_w_s[1, 2, 1, ]))

for (w in 2:num_freqs) {
    lines(density(Lambdak_w_s[1, 2, w, ]), col = w)
}

abline(v = Lambdak0_w[1, 2, ])

plot(density(Lambdak_w_s[2, 2, 1, ]))

for (w in 2:num_freqs) {
    lines(density(Lambdak_w_s[2, 2, w, ]), col = w)
}

abline(v = Lambdak0_w[2, 2, ])

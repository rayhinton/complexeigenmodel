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

# True parameters ---------------------------------------------------------

# dataseed <- 21092025
# dataseed <- 22092025
# dataseed <- 10102025
dataseed <- 17102025
parseed <- 314

P <- 4
d <- 2
K <- 2
Tt <- 2048 # length of time series
LL <- round(sqrt(Tt))

num_freqs <- Tt/2 - 2

omegaw <- seq(1/Tt, by = 1/Tt, length.out = num_freqs)

gibbsIts <- 2000

tau2_Lambda <- 0.5
# hyperparameters to the prior for tau2
tau2_a <- 1
tau2_b <- 1

set.seed(parseed)
sigmak02 <- c(5, 10)
Sigma0 <- rFTCW(diag(P), P+1, P)

U_k0 <- array(NA, c(P, d, K))
for (k in 1:K) {
    U <- qr.Q(qr(matrix(rnorm(P*d), ncol = d)))
    U_k0[, , k] <- U
}

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

k <- 2
curve(sdf_ar2(x, phis[, 1, k]), 
      from = 0, to = 1/2, n = 501, 
      ylab = "spec. density", xlab = "frequency",
      ylim = c(0, 16))
curve(sdf_ar2(x, phis[, 2, k]), 
      from = 0, to = 1/2, add = TRUE, n = 501, col = "red")

Lambdak0_w[1:2, , 1]

# simulate time series ----------------------------------------------------
TS_data <- list()
SDMests <- list()
kernel_SDMests <- list()

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

    kernel_SDMests[[k]] <- 
        # astsa::mvspec(Yts, kernel("modified.daniell", c(8,8,8)), 
        astsa::mvspec(Yts, kernel("fejer", 100, 6), 
                      log = 'no', taper = .5, plot = FALSE)
}

plot(pi*Re(kernel_SDMests[[1]]$fxx[1, 1, 1 + 2:num_freqs]), type = "l", 
     ylim = c(0, 60))
lines(pi*Re(kernel_SDMests[[1]]$fxx[2, 2, 1 + 2:num_freqs]), col = 2)
lines(pi*Re(kernel_SDMests[[1]]$fxx[3, 3, 1 + 2:num_freqs]), col = 3)
lines(pi*Re(kernel_SDMests[[1]]$fxx[4, 4, 1 + 2:num_freqs]), col = 4)

plot(pi*Re(kernel_SDMests[[2]]$fxx[1, 1, 1 + 2:num_freqs]), type = "l", 
     ylim = c(0, 80))
lines(pi*Re(kernel_SDMests[[2]]$fxx[2, 2, 1 + 2:num_freqs]), col = 2)
lines(pi*Re(kernel_SDMests[[2]]$fxx[3, 3, 1 + 2:num_freqs]), col = 3)
lines(pi*Re(kernel_SDMests[[2]]$fxx[4, 4, 1 + 2:num_freqs]), col = 4)

# comparing scales --------------------------------------------------------

kernel_SDMests[[1]]$freq[1]
kernel_SDMests[[1]]$freq[Tt/2]

# Data (exact SDMs, for now) ----------------------------------------------

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
        data_list[[k]] <- LL * SDMests[[k]][, , w+1]
    }
    
    data_list_w[[w]] <- data_list
}


# plot estimated 4x4 SDMs -------------------------------------------------

dim(SDMests)
length(SDMests)
dim(SDMests[[1]])

plot(Re(SDMests[[1]][1,1, ]), type = "l", ylim = c(0, 55))
lines(Re(SDMests[[1]][2,2, ]), col = 2)
lines(Re(SDMests[[1]][3,3, ]), col = 3)
lines(Re(SDMests[[1]][4,4, ]), col = 4)

lines(Re(trueSDMs[1, 1, 1, ]), col = 1, lty = 2)
lines(Re(trueSDMs[1, 2, 2, ]), col = 2, lty = 2)
lines(Re(trueSDMs[1, 3, 3, ]), col = 3, lty = 2)
lines(Re(trueSDMs[1, 4, 4, ]), col = 4, lty = 2)

plot(Re(SDMests[[2]][1,1, ]), type = "l", ylim = c(0, 55))
lines(Re(SDMests[[2]][2,2, ]), col = 2)
lines(Re(SDMests[[2]][3,3, ]), col = 3)
lines(Re(SDMests[[2]][4,4, ]), col = 4)

lines(Re(trueSDMs[2, 1, 1, ]), col = 1, lty = 2)
lines(Re(trueSDMs[2, 2, 2, ]), col = 2, lty = 2)
lines(Re(trueSDMs[2, 3, 3, ]), col = 3, lty = 2)
lines(Re(trueSDMs[2, 4, 4, ]), col = 4, lty = 2)

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
Lambdak_w_s[, , , 1] <- Lambdak0_w
# tau2_s <- rep(NA, gibbsIts)
# tau2_s[1] <- 1/rgamma(1, tau2_a, rate = tau2_b)

taujk2_s <- array(NA, c(d, K, gibbsIts))
taujk2_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)

# sampling loop -----------------------------------------------------------
library(Runuran)
Runuran::Runuran.options(error.level="none")
Runuran::Runuran.options("error.level")

# sampling for lambda, omega_0: can use current FCD sampler, with ordering TRUE

print(Sys.time())
# profvis::profvis({
    
for (s in 2:gibbsIts) {
    
    if (s %% 100 == 0) cat(paste0(s, ": ", Sys.time(), "\n"))
    
    #param_list_w0$Lambda_ks <-
    #    Lambdak_gibbs_densCovar(data_list_w[[1]], param_list_w0, 
    #                            doOrdered = TRUE)

    # sampling for lambda, omega_1: need new code, with a Gamma distribution
    # and in fact, we're sampling xijk, for omega_1, so be sure to do the transformation
    
    result_Lambdas <- Lambdak_w_s[, , , s-1]
    result_Lambdas[, , c(1:2, num_freqs)] <- Lambdak0_w[, , c(1:2, num_freqs)]
    
    for (w in sample(3:(num_freqs-1)) ) {
        for (k in 1:K) {
            # TODO need to change this to be the actual updated values for a particular frequency.
            # for now, I am not sampling them, so they are constant and do not need proper indexing.
            # Uw1 <- param_list_w1$U_ks[, , k]
            Uw1 <- U_k0[, , k]
            Dkw1 <- data_list_w[[w]][[k]]
            
            for (j in 1:d) {
                ajk <- t(Conj(Uw1[, j])) %*% Dkw1 %*% Uw1[, j] / 
                    sigmak02[k]
                # should be strictly real - this is done in the existing Lambda
                # FCD sampler.
                ajk <- Re(ajk)
                
                # 2nd order RW 
                lp <- result_Lambdas[j, k, w-1]
                ln <- result_Lambdas[j, k, w+1]
                mu <- (lp + ln)/2

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

    # save the sampled Lambdas    
    Lambdak_w_s[, , , s] <- result_Lambdas
        
    # sample the separate taujk2 parameters
    for (j in 1:d) {
        for (k in 1:K) {
            
            # 2nd order RW, next and previous neighbors
            sumalljk <- sum(
                (result_Lambdas[j, k, 2:(num_freqs-1)] -
                     .5*result_Lambdas[j, k, 1:(num_freqs-2)] -
                     .5*result_Lambdas[j, k, 3:(num_freqs)])^2)

            taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-2),
                                          rate = tau2_b + .5*sumalljk)
        
            # 1st order RW
            # sumalljk <- sum(
            #     (result_Lambdas[j, k, 2:(num_freqs-1)] -
            #          .5*result_Lambdas[j, k, 1:(num_freqs-2)])^2)
            # 
            # taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-1),
            #                               rate = tau2_b + .5*sumalljk)
        }
    }
} # end of Gibbs sampling, over s

# }) # end of profvis

# timing ------------------------------------------------------------------

microbenchmark::microbenchmark(
    sumall <- sum( (result_Lambdas[, , 2:num_freqs] - 
                        result_Lambdas[, , 1:(num_freqs-1)])^2 ),
    times = 1000
)

# assess ------------------------------------------------------------------

# gibbsPostBurn <- 2:gibbsIts
gibbsPostBurn <- (gibbsIts/2):gibbsIts
w <- 1800

# show the true values
# Lambdak0_w
Lambdak0_w[, , w]

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
# Lambda_means <- apply(Lambdak_w_s, c(1, 2, 3), mean)
Lambda_means <- apply(Lambdak_w_s[, , , gibbsPostBurn], c(1, 2, 3), mean)

k <- 2

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, 
     # xlim = c(400, 510),
     ylim = c(0, 11)
     )
lines(lower_q[1, k, ], type = "l", lty = 2)
lines(Lambda_means[1, k, ])
lines(Lambdak0_w[1, k, ], lty = 3)

lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
lines(Lambda_means[2, k, ], , col = "red")
lines(Lambdak0_w[2, k, ], lty = 3, , col = "red")

abline(v = w, lty = 4)

mean(Lambdak_w_s[2, 2, 32, ])

# the other k -------------------------------------------------------------

k <- 1

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, ylim = c(0, 17))
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

# assess tau2 samples -----------------------------------------------------

apply(taujk2_s, c(1, 2), median)

apply(taujk2_s, c(1, 2), quantile, probs = 0.025)
apply(taujk2_s, c(1, 2), quantile, probs = 0.975)

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
unnorm_logPDF <- function(x, tau2, ajk, l1, l2, N, logscale = FALSE) {
    logf1 <- -N * log(1 + x)
    # logf2 <- a * x / (1+x) - (.5/tau2) * (x^2 + 2*x*(l1 - 2*l2))
    logf2 <- -ajk / (1+x) - (.5/tau2) * (x - (2*l2 - l1))^2
    
    result <- ifelse(x <= 0, -Inf, logf1 + logf2)
    
    if (logscale) {
        return(result)
    } else {
        return(exp(result))
    }
}

d_unnorm_logPDF <- function(x, tau2, a, l1, l2, N, logscale = TRUE) {
    f1 <- -N/(1 + x) 
    f2 <- a*(1 / (1+x) -x/(1 + x)^2)
    f3 <- -(.5/tau2) * (2*x + 2*(l1 - 2*l2))
    
    return(f1 + f2 + f3)
}

# True parameters ---------------------------------------------------------

dataseed <- 21092025
# dataseed <- 22092025
parseed <- 314

P <- 4
d <- 2
K <- 2
Tt <- 1024 # length of time series
LL <- round(sqrt(Tt))

num_freqs <- 384

omegaw <- seq(1/Tt, by = 1/Tt, length.out = num_freqs)

gibbsIts <- 100

tau2_Lambda <- 0.5

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
curve(sdf_ar2(x, phis[, 1, 2]), 
      from = 0, to = 1/2, n = 501, ylab = "spec. density", xlab = "frequency")
curve(sdf_ar2(x, phis[, 2, 2]), 
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
library(Runuran)
Runuran::Runuran.options(error.level="none")
Runuran::Runuran.options("error.level")

# sampling for lambda, omega_0: can use current FCD sampler, with ordering TRUE

print(Sys.time())
profvis::profvis({
    
for (s in 1:gibbsIts) {
    
    if (s %% 20 == 0) cat(paste0(s, ": ", Sys.time(), "\n"))
    
    #param_list_w0$Lambda_ks <-
    #    Lambdak_gibbs_densCovar(data_list_w[[1]], param_list_w0, 
    #                            doOrdered = TRUE)
    #Lambdak_w_s[, , 1, s] <- param_list_w0$Lambda_ks
    
    Lambdak_w_s[, , 1:2, s] <- Lambdak0_w[, , 1:2]
    
    # sampling for lambda, omega_1: need new code, with a Gamma distribution
    # and in fact, we're sampling xijk, for omega_1, so be sure to do the transformation
    
    result_Lambdas <- Lambdak_w_s[, , , s]
    
    for (w in 3:num_freqs) {
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
                
                # TODO I had to change this to result_lambdas instead of the
                # larger array, Lambdak_w_s, because values were missing
                l2 <- result_Lambdas[j, k, w-1]
                l1 <- result_Lambdas[j, k, w-2]
                
                # calculate the mode of the FCD, needed for Runuran sampler.
                # this is accomplished by solving a cubic equation
                mu <- 2*l2 - l1
                
                a <- 1
                b <- 2 - mu
                c <- LL*tau2_Lambda + 1 - 2*mu
                dd <- LL*tau2_Lambda - ajk*tau2_Lambda - mu
                
                xroots <- polyroot(c(dd, c, b, a))
                
                DISC <- 18*a*b*c*dd - 4*b^3*dd + b^2*c^2 - 4*a*c^3 - 27*a^2*dd^2
                
                if (DISC < 0) {
                    # TODO this is not ideal, but helps Runuran tabl sampler run
                    # xmode <- max(0, Re(xroots[zapsmall(Im(xroots)) == 0]))
                    xmode <- max(2e-8, Re(xroots[zapsmall(Im(xroots)) == 0]))
                } else if (DISC > 0) {
                    numaccroots <- sum(Re(xroots) >= 0)
                    if (numaccroots <= 1) {
                        xmode <- max(2e-8, Re(xroots))
                    } else {
                        stop("There are multiple acceptable real roots.")
                    }
                } else {
                    # TODO just report an error for now
                    # but, there could be ways to handle if there are "multiple roots;
                    # - good: it could that they are 3 repeated real roots
                    # - maybe good: 1 and 2 mult. roots; maybe one is <= 0, outside domain?
                    stop(paste0("There are multiple real roots. DISC = ", DISC))
                }
                
                # create the generator object
                # gen <- 
                #     Runuran::pinv.new(unnorm_logPDF, lb = 0, ub = Inf,
                #                       center = xmode, islog = TRUE,
                #                       uresolution = 1.0e-9, # 1.0e-10 is default
                #                       smooth = FALSE, # FALSE is default
                #                       tau2 = tau2_Lambda, a = ajk, l1 = l1, l2 = l2, N = LL, 
                #                       logscale = TRUE)

                gen <- 
                Runuran::tabl.new(unnorm_logPDF, 
                                  lb = 0, ub = Inf, mode = xmode, 
                                  islog = TRUE,
                                  tau2 = tau2_Lambda, 
                                  ajk = ajk, l1 = l1, l2 = l2, N = LL, 
                                  logscale = TRUE)
                # draw a random sample
                result_Lambdas[j, k, w] <- Runuran::ur(gen, 1)
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
}) # end of profvis

# assess ------------------------------------------------------------------

w <- 250

# show the true values
# Lambdak0_w
Lambdak0_w[, , w]

plot(Lambdak_w_s[1, 2, w, ], type = "l")
plot(Lambdak_w_s[2, 2, w, ], type = "l")

mean(Lambdak_w_s[1, 2, w, ])
quantile(Lambdak_w_s[1, 2, w, ], c(.025, .5, .975))

mean(Lambdak_w_s[2, 2, w, ])
quantile(Lambdak_w_s[2, 2, w, ], c(.025, .5, .975))

# plot density - not multi-modal (good), but overlapping (understandable)
plot(density(Lambdak_w_s[1, 2, w, ]))
lines(density(Lambdak_w_s[2, 2, w, ]), col = "red")
abline(v = Lambdak0_w[1:2, 2, w], col = 1:2)

# apply(Lambdak_w_s, c(1, 2, 3), mean)

upper_q <- apply(Lambdak_w_s, c(1, 2, 3), quantile, probs = 0.975)
lower_q <- apply(Lambdak_w_s, c(1, 2, 3), quantile, probs = 0.025)
# Lambda_means <- apply(Lambdak_w_s, c(1, 2, 3), mean)
Lambda_means <- apply(Lambdak_w_s, c(1, 2, 3), mean)

k <- 2

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, ylim = c(0, 20))
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
plot(upper_q[1, k, ], type = "l", lty = 2, ylim = c(0, 50))
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

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

# Generate a smooth positive random function using splines
generate_smooth_positive_function <- function(N = 100, n_knots = 7, 
                                              log_mean = 0, log_sd = 1) {
    # Define evaluation points
    x <- 1:N
    
    # Sample knot locations (including boundaries)
    knot_locs <- sort(c(1, sample(2:(N-1), n_knots - 2, replace = FALSE), N))
    
    # Sample log-values at knots
    log_y_knots <- rnorm(n_knots, mean = log_mean, sd = log_sd)
    
    # Fit smooth spline through knots
    spline_fit <- splinefun(knot_locs, log_y_knots, method = "natural")
    
    # Evaluate at all points and exponentiate
    log_f <- spline_fit(x)
    f <- exp(log_f)
    
    return(list(x = x, f = f, knot_locs = knot_locs, log_y_knots = log_y_knots))
}

generate_VAR1_coef <- function(P, max_eigenvalue = 0.9) {
    # Generate random matrix and scale to control eigenvalues
    A1 <- matrix(rnorm(P * P), P, P)
    
    # Eigen decomposition
    eig <- eigen(A1)
    
    # Scale eigenvalues to be inside unit circle
    eig$values <- eig$values / max(Mod(eig$values)) * max_eigenvalue
    
    # Reconstruct matrix with scaled eigenvalues
    A1 <- Re(eig$vectors %*% diag(eig$values) %*% solve(eig$vectors))
    
    return(A1)
}

generate_AR1_covariance <- function(P, sigma2 = 1, rho = 0.5) {
    # Check validity
    if (abs(rho) >= 1) stop("rho must be in (-1, 1)")
    
    # Create AR(1) covariance matrix
    Sigma <- sigma2 * rho^abs(outer(1:P, 1:P, "-"))
    
    return(Sigma)
}

# setup -------------------------------------------------------------------

# dataseed <- 21092025
# dataseed <- 22092025
# dataseed <- 10102025
dataseed <- 17102025
# dataseed <- 28112025

# parseed <- 314
parseed <- 3141
# parseed <- 963456789

P <- 4
d <- 2
K <- 2
Tt <- 1024 # length of time series
LL <- round(sqrt(Tt))
# LL <- round(5)

num_freqs <- Tt/2 - 1

# time series generation parameters
# number of knots in the function that generates Lambda curves
n_knots <- 4

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

# generate true parameters ------------------------------------------------

set.seed(parseed)

VARpars <- array(NA, c(P, P, K))
for (k in 1:K) {
    VARpars[, , k] <- generate_VAR1_coef(P, 0.8)
}
noiseSigma <- generate_AR1_covariance(P, sigma2 = 1, rho = 0.5)

# sigmakl02, scale parameter
# sigmakl02 <- matrix(rgamma(K*(num_freqs + 2), 1, 1), 
#                     K, num_freqs + 2)
# this makes it constant across all frequencies for each k
sigmakl02 <- matrix(rgamma(K, 1, 1),
                    K, Tt)

fkTR <- array(NA, c(P, P, K, Tt))
U_kl0 <- array(NA, c(P, d, K, Tt))
Lambdakl0 <- array(NA, c(d, K, Tt))

for (k in 1:K) {
    for (t in 1:Tt) {
        Hz <- solve( diag(P) - exp(-1i*2*pi * t/Tt) * VARpars[, , k] )
        fomega <- Hz %*% noiseSigma %*% t(Conj(Hz))
        
        f_evd <- eigen(fomega)
        
        thisU <- f_evd$vectors[, 1:d]
        thisLambda <- f_evd$values[1:d]
        
        fkTR[, , k, t] <- sigmakl02[k, t] * ( thisU %*% diag(thisLambda) %*% 
                                          t(Conj(thisU)) + diag(P) )
        U_kl0[, , k, t] <- thisU
        Lambdakl0[, k, t] <- thisLambda
    }
}

# calculate a reasonable Sigmal0 parameter
# # Sigmal0, prior parameter to Ukl0
Sigmal0 <- array(NA, c(P, P, num_freqs))
invSigmal0 <- array(NA, c(P, P, num_freqs))

for (l in 1:num_freqs) {
    avg_Uk0 <- apply(U_kl0[, , , l], c(1, 2), mean)
    avg_Uk0 <- avg_Uk0 %*%
        solve( EigenR::Eigen_sqrt( t(Conj(avg_Uk0)) %*% avg_Uk0 ) )
    
    avg_Uk0_perp <- (qr(avg_Uk0, complete = TRUE) |>
                         qr.Q(complete = TRUE))[, (d+1):P]

    V_Sigma0 <- cbind(avg_Uk0, avg_Uk0_perp)
    Sigmal0[, , l] <- V_Sigma0 %*% diag(P*(P:1) / sum(P:1)) %*% t(Conj(V_Sigma0))
    invSigmal0[, , l] <- solve(Sigmal0[, , l])
}

# check Lambdakl0 ---------------------------------------------------------

k <- 1
plot(Lambdakl0[1, k, 1:num_freqs], type = "l", ylab = "lambda",
     main = paste0("k = ", k))
lines(Lambdakl0[2, k, 1:num_freqs], col = 2)

k <- 2
plot(Lambdakl0[1, k, 1:num_freqs], type = "l", ylab = "lambda",
     main = paste0("k = ", k))
lines(Lambdakl0[2, k, 1:num_freqs], col = 2)

# calculate Cholesky decompositions ---------------------------------------

Rfs <- array(NA, c(P, P, K, Tt))

for (k in 1:K) {
    for (t in 1:Tt) {
        # need lower triangular part, take t(Conj())
        # AFAICT, this is the unique Chol. with positive real diagonals
        Rfs[, , k, t] <- t(Conj(EigenR::Eigen_chol(fkTR[, , k, t])))
    }
}

# generate random vectors -------------------------------------------------

set.seed(dataseed)
Zs <- array(NA, c(P, K, Tt))

for (k in 1:K) {
    Zs[, k, 1:(Tt/2 - 1)] <- rscnorm(P * (Tt/2 - 1)) / sqrt(Tt)
    Zs[, k, c(Tt/2, Tt)] <- rnorm(P * 2) / sqrt(Tt)
    
    Zs[, k, (Tt/2 + 1):(Tt - 1)] <- Conj(Zs[, k, (Tt/2 - 1):1])
}

# generate the time series ------------------------------------------------

Yts <- array(0, c(P, K, Tt))

for (k in 1:K) {
    for (t in 1:Tt) {
        for (l in 1:Tt) {
            Yts[, k, t] <- Yts[, k, t] + 
                Rfs[, , k, l] %*% Zs[, k, l] * exp(2 * pi * 1i * l/Tt * t)
        }
    }
}

# verify the imaginary parts are small (i.e. to numerical precision)
max(abs(Im(Yts)))

# set Y to the real components, only
Yts <- Re(Yts)

plot(Yts[1, 1, ], type = "l")
lines(Yts[2, 1, ], col = 2)
lines(Yts[3, 1, ], col = 3)
lines(Yts[4, 1, ], col = 4)

plot(Yts[1, 2, ], type = "l")
lines(Yts[2, 2, ], col = 2)
lines(Yts[3, 2, ], col = 3)
lines(Yts[4, 2, ], col = 4)

# estimate SDMs -----------------------------------------------------------

SDMests <- list()
TS_data <- list()

len_freq <- Tt/2

Utp <- waveslim::sine.taper(Tt, LL)

for (k in 1:K) {
    thisYts <- ts(t(Yts[, k, ]), start = c(1, 1), freq = 1)
    
    TS_data[[k]] <- thisYts
    
    # estimate SDMs with multitaper
    Y_tp <- apply(Utp, MARGIN = 2, function(u) u*thisYts, simplify = FALSE)
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

# calculate true SDMs as data
data_list_w <- list()

for (w in 1:num_freqs){
    data_list <- list()
    
    for (k in 1:K) {
        # use the estimated SDM
        data_list[[k]] <- LL * SDMests[[k]][, , w]
    }
    
    data_list_w[[w]] <- data_list
}

# compare SDM ests and true SDMs ------------------------------------------

plot(Re(SDMests[[1]][1,1, ]), type = "l", ylim = c(0, 12), 
     ylab = "spectral density", main = paste0("k = ", 1))
lines(Re(SDMests[[1]][2,2, ]), col = 2)
lines(Re(SDMests[[1]][3,3, ]), col = 3)
lines(Re(SDMests[[1]][4,4, ]), col = 4)

lines(Re(fkTR[1, 1, 1, 1:num_freqs]), col = 1, lty = 2)
lines(Re(fkTR[2, 2, 1, 1:num_freqs]), col = 2, lty = 2)
lines(Re(fkTR[3, 3, 1, 1:num_freqs]), col = 3, lty = 2)
lines(Re(fkTR[4, 4, 1, 1:num_freqs]), col = 4, lty = 2)

plot(Re(SDMests[[2]][1,1, ]), type = "l", ylim = c(0, 20),
     ylab = "spectral density", main = paste0("k = ", 2))
lines(Re(SDMests[[2]][2,2, ]), col = 2)
lines(Re(SDMests[[2]][3,3, ]), col = 3)
lines(Re(SDMests[[2]][4,4, ]), col = 4)

lines(Re(fkTR[1, 1, 2, 1:num_freqs]), col = 1, lty = 2)
lines(Re(fkTR[2, 2, 2, 1:num_freqs]), col = 2, lty = 2)
lines(Re(fkTR[3, 3, 2, 1:num_freqs]), col = 3, lty = 2)
lines(Re(fkTR[4, 4, 2, 1:num_freqs]), col = 4, lty = 2)

# initialize arrays -------------------------------------------------------

# U_ks_all <- array(NA, c(P, d, K, gibbsIts))
# U_ks_all[, , , 1] <- U_k0
U_kls_all <- array(NA, c(P, d, K, num_freqs, gibbsIts))
U_kls_all[, , , , 1] <- U_kl0[, , , 1:num_freqs]

Lambdak_w_s <- array(NA, c(d, K, num_freqs, gibbsIts))
Lambdak_w_s[, , , 1] <- Lambdakl0[, , 1:num_freqs]

taujk2_s <- array(NA, c(d, K, gibbsIts))
taujk2_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)

zetajk2_s <- array(NA, c(d, K, gibbsIts))
zetajk2_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)

accCount_s <- array(NA, c(K, num_freqs, gibbsIts))
accCount_s[, , 1] <- TRUE

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
                tjk <- Re( t(Conj(Ukl[, j])) %*% LSkw %*% Ukl[, j] ) / 
                    sigmakl02[k, 1]
                
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
                        sigmakl02[k, w]
                    # should be strictly real
                    ajk <- Re(ajk)
                    
                    # 2nd order RW, previous and next neighbors
                    # lp <- result_Lambdas[j, k, w-1]
                    # ln <- result_Lambdas[j, k, w+1]
                    # mu <- (lp + ln)/2
                    
                    # 2nd order RW, previous two
                    # lm1 <- result_Lambdas[j, k, w-1]
                    # lm2 <- result_Lambdas[j, k, w-2]
                    # mu <- 2*lm1 - lm2
                    
                    # 1st order RW
                    mu <- result_Lambdas[j, k, w-1]
                    
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
                    sigmakl02[k, num_freqs]
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
                # sumalljk <- sum(
                #     (result_Lambdas[j, k, 2:(num_freqs-1)] -
                #          .5*result_Lambdas[j, k, 1:(num_freqs-2)] -
                #          .5*result_Lambdas[j, k, 3:(num_freqs)])^2)
                # 
                # taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-2),
                #                               rate = tau2_b + .5*sumalljk)
                
                # 2nd order RW, previous two
                # sumalljk <- sum(
                #     (result_Lambdas[j, k, 3:num_freqs] -
                #          result_Lambdas[j, k, 2:(num_freqs-1)] +
                #          .5*result_Lambdas[j, k, 1:(num_freqs-2)])^2)
                # 
                # taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-2),
                #                               rate = tau2_b + .5*sumalljk)
                
                # 1st order RW
                sumalljk <- sum(
                    (result_Lambdas[j, k, 2:(num_freqs-1)] -
                        result_Lambdas[j, k, 1:(num_freqs-2)])^2)
                taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-2),
                                            rate = tau2_b + .5*sumalljk)
                
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
            
            # newU_k <- array(NA, c(P, d, num_freqs))
            # accCount_k <- rep(TRUE, num_freqs)
            
            # calculate traces and other terms, for the MH acceptance ratio
            for (l in 1:num_freqs) {
                sigma_k2 <- sigmakl02[k, l]
                # Sigmas <- Sigmal0[, , l]
                invSigmas <- invSigmal0[, , l]
                
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
        
        #####
        # Sigmal sampling
        #####
        
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

print(endtime - starttime)

# check the acceptance rates
gibbsPostBurn <- round(gibbsIts * burnin):gibbsIts

apply(accCount_s[, , 2:gibbsIts], c(1, 2), mean) |> t() |> summary()
apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) |> t() |> summary()

apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) |>
    as.vector() |> 
    quantile(probs = c(0, 0.025, .25, .5, .75, .975, 1))

# evaluate Lambdas --------------------------------------------------------

upper_q <- apply(Lambdak_w_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.975)
lower_q <- apply(Lambdak_w_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.025)
# Lambda_means <- apply(Lambdak_w_s, c(1, 2, 3), mean)
Lambda_means <- apply(Lambdak_w_s[, , , gibbsPostBurn], c(1, 2, 3), median)

k <- 1

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, ylim = c(0, 50),
     main = paste0("k = ", k))
lines(lower_q[1, k, ], type = "l", lty = 2)
lines(Lambda_means[1, k, ])
lines(Lambdakl0[1, k, 1:num_freqs], lty = 3)

lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
lines(Lambda_means[2, k, ], , col = "red")
lines(Lambdakl0[2, k, 1:num_freqs], lty = 3, , col = "red")

k <- 2

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, 
     # xlim = c(400, 510),
     ylim = c(0, 60),
     main = paste0("k = ", k)
)
lines(lower_q[1, k, ], type = "l", lty = 2)
lines(Lambda_means[1, k, ])
lines(Lambdakl0[1, k, 1:num_freqs], lty = 3)

lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
lines(Lambda_means[2, k, ], , col = "red")
lines(Lambdakl0[2, k, 1:num_freqs], lty = 3, , col = "red")

# distances to true Ukl0 for all Ukl --------------------------------------

ds_to_true <- array(NA, c(K, num_freqs))

for (k in 1:K) {
    for (l in 1:num_freqs) {
        avgUkl <- apply(U_kls_all[, , k, l, gibbsPostBurn],
                        c(1, 2), mean)
        avgUkl <- avgUkl %*% 
            solve( EigenR::Eigen_sqrt( t(Conj(avgUkl)) %*% avgUkl ) )
        
        ds_to_true[k, l] <- fast_evec_Frob_stat(U_kl0[, , k, l], avgUkl)
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
l <- 500

apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean)[k, l]

avgUkl <- apply(U_kls_all[, , k, l, gibbsPostBurn],
                c(1, 2), mean)
avgUkl <- avgUkl %*% solve( EigenR::Eigen_sqrt( t(Conj(avgUkl)) %*% avgUkl ) )
print(fast_evec_Frob_stat(U_kl0[, , k, l], avgUkl))

d_to_avgUkl <- rep(NA, gibbsIts)
for (s in 1:gibbsIts) {
    # d_to_avgUkl[s] <- fast_evec_Frob_stat(U_kls_all[, , k, l, s], avgUkl)
    d_to_avgUkl[s] <- fast_evec_Frob_stat(U_kls_all[, , k, l, s], U_kl0[, , k, l])
}

plot(d_to_avgUkl, type = "l",
     main = paste0("k = ", k, ", l = ", l))

evec_Frob_stat(U_kl0[, , k, l], avgUkl, returnMats = TRUE)
evec_Frob_stat(avgUkl, U_kl0[, , k, l], returnMats = TRUE)

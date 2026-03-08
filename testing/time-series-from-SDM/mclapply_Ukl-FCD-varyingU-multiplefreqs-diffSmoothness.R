library(splines)
library(expm)
library(ggplot2)
library(foreach)
library(doParallel)

if (!interactive()) {
    pdf(tempfile(fileext = ".pdf"))
    dev.control("enable")
}

source("testing/time-series-from-SDM/model-simulation-parameters.R")

source("functions/utility.R")
source("functions/FCD_Uk_CMACG.R")
source("functions/rcomplex_wishart.R")
source("functions/matrix_distances.R")
source("testing/first-time-series/geoSS.R")

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
    logf1 <- -N * log(1 + x) - ajk/(1 + x)
    logf2 <-  -(.5/tau2) * (x - mu)^2
    
    result <- ifelse(x <= 0, -Inf, logf1 + logf2)
    
    if (logscale) {
        return(result)
    } else {
        return(exp(result))
    }
}

logd_Uk <- function(Us, lw, invSigmas, sigma_k2, result_Lambdas, data_list_w, 
                    num_freqs, k) {
    # calculate traces and other terms, for the MH acceptance ratio
    Skl <- data_list_w[[lw]][[k]]
    lambdakl <- result_Lambdas[, k, lw]
    Omegakl <- diag(lambdakl / (1 + lambdakl))
    
    traces <- Re(sum(Conj(Us) * (Skl %*% Us %*% Omegakl)))
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

logdet <- function(X) {
    return( log(EigenR::Eigen_det(X)) )
}

save_plot_pdf <- function(plot_path, width = 6, height = 4) {
    p <- recordPlot()
    pdf(plot_path, width = width, height = height)
    replayPlot(p)
    dev.off()    
}

log_and_print_obj <- function(x) {
    print(x)
    sink(file.path(result_dir, "output.txt"), append = TRUE)
    print(x)
    sink()
}

# check simulation parameters ---------------------------------------------

if (gibbsIts * burnin != round(gibbsIts * burnin)) {
    stop(paste0("Burn-in proportion (", burnin, 
                ") must evenly divide the number of iterations (",
                gibbsIts, ")."))
}

if ((gibbsIts * burnin) %% num_tau_check != 0) {
    stop(paste0("Burn-in length (", gibbsIts * burnin,
                ") must be a multiple of the tau adaptation interval (",
                num_tau_check, ")."))
}

if ((gibbsIts * burnin) %% t_thin != 0) {
    stop(paste0("Burn-in length (", gibbsIts * burnin,
                ") must be a multiple of the thinning interval (",
                t_thin, ")."))
}

# setup -------------------------------------------------------------------

# the parameters file is run above

job_id <- Sys.getenv("SLURM_JOB_ID")

# Create timestamped results directory
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")

result_dir <- file.path("results",
                        paste("full-SDM-diffSmoothness", job_id, timestamp, 
                              sep = "_"))
result_dir <- gsub("__", "_", result_dir)

dir.create(result_dir, recursive = TRUE)
dir.create(file.path(result_dir, "post-SDM-est-dist-density"), 
           recursive = TRUE)
dir.create(file.path(result_dir, "compare-SDM-powers"), 
           recursive = TRUE)
dir.create(file.path(result_dir, "post-Lambda-and-true"), 
           recursive = TRUE)
dir.create(file.path(result_dir, "SDM-est-coherence-and-phase"), 
           recursive = TRUE)

# save the parameters file in the results directory, for later reference
file.copy("testing/time-series-from-SDM/model-simulation-parameters.R", 
    file.path(result_dir, "model-simulation-parameters.R"))

num_samp <- gibbsIts / t_thin

# generate true parameters ------------------------------------------------

set.seed(parseed)

# sigmak02, scale parameter
sigmak02 <- rgamma(K, 1, 1)

fkTR <- array(NA, c(P, P, K, Tt))
U_kl0 <- array(NA, c(P, d, K, Tt))
Lambdakl0 <- array(NA, c(d, K, Tt))

if (TS_par_gen_method == "truncVAR") {
    source("testing/time-series-from-SDM/sim-setup-generate-truncVAR-pars.R")
} else if (TS_par_gen_method == "smoothly-similar-Ukl") { 
    source("testing/time-series-from-SDM/sim-setup-generate-smoothly-similar-Ukl-pars.R")
} else if (TS_par_gen_method == "baseVAR") {
    source("testing/time-series-from-SDM/sim-setup-generate-baseVAR-pars.R")
}

# calculate a reasonable Sigmal0 parameter
# # Sigmal0, prior parameter to Ukl0
Sigmal0 <- array(NA, c(P, P, num_freqs))
invSigmal0 <- array(NA, c(P, P, num_freqs))

d_Ukl0_to_avg_Uk0 <- array(NA, c(K, num_freqs))

for (l in 1:num_freqs) {
    avg_Uk0 <- apply(U_kl0[, , , l], c(1, 2), mean)
    avg_Uk0 <- avg_Uk0 %*%
        solve( EigenR::Eigen_sqrt( t(Conj(avg_Uk0)) %*% avg_Uk0 ) )
    
    for (k in 1:K) {
        d_Ukl0_to_avg_Uk0[k, l] <- evec_Frob_stat(avg_Uk0, U_kl0[, , k, l])
    }
    
    avg_Uk0_perp <- (qr(avg_Uk0, complete = TRUE) |>
                         qr.Q(complete = TRUE))[, (d+1):P]

    V_Sigma0 <- cbind(avg_Uk0, avg_Uk0_perp)
    Sigmal0[, , l] <- V_Sigma0 %*% diag(P*(P:1) / sum(P:1)) %*% t(Conj(V_Sigma0))
    invSigmal0[, , l] <- solve(Sigmal0[, , l])
}

# check Lambdakl0 ---------------------------------------------------------

k <- 1
plot(Lambdakl0[1, k, 1:num_freqs], type = "l", ylab = "lambda",
     main = paste0("True Lambda, k = ", k), 
     ylim = c(0, max(Lambdakl0[, k, ])))
lines(Lambdakl0[2, k, 1:num_freqs], col = 2)
save_plot_pdf(file.path(result_dir, "trueLambda_1.pdf"))

k <- 2
plot(Lambdakl0[1, k, 1:num_freqs], type = "l", ylab = "lambda",
     main = paste0("True Lambda, k = ", k), 
     ylim = c(0, max(Lambdakl0[, k, ])))
lines(Lambdakl0[2, k, 1:num_freqs], col = 2)
save_plot_pdf(file.path(result_dir, "trueLambda_2.pdf"))

# check Ukl0 distances to means -------------------------------------------

avg_d_Ukl0 <- colMeans(d_Ukl0_to_avg_Uk0)
sd_d_Ukl0 <- apply(d_Ukl0_to_avg_Uk0, 2, sd)

plot(avg_d_Ukl0, type = "l", 
     main = "avg. of axis Frobenius distance of Ukl0 to avg. Ukl0")
save_plot_pdf(file.path(result_dir, "Ukl0-avg-dist-to-avg.pdf"))

plot(sd_d_Ukl0, type = "l", 
     main = "st. dev. of axis Frobenius distance of Ukl0 to avg. Ukl0")
save_plot_pdf(file.path(result_dir, "Ukl0-SD-dist-to-avg.pdf"))

# calculate Cholesky decompositions ---------------------------------------

Rfs <- array(NA, c(P, P, K, Tt))

for (k in 1:K) {
    for (t in 1:Tt) {
        Rkt <- EigenR::Eigen_sqrt(fkTR[, , k, t])
        diag(Rkt) <- Re(diag(Rkt))
        Rkt <- ( Rkt + t(Conj(Rkt)) )/2
        Rfs[, , k, t] <- Rkt
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
    W <- sapply(1:Tt, function(l) Rfs[, , k, l] %*% Zs[, k, l])
    
    phase_W <- exp(2i * pi * (1:Tt) / Tt)
    phase_t <- exp(2i * pi * (0:(Tt-1)) / Tt)
    
    for (j in 1:P) {
        Yts[j, k, ] <- fft(W[j, ] * phase_W, inverse = TRUE) * phase_t
    }
}

# verify the imaginary parts are small (i.e. to numerical precision)
max(abs(Im(Yts)))

# set Y to the real components, only
Yts <- Re(Yts)

plot(Yts[1, 1, ], type = "l", main = "Observed time series 1",
     ylab = "Y", xlab = "t")
lines(Yts[2, 1, ], col = 2)
lines(Yts[3, 1, ], col = 3)
lines(Yts[4, 1, ], col = 4)
save_plot_pdf(file.path(result_dir, "observed-TS-1.pdf"))

plot(Yts[1, 2, ], type = "l", main = "Observed time series 2",
     ylab = "Y", xlab = "t")
lines(Yts[2, 2, ], col = 2)
lines(Yts[3, 2, ], col = 3)
lines(Yts[4, 2, ], col = 4)
save_plot_pdf(file.path(result_dir, "observed-TS-2.pdf"))

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

for (l in 1:num_freqs){
    data_list <- list()
    
    for (k in 1:K) {
        if (use_true_SDMs) {
            # use the true SDM
            data_list[[k]] <- LL * fkTR[, , k, l]
        } else {
            data_list[[k]] <- LL * SDMests[[k]][, , l]
        }
    }
    
    data_list_w[[l]] <- data_list
}

# compare SDM ests and true SDMs ------------------------------------------

plot(Re(SDMests[[1]][1,1, ]), type = "l", ylim = c(0, max(Re(SDMests[[1]]))), 
     ylab = "spectral density", 
     main = paste0("Diag. entries of multitaper est., k = ", 1))
lines(Re(SDMests[[1]][2,2, ]), col = 2)
lines(Re(SDMests[[1]][3,3, ]), col = 3)
lines(Re(SDMests[[1]][4,4, ]), col = 4)

lines(Re(fkTR[1, 1, 1, 1:num_freqs]), col = 1, lty = 2)
lines(Re(fkTR[2, 2, 1, 1:num_freqs]), col = 2, lty = 2)
lines(Re(fkTR[3, 3, 1, 1:num_freqs]), col = 3, lty = 2)
lines(Re(fkTR[4, 4, 1, 1:num_freqs]), col = 4, lty = 2)

save_plot_pdf(file.path(result_dir, "SDM-est-and-true-1.pdf"))

plot(Re(SDMests[[2]][1,1, ]), type = "l", ylim = c(0, max(Re(SDMests[[2]]))),
     ylab = "spectral density", 
     main = paste0("Diag. entries of multitaper est., k = ", 2))
lines(Re(SDMests[[2]][2,2, ]), col = 2)
lines(Re(SDMests[[2]][3,3, ]), col = 3)
lines(Re(SDMests[[2]][4,4, ]), col = 4)

lines(Re(fkTR[1, 1, 2, 1:num_freqs]), col = 1, lty = 2)
lines(Re(fkTR[2, 2, 2, 1:num_freqs]), col = 2, lty = 2)
lines(Re(fkTR[3, 3, 2, 1:num_freqs]), col = 3, lty = 2)
lines(Re(fkTR[4, 4, 2, 1:num_freqs]), col = 4, lty = 2)

save_plot_pdf(file.path(result_dir, "SDM-est-and-true-2.pdf"))

# initialize arrays -------------------------------------------------------

U_kls_all <- array(NA, c(P, d, K, num_freqs, num_samp))
U_kls_all[, , , , 1] <- U_kl0[, , , 1:num_freqs]

Lambdak_l_s <- array(NA, c(d, K, num_freqs, num_samp))
Lambdak_l_s[, , , 1] <- Lambdakl0[, , 1:num_freqs]

taujkl2_s <- array(NA, c(d, K, num_freqs-2, num_samp))
taujkl2_s[, , , 1] <- 1/rgamma(d*k*(num_freqs-2), tau2_a, rate = tau2_b)

zetajk2_s <- array(NA, c(d, K, num_samp))
zetajk2_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)

sigmak2_s <- array(NA, c(K, num_samp))
sigmak2_s[, 1] <- sigmak02

Sigmal_s <- array(NA_complex_, c(P, P, num_freqs, num_samp))

if (use_Id_Sigmal_init) {
    Sigmal_s[, , , 1] <- array(diag(P), c(P, P, num_freqs))
    result_Sigmals <- array(diag(P), c(P, P, num_freqs))    
} else {
    Sigmal_s[, , , 1] <- Sigmal0
    result_Sigmals <- Sigmal0
}

result_invSigmals <- array(NA_complex_, c(P, P, num_freqs))
for (l in 1:num_freqs) {
    result_invSigmals[, , l] <- solve(result_Sigmals[, , l])
}

accCount_s <- array(NA, c(K, num_freqs, gibbsIts))
accCount_s[, , 1] <- TRUE

accCount_Sigma_s <- array(NA, c(num_freqs, gibbsIts))
accCount_Sigma_s[, 1] <- TRUE

# temporary matrices with the first samples
U_kls <- U_kls_all[, , , , 1]
result_Lambdas <- Lambdak_l_s[, , , 1]
result_taujkl2 <- taujkl2_s[, , , 1]
result_zetajk2 <- zetajk2_s[, , 1]
result_sigmak2 <- sigmak2_s[, 1]
# result_Sigmals is set above

# parallel ----------------------------------------------------------------

if (useMclapply) {
    cluster <- makeCluster(n_cores)
    registerDoParallel(cluster)
}

# function for mclapply ---------------------------------------------------

ksampler <- function(k) {
    thisLambda <- result_Lambdas[, k, ]
    thisUkl <- U_kls[, , k, ]
    thisAccCount <- rep(TRUE, num_freqs)

    thisTaujkl2 <- result_taujkl2[, k, ]
    thisZetajk2 <- result_zetajk2[, k]
    thisSigmak2 <- result_sigmak2[k]

    ##### Lambda sampling
    
    ### Lambda: 1st frequency
    LSkw <- data_list_w[[1]][[k]]
    Ukl <- thisUkl[, , 1]
    
    # for each Lambdak value, j
    for (j in sample(d)) {
        # xi_jk has a truncated Gamma(nk, tjk) distribution (shape, rate), 
        # and tjk is a temporary parameter defined as follows.
        # (take the Real part, since the quadratic form should be real, 
        # but may have a small complex part due to numerical issues.)
        tjk <- Re( t(Conj(Ukl[, j])) %*% LSkw %*% Ukl[, j] ) / 
            thisSigmak2
        
        ##### sample from a truncated Gamma distribution
        # lower and upper bounds for truncating the Gamma distribution
        if (j == 1) lb <- 0 else lb <- 1 / (thisLambda[j-1, 1] + 1)
        if (j == d) ub <- 1 else ub <- 1 / (thisLambda[j+1, 1] + 1)
        
        # by slice sampling, instead
        xi_jk <- uni.slice(1/(1 + thisLambda[j, 1]),
                           dgamma,
                           w = 1, m = Inf, lower = lb, upper = ub,
                           shape = LL+1, scale = 1/tjk,
                           log = TRUE)
        
        # convert xi to Lambda
        thisLambda[j, 1] <- 1/xi_jk - 1
    } # end of l = 1 Lambda sampling
    
    ### Lambda: frequencies in 2 to num_freqs-1
    if (Lambda_prior %in% c("1RW", "2RWPN")) {
        for ( l in sample(2:(num_freqs-1)) ) {
            Ukw <- thisUkl[, , l]
            Dkw1 <- data_list_w[[l]][[k]]
            
            for (j in sample(d)) {
                ajk <- t(Conj(Ukw[, j])) %*% Dkw1 %*% Ukw[, j] / 
                    thisSigmak2
                # should be strictly real
                ajk <- Re(ajk)
                
                if (Lambda_prior == "2RWPN") {
                    # 2nd order RW, previous and next neighbors
                    lp <- thisLambda[j, l-1]
                    ln <- thisLambda[j, l+1]
                    mu <- (lp + ln)/2
                }
                
                else if (Lambda_prior == "1RW") {
                    # 1st order RW
                    mu <- thisLambda[j, l-1]
                }
                
                # for single smoothing parameter
                varpar_jkl <- thisTaujkl2[j, l-1]
                
                # by slice sampling, instead
                newdraw <- uni.slice(thisLambda[j, l],
                                     unnorm_logPDF,
                                     w = 1, m = Inf, lower = 0, 
                                     upper = Inf, 
                                     tau2 = varpar_jkl,
                                     ajk = ajk, mu = mu, N = LL, 
                                     logscale = TRUE)    
                # result_Lambdas[j, k, l] <- newdraw
                thisLambda[j, l] <- newdraw
                
            } # end of sampling j in 1:d
        } # end of sampling over frequencies
    } # end of the conditional to check the type of prior for Lambda
    
    ### Lambda: last frequency
    Uw1 <- thisUkl[, , num_freqs]
    Dkw1 <- data_list_w[[num_freqs]][[k]]
    
    for (j in sample(d)) {
        ajk <- t(Conj(Uw1[, j])) %*% Dkw1 %*% Uw1[, j] / 
            thisSigmak2
        # should be strictly real - this is done in the existing Lambda
        # FCD sampler.
        ajk <- Re(ajk)
        
        # 1st order random walk
        mu <- thisLambda[j, num_freqs - 1]
        
        # by slice sampling, instead
        newdraw <- uni.slice(thisLambda[j, num_freqs], 
                             unnorm_logPDF,
                             w = 1, m = Inf, lower = 0, upper = Inf, 
                             tau2 = thisZetajk2[j],
                             ajk = ajk, mu = mu, N = LL, 
                             logscale = TRUE)    
        thisLambda[j, num_freqs] <- newdraw
        
    } # end of sampling j in 1:d
    
    ##### end of Lambda sampling
    
    ##### taujkl2 and zetajk2 sampling
    
    # sample the separate taujk2 and zetajk2 parameters
    for (j in 1:d) {
        if (Lambda_prior == "2RWPN") {
            # 2nd order RW, next and previous neighbors
            ### BEGIN STANDARD SAMPLING
            sumalljk <- (
                (thisLambda[j, 2:(num_freqs-1)] -
                     .5*thisLambda[j, 1:(num_freqs-2)] -
                     .5*thisLambda[j, 3:(num_freqs)])^2)
            
            thisTaujkl2[j, ] <-
                1/rgamma(length(sumalljk), 
                         tau2_a + .5,
                         rate = tau2_b + .5*sumalljk)
            ### END STANDARD SAMPLING
        } else {
            stop(paste0("No sampler implemented for taujkl when Lambda_prior is ", Lambda_prior))
        }
        
        # zetajk2 sampling
        # based on 1st order random walk for the last Lambda
        sum_zeta <- (thisLambda[j, num_freqs] - 
                         thisLambda[j, num_freqs-1])^2
        thisZetajk2[j] <- 1/rgamma(1, tau2_a + .5,
                                       rate = tau2_b + .5*sum_zeta)
    }
    
    ##### Ukl sampling    
                
    if (sample_true_Ukl0) {
        thisUkl <- U_kl0[, , k, 1:num_freqs]
    } else {
    
    # calculate traces and other terms, for the MH acceptance ratio
    for (l in 1:num_freqs) {
        invSigmas <- result_invSigmals[, , l]
        
        Us <- thisUkl[, , l]
        
        ### MH sampling
        # propose U by small rotation of Us
        Up <- Uk_MH_Cayley(Us, tau_Ukl[k, l], doCayleyZeros, CayleyZeroProb)
        # let Skl be the scaled data matrix (i.e. SDM est. at freq. index l)
        Skl <- data_list_w[[l]][[k]]
        lambdakl <- thisLambda[, l]
        Omegakl <- diag(lambdakl / (1 + lambdakl))
        
        tracep <- Re(sum(Conj(Up) * (Skl %*% Up %*% Omegakl)))
        traces <- Re(sum(Conj(Us) * (Skl %*% Us %*% Omegakl)))
        
        # remaining terms
        logdetp <- log(Re(EigenR::Eigen_det( t(Conj(Up)) %*%
                                                 invSigmas %*% Up )))
        logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*%
                                                 invSigmas %*% Us )))
        
        # calculate acceptance ratio
        logr <- tracep/thisSigmak2 - P*logdetp - traces/thisSigmak2 + P*logdets
        
        # accept or reject
        if (log(runif(1)) <= logr) {
            thisUkl[, , l] <- Up
        } else {
            thisUkl[, , l] <- Us
            thisAccCount[l] <- FALSE
        }
                        
    } # end of sampling over frequencies
    
    } # end of conditional to sample random or true values
    
    ### sigmakl2 sampling
    par1 <- num_freqs*P*LL
    par2 <- 0
    
    for (l in 1:num_freqs) {
        LSkw <- data_list_w[[l]][[k]]

        Ukl <- thisUkl[, , l]
        Lambdakl <- thisLambda[, l]

        mat1 <- diag(P) - Ukl %*% diag(1/(1/Lambdakl + 1)) %*% t(Conj(Ukl))
        par2 <- par2 + Re(sum(t(mat1) * LSkw))
    }
    thisSigmak2 <- 1/rgamma(1, par1, rate = par2)
    
    ### end of sigmakl2 sampling
    
    # return this value for the dopar for loop
    list(Lambdakl = thisLambda, 
         taujkl2 = thisTaujkl2, zetajk2 = thisZetajk2,
         Ukl = thisUkl, accCount = thisAccCount,
         sigmak2 = thisSigmak2)
} # end of ksampler function

# sampling ----------------------------------------------------------------
{
    starttime <- Sys.time()
    for (s in 2:gibbsIts) {
        if (s %% gibbsPrint == 0) cat(paste0(s, ": ", Sys.time(), "\n"))
        
        # before sampling
        accCount <- array(TRUE, c(K, num_freqs))
        
        ##### Parallel sampling
        if (useMclapply) {
            ksamples <- mclapply(1:K, ksampler, mc.cores = n_cores)
        } else {
            ksamples <- lapply(1:K, ksampler)
        }
        
        # save sampled values from the dopar list into iteration-level arrays
        for (k in 1:K) {
            if (inherits(ksamples[[k]], "try-error") || 
                is.character(ksamples[[k]])) {
                stop(paste0("ksampler failed for k=", k, ": ", ksamples[[k]]))
            }
            
            U_kls[, , k, ] <- ksamples[[k]]$Ukl
            accCount[k, ] <- ksamples[[k]]$accCount
            result_Lambdas[, k, ] <- ksamples[[k]]$Lambdakl
            
            result_taujkl2[, k, ] <- ksamples[[k]]$taujkl2
            result_zetajk2[, k] <- ksamples[[k]]$zetajk2
            result_sigmak2[k] <- ksamples[[k]]$sigmak2
        }
        
        accCount_s[, , s] <- accCount    
        
        ###
        # Sigmal sampling
        ###
        
        accCount_Sigma <- rep(TRUE, num_freqs)
        
        for (l in 1:num_freqs) {
        
            Sigmals <- result_Sigmals[, , l]
            invSigmals <- result_invSigmals[, , l]
            
            prop_par_s <- (Sigmals + Sigma_add*diag(P)) / (1 + Sigma_add)
            
            Sigmap <- tryCatch(
                rFTCW(prop_par_s, n_Sig[l], P, useEigenR, byCholesky),
                error = function(e) {
                    cat("rFTCW failed at s =", s, ", l =", l, "\n")
                    cat("n_Sig[l] =", n_Sig[l], "\n")
                    cat("eigenvalues of prop_par_s:\n")
                    print(eigen(prop_par_s)$values)
                    stop(e)
                }
            )
            
            prop_par_p <- (Sigmap + Sigma_add*diag(P)) / (1 + Sigma_add)
            
            invSigmap <- tryCatch(
                solve(Sigmap),
                error = function(e) {
                    cat("Singularity at l =", l, ", n_Sig =", n_Sig[l], "\n")
                    # eigenvalues of proposed Sigma
                    print(eigen(Sigmap)$values)
                    # eigenvalues of the previous Sigma
                    print(eigen(result_Sigmals[, , l])$values)
                    stop(e)
                    })

            # sumlog_p
            sumlog_p <- 0
            sumlog_s <- 0
            for (k in 1:K) {
                sumlog_p <- sumlog_p + 
                    logdet(t(Conj(U_kls[, , k, l])) %*% invSigmap %*% 
                                U_kls[, , k, l])
                sumlog_s <- sumlog_s + 
                    logdet(t(Conj(U_kls[, , k, l])) %*% invSigmals %*% 
                               U_kls[, , k, l])
            }
                            
            logdens_num <- -d*K* logdet(Sigmap) - P * sumlog_p - 
                n_Sig[l] * logdet(prop_par_p) + # Sigmap is the parameter
                (n_Sig[l] - P) * logdet(Sigmals) - 
                P*n_Sig[l] * log( Re(sum( t(solve(prop_par_p)) * Sigmals )) ) # Sigmap is the parameter
            
            logdens_den <- -d*K* logdet(Sigmals) - P * sumlog_s - 
                n_Sig[l] * logdet(prop_par_s) + # Sigmals is the parameter
                (n_Sig[l] - P) * logdet(Sigmap) - 
                P*n_Sig[l] * log( Re(sum( t(solve(prop_par_s)) * Sigmap )) ) # Sigmals is the parameter
            
            logr <- Re(logdens_num - logdens_den)
            
            if (log(runif(1)) <= logr) {
                result_Sigmals[, , l] <- Sigmap
                result_invSigmals[, , l] <- invSigmap

            } else {
                # TODO this is technically redundant, I think
                result_Sigmals[, , l] <- Sigmals
                result_invSigmals[, , l] <- invSigmals
                accCount_Sigma[l] <- FALSE

            }
        } # end of frequencies for Sigmal sampling
        
        accCount_Sigma_s[, s] <- accCount_Sigma
        
        ### end of Sigmal sampling
        
        ### save thinned samples
        if ((s-1) %% t_thin == 0) {    
            # save this iteration into sampler-level arrays
            s_thin <- ceiling(s / t_thin)
            U_kls_all[, , , , s_thin] <- U_kls
            Lambdak_l_s[, , , s_thin] <- result_Lambdas
            taujkl2_s[, , , s_thin] <- result_taujkl2
            zetajk2_s[, , s_thin] <- result_zetajk2
            sigmak2_s[, s_thin] <- result_sigmak2
            Sigmal_s[, , , s_thin] <- result_Sigmals
        }
        
        ### do adaptation of the Ukl MH tuning parameter
        if (s %in% tau_s_check) {
            cat(paste0("\n", s, ": Check for tau_Ukl adaptation\n"))
            
            curr_Ukl_acc_rate <-
                apply(accCount_s[, , (s - tau_numin + 1):s], c(1, 2), mean)
            
            if (show_tau_tune_summ) {
                print(which(curr_Ukl_acc_rate == 0))
            }
            
            tau_Ukl[curr_Ukl_acc_rate >= .45] <-
                tau_Ukl[curr_Ukl_acc_rate >= .45] * 2
            tau_Ukl[curr_Ukl_acc_rate <= .15] <-
                tau_Ukl[curr_Ukl_acc_rate <= .15] / 4
            
            if (show_tau_tune_summ) {
                apply(accCount_s[, , (s - tau_numin + 1):s], c(1, 2), mean) |>
                    t() |>
                    summary() |>
                    print()
            }
            
            ### Sigmal adaptation
            curr_Sigmal_acc_rate <- 
                rowMeans(accCount_Sigma_s[, (s - tau_numin + 1):s])
            
            n_Sig[curr_Sigmal_acc_rate >= .45] <- 
                pmax(round(n_Sig[curr_Sigmal_acc_rate >= .45] / 4), P)
            n_Sig[curr_Sigmal_acc_rate <= .15] <- 
                n_Sig[curr_Sigmal_acc_rate <= .15] * 2
            
            if (show_n_Sig_summary) {
                print("Summary of n_Sig tuning par.:")
                print(summary(n_Sig))
                print("Summary of recent Sigmal acc. rates:")
                print(summary(curr_Sigmal_acc_rate))
            }
            
        } 
        
    } # end of s, Gibbs sampling
    endtime <- Sys.time()
}

# analysis ----------------------------------------------------------------

source("testing/time-series-from-SDM/post-sampling-analysis.R", 
       print.eval = TRUE)

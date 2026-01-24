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

logdet <- function(X) {
    return( log(EigenR::Eigen_det(X)) )
}

# save_plot_png <- function(plot_path, width = 1600, height = 900, res = 300) {
#     p <- recordPlot()
#     png(plot_path, width = width, height = height, res = res, type = "Xlib")
#     replayPlot(p)
#     dev.off()    
# }
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

# setup -------------------------------------------------------------------

# the parameters file is run above

# Create timestamped results directory
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
result_dir <- file.path("results", paste0("full-SDM-diffSmoothness_", timestamp))
dir.create(result_dir, recursive = TRUE)

# save the parameters file in the results directory, for later reference
file.copy("testing/time-series-from-SDM/model-simulation-parameters.R", 
    file.path(result_dir, "model-simulation-parameters.R"))

# generate true parameters ------------------------------------------------

set.seed(parseed)

VARpars <- array(NA, c(P, P, K))
# for (k in 1:K) {
#     VARpars[, , k] <- generate_VAR1_coef(P, 0.8)
# }
VARpars[, , 1:K] <- generate_VAR1_coef(P, 0.8)
noiseSigma <- generate_AR1_covariance(P, sigma2 = 1, rho = 0.5)

# sigmak02, scale parameter
sigmak02 <- rgamma(K, 1, 1)

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
        
        fkTR[, , k, t] <- sigmak02[k] * ( thisU %*% diag(thisLambda) %*% 
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

for (w in 1:num_freqs){
    data_list <- list()
    
    for (k in 1:K) {
        # use the estimated SDM
        data_list[[k]] <- LL * SDMests[[k]][, , w]
        
        # use the true SDM
        # data_list[[k]] <- LL * fkTR[, , k, w]
        
    }
    
    data_list_w[[w]] <- data_list
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

U_kls_all <- array(NA, c(P, d, K, num_freqs, gibbsIts))
U_kls_all[, , , , 1] <- U_kl0[, , , 1:num_freqs]

Lambdak_w_s <- array(NA, c(d, K, num_freqs, gibbsIts))
Lambdak_w_s[, , , 1] <- Lambdakl0[, , 1:num_freqs]

#taujk2_s <- array(NA, c(d, K, gibbsIts))
#taujk2_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)
taujkl2_s <- array(NA, c(d, K, num_freqs-2, gibbsIts))
taujkl2_s[, , , 1] <- 1/rgamma(d*k*(num_freqs-2), tau2_a, rate = tau2_b)

# taujk2B_s <- array(NA, c(d, K, gibbsIts))
# taujk2B_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)

zetajk2_s <- array(NA, c(d, K, gibbsIts))
zetajk2_s[, , 1] <- 1/rgamma(d*K, tau2_a, rate = tau2_b)

sigmak2_s <- array(NA, c(K, gibbsIts))
sigmak2_s[, 1] <- sigmak02

Sigmal_s <- array(NA, c(P, P, num_freqs, gibbsIts))
Sigmal_s[, , , 1] <- Sigmal0
result_Sigmals <- Sigmal0
result_invSigmals <- array(NA, c(P, P, num_freqs))
for (l in 1:num_freqs) {
    result_invSigmals[, , l] <- solve(result_Sigmals[, , l])
}

accCount_s <- array(NA, c(K, num_freqs, gibbsIts))
accCount_s[, , 1] <- TRUE

accCount_Sigma_s <- array(NA, c(num_freqs, gibbsIts))
accCount_Sigma_s[, 1] <- TRUE

# parallel ----------------------------------------------------------------

# library(foreach)
# library(doParallel)
# 
# n_cores <- 12
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
                tjk <- Re( t(Conj(Ukl[, j])) %*% LSkw %*% Ukl[, j] ) / 
                    #sigmakl2_s[k, 1, s-1]
                    sigmak2_s[k, s-1]
                
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
        
        if (Lambda_prior %in% c("1RW", "2RWPN")) {
            
            for ( w in sample(2:(num_freqs-1)) ) {
                for (k in 1:K) {
                    Ukw <- U_kls[, , k, w]
                    Dkw1 <- data_list_w[[w]][[k]]
                    
                    for (j in sample(d)) {
                        ajk <- t(Conj(Ukw[, j])) %*% Dkw1 %*% Ukw[, j] / 
                            # sigmakl2_s[k, w, s-1]
                            sigmak2_s[k, s-1]
                        # should be strictly real
                        ajk <- Re(ajk)
                        
                        if (Lambda_prior == "2RWPN") {
                            # 2nd order RW, previous and next neighbors
                            lp <- result_Lambdas[j, k, w-1]
                            ln <- result_Lambdas[j, k, w+1]
                            mu <- (lp + ln)/2
                        }
                        
                        # SPECIAL for different smoothing parameters
                        # if (w < 300) {
                        #     varpar_jkl <- taujk2_s[j, k, s-1]
                        # } else {
                        #     varpar_jkl <- taujk2B_s[j, k, s-1]
                        # }
                        
                        # 2nd order RW, previous two
                        # lm1 <- result_Lambdas[j, k, w-1]
                        # lm2 <- result_Lambdas[j, k, w-2]
                        # mu <- 2*lm1 - lm2
                        
                        else if (Lambda_prior == "1RW") {
                        # 1st order RW
                            mu <- result_Lambdas[j, k, w-1]
                        }
                        
                        
                        # for single smoothing parameter
                        # varpar_jkl <- taujk2_s[j, k, s-1]
                        varpar_jkl <- taujkl2_s[j, k, w-1, s-1]
                        
                        # by slice sampling, instead
                        newdraw <- uni.slice(result_Lambdas[j, k, w], unnorm_logPDF,
                                             w = 1, m = Inf, lower = 0, upper = Inf, 
                                             #tau2 = taujk2_s[j, k, s-1], 
                                             tau2 = varpar_jkl,
                                             ajk = ajk, mu = mu, N = LL, 
                                             logscale = TRUE)    
                        result_Lambdas[j, k, w] <- newdraw
                        
                    } # end of sampling j in 1:d
                } # end of sampling k in 1:K
            } # end of sampling over frequencies
        }
        
        ##### for the last frequency, l = num_freq
        # 1st order RW
        for (k in 1:K) {
            Uw1 <- U_kls[, , k, num_freqs]
            Dkw1 <- data_list_w[[num_freqs]][[k]]
            
            for (j in sample(d)) {
                ajk <- t(Conj(Uw1[, j])) %*% Dkw1 %*% Uw1[, j] / 
                    # sigmakl2_s[k, num_freqs, s-1]
                    sigmak2_s[k, s-1]
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
                
                if (Lambda_prior == "2RWPN") {
                    # 2nd order RW, next and previous neighbors
                    ### BEGIN STANDARD SAMPLING
                    sumalljk <- (
                        (result_Lambdas[j, k, 2:(num_freqs-1)] -
                             .5*result_Lambdas[j, k, 1:(num_freqs-2)] -
                             .5*result_Lambdas[j, k, 3:(num_freqs)])^2)
    
                    #taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5,
                    #                              rate = tau2_b + .5*sumalljk)
                    taujkl2_s[j, k, , s] <- 
                        1/rgamma(length(sumalljk), 
                                 tau2_a + .5,
                                 rate = tau2_b + .5*sumalljk)
                    ### END STANDARD SAMPLING
                }
                
                # 2nd order RW, next and previous; DIFFERENT FREQUENCIES
                # sumalljk <- sum(
                #     (result_Lambdas[j, k, 2:(299)] -
                #          .5*result_Lambdas[j, k, 1:(299-1)] -
                #          .5*result_Lambdas[j, k, 3:(300)])^2)
                # 
                # sumalljkB <- sum(
                #     (result_Lambdas[j, k, 300:(num_freqs-1)] -
                #          .5*result_Lambdas[j, k, 299:(num_freqs-2)] -
                #          .5*result_Lambdas[j, k, 301:(num_freqs)])^2)
                # 
                # # this is for portion A
                # taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(298),
                #                               rate = tau2_b + .5*sumalljk)
                # # taujk2_s[j, k, s] <- 10000
                # 
                # # this is for portion B
                # taujk2B_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs - 300),
                #                               rate = tau2_b + .5*sumalljkB)
                ### END different frequencies
                
                # 2nd order RW, previous two
                # sumalljk <- sum(
                #     (result_Lambdas[j, k, 3:num_freqs] -
                #          result_Lambdas[j, k, 2:(num_freqs-1)] +
                #          .5*result_Lambdas[j, k, 1:(num_freqs-2)])^2)
                # 
                # taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-2),
                #                               rate = tau2_b + .5*sumalljk)
                
                #else if (Lambda_prior == "1RW") {
                #    # 1st order RW
                #    sumalljk <- sum(
                #        (result_Lambdas[j, k, 2:(num_freqs-1)] -
                #            result_Lambdas[j, k, 1:(num_freqs-2)])^2)
                #    taujk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5*(num_freqs-2),
                #                                rate = tau2_b + .5*sumalljk)
                #}
                
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
                # sigma_k2 <- sigmakl2_s[k, l, s-1]
                sigma_k2 <- sigmak2_s[k, s-1]
                # Sigmas <- Sigmal0[, , l]
                # invSigmas <- invSigmal0[, , l]
                invSigmas <- result_invSigmals[, , l]
                
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
        
        accCount_Sigma <- rep(TRUE, num_freqs)
        
        for (l in 1:num_freqs) {
        # this_Sigma_result <- foreach(l = 1:num_freqs) %dopar% {
            Sigmap <- rFTCW(result_Sigmals[, , l], n_Sig[l], P, TRUE, TRUE)
            # invSigmap <- solve(Sigmap)
            
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
            
            Sigmals <- result_Sigmals[, , l]
            invSigmals <- result_invSigmals[, , l]
            
            # sumlog_p
            sumlog_p <- 0
            sumlog_s <- 0
            for (k in 1:K) {
                sumlog_p <- sumlog_p + 
                    logdet(t(Conj(newU_kls[, , k, l])) %*% invSigmap %*% 
                                newU_kls[, , k, l])
                sumlog_s <- sumlog_s + 
                    logdet(t(Conj(newU_kls[, , k, l])) %*% invSigmals %*% 
                               newU_kls[, , k, l])
            }
            
            logdens_num <- -d*K* logdet(Sigmap) - P * sumlog_p - 
                n_Sig[l] * logdet(Sigmap) +
                (n_Sig[l] - P) * logdet(Sigmals) - 
                # P*n_Sig[l] * log( Re(sum(diag( invSigmap %*% Sigmals))) ) 
                P*n_Sig[l] * log( Re(sum(t(invSigmap) * Sigmals)) ) 
            
            logdens_den <- -d*K* logdet(Sigmals) - P * sumlog_s - 
                n_Sig[l] * logdet(Sigmals) +
                (n_Sig[l] - P) * logdet(Sigmap) - 
                # P*n_Sig[l] * log( Re(sum(diag( invSigmals %*% Sigmap ))) ) 
                P*n_Sig[l] * log( Re(sum(t(invSigmals) * Sigmap )) ) 
            
            logr <- Re(logdens_num - logdens_den)
            
            if (log(runif(1)) <= logr) {
                result_Sigmals[, , l] <- Sigmap
                result_invSigmals[, , l] <- invSigmap

                # Sigmals <- Sigmap
                # invSigmals <- invSigmap
                # acc_this_Sigma <- TRUE
            } else {
                # TODO this is technically redundant, I think
                result_Sigmals[, , l] <- Sigmals
                result_invSigmals[, , l] <- invSigmals
                accCount_Sigma[l] <- FALSE

                # acc_this_Sigma <- FALSE
            }
            
            # list(Sigmals, invSigmals, acc_this_Sigma)
            # list(Sigmal = Sigmals, invSigmal = invSigmals, acc = acc_this_Sigma)
        } # end of frequencies for Sigmal sampling
        
        # for (l in 1:num_freqs) {
        #     result_Sigmals[, , l] <- this_Sigma_result[[l]][[1]]
        #     result_invSigmals[, , l] <- this_Sigma_result[[l]][[2]]
        #     accCount_Sigma[l] <- this_Sigma_result[[l]][[3]]
        # }
        
        # Vectorized extraction (faster than for loop)
        # result_Sigmals <- simplify2array(lapply(this_Sigma_result, `[[`, "Sigmal"))
        # result_invSigmals <- simplify2array(lapply(this_Sigma_result, `[[`, "invSigmal"))
        # accCount_Sigma <- vapply(this_Sigma_result, `[[`, logical(1), "acc")
        
        Sigmal_s[, , , s] <- result_Sigmals
        accCount_Sigma_s[, s] <- accCount_Sigma
        
        ### end of Sigmal sampling
        
        ### sigmakl2 sampling
        
        for (k in 1:K) {
            par1 <- num_freqs*P*LL
            par2 <- 0
            
            for (l in 1:num_freqs) {
                # TODO rest of the calculations
                # zetajk2_s[j, k, s] <- 1/rgamma(1, tau2_a + .5,
                #                                rate = tau2_b + .5*sum_zeta)
                
                # Re(sum(t(invSigmap) * Sigmals))
                
                LSkw <- data_list_w[[l]][[k]]
                Ukl <- newU_kls[, , k, l]
                Lambdakl <- diag(result_Lambdas[, k, l])
                
                mat1 <- diag(P) - Ukl %*% 
                    diag(1/(1/result_Lambdas[, k, l] + 1)) %*% t(Conj(Ukl))
                # par2 <- 
                # t(solve(newU_kls[, , k, l] %*% diag(result_Lambdas[, k, l]) %*% 
                    # t(Conj(newU_kls[, , k, l])) + diag(P))) * data
                # par2 <- Re(sum(
                    # solve(Ukl %*% Lambdakl %*% t(Conj(Ukl)) + diag(P)) * LSkw))
                
                par2 <- par2 + Re(sum(t(mat1) * LSkw))
                
                # sigmakl2_s[k, l, s] <- 1/rgamma(1, par1, rate = par2)
            }
            
            sigmak2_s[k, s] <- 1/rgamma(1, par1, rate = par2)
        }
        
        ### end of sigmakl2 sampling
        
        ### do adaptation of the Ukl MH tuning parameter
        if (s %in% tau_s_check) {
            cat(paste0("\n", s, ": Check for tau_Ukl adaptation\n"))
            
            curr_Ukl_acc_rate <-
                apply(accCount_s[, , (s - tau_numin + 1):s], c(1, 2), mean)
            
            # print(tau_Ukl[curr_Ukl_acc_rate == 0])
            if (show_tau_tune_summ) {
                print(which(curr_Ukl_acc_rate == 0))
            }
            
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
            
            ### Sigmal adaptation
            curr_Sigmal_acc_rate <- 
                rowMeans(accCount_Sigma_s[, (s - tau_numin + 1):s])
            
            n_Sig[curr_Sigmal_acc_rate >= .45] <- 
                pmax(round(n_Sig[curr_Sigmal_acc_rate >= .45] / 4), P+1)
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

# evaluate ----------------------------------------------------------------
# stopCluster(cl = cluster)

print(endtime - starttime)

# check the acceptance rates
gibbsPostBurn <- round(gibbsIts * burnin):gibbsIts

apply(accCount_s[, , 2:gibbsIts], c(1, 2), mean) |> t() |> summary()
apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) |> t() |> summary()

apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean) |>
    as.vector() |> 
    quantile(probs = c(0, 0.025, .25, .5, .75, .975, 1))

# check Sigmal acceptance rates
rowMeans(accCount_Sigma_s) |> 
    quantile(probs = c(0, 0.025, .25, .5, .75, .975, 1))
rowMeans(accCount_Sigma_s[, gibbsPostBurn]) |> 
    quantile(probs = c(0, 0.025, .25, .5, .75, .975, 1))

# evaluate Lambdas --------------------------------------------------------

upper_q <- apply(Lambdak_w_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.975)
lower_q <- apply(Lambdak_w_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.025)
# Lambda_means <- apply(Lambdak_w_s, c(1, 2, 3), mean)
Lambda_means <- apply(Lambdak_w_s[, , , gibbsPostBurn], c(1, 2, 3), mean)

k <- 1

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, 
     ylim = c(0, max(Lambdakl0[1, k, 1:num_freqs])),
     main = paste0("post. mean and true Lambda, k = ", k),
     ylab = "Lambda")
lines(lower_q[1, k, ], type = "l", lty = 2)
lines(Lambda_means[1, k, ])
lines(Lambdakl0[1, k, 1:num_freqs], lty = 3)

lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
lines(Lambda_means[2, k, ], , col = "red")
lines(Lambdakl0[2, k, 1:num_freqs], lty = 3, , col = "red")

save_plot_pdf(file.path(result_dir, "post-Lambda-and-true-1.pdf"))

k <- 2

# compare quantiles and means for all frequencies of k = 2
plot(upper_q[1, k, ], type = "l", lty = 2, 
     # xlim = c(400, 510),
     ylim = c(0, max(Lambdakl0[1, k, 1:num_freqs])),
     main = paste0("post. mean and true Lambda, k = ", k),
     ylab = "Lambda")
lines(lower_q[1, k, ], type = "l", lty = 2)
lines(Lambda_means[1, k, ])
lines(Lambdakl0[1, k, 1:num_freqs], lty = 3)

lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
lines(Lambda_means[2, k, ], , col = "red")
lines(Lambdakl0[2, k, 1:num_freqs], lty = 3, , col = "red")

save_plot_pdf(file.path(result_dir, "post-Lambda-and-true-2.pdf"))

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

k <- 2
l <- 396

apply(accCount_s[, , gibbsPostBurn], c(1, 2), mean)[k, l]

avgUkl <- apply(U_kls_all[, , k, l, gibbsPostBurn],
                c(1, 2), mean)
avgUkl <- avgUkl %*% solve( EigenR::Eigen_sqrt( t(Conj(avgUkl)) %*% avgUkl ) )
print(fast_evec_Frob_stat(U_kl0[, , k, l], avgUkl))

d_to_avgUkl <- rep(NA, gibbsIts)
for (s in 1:gibbsIts) {
    d_to_avgUkl[s] <- fast_evec_Frob_stat(U_kls_all[, , k, l, s], avgUkl)
    # d_to_avgUkl[s] <- fast_evec_Frob_stat(U_kls_all[, , k, l, s], U_kl0[, , k, l])
}

plot(d_to_avgUkl, type = "l",
     main = paste0("Trace plot, Ukl axis dist. to mean, k = ", k, ", l = ", l),
     ylab = "axis Frobenius dist.")
save_plot_pdf(file.path(result_dir, 
    paste0("post-Ukl-trace-axis-dist-k-", k, "-l-", l, ".pdf")))

evec_Frob_stat(U_kl0[, , k, l], avgUkl, returnMats = TRUE)
evec_Frob_stat(avgUkl, U_kl0[, , k, l], returnMats = TRUE)

# evaluate Sigmal ---------------------------------------------------------

quantile(n_Sig, c(0, .025, .25, .5, .75, .975, 1))

l <- 396

# check Sigmal acceptance rates
mean(accCount_Sigma_s[l, ])
mean(accCount_Sigma_s[l, gibbsPostBurn])

avgSigmal <- apply(Sigmal_s[, , l, gibbsPostBurn], c(1, 2), mean)

d_to_avgSigmal <- rep(NA, gibbsIts)
for (s in 1:gibbsIts) {
    d_to_avgSigmal[s] <- frob_dist(Sigmal_s[, , l, s], avgSigmal)
}

plot(d_to_avgSigmal, type = "l",
     main = paste0("Trace plot, Sigmal dist. to mean, l = ", l),
     ylab = "Frobenius dist.")
save_plot_pdf(file.path(result_dir, 
    paste0("post-Sigmal-trace-l-", l, ".pdf")))

par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
for(j in 1:4) {
    plot(Re(Sigmal_s[j, j, l, ]), type = "l", 
         main = paste0("Trace plot, Sigmal[j,j], l = ", l, ", j = ", j),
         xlab = "", ylab = "")
}
save_plot_pdf(file.path(result_dir, 
    paste0("post-Sigmal-diag-trace-l-", l, "-j-", j, ".pdf")))

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

# smoothing parameter summaries -------------------------------------------

dim(taujkl2_s)
dim(taujkl2_s[1, 1, 1, gibbsPostBurn])
length(taujkl2_s[1, 1, 1, gibbsPostBurn])
length(taujkl2_s[1, 1, round(num_freqs/2), gibbsPostBurn])

plot(taujkl2_s[1, 1, 1, gibbsPostBurn], type = "l")

median(taujkl2_s[1, 1, 1, gibbsPostBurn])

median(taujkl2_s[1, 1, 150, gibbsPostBurn])
median(taujkl2_s[2, 1, 150, gibbsPostBurn])

median(taujkl2_s[1, 1, 256, gibbsPostBurn])
median(taujkl2_s[1, 1, 350, gibbsPostBurn])

# evaluate sigmak2 -------------------------------------------------------

sigmak02[1]
sigmak02[2]

k <- 1
sigmak02[k]
quantile(sigmak2_s[k, gibbsPostBurn], c(.025, .5, .975))

plot(sigmak2_s[k, ], type = "l",
     main = paste0("Trace plot, sigmak2, k = ", k),
     ylab = "sigmak2")
abline(h = sigmak02[k])
save_plot_pdf(file.path(result_dir, paste0("post-sigmak2-k-", k, ".pdf")))

k <- 2
sigmak02[k]
quantile(sigmak2_s[k, gibbsPostBurn], c(.025, .5, .975))

plot(sigmak2_s[k, ], type = "l",
     main = paste0("Trace plot, sigmak2, k = ", k),
     ylab = "sigmak2")
abline(h = sigmak02[k])
save_plot_pdf(file.path(result_dir, paste0("post-sigmak2-k-", k, ".pdf")))

# evaluate full SDM estimates ---------------------------------------------

# sigmakl2 * (Ukl %*% Lambdakl %*% t(Conj(Ukl)) + I_P)
posterior_dists <- array(NA, c(K, num_freqs))
multitaper_dists <- array(NA, c(K, num_freqs))

for (k in 1:K) {
    for (l in 1:num_freqs) {
        thisSDM <- matrix(0 + 0i, P, P)
        for (s in gibbsPostBurn) {
            thisSDM <- thisSDM + sigmak2_s[k, s] * 
                (U_kls_all[, , k, l, s] %*% diag(Lambdak_w_s[, k, l, s])
                 %*% t(Conj(U_kls_all[, , k, l, s])) + diag(P))
        }
        
        thisSDM <- thisSDM / length(gibbsPostBurn)
        
        # compare distances
        posterior_dists[k, l] <- frob_dist(thisSDM, fkTR[, , k, l])
        multitaper_dists[k, l] <- frob_dist(data_list_w[[l]][[k]] / LL, 
                                      fkTR[, , k, l])
    }
}

summary(as.vector(posterior_dists))
summary(as.vector(multitaper_dists))

k <- 1
quantile(as.vector(posterior_dists[k, ]), 
         probs = c(0, .025, .5, .975, 1))
quantile(as.vector(multitaper_dists[k, ]), 
         probs = c(0, .025, .5, .975, 1))

plot(density(as.vector(posterior_dists[k, ]), from = 0), col = 1,
     main = "densities of SDM estimate distances",
     xlab = "Frob. distance", ylab = "density")
lines(density(as.vector(multitaper_dists[k, ]), from = 0), col = 2)
legend(x = "topright", legend = c("posterior", "multitaper"),
       col = c(1, 2), lwd = 2)
save_plot_pdf(file.path(result_dir, 
    paste0("post-SDM-est-dist-density-k-", k, ".pdf")))
    
k <- 2
quantile(as.vector(posterior_dists[k, ]), 
         probs = c(0, .025, .5, .975, 1))
quantile(as.vector(multitaper_dists[k, ]), 
         probs = c(0, .025, .5, .975, 1))

plot(density(as.vector(posterior_dists[k, ]), from = 0), col = 1,
     main = "densities of SDM estimate distances",
     xlab = "Frob. distance", ylab = "density")
lines(density(as.vector(multitaper_dists[k, ]), from = 0), col = 2)
legend(x = "topright", legend = c("posterior", "multitaper"),
       col = c(1, 2), lwd = 2)
save_plot_pdf(file.path(result_dir, 
    paste0("post-SDM-est-dist-density-k-", k, ".pdf")))

k <- 1
l <- 75
thisSDM <- matrix(0 + 0i, P, P)
for (s in gibbsPostBurn) {
    thisSDM <- thisSDM + sigmak2_s[k, s] * 
        (U_kls_all[, , k, l, s] %*% diag(Lambdak_w_s[, k, l, s])
         %*% t(Conj(U_kls_all[, , k, l, s])) + diag(P))
}

thisSDM <- thisSDM / length(gibbsPostBurn)

# true SDM
Re(diag(fkTR[, , k, l]))
# posterior mean
Re(diag(thisSDM))
# multitaper estimate
Re(diag(data_list_w[[l]][[k]]) / LL)

# compare distances
frob_dist(thisSDM, fkTR[, , k, l])
frob_dist(data_list_w[[l]][[k]] / LL, fkTR[, , k, l])


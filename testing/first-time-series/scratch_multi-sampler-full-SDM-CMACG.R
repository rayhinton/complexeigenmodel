# multi-sampler for covariance model with Uk ~ CMACG prior 
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcomplex_wishart.R")

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Lambdak.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_sigmak2.R")

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Uk_CMACG.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Sigma.R")


# functions ---------------------------------------------------------------

frob_dist <- function(A, B, returnDists = FALSE) {
    diffAB <- A - B
    
    if (returnDists) {
        sqDists <- Re(diag( t(Conj(diffAB)) %*% diffAB ))
        Fdist <- Re(sum(sqDists))
        return(list(Fdist = sqrt(Fdist),
                    sqDists = sqDists))
    } else {
        return(norm(diffAB, "F"))
    }
}

fast_evec_Frob_stat <- function(X, Y) {
    k <- ncol(X)
    return(
        sqrt(2*k - 2*sum(Mod( diag(t(Conj(X)) %*% Y)) ))
    )
}

evec_Frob_stat <- function(X, Y, returnDists = FALSE, returnMats = FALSE) {
    Rmat <- diag(complex(modulus = 1, 
                         argument = -Arg(diag(t(Conj(X)) %*% Y))))  
    dist_list <- frob_dist(X, Y %*% Rmat, returnDists = returnDists)
    
    if (returnMats) {
        return(list(dist_obj = dist_list,
                    Xopt = X,
                    Yopt = Y %*% Rmat))
    } else {
        return(dist_list)
    }
}

# need spectral density function for AR(2)
sdf_ar2 <- function(phis, omega, sigma2 = 1) {
    phi1 <- phis[1]
    phi2 <- phis[2]
    
    denom <- 1 + phi1^2 + phi2^2 + 2*phi1*(phi2 - 1) * cos(2*pi*omega) - 
        2*phi2*cos(4*pi*omega)
    
    return(sigma2 / denom)
}

# set up ------------------------------------------------------------------

P <- 4
d <- 2
n_k_approx <- 1024
K <- 2
# n_k_approx <- 5000
L_numtprs <- 32

sampleUk <- TRUE
sampleLambdak <- TRUE
samplesigmak2 <- TRUE
sampleSigma <- TRUE

data_seed <- 10072025
# data_seed <- 12072025
parameter_seed <- 314

randomInit <- TRUE
U_ks_est_init <- TRUE
Lambda_ks_est_init <- TRUE
sigma_k2s_est_init <- TRUE

doLambdaOrdered <- FALSE

gibbsIts <- 1e5
thinBy <- 10
numSave <- ceiling(gibbsIts/thinBy)
burnin <- 0.5
gibbsPrint <- gibbsIts/10

### tuning parameters
# Uk sampler
num_tau_check <- 10

# tau_U <- 0.005
# tau_Uk <- c(.01, .0025, .0025, .01)
tau_Uk <- rep(.05, K) # good starting point for P = 8
# tau_Uk <- rep(.0005, K) # good starting point for P = 64
doCayleyZeros <- FALSE
CayleyZeroProb <- 0.5

# Sigma sampler
n_Sig <- P + 250 # works well for P = 8, d = 4, K = 4
# n_Sig <- P + 5000 # works well for P = 8, d = 4, K = 100
# n_Sig <- P + 2500 # works well for P = 8, d = 4, K = 4
useEigenR <- TRUE
byCholesky <- TRUE

# adapting tuning parameters
burninS <- floor(gibbsIts * burnin)
tau_numin <- floor(burninS / num_tau_check)
tau_s_check <- seq(1+tau_numin, burninS, tau_numin)

# true parameters ---------------------------------------------------------

##### TODO all of this is written assuming K = 1, just as a proof of concept

# frequency to analyze
# naively, choose pi/2, or in terms of non-pi frequencies, 1/4
# omegaish <- 3/32
omegaish <- 3/16
k_freq <- round(n_k_approx * omegaish)  # for example
omega0 <- k_freq / n_k_approx

# AR parameters
# AR parameters
# phis1 has the AR(2) parameters for the 1st factor
# phis2 has the AR(2) parameters for the 2nd factor
# each row is for one time series; they have K rows each

# phis1 <- matrix(c(0.5, -0.5, 1, -0.9), nrow = K, byrow = TRUE)
# phis2 <- matrix(c(-0.5, 0.25, -0.25, -0.75), nrow = K, byrow = TRUE)

# phis1 <- matrix(c(0.5, -0.5, -0.25, -0.75), nrow = K, byrow = TRUE)
# phis2 <- matrix(c(-0.5, 0.25, 1, -0.9), nrow = K, byrow = TRUE)

# these parameters are well separated at omega0 = 3/16
# phis1 <- matrix(c(1, -0.85, .75, -0.65), nrow = K, byrow = TRUE)
# phis2 <- matrix(c(1, -.5, .75, -.25), nrow = K, byrow = TRUE)

# these parameters are also well separated at omega0 = 3/16, but one of them is much closer to 0
phis1 <- matrix(c(1, -0.85, .75, -0.65), nrow = K, byrow = TRUE)
phis2 <- matrix(c(1, -.5, -1, -.05), nrow = K, byrow = TRUE)

# scale
sigma02 <- c(5, 10)

# U, the matrix that creates the observed time series
# U <- qr.Q(qr(matrix(rnorm(P*d), ncol = d)))
# confirm they are orthogonal, the result should be 2x2 I
# crossprod(U)

# generate true parameters and data ---------------------------------------

# n_k <- rep(n_k_approx, K) + sample(-10:10, K, replace = TRUE)
n_k <- rep(n_k_approx, K)
L_k <- round(sqrt(n_k))

Sigma0 <- rFTCW(diag(P), P+1, P)

Lambda_k0 <- array(NA, c(d, K))
sigma_k20 <- rep(NA, K)
U_k0 <- array(NA, c(P, d, K))
data_list <- list()

# parameter and data generation

set.seed(data_seed)
for (k in 1:K) {
    # generate parameters
    Lambda_k0[, k] <- c(sdf_ar2(phis1[k, ], omega0), 
                        sdf_ar2(phis2[k, ], omega0))
    sigma_k20[k] <- sigma02[k]
    U <- qr.Q(qr(matrix(rnorm(P*d), ncol = d)))
    U_k0[, , k] <- U
    
    # generate the d-dim factor time series
    arcomp1 <- list(ar = phis1[k, ])
    arcomp2 <- list(ar = phis2[k, ])
    f1 <- arima.sim(arcomp1, n_k[k])
    f2 <- arima.sim(arcomp2, n_k[k])

    # call X the multivariate time series, which is a 2xT matrix
    X <- ts(cbind(f1, f2), start = start(f1), frequency = frequency(f1))
    # create final time series, Y = sigma * (UX + epsilon), NOT sigma^2!!!
    # generate more unit variance noise, PxT, call it epsilon
    Ymat <- sqrt(sigma02[k]) * (U %*% t(X) + 
                                    matrix(rnorm(P*n_k[k]), ncol = n_k[k]))
    Yts <- ts(t(Ymat), start = start(f1), frequency = frequency(f1))

    # these are the beginning indices of the L intervals, plus a final entry for convenience
    starts <- c(seq(1, n_k[k], by = L_k[k]), n_k[k]+1)
    # calculate the averaged periodogram, which will be the data to the model
    Shatomega <- matrix(0, P, P)

    for (l in 1:L_k[k]) {
        tband <- starts[l] : (starts[l+1] - 1)
        nl <- length(tband)
        ytilde <- Yts[tband, ]
        dlomega <- matrix(0, nrow = nl, ncol = P)
        for (t in 1:nl) {
            dlomega[t, ] <- ytilde[t, ] * exp(-omega0*1i*t)
        }
        dlomega <- colSums(dlomega)

        Shatomega <- Shatomega + 
            1/(2*pi*nl*L_k[k]) * dlomega %*% t(Conj(dlomega))
    }
    
    # data_list[[k]] <- L_k[k]*Shatomega
    
    fXomega0 <- diag(c(sdf_ar2(phis1[k, ], omega0),
                       sdf_ar2(phis2[k, ], omega0)))
    ##### calculate the exact SDM for the P-dim time series
    # need to confirm the SDM relationship for this sort of factor model

    fYomega0 <- sigma02[k] * ( U %*% fXomega0 %*% t(Conj(U)) + diag(P) )

    data_list[[k]] <- L_numtprs*fYomega0
    
}

# estimates of sigmak2 ----------------------------------------------------

eigen(data_list[[1]] / n_k[1])$values
eigen(data_list[[1]])$values
eigen(data_list[[1]])$values/L_numtprs
sigma_k20

# initialize arrays to store posterior samples ----------------------------

accCount_U_k <- array(NA, c(K, gibbsIts))
accCount_Sigma <- rep(NA, gibbsIts)

accCount_U_k[, 1] <- TRUE
accCount_Sigma[1] <- TRUE

U_ks_init <- array(NA, c(P, d, K))
Lambda_ks_init <- array(NA, c(d, K))
sigma_k2s_init <- rep(NA, K)

for (k in 1:K) {
    # Uk initialization
    if (U_ks_est_init) {
        U_ks_init[, , k] <- eigen(data_list[[k]])$vectors[, 1:d]
    } else {
        U_ks_init[, , k] <- runif_stiefel(P, d, 2)
    }
    
    # Lambdak initialization
    if (Lambda_ks_est_init) {
        sigmaK2_hat <- mean(eigen(data_list[[k]]/L_numtprs)$values[(d+1):P])
        Lambda_ks_init[, k] <- 
            eigen(data_list[[k]] / L_numtprs / sigmaK2_hat)$values[1:d] - 1
    } else {
        omegak <- runif(d)
        if (doLambdaOrdered) omegak <- sort(omegak, decreasing = TRUE)
        # omegak <- sort(omegak, decreasing = TRUE)
        Lambda_ks_init[, k] <- omegak / (1 - omegak)
    }
    
    # sigmak2 initialization
    # use sigmaK2_hat calculated above
    if (sigma_k2s_est_init) {
        sigma_k2s_init[k] <- sigmaK2_hat
    } else {
        sigma_k2s_init[k] <- rgamma(1, 1, scale = 100)
    }
}

U_kS <- array(NA, c(P, d, K, numSave))
Lambda_kS <- array(NA, c(d, K, numSave))
sigma_k2S <- array(NA, c(K, numSave))
SigmaS <- array(NA, c(P, P, numSave)) 

if (randomInit) {
    U_kS[, , , 1] <- U_ks_init
    Lambda_kS[, , 1] <- Lambda_ks_init
    # sigma_k2S[, 1] <- rgamma(K, 1, 1/100)
    sigma_k2S[, 1] <- sigma_k2s_init
    SigmaS[, , 1] <- rFTCW(diag(P), P+1, P)
} else {
    U_kS[, , , 1] <- U_k0
    Lambda_kS[, , 1] <- Lambda_k0
    sigma_k2S[, 1] <- sigma_k20
    SigmaS[, , 1] <- Sigma0
}

# TODO I've input L for n_k, since that is the deg. of freedom in time series case
param_list <- list(
    P = P, d = d, K = K, n_k = rep(L_numtprs, K),
    U_ks = U_kS[, , , 1], Lambda_ks = Lambda_kS[, , 1], 
    sigma_k2s = sigma_k2S[, 1], Sigmas = SigmaS[, , 1]
)

# estimates ---------------------------------------------------------------

Lambda_k0
Lambda_kS[, , 1]

sigma_k20
sigma_k2S[, 1]

# sampling ----------------------------------------------------------------

for (s in 2:gibbsIts) {
    if (s %% gibbsPrint == 0) {
        print(paste0("s = ", s, ": ", Sys.time()))
    }
    
    # sample Uk
    if (sampleUk) {
        Uk_out <-
            Uk_CMACG_gibbs_densCovar(data_list, param_list, tau_Uk, 
                                     doCayleyZeros, CayleyZeroProb)
        accCount_U_k[, s] <- Uk_out$accCount_U_k
        param_list$U_ks <- Uk_out$U_ks
    } else {
        param_list$U_ks <- U_k0
    }
    
    # adjust Uk tuning parameter tau during burnin
    if (s %in% tau_s_check & sampleUk) {
        curr_Uk_acc_rate <- rowMeans(accCount_U_k[, (s - tau_numin + 1):s])
        
        # tau_Uk <- ifelse(curr_Uk_acc_rate >= .3, tau_Uk * 10, tau_Uk / 4)
        
        tau_Uk[curr_Uk_acc_rate >= .45] <- tau_Uk[curr_Uk_acc_rate >= .45] * 5
        tau_Uk[curr_Uk_acc_rate <= .15] <- tau_Uk[curr_Uk_acc_rate <= .15] / 10
        
        print(curr_Uk_acc_rate)
        print(tau_Uk)
    }

    # sample Lambdak
    if (sampleLambdak) {
        param_list$Lambda_ks <-
            Lambdak_gibbs_densCovar(data_list, param_list, 
                                    doOrdered = doLambdaOrdered) 
    } else {
        param_list$Lambda_ks <- Lambda_k0
    }
    
    # sample sigmak2
    if (samplesigmak2) {
        param_list$sigma_k2s <- sigmak2_gibbs_densCovar(data_list, param_list)
    } else {
        param_list$sigma_k2s <- sigma_k20
    }
    
    # sample Sigma
    if (sampleSigma) {
        Sigma_out <-
            Sigma_gibbs_densCovar(data_list, param_list, n_Sig,
                                  useEigenR, byCholesky)
        accCount_Sigma[s] <- Sigma_out$accCount_Sigma
        param_list$Sigmas <- Sigma_out$Sigmas
    } else {
        param_list$Sigmas <- Sigma0
    }
    
    if (s %in% tau_s_check & sampleSigma) {
        curr_Sigma_acc_rate <- mean(accCount_Sigma[(s - tau_numin + 1):s])
        if (curr_Sigma_acc_rate >= .45) {
            n_Sig <- P + (n_Sig - P)/5
        } else if ( curr_Sigma_acc_rate <= .15) {
            n_Sig <- P + (n_Sig - P)*10
        }
        print(c(curr_Sigma_acc_rate, n_Sig))
    }
    
    # save thinned samples
    if ((s-1) %% thinBy == 0 ) {
        saveInd <- ceiling(s/thinBy)
        
        U_kS[, , , saveInd] <- param_list$U_ks
        Lambda_kS[, , saveInd] <- param_list$Lambda_ks
        sigma_k2S[, saveInd] <- param_list$sigma_k2s
        SigmaS[, , saveInd] <- param_list$Sigmas
    }
}

# acceptance rates --------------------------------------------------------

gibbsKeepInds <- seq(1, gibbsIts, by = thinBy)
gibbsPostBurn <- gibbsKeepInds > gibbsIts*burnin
# TODO these are the indices to keep after burn in, if burnin < 1
# TODO make a case to handle if burnin > 1, i.e. burnin is explicit
# gibbsKeepInds > gibbsIts*burnin

rowMeans(accCount_U_k, na.rm = TRUE)
mean(accCount_Sigma, na.rm = TRUE)

rowMeans(accCount_U_k[, burninS:gibbsIts], na.rm = TRUE)
tau_Uk

mean(accCount_Sigma[burninS:gibbsIts])
n_Sig

# assess Lambdak ----------------------------------------------------------

Lambda_kS[, , 1]

apply(Lambda_kS[ , , gibbsPostBurn], c(1, 2), mean)
Lambda_k0

Lj <- 1
Lk <- 1

dev.new()
par(mar = c(5.1, 6, 4.1, 2.1))
par(mfrow = c(2, 2))
for (Lk in 1:K) {
    for (Lj in 1:d) {
        plot(Lambda_kS[Lj, Lk, ], type = "l",
             main = bquote(j == .(Lj) ~","~ k == .(Lk)),
             ylab = bquote(lambda[j *","* k]),
             cex.lab = 2,
             cex.main = 2)
        abline(h = Lambda_k0[Lj, Lk])
    }
}

par(mar = c(5.1, 4.1, 4.1, 2.1))
par(mfrow = c(1, 1))

Lambda_k0[Lj, Lk]

quantile(Lambda_kS[Lj, Lk, gibbsPostBurn],
         probs = c(.025, .5, .975))

plot(density(Lambda_kS[Lj, Lk, gibbsPostBurn]))

Lambdak_df <- NULL

for (Lk in 1:K) {
    for (Lj in 1:d) {
        CIjk <- quantile(Lambda_kS[Lj, Lk, gibbsPostBurn],
                         probs = c(.025, 0.5, .975))

        summ_tosave <- c(
            Lj, Lk, Lambda_k0[Lj, Lk], unname(CIjk), 
            mean(Lambda_kS[Lj , Lk, gibbsPostBurn]), Lambda_kS[Lj, Lk, 1],
            unname(CIjk[1] <= Lambda_k0[Lj, Lk] & Lambda_k0[Lj, Lk] <= CIjk[3])
        )
        Lambdak_df <- rbind(Lambdak_df, unname(summ_tosave))
    }
}

colnames(Lambdak_df) <- c("j", "k", "true", ".025", "0.5", ".975", "mean", "start", "in_CI")
View(round(Lambdak_df, 3))

# assess sigmak2 ----------------------------------------------------------

sigma_k2_df <- 
    cbind(1:K,
          unname(sigma_k20),
          t(unname(apply(sigma_k2S[, gibbsPostBurn], 1, quantile, 
                         probs = c(.025, 0.5, .975)))),
          apply(sigma_k2S[, gibbsPostBurn], 1, mean),
          sigma_k2S[, 1]
      )

colnames(sigma_k2_df) <- c("k", "true", ".025", "0.5", ".975", "mean", "start")
View(sigma_k2_df)

# trace plot
sigmakplot <- 1
dev.new()
par(mar = c(5.1, 6, 4.1, 2.1))
plot(sigma_k2S[sigmakplot, ], type = "l",
     main = bquote(k == .(sigmakplot)),
     ylab = bquote(sigma[k]^2),
     cex.lab = 2,
     cex.main = 2)
par(mar = c(5.1, 4.1, 4.1, 2.1))
abline(h = sigma_k20[sigmakplot])

# assess Sigma ------------------------------------------------------------

avgSigma <- apply(SigmaS[, , gibbsPostBurn],
                  c(1, 2), mean)
norm(Sigma0 - avgSigma, "F")

cbind(avgSigma[, 1], Sigma0[, 1])

d_to_avgSigma <- rep(NA, length(gibbsKeepInds))
for (s in 1:length(gibbsKeepInds)) {
    d_to_avgSigma[s] <- norm(avgSigma - SigmaS[, , s], "F")
}

dev.new()
plot(d_to_avgSigma, type = "l", ylab = "Frob. dist.",
     main = "Frobenius distance of samples to mean Sigma",
     cex.main = 2, cex.lab = 2)

quantile(d_to_avgSigma[gibbsPostBurn], c(.025, .5, .975))

diagii <- 4
plot(Re(SigmaS[diagii, diagii, ]), type = "l",
     main = bquote(Sigma[.(diagii) *","* .(diagii)]),
     cex.main = 3)
abline(h = Re(Sigma0[diagii, diagii]))

# assess Uk ---------------------------------------------------------------

K <- 2
# K <- param_list$K

for (k in 1:K) {
    avgUk <- apply(U_kS[, , k, gibbsPostBurn],
                   c(1, 2), mean)
    avgUk <- avgUk %*% solve( EigenR::Eigen_sqrt( t(Conj(avgUk)) %*% avgUk ) )
    print(fast_evec_Frob_stat(U_k0[, , k], avgUk))
}

# recall that the usual Frobenius norm is not good for comparing axes
norm(U_k0[, , k] - avgUk, "F")
fast_evec_Frob_stat(U_k0[, , k], avgUk)
Uk_dist <- evec_Frob_stat(U_k0[, , k], avgUk, returnMats = TRUE)

# the average axes are quite close, once you adjust them
cbind(U_k0[, 1, k], avgUk[, 1])
cbind(Uk_dist$Yopt[, 1], Uk_dist$Xopt[, 1])

Uk_dist |> lapply(round, digits = 4)

Uk_df <- NULL

dev.new()
par(mfrow = c(2, 1))
for (k in 1:2) {
    avgUk <- apply(U_kS[, , k, gibbsPostBurn],
                   c(1, 2), mean)
    avgUk <- avgUk %*% solve( EigenR::Eigen_sqrt( t(Conj(avgUk)) %*% avgUk ) )
    
    d_to_avgUk <- rep(NA, length(gibbsKeepInds))
    for (s in 1:length(gibbsKeepInds)) {
        d_to_avgUk[s] <- fast_evec_Frob_stat(U_kS[, , k, s], avgUk)
    }
    
    plot(d_to_avgUk, type = "l", 
         main = bquote(k == .(k)),
         cex.main = 3, cex.lab = 2)
    
    Uk_df <- rbind(Uk_df,
                   c(k, 
                     unname(quantile(d_to_avgUk[gibbsPostBurn], probs = c(.025, .975))),
                     fast_evec_Frob_stat(avgUk, U_k0[, , k]))
                     )
}
par(mfrow = c(1, 1))

colnames(Uk_df) <- c("k", ".025", ".975", "d_avg_to_true")
Uk_df


# compare posterior full SDMs ---------------------------------------------

# if we compute the full SDM with the posterior samples, even if there is
# bi-modality, is the average SDM close to the true SDM?

# SDM is sigmak2 * (U Lambda U^H + I_P)
length(gibbsPostBurn)

k <- 2
whichit <- which(gibbsPostBurn)[1]

SDM0 <- sigma_k20[k] * (U_k0[, , k] %*% diag(Lambda_k0[, k]) %*% 
                            t(Conj(U_k0[, , k])) + diag(P))

SDMavg <- matrix(0, P, P)
proj_sum <- matrix(0, P, P)

for (s in which(gibbsPostBurn)) {
    SDMavg <- SDMavg + (1 / sum(gibbsPostBurn)) *
        sigma_k2S[k, s] * (U_kS[, , k, s] %*% diag(Lambda_kS[, k, s]) %*%
                               t(Conj(U_kS[, , k, s])) + diag(P))
    
    proj_sum <- proj_sum + 
        U_kS[, , k, s] %*% t(Conj(U_kS[, , k, s]))
}

eigen_proj_sum <- eigen(proj_sum)$vectors[, 1:d]
proj_avg <- eigen_proj_sum %*% t(Conj(eigen_proj_sum))

true_proj <- U_k0[, , k] %*% t(Conj(U_k0[, , k]))

norm(SDM0 - SDMavg, type = "F")
SDM0 |> round(3)
SDMavg |> round(3)

cbind(SDM0[, 1], SDMavg[, 1])

# swapping posterior samples ----------------------------------------------

Lambda_k0

# if Lambda1 > Lambda2, keep the order of Lambda and columns of Uk
# else, switch the order of Lambda and the columns of U

newUks <- array(NA, c(P, d, sum(gibbsPostBurn)))
newLambdaks <- array(NA, c(d, sum(gibbsPostBurn)))

new_proj_sum <- matrix(0, P, P)

# thisi is just an index counter used in the loop
thisi <- 1
for (s in which(gibbsPostBurn)) {
    if (Lambda_kS[1, k, s] > Lambda_kS[2, k, s]) {
        newLambdaks[, thisi] <- Lambda_kS[, k, s]
        newUks[, , thisi] <- U_kS[, , k, s]
    } else {
        newLambdaks[, thisi] <- Lambda_kS[2:1, k, s]
        newUks[, , thisi] <- U_kS[, 2:1, k, s]
    }
    
    new_proj_sum <- new_proj_sum + 
        U_kS[, , k, s] %*% t(Conj(U_kS[, , k, s]))
        
    thisi <- thisi + 1
}

newavgUk <- apply(newUks, c(1, 2), mean)
newavgUk <- newavgUk %*% 
    solve( EigenR::Eigen_sqrt( t(Conj(newavgUk)) %*% newavgUk ) )
evec_Frob_stat(U_k0[, , k], newavgUk, returnMats = TRUE) |> lapply(round, digits = 4)

evec_Frob_stat(U_k0[, , k], avgUk, returnMats = TRUE)

# eigen_proj_sum <- eigen(proj_sum)$vectors[, 1:d]
# proj_avg <- eigen_proj_sum %*% t(Conj(eigen_proj_sum))

eigen_new_proj_sum <- eigen(new_proj_sum)$vectors[, 1:d]
new_proj_avg <- eigen_new_proj_sum %*% t(Conj(eigen_new_proj_sum))



# comparing average projection matrix -------------------------------------

norm(true_proj - proj_avg, type = "F")
true_proj |> round(4)
proj_avg |> round(3)

norm(true_proj - new_proj_avg, type = "F")
true_proj |> round(4)
new_proj_avg |> round(3)


# only selecting samples from certain clusters ----------------------------


cluster_avgUk <- 
    apply(U_kS[, , k, (Lambda_kS[1, k, gibbsPostBurn] > 
                           Lambda_kS[2, k, gibbsPostBurn])],
      c(1, 2), mean)
cluster_avgUk <- cluster_avgUk %*% 
    solve( EigenR::Eigen_sqrt( t(Conj(cluster_avgUk)) %*% cluster_avgUk ) )

evec_Frob_stat(U_k0[, , k], cluster_avgUk, returnMats = TRUE) |> lapply(round, digits = 4)

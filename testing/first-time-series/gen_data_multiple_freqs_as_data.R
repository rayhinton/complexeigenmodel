# generate data for the model using multiple frequencies as the data

# that is, we actually only are modeling one time series, but in the model, we
# use K = T, the number of time points collected.
# for each data value, we use an estimate of the SDM at a Fourier frequency

# Need to:
# - generate true parameters for this model
# - generate data for this model

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
n_k_approx <- 64
K <- n_k_approx - 2
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
# omegaish <- 3/16
# k_freq <- round(n_k_approx * omegaish)  # for example
# omega0 <- k_freq / n_k_approx

# AR parameters
phis1 <- c(1, -.25)
phis2 <- c(-.75, -0.5)

# scale
sigma02 <- c(5)

# U, the matrix that creates the observed time series
U <- qr.Q(qr(matrix(rnorm(P*d), ncol = d)))
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

# n_k_approx is the length of the time series
omegas <- seq(1/n_k_approx, 1 - 1/n_k_approx, by = 1/n_k_approx) 
# K is the number of SDMs to estimate, which we will do n_k - 2

set.seed(data_seed)
for (k in 1:K) {
    omega0 <- omegas[k]
    # generate parameters
    Lambda_k0[, k] <- c(sdf_ar2(phis1, omega0), 
                        sdf_ar2(phis2, omega0))
    sigma_k20[k] <- sigma02

    U_k0[, , k] <- U
    
    fXomega0 <- diag(c(sdf_ar2(phis1, omega0),
                       sdf_ar2(phis2, omega0)))
    ##### calculate the exact SDM for the P-dim time series
    # need to confirm the SDM relationship for this sort of factor model
    
    fYomega0 <- sigma02 * ( U %*% fXomega0 %*% t(Conj(U)) + diag(P) )
    
    data_list[[k]] <- L_numtprs*fYomega0
    
}

data_list[[1]]
data_list[[2]]
data_list[[3]]

# estimates of sigmak2 ----------------------------------------------------

# eigen(data_list[[1]] / n_k[1])$values
# eigen(data_list[[1]])$values
# eigen(data_list[[1]])$values/L_numtprs
# sigma_k20

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

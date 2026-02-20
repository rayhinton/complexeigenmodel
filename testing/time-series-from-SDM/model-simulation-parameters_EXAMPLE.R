
# dataseed <- 21092025
# dataseed <- 22092025
# dataseed <- 10102025
# dataseed <- 17102025
dataseed <- 28112025

# parseed <- 314
# parseed <- 3141
parseed <- 3141
# parseed <- 963456789

P <- 4
d <- 2
K <- 2
Tt <- 1024 # length of time series
LL <- round(sqrt(Tt))
# LL <- round(P+1)

useMclapply <- FALSE
n_cores <- 2

use_true_SDMs <- FALSE
use_Id_Sigmal_init <- FALSE
sample_true_Ukl0 <- FALSE
all_same_VAR_pars <- TRUE

# options include: 1RW, 2RWPN
Lambda_prior <- "2RWPN"

gibbsIts <- 5000
burnin <- 0.5
gibbsPrint <- 100

num_freqs <- Tt/2 - 1

# time series generation parameters
# possible par. gen. methods: "smoothly-similar-Ukl", "truncVAR", "baseVAR"
TS_par_gen_method <- "baseVAR"
# number of knots in the function that generates Lambda curves
n_knots <- 4
U_k_n_basis <- 3
U_k_scale_base <- 0.5
U_k_scale_k <- 0.1

# Geodesic slice sampling parameters
w <- 10
m <- 2

### Ukl MH tuning parameters
# tau_Uk <- rep(.1, K)
tau_Ukl <- array(0.1, c(K, num_freqs))
num_tau_check <- 20
show_tau_tune_summ <- FALSE
show_n_Sig_summary <- FALSE
doCayleyZeros <- FALSE
CayleyZeroProb <- 0.5

### Sigmal MH tuning parameters
n_Sig <- rep(max(50, P*2), num_freqs)
Sigma_add <- 0.001
useEigenR <- TRUE 
byCholesky <- TRUE

# grid of frequencies to calculate over
omegaw <- seq(1/Tt, by = 1/Tt, length.out = num_freqs)

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

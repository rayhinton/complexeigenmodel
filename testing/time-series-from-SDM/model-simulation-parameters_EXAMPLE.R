
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

# options include: 1RW, 2RWPN
Lambda_prior <- "2RWPN"

gibbsIts <- 5000
burnin <- 0.5
gibbsPrint <- 100

num_freqs <- Tt/2 - 1

# time series generation parameters
# number of knots in the function that generates Lambda curves
n_knots <- 4

# Geodesic slice sampling parameters
w <- 10
m <- 2

### Ukl MH tuning parameters
# tau_Uk <- rep(.1, K)
tau_Ukl <- array(0.1, c(K, num_freqs))
num_tau_check <- 20
show_tau_tune_summ <- FALSE
doCayleyZeros <- FALSE
CayleyZeroProb <- 0.5

### Sigmal MH tuning parameters
n_Sig <- rep(50, num_freqs)

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

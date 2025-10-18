# generating simple time series

##### setup

P <- 4
d <- 2
K <- 1
# length of each time series, for k = 1, ..., K
n_k <- c(1024)

##### choose a frequency
# naively, choose pi/2, or in terms of non-pi frequencies, 1/4
# k_freq <- 1024  # for example
k_freq <- round(n_k[1] * .1)  # for example
omega0 <- k_freq / n_k[1]

# omega0 <- .0625


##### simulate 2 AR(2) processes
# arima.sim
# need to choose the 2 AR parameters for each, 4 parameters total
# and need to make sure to choose them so that the TS are stationary

phis1 <- c(.5, -0.5)
phis2 <- c(-0.5, 0.25)

arcomp1 <- list(ar = phis1)
arcomp2 <- list(ar = phis2)
f1 <- arima.sim(arcomp1, n_k[1])
f2 <- arima.sim(arcomp2, n_k[1])

# call X the multivariate time series, which is a 2xT matrix
X <- ts(cbind(f1, f2), start = start(f1), frequency = frequency(f1))
plot(X)

##### generate sigma, scale parameter
# generate from the prior
sigma02 <- 5

##### calculate observed P-dim time series
# create U matrix, Px2, with orthogonal columns
# make sure the time series is arranged as 2xT, where T is the number of time points
# create new time series as UX, a PxT matrix

U <- qr.Q(qr(matrix(rnorm(P*2), ncol = 2)))
# confirm they are orthogonal, the result should be 2x2 I
crossprod(U)

# generate more unit variance noise, PxT, call it epsilon
# create final time series, Y = sigma * (UX + epsilon), NOT sigma^2!!!
Ymat <- sqrt(sigma02) * (U %*% t(X) + matrix(rnorm(P*n_k[1]), ncol = n_k[1]))
Yts <- ts(t(Ymat), start = start(f1), frequency = frequency(f1))

##### calculate the exact SDM for the 2 factors
# need spectral density function for AR(2)
sdf_ar2 <- function(phis, omega, sigma2 = 1) {
    phi1 <- phis[1]
    phi2 <- phis[2]
    
    denom <- 1 + phi1^2 + phi2^2 + 2*phi1*(phi2 - 1) * cos(2*pi*omega) - 
        2*phi2*cos(4*pi*omega)
    
    return(sigma2 / denom)
}

fXomega0 <- diag(c(sdf_ar2(phis1, omega0), sdf_ar2(phis2, omega0)))
##### calculate the exact SDM for the P-dim time series
# need to confirm the SDM relationship for this sort of factor model

fYomega0 <- sigma02 * ( U %*% fXomega0 %*% t(Conj(U)) + diag(P) )

##### multitaper estimator: Sine tapers

num_tprs <- round(sqrt(n_k))
cat(paste("num_tprs =", num_tprs))
Utp <- waveslim::sine.taper(n_k, num_tprs)
Y_tp <- apply(Utp, MARGIN = 2, function(u) u*Yts, simplify = FALSE)
X_tp <- apply(Utp, MARGIN = 2, function(u) u*X, simplify = FALSE)

# R recycling

fr_exp <- exp(-1i*2*pi*omega0*(1:n_k[1]))

Shatf_Y <- matrix(0, P, P)
Shatf_X <- matrix(0, d, d)
for (ell in 1:num_tprs) {
    JY_ell <- colSums(fr_exp * Y_tp[[ell]])
    JX_ell <- colSums(fr_exp * X_tp[[ell]])
    
    Shatf_Y <- Shatf_Y + 1/num_tprs * JY_ell %*% t(Conj(JY_ell))
    Shatf_X <- Shatf_X + 1/num_tprs * JX_ell %*% t(Conj(JX_ell))
}

diag(Shatf_Y)
diag(fYomega0)

diag(Shatf_X)
diag(fXomega0)

##### alternative multitaper estimator method, from Dr. Namdari

## Multitaper estimate of the spectral density matrix

F_tp_list <- lapply(Y_tp, FUN = function(Y) astsa::mvspec(Y,plot = FALSE))

len_freq <- n_k[1]/2
F_tp1 <- array(0, c(P, P, len_freq))

for (ell in 1:len_freq) {
    for(j in 1:length(F_tp_list)){
        F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
    }
    F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
}

f_xx1 <- F_tp1*n_k[1]

# thefreq <- which(F_tp_list[[1]]$freq >= omega0)[1]
# F_tp_list[[1]]$freq[thefreq]
thefreq <- k_freq  # mvspec uses 1-based indexing starting from freq 0
F_tp_list[[1]]$freq[thefreq]

diag(F_tp1[, , thefreq])
diag(F_tp1[, , thefreq]*n_k[1])

diag(fYomega0)

##### are they giving similar estimates?
# Calculate relative differences
rel_diff <- abs(diag(Shatf_Y) - diag(F_tp1[,,thefreq]*n_k[1])) / 
    abs(diag(F_tp1[,,thefreq]*n_k[1]))
print(rel_diff)
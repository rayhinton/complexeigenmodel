# simulate data
library(astsa)

# dd <- astsa::sarima.sim()

source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcomplex_wishart.R")

N <- 2000
P <- 2

f0 <- array(NA, c(P, P, N))

# spectral density matrices are Hermitian p.d. - generate from CW dist.
f0[, , N/2] <- rWishart(1, P+1, diag(P))[, , 1]
f0[, , N] <- rWishart(1, P+1, diag(P))[, , 1]
for (k in 1:(N/2-1)) {
    f0[, , k] <- rcomplex_wishart(P+1, diag(P))
    f0[, , N-k] <- t(f0[, , k])
}

# Implement Chambers 1995

# calculate Zk
Zk <- matrix(NA, P, N)

# Case 1: k = N/2 and k = N
Zk[, N/2] <- rnorm(P, sd = sqrt(2))
Zk[, N] <- rnorm(P, sd = sqrt(2))

for (k in 1:(N/2-1)) {
    ReZk <- rnorm(P)
    ImZk <- rnorm(P)

    # Case 2: 1 <= k <= N/2 - 1, generate Zk as:
    Zk[, k] <- ReZk + 1i*ImZk
    # Case 3: set Z[n-k] in terms of Zk
    Zk[, N-k] <- ReZk - 1i*ImZk
}

t(Zk[, 1:5])

# Intermediate step
Vk <- matrix(NA, P, N)

for (k in c(1:(N/2), N)) {
    Uk <- eigen(f0[, , k])$vectors
    Mksqrt <- diag(sqrt(eigen(f0[, , k])$values))
    Vk[, k] <- Uk %*% Mksqrt %*% Zk[, k]
    
    if (k <= (N/2-1)) {
        Vk[, N-k] <- Conj(Vk[, k])
    }
}

lambdak <- 2*pi*(1:N)/N

Xt <- matrix(0, P, N)

for (t in 1:N) {
    for (k in 1:N) {
        Xt[, t] <- Xt[, t] + sqrt(pi/N) * Vk[, k] * exp(1i * t * lambdak[k])
    }
}
t(Xt[, 1:5])
max(abs(Im(Xt)))

matplot(t(Re(Xt[, 1:100])), type = "l")

fxxspec <- astsa::mvspec(t(Re(Xt)), spans = c(25, 25), taper = .1)

freq_i <- 903
fxxspec$fxx[, , freq_i]
f0[, , freq_i]

fxxspec$fxx[, , freq_i] / f0[, , freq_i]

# What if I sample something with a known SDM, like 2 independent AR(1) processes?
ts1 <- astsa::sarima.sim(ar = 0.5, n = 5000)
ts2 <- astsa::sarima.sim(ar = 0.25, n = 5000)

ts12 <- cbind(ts1, ts2)

spec12 <- astsa::mvspec(ts12, spans = c(51, 51), taper = .5)
spec12$fxx[, , 1:5]
1 / Mod(1 - 0.25 * exp(-1i*(1:5)/5000))^2

plot(1 / Mod(1 - 0.5 * exp(-2*pi*1i*(0:4999)/5000))^2)



# Claude example ----------------------------------------------------------

# Set parameters for two independent AR(1) processes
set.seed(123)
n <- 500  # time series length
phi1 <- 0.7  # AR parameter for first process
phi2 <- -0.4  # AR parameter for second process
sigma1 <- 1.0  # innovation SD for first process
sigma2 <- 1.5  # innovation SD for second process

# Generate two independent AR(1) time series
ts1 <- arima.sim(model = list(ar = phi1), n = n, sd = sigma1)
ts2 <- arima.sim(model = list(ar = phi2), n = n, sd = sigma2)

# Combine into bivariate time series
Y <- cbind(ts1, ts2)

# Estimate spectral density matrix using mvspec from astsa
# This uses smoothed periodogram with default settings
spec_est <- mvspec(Y, spans = c(17, 17), taper = 0.1, plot = FALSE)

# Extract frequency grid and estimated spectral density matrices
freqs <- spec_est$freq
spec_matrices <- spec_est$fxx  # array of spectral density matrices

# Function to compute true spectral density for AR(1) process
true_ar1_spec <- function(phi, sigma2, omega) {
    sigma2 / abs(1 - phi * exp(-1i * omega))^2
}

# Compute true spectral density matrices
n_freqs <- length(freqs)
true_spec_matrices <- array(0, dim = c(2, 2, n_freqs))

for (k in 1:n_freqs) {
    omega <- 2 * pi * freqs[k]
    
    # True SDM is block diagonal for independent processes
    true_spec_matrices[1, 1, k] <- true_ar1_spec(phi1, sigma1^2, omega)
    true_spec_matrices[2, 2, k] <- true_ar1_spec(phi2, sigma2^2, omega)
    true_spec_matrices[1, 2, k] <- 0  # off-diagonal should be zero
    true_spec_matrices[2, 1, k] <- 0
}

# Comparison metrics
# 1. Frobenius norm of difference at each frequency
frobenius_errors <- numeric(n_freqs)
for (k in 1:n_freqs) {
    diff_matrix <- spec_matrices[,,k] - true_spec_matrices[,,k]
    frobenius_errors[k] <- norm(diff_matrix, type = "F")
}

spec_matrices[,,100]
true_spec_matrices[,,100]

# Plotting comparisons
par(mfrow = c(2, 2))

# Plot 1: Diagonal element (1,1) comparison
plot(freqs, Re(true_spec_matrices[1,1,]), type = "l", col = "red", lty = 2, 
     xlab = "Frequency", ylab = "Spectral Density", 
     main = "Diagonal (1,1): Estimated vs True")
lines(freqs, Re(spec_matrices[1,1,]), col = "blue")
legend("topright", c("Estimated", "True"), col = c("blue", "red"), lty = c(1, 2))

# Plot 2: Diagonal element (2,2) comparison
plot(freqs, Re(true_spec_matrices[2,2,]), type = "l", col = "red", lty = 2,
     xlab = "Frequency", ylab = "Spectral Density",
     main = "Diagonal (2,2): Estimated vs True")
lines(freqs, Re(spec_matrices[2,2,]), col = "blue")
legend("topright", c("Estimated", "True"), col = c("blue", "red"), lty = c(1, 2))

# Plot 3: Off-diagonal magnitude (should be near zero)
plot(freqs, abs(spec_matrices[1,2,]), type = "l", col = "blue",
     xlab = "Frequency", ylab = "Magnitude",
     main = "Off-diagonal (1,2) Magnitude")
abline(h = 0, col = "red", lty = 2)

# Plot 4: Frobenius error across frequencies
plot(freqs, frobenius_errors, type = "l", col = "purple",
     xlab = "Frequency", ylab = "Frobenius Error",
     main = "Frobenius Norm Error")

par(mfrow = c(1, 1))

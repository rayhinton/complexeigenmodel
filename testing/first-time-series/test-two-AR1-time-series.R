# Load required library
library(astsa)

# Set parameters for two independent AR(1) processes
set.seed(123)  # for reproducibility
n <- 500       # sample size
phi1 <- 0.7    # AR coefficient for series 1
phi2 <- -0.5   # AR coefficient for series 2
sigma1 <- 1    # innovation sd for series 1
sigma2 <- 1.5  # innovation sd for series 2

# Generate the two independent AR(1) series
x1 <- arima.sim(n = n, list(ar = phi1), sd = sigma1)
x2 <- arima.sim(n = n, list(ar = phi2), sd = sigma2)

# Combine into multivariate series
X <- cbind(x1, x2)

# Choose a frequency to evaluate (e.g., 0.1 cycles per unit)
freq_eval <- 0.1
omega <- 2 * pi * freq_eval  # convert to radians

# Calculate TRUE spectral density matrix at this frequency
# For AR(1): f(omega) = sigma^2 / |1 - phi * exp(-i*omega)|^2
#                     = sigma^2 / ((1 - phi*cos(omega))^2 + (phi*sin(omega))^2)

denom1 <- (1 - phi1*cos(omega))^2 + (phi1*sin(omega))^2
denom2 <- (1 - phi2*cos(omega))^2 + (phi2*sin(omega))^2

f11_true <- sigma1^2 / denom1
f22_true <- sigma2^2 / denom2
f12_true <- 0 + 0i  # zero cross-spectrum since independent (complex number)
f21_true <- 0 + 0i  # conjugate of f12

# True SDM (Hermitian matrix)
SDM_true <- matrix(c(f11_true, f21_true, f12_true, f22_true), 2, 2)

cat("True Spectral Density Matrix at frequency", freq_eval, ":\n")
print(SDM_true)
cat("\n")

# Estimate using mvspec
spec_est <- mvspec(X, plot = FALSE, detrend = TRUE)

# Find the closest frequency in the estimate to our evaluation frequency
freq_idx <- which.min(abs(spec_est$freq - freq_eval))
actual_freq <- spec_est$freq[freq_idx]

# Extract estimated SDM at this frequency
# spec_est$fxx is a 3D array: [series1, series2, frequency]
SDM_est <- spec_est$fxx[, , freq_idx]

cat("Estimated Spectral Density Matrix at frequency", actual_freq, ":\n")
print(SDM_est)
cat("\n")

# Verify Hermitian property (should be TRUE or very close to TRUE)
cat("Is estimated SDM Hermitian? Max difference from conjugate transpose:\n")
hermitian_diff <- max(abs(SDM_est - Conj(t(SDM_est))))
cat(hermitian_diff, "\n\n")

# Calculate differences (using absolute values since SDM_est might be complex)
abs_diff <- abs(SDM_est - SDM_true)
cat("Absolute differences:\n")
print(abs_diff)
cat("\n")

# Calculate relative differences for diagonal elements (which are real)
rel_diff_diag <- abs(Re(diag(SDM_est)) - diag(SDM_true)) / diag(SDM_true)
cat("Relative differences for diagonal elements (%):\n")
cat("f11:", round(rel_diff_diag[1] * 100, 2), "%\n")
cat("f22:", round(rel_diff_diag[2] * 100, 2), "%\n\n")

# Check coherence (should be near 0 for independent series)
coherence <- abs(SDM_est[1,2])^2 / (SDM_est[1,1] * SDM_est[2,2])
cat("Estimated squared coherence at frequency", actual_freq, ":", coherence, "\n")
cat("(Should be near 0 for independent series)\n\n")

# Additional diagnostics
cat("Some potential sources of discrepancy:\n")
cat("1. Sample size used:", n, "\n")
cat("2. Periodogram ordinates available:", length(spec_est$freq), "\n")
cat("3. Frequency resolution:", spec_est$freq[2] - spec_est$freq[1], "\n")
cat("4. Taper used:", spec_est$taper, "\n")
cat("5. Bandwidth:", spec_est$bandwidth, "\n")
cat("6. Dimensions of fxx array:", paste(dim(spec_est$fxx), collapse = " x "), "\n")
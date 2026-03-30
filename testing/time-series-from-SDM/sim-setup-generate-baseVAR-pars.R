# sim-setup

# sigmak02 <- rgamma(K, 1, 1)
# 
# fkTR <- array(NA, c(P, P, K, Tt))
# U_kl0 <- array(NA, c(P, d, K, Tt))
# Lambdakl0 <- array(NA, c(d, K, Tt))

VAR_pars <- generate_VAR1_coef(P, 0.9)
noiseSigma <- generate_AR1_covariance(P, sigma2 = 1, rho = 0.5)

U_kl0 <- array(NA, c(P, d, K, Tt))
base_VAR1_SDM <- array(NA, c(P, P, Tt))
base_trunc_SDM <- array(NA, c(P, P, Tt))
fkTR <- array(NA, c(P, P, K, Tt))
Q_k <- array(NA, c(P, P, K))

for (k in 1:K) {
    temp <- matrix(rnorm(P^2), P, P) * U_k_scale_k
    A_k <- (temp - t(temp)) / 2
    Q_k[, , k] <- expm::expm(A_k)
}

angles_to_Ukl0 <- array(NA, c(d, K, num_freqs))
princ_angles_to_Ukl0 <- array(NA, c(d, K, num_freqs))
for (t in 1:Tt) {
    Hz <- solve( diag(P) - exp(-1i*2*pi * t/Tt) * VAR_pars )
    fomega <- Hz %*% noiseSigma %*% t(Conj(Hz))
    f_evd <- eigen(fomega)
    
    thisU <- f_evd$vectors[, 1:d]
    thisLambda <- f_evd$values[1:d]
    
    base_VAR1_SDM[, , t] <- fomega
    base_trunc_SDM[, , t] <- thisU %*% diag(thisLambda) %*% t(Conj(thisU)) +
        diag(P)
    
    for (k in 1:K) {
        U_kl0[, , k, t] <- Q_k[, , k] %*% thisU
        fkTR[, , k, t] <- 
            sigmak02[k] * (U_kl0[, , k, t] %*% diag(thisLambda) %*% 
                               t(Conj(U_kl0[, , k, t])) + diag(P))
        
        Lambdakl0[, k, t] <- thisLambda
        
        if (t <= num_freqs) {
            angles_to_Ukl0[, k, t] <- 
                acos(Mod( diag(t(Conj(U_kl0[, , k, t])) %*% thisU) ))
            
            princ_angles_to_Ukl0[, k, t] <-
                acos(svd( t(Conj(U_kl0[, , k, t])) %*% thisU )$d)
        }
    }
}

plot(angles_to_Ukl0[1, 1, ], type = "l")
plot(princ_angles_to_Ukl0[1, 1, ], type = "l")

# find mean angle at each frequency and convert to data frame for ggplot
mean_angles_to_Ukl0 <- apply(angles_to_Ukl0, c(1, 3), mean)
mean_angles_to_Ukl0 <- as.data.frame(t(mean_angles_to_Ukl0))
colnames(mean_angles_to_Ukl0) <- paste0("j", 1:d)
mean_angles_to_Ukl0$freq_index <- 1:num_freqs
mean_angles_to_Ukl0 <- 
    tidyr::pivot_longer(mean_angles_to_Ukl0, 
                        -freq_index, names_to = "j", values_to = "angle")

# plot
plotp <- ggplot(mean_angles_to_Ukl0, aes(freq_index, angle, color = j)) + 
    geom_line() +
    labs(title = "mean of vector angles from Ukl0 to base VAR Ul0",
         x = "freq. index",
         y = "angle (radians)")
# save PDF
print(plotp)
save_plot_pdf(file.path(result_dir, "Ukl0-vector-angles-to-Ul0.pdf"))

# calculate mean principal angle at each freq.; save as data frame for ggplot 
mean_princ_angles_to_Ukl0 <- apply(princ_angles_to_Ukl0, c(1, 3), mean)

mean_princ_angles_to_Ukl0 <- as.data.frame(t(mean_princ_angles_to_Ukl0))
colnames(mean_princ_angles_to_Ukl0) <- paste0("j", 1:d)
mean_princ_angles_to_Ukl0$freq_index <- 1:num_freqs
mean_princ_angles_to_Ukl0 <- 
    tidyr::pivot_longer(mean_princ_angles_to_Ukl0, 
                        -freq_index, names_to = "j", values_to = "angle")
# plot
plotp <- ggplot(mean_princ_angles_to_Ukl0, aes(freq_index, angle, color = j)) + 
    geom_line() +
    labs(title = "mean of principal angles from Ukl0 to base VAR Ul0",
         x = "freq. index",
         y = "angle (radians)")
# save PDF
print(plotp)
save_plot_pdf(file.path(result_dir, "Ukl0-principal-angles-to-Ul0.pdf"))

generate_smooth_Uk <- function(P, d, K, Tt, n_basis = 3, scale_base = 0.5, scale_k = 0.1) {
    
    Uk <- array(NA_complex_, c(P, d, K, Tt))
    
    # Frequencies from 0 to π
    n_freq <- Tt/2 + 1
    omega <- 2 * pi * (0:(Tt/2)) / Tt
    
    # Random real orthonormal starting point
    U_0 <- qr.Q(qr(matrix(rnorm(P * d), P, d)))
    
    # Coefficient matrices for smooth base path
    C_R <- array(0, c(P, P, n_basis))  # skew-symmetric
    C_S <- array(0, c(P, P, n_basis))  # symmetric
    
    for (j in 1:n_basis) {
        temp <- matrix(rnorm(P^2), P, P) * scale_base
        C_R[,,j] <- (temp - t(temp)) / 2
        temp <- matrix(rnorm(P^2), P, P) * scale_base
        C_S[,,j] <- (temp + t(temp)) / 2
    }
    
    # Compute base path U*(ω)
    U_star <- array(NA_complex_, c(P, d, n_freq))
    
    for (idx in 1:n_freq) {
        om <- omega[idx]
        A <- matrix(0 + 0i, P, P)
        for (j in 1:n_basis) {
            A <- A + cos(j * om) * C_R[,,j]
            A <- A + 1i * sin(j * om) * C_S[,,j]  # vanishes at 0 and π
        }
        U_star[,,idx] <- expm::expm(A) %*% U_0
    }
    
    # Series-specific rotations (real skew-symmetric preserves real-ness at boundaries)
    for (k in 1:K) {
        temp <- matrix(rnorm(P^2), P, P) * scale_k
        A_k <- (temp - t(temp)) / 2
        Q_k <- expm::expm(A_k)
        
        for (idx in 1:n_freq) {
            U_k_omega <- Q_k %*% U_star[,,idx]
            
            if (idx == 1) {
                Uk[,,k,Tt] <- Re(U_k_omega)
            } else if (idx == n_freq) {
                Uk[,,k,Tt/2] <- Re(U_k_omega)
            } else {
                t_idx <- idx - 1
                Uk[,,k,t_idx] <- U_k_omega
                Uk[,,k,Tt - t_idx] <- Conj(U_k_omega)
            }
        }
    }
    
    return(Uk)
}

U_kl0 <- generate_smooth_Uk(P, d, K, Tt,
                            n_basis = U_k_n_basis, scale_base = U_k_scale_base,
                            scale_k = U_k_scale_k)
for (k in 1:K) {
    for(j in 1:d) {
        Lambdakl0[j, k, c(Tt, 1:(Tt/2))] <- 
            generate_smooth_positive_function(Tt/2 + 1, n_knots = n_knots, 
                                              log_mean = 1, log_sd = 0.25)$f
        Lambdakl0[j, k, (Tt/2 + 1) : (Tt - 1)] <- Lambdakl0[j, k, (Tt/2 - 1) : 1]
    }
    
    Lambdakl0[, k, ] <- Lambdakl0[order(Lambdakl0[, k, 1], decreasing = TRUE), k, ]
}

for (k in 1:K) {
    for (t in 1:Tt) {
        thisU <- U_kl0[, , k, t]
        thisLambda <- Lambdakl0[, k, t]

        fkTR[, , k, t] <- sigmak02[k] * ( thisU %*% diag(thisLambda) %*%
                                              t(Conj(thisU)) + diag(P) )
    }
}

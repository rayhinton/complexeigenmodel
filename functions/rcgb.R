# sample from complex generalized Bingham distribution

# sampling from complex matrix-variate generalized Bingham Langevin dist.
# as in Abdallah et al. 2020

library(rstiefel)
library(cmvnorm)
library(mvtnorm)

fr <- function(b, betas) {
    return(abs(sum(1 / (b + 2*betas)) - 1))
}

# rvB - sample from the real vector Bingham distribution -----------------

# input:
# A, NxN symmetric matrix parameter
# 
# output:
# u, Nx1 real-valued unit vector
rvB <- function(A, opt_upper = 1e5) {
    # based on Algorithm 4 of Abdallah 2020
    
    mat1 <- A
    N <- nrow(A)
    
    # Line 1: compute gamma
    gamm <- max(eigen(mat1)$values)
    
    # Preparation for Line 2
    mat2 <- diag(gamm, nrow = N) - mat1 
    betas <- eigen(mat2)$values
    
    # TODO can I solve this in another way?
    # TODO I should perform some sort of check that determines where I need to use a larger "upper" value
    # Line 2: solve for b
    b <- optim(0.1, fr, method = "Brent", lower = 0, upper = opt_upper, 
               betas = betas)$par
    
    # Line 4: compute Omega
    Om <- diag(N) + (2/b) * mat2
    # Line 3: compute Mstar, as a function of b
    Mst <- exp(.5 + gamm) * exp(-(N - b)/2) * (N/b)^(N/2) * det(Om)^-.5
    
    # Line 5: sample u by acceptance-rejection scheme    
    acc_prob <- 0
    u_draw <- 1
    while(u_draw >= acc_prob) {
        # Line 6: sample y from a Normal dist.
        y <- rmvnorm(1, sigma = solve(Om))
        # Line 7: normalize y
        u <- t(y) / sqrt(sum(y^2))
        # Line 8: compute the ACG density value at u
        f_ACG <- det(Om)^.5 * (t(u) %*% Om %*% u)^(-N/2)
        # Line 9: compute the vector Bingham density value at u
        f_vBL <- exp(t(u) %*% A %*% u)
        
        # Lines 10-11: accept u with the calculated probability
        u_draw <- runif(1)
        acc_prob <- f_vBL / (Mst * f_ACG)
    }
    
    # Line 12: return u, after the acceptance criteria is met
    return(u)
}

# rcvb - using the vBL, sample from complex VGBL --------------------------

# AA <- rbind(cbind(Re(A), -Im(A)),
#             cbind(Im(A), Re(A)))
# cc <- matrix(c(Re(c), Im(c)), ncol = 1)
# 
# uu <- rvBL(cc, AA) |> matrix(ncol = 2)
# complex(real = uu[, 1], imaginary = uu[, 2])

rcvb <- function(A, opt_upper = 1e5) {
    AA <- rbind(cbind(Re(A), -Im(A)),
                cbind(Im(A), Re(A)))
    print(dim(AA))
    # cc <- matrix(c(Re(c), Im(c)), ncol = 1)
    
    uu <- rvB(AA, opt_upper = opt_upper) |> matrix(ncol = 2)
    
    return(complex(real = uu[, 1], imaginary = uu[, 2]))
}


# rcvb_LN - a function to sample new column, orthogonal to others ---------

# input:
# U, an NxP matrix with orthonormal columns, for which to sample a column
# p, a scalar, the index of the column to sample
# A, NxN psd Hermitian matrix, the original matrix parameter
# B, PxP diagonal matrix, or a Px1 vector representing the diagonal entries of B

# output:
# u, Nx1 vector, to be inserted into matrix U

# Description

# This function samples a single column for U, so that the sampled column is
# orthogonal to the other columns in U. The sample column, u_p, has a vector
# Bingham distribution with parameters (A, B), restricted to the space that is
# orthogonal to the columns U[, -p]. A transformed version of the column is
# sampled, and then this result is transformed back to the originality
# dimensionality, using a basis of the left null space of U[, -p], ensuring u_p
# is orthogonal to the rest of the columns in U.
rcvb_LN <- function(U, p, A, B) {
    
    # if B is provided as a matrix, instead, use the diagonal of B
    if (is.matrix(B)) {
        B <- diag(B)
    }
    
    LN <- NullC(U[, -p])
    # this value is actually not used later
    u <- Conj(t(LN)) %*% U[, -p]
    
    # cbar <- Conj(t(LN)) %*% C[, p] 
    # Abar <- Conj(t(LN)) %*% AP[, , p] %*% LN
    Abar <- Conj(t(LN)) %*% (B[p] * A) %*% LN
    
    u <- LN %*% rcvb(Abar)
    return(u)
}

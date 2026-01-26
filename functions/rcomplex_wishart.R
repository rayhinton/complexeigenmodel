rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

#' Generate a random complex Wishart matrix
#'
#' This function generates a random complex Wishart matrix with specified
#' degrees of freedom and scale matrix.
#'
#' @param n Integer, degrees of freedom (must be >= p)
#' @param p Integer, dimension of the matrix
#' @param Sigma Complex matrix of dimension p x p, the scale matrix (must be Hermitian positive-definite)
#' @return A complex p x p Wishart distributed random matrix
#'
#' @examples
#' # Create a 3x3 identity matrix as the scale matrix
#' Sigma <- diag(3)
#' # Generate a 3x3 complex Wishart matrix with 5 degrees of freedom
#' W <- rcomplex_wishart(5, 3, Sigma)
#'
rcomplex_wishart <- function(df, Sigma, useEigenR = FALSE, byCholesky = FALSE, Htol = 1e-10) {
    
    p <- nrow(Sigma)
    
    # Parameter validation
    if (!is.numeric(df) || df < p || df != round(df)) {
        stop("df must be an integer and df >= p")
    }
    
    if (!is.numeric(p) || p <= 0 || p != round(p)) {
        stop("p must be a positive integer")
    }
    
    # Validate scale matrix
    if (!is.matrix(Sigma) || nrow(Sigma) != p || ncol(Sigma) != p) {
        stop("Sigma must be a p x p matrix")
    }
    
    # Check if Sigma is Hermitian
    if (!all(abs(Sigma - Conj(t(Sigma))) < Htol)) {
        stop("Sigma must be Hermitian (equal to its conjugate transpose)")
    }

    # calculate square root of Sigma
    if (useEigenR) {
        # C <- EigenR::Eigen_sqrt(Sigma)
        C <- EigenR::Eigen_Chol(Sigma)
    } else {
        Sigmaevd <- eigen(Sigma)
        
        if (any(Sigmaevd$values <= 0)) {
            stop("Sigma must be positive definite")
        }
        
        C <- Sigmaevd$vectors %*% diag(sqrt(Sigmaevd$values)) %*%
            t(Conj(Sigmaevd$vectors))
    }
    
    # generate W by Cholesky decomposition or matrix of CN vectors
    if (byCholesky) {
        Tt <- matrix(0 + 0i, p, p)
        diag(Tt) <- sqrt(rgamma(p, df - p + 1:p))
        Tt[upper.tri(Tt)] <- rscnorm(choose(p, 2))
        WI <- Tt %*% t(Conj(Tt))
        
        # TODO get rid of t Conj, since C is Hermitian?
        # for small matrices, it can be faster
        W <- t(Conj(C)) %*% WI %*% C
        
    } else {
        # Generate standard complex normal matrix (real and imaginary parts are independent normal)
        Z <- matrix(rnorm(df * p, sd = sqrt(0.5)), nrow = df, ncol = p) + 
            1i * matrix(rnorm(df * p, sd = sqrt(0.5)), nrow = df, ncol = p)
        
        # Transform Z using the square root of Sigma
        X <- Z %*% C
        
        # Calculate Wishart matrix: X^H * X
        W <- Conj(t(X)) %*% X
    }

    W <- ( W + t(Conj(W)) )/2
    diag(W) <- Re(diag(W))
    return(W)
}

#' Test if a matrix is Hermitian positive definite
#'
#' @param M Complex matrix to test
#' @return Logical value indicating whether M is Hermitian positive definite
#'
is_hermitian_positive_definite <- function(M) {
    # Check if M is Hermitian
    is_hermitian <- all(abs(M - Conj(t(M))) < 1e-10)
    
    # Check if M is positive definite
    eigen_values <- eigen(M, only.values = TRUE)$values
    is_positive_definite <- all(eigen_values > 0)
    
    return(is_hermitian && is_positive_definite)
}

# Example usage:
#
# # Set seed for reproducibility
# set.seed(123)
#
# # Create a 3x3 Hermitian positive definite scale matrix
# A <- matrix(rnorm(9) + 1i*rnorm(9), 3, 3)
# Sigma <- A %*% Conj(t(A)) + 3*diag(3)  # Ensure positive definiteness
#
# # Verify Sigma is Hermitian positive definite
# print(is_hermitian_positive_definite(Sigma))
#
# # Generate a complex Wishart matrix with 5 degrees of freedom
# W <- rcomplex_wishart(5, 3, Sigma)
#
# # Verify the result is Hermitian
# print(all(abs(W - Conj(t(W))) < 1e-10))

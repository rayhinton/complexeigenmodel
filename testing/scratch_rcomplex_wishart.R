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
rcomplex_wishart <- function(n, p, Sigma = NULL) {
    # Parameter validation
    if (!is.numeric(n) || n < p || n != round(n)) {
        stop("n must be an integer and n >= p")
    }
    
    if (!is.numeric(p) || p <= 0 || p != round(p)) {
        stop("p must be a positive integer")
    }
    
    # Default scale matrix is identity
    if (is.null(Sigma)) {
        Sigma <- diag(p)
    }
    
    # Validate scale matrix
    if (!is.matrix(Sigma) || nrow(Sigma) != p || ncol(Sigma) != p) {
        stop("Sigma must be a p x p matrix")
    }
    
    # Check if Sigma is Hermitian
    if (!all(abs(Sigma - Conj(t(Sigma))) < 1e-10)) {
        stop("Sigma must be Hermitian (equal to its conjugate transpose)")
    }
    
    # Use eigendecomposition approach since chol() doesn't support complex matrices
    eigen_decomp <- eigen(Sigma)
    eigen_values <- eigen_decomp$values
    eigen_vectors <- eigen_decomp$vectors
    
    # Check if Sigma is positive definite
    if (any(eigen_values <= 0)) {
        stop("Sigma must be positive definite")
    }
    
    # Compute square root of Sigma using eigendecomposition: Sigma^(1/2)
    # Sigma^(1/2) = U * D^(1/2) * U^H, where U contains eigenvectors and D is diagonal with eigenvalues
    Sigma_sqrt <- eigen_vectors %*% diag(sqrt(eigen_values)) %*% Conj(t(eigen_vectors))
    
    # Generate standard complex normal matrix (real and imaginary parts are independent normal)
    Z <- matrix(rnorm(n * p, sd = sqrt(0.5)), nrow = n, ncol = p) + 
        1i * matrix(rnorm(n * p, sd = sqrt(0.5)), nrow = n, ncol = p)
    
    # Transform Z using the square root of Sigma
    X <- Z %*% Sigma_sqrt
    
    # Calculate Wishart matrix: X^H * X
    W <- Conj(t(X)) %*% X
    
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
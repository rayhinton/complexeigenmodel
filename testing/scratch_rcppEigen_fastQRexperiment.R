# QR decomposition

# Load required packages
library(Rcpp)
library(RcppEigen)

# Create the C++ function for QR decomposition
cppFunction(depends = "RcppEigen", '
List eigenQR(const Eigen::MatrixXcd& A) {
  // Use Eigen\'s HouseholderQR for the decomposition
  Eigen::HouseholderQR<Eigen::MatrixXcd> qr(A);
  
  // Extract the Q matrix
  Eigen::MatrixXcd Q = qr.householderQ();
  
  // Extract the R matrix (upper triangular part of the triangularView)
  Eigen::MatrixXcd R = qr.matrixQR().triangularView<Eigen::Upper>();
  
  // Return the Q and R matrices as a list
  return List::create(Named("Q") = Q, Named("R") = R);
}')

# Create a wrapper R function
fast_qr <- function(A) {
    # Convert the R matrix to a complex matrix if it's not already
    if (!is.complex(A)) {
        A <- matrix(complex(real = A, imaginary = 0), nrow = nrow(A), ncol = ncol(A))
    }
    
    # Call the C++ function
    result <- eigenQR(A)
    
    # Return the result
    return(result)
}


# Create a random complex matrix
set.seed(123)
n <- 8
m <- 4
real_part <- matrix(rnorm(n*m), nrow = n)
imag_part <- matrix(rnorm(n*m), nrow = n)
A <- matrix(complex(real = real_part, imaginary = imag_part), nrow = n)

# Use the fast QR decomposition
system.time({
    qr_result <- fast_qr(A)
    Q <- qr_result$Q
    R <- qr_result$R
})

microbenchmark::microbenchmark(
    {
        qr_result <- fast_qr(A)
        Q <- qr_result$Q
        R <- qr_result$R
    },
    times = 1000
)

# Compare with R's built-in QR
system.time({
    qr_r <- qr(A)
    Q_r <- qr.Q(qr_r)
    R_r <- qr.R(qr_r)
})

microbenchmark::microbenchmark(
    {
        qr_r <- qr(A)
        Q_r <- qr.Q(qr_r)
        R_r <- qr.R(qr_r)
    },
    times = 1000
)

# Verify correctness
max(abs(Q %*% R - A)) # Should be very close to zero

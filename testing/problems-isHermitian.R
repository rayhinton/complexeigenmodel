# problems with isHermitian

library(cmvnorm)

# Create a test matrix
create_test_matrix <- function(n = 4, deviation = 1e-10) {
    # Start with a truly Hermitian matrix
    A <- matrix(complex(real = rnorm(n*n), imaginary = rnorm(n*n)), nrow = n)
    A <- (A + Conj(t(A)))/2  # Make it perfectly Hermitian
    
    # Now introduce a specific deviation in the imaginary part
    # Add a small value to an off-diagonal element
    i <- 1
    j <- 2
    
    # This is the critical part: we'll add a deviation that might be handled
    # differently by the two functions
    A[i, j] <- A[i, j] + complex(real = 0, imaginary = deviation)
    
    # Important: don't adjust the symmetric counterpart
    # This will break the skew-symmetry of the imaginary part
    # while keeping the overall matrix close to Hermitian
    
    return(A)
}

set.seed(6032025)
# Create a test matrix with a carefully chosen deviation
test_mat <- create_test_matrix(n = 4, deviation = 3e-14)

# Test both functions
(is_sym <- isSymmetric(test_mat))
(is_herm <- isHermitian(test_mat))

# Print results
cat("isSymmetric result:", is_sym, "\n")
cat("isHermitian result:", is_herm, "\n")

# If needed, examine the matrix
cat("\nReal part:\n")
print(Re(test_mat))
cat("\nImaginary part:\n")
print(Im(test_mat))
cat("\nImaginary part skew-symmetry check (should be near zero for Hermitian):\n")
print(Im(test_mat) + t(Im(test_mat)))

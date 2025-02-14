# scratch and testing for FCDs

library(cmvnorm)

# generate testing data ---------------------------------------------------

# number of groups
K <- 3

# number rows per observation matrix
nk <- c(20, 25, 30)
stopifnot(length(nk) == K)

# dimension of the observed vectors
P <- 6
# reduced dimension of interest
d <- 2

# generate true underlying covariance matrices
# function of sigma_k, U_k, Lambda_k
1/rgamma(K, 1, 1)
# generating U_k, Lambda_k - these are just the eigenvectors, eigenvalues of some complex symmetric matrix

# generate true mean vectors, mu_k

# List of observed data matrices, Y_k
Y_k <- list()

########## Begin, Claude idea for generating a random rectangular matrix with orthonormal columns
n_rows <- 5  # Number of rows (can be any number)
n_cols <- 3  # Number of columns (must be less than or equal to n_rows)

# Generate a random complex matrix
random_matrix <- matrix(complex(real = rnorm(n_rows * n_cols), 
                                imaginary = rnorm(n_rows * n_cols)), 
                        nrow = n_rows)

# Perform QR decomposition
qr_result <- qr(random_matrix)

# Extract Q matrix
Q <- qr.Q(qr_result)

# Take the first n_cols columns of Q
orthonormal_matrix <- Q[, 1:n_cols]

# verify orthonormal
Conj(t(orthonormal_matrix)) %*% orthonormal_matrix

########## end Claude idea

# sigma_k^2 ---------------------------------------------------------------

# sigma_k^2: an Inverse Gamma distribution

# depends on parameters U_k, Lambda_k, and observed P_k
# also depends on the complex Wishart degrees of freedom
# and the dimension the observation vectors, i.e. Px1

# TODO need to generate data, and the corresponding sum of squares matrix, P_k

# data should be a list, since it could be different dimensions (numbers of rows), using an array would not be appropriate


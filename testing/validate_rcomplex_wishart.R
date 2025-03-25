# validate rcomplex_wishart

# In this file, I validate the rcomplex_wishart function and show that the
# samples it generates match some known characteristics of the Complex Wishart
# distribution.

# In particular, the mean trace is compared to the expected trace, as a simple
# scalar comparison. Also, the sample mean is compared to the true mean.
# Finally, distributions of elements of a transformation of the samples are
# compared to known distributions from Nagar and Gupta 2011.

# In all cases, the generated samples have values or distributions close to the
# known true values, even if the eigenvalues of the CW matrix parameter are
# large or of different orders of magnitude.

library(EigenR)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")

# df of the CW distribution to test
n <- 100
# dimension of the CW dist. to test
P <- 4

# number of samples to draw
its <- 2000

# arrays to store results
all_my_Ts <- array(NA, c(P, P, its))
all_my_trs <- rep(NA, its)
all_my_As <- array(NA, c(P, P, its))

#####
# Generate a matrix parameter
#####

set.seed(22032025)
S_vecs <- runitary(P, P)
S_vals <- diag(seq(1000, 2, length.out = P))
Sigma_0 <- S_vecs %*% S_vals %*% t(Conj(S_vecs))

# Cholesky decomposition of inverse is used to "standardize" the samples
C <- Eigen_chol(solve(Sigma_0))

# Confirm Sigma_0 is Hermitian and positive definite (positive eigenvalues)
isSymmetric(Sigma_0)
eigen(Sigma_0)$values

# Confirm the Cholesky decomposition still returns inverse of Sigma_0
all.equal(t(Conj(C)) %*% C, solve(Sigma_0))

# generate samples. store observations and certain statistics
set.seed(21032025)
for (i in 1:its) {
    A_my_unt <- rcomplex_wishart(n, P, Sigma_0)
    
    # standardize the observation and store its Cholesky decomposition
    A_my <- C %*% A_my_unt %*% t(Conj(C))
    all_my_Ts[, , i] <- t(Conj(Eigen_chol( A_my )))
    
    # store the trace and observation
    all_my_trs[i] <- Re(sum(diag(A_my_unt)))
    all_my_As[, , i] <- A_my_unt
}

### compare Traces
n*sum(diag(Sigma_0))
mean(all_my_trs)

# plot observed trace densities
plot(density(all_my_trs), 
     main = "observed density of traces")
abline(v = Re(n*sum(diag(Sigma_0))), col = "red")

### compare mean of the observed samples
n*Sigma_0
apply(all_my_As, c(1, 2), mean)

### compare distribution of diagonal elements of Cholesky decomposition
{
    # which diagonal element
    ii <- 4
    
    # get limits to use for plotting
    truemin <- min(Re(all_my_Ts[ii, ii, ])^2)
    truemax <- max(Re(all_my_Ts[ii, ii, ])^2)
    
    # observed density from rcomplex_wishart and true Gamma density
    plot(density( Re(all_my_Ts[ii, ii, ])^2 ), lty = 2,
         main = "rcomplex_wishart and true Gamma")
    curve(dgamma(x, shape = n-ii + 1, 
                 scale = 1),
          from = truemin, to = truemax, add = TRUE,
          col = "red")
    legend("topright", legend = c("sampled", "true"), col = c(1, 2), lty = c(2, 1))
}


# References --------------------------------------------------------------

#' @article{nagar2011expectations,
#'     title={Expectations of functions of complex Wishart matrix},
#'     author={Nagar, Daya K and Gupta, Arjun K},
#'     journal={Acta applicandae mathematicae},
#'     volume={113},
#'     pages={265--288},
#'     year={2011},
#'     publisher={Springer}
#' }

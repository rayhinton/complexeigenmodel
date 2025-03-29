# sampling for V matrices

# Assuming that V has positive dimensions

# library(cmvnorm)
library(rstiefel)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

P <- 4

set.seed(9032025)
A <- rcomplex_wishart(P, P, diag(P))
B <- runif(P, 0, 10) |> sort(decreasing = TRUE)

# generate a random initial matrix, X, which should be dist. as Complex Bing.
X <- (rcomplex_wishart(P, P, diag(P)) |> eigen())$vector
# the matrix has orthonormal columns
t(Conj(X)) %*% X

# sample matrix using the function
newX <- rcBingUP_gibbs(X, A, B)
    
# the matrix is orthonormal
zapsmall(t(Conj(newX)) %*% newX)
zapsmall(newX %*% t(Conj(newX)))

# sample an X many times

gibbsIts <- 2000
Xs <- array(NA, c(P, P, gibbsIts))
# array to hold some projection matrices
projs <- array(NA, c(P, P, gibbsIts))

set.seed(25032025)
# generate a random initial matrix, X, which should be dist. as Complex Bing.
Xs[, , 1] <- (rcomplex_wishart(P, P, diag(P)) |> eigen())$vector
for (i in 2:gibbsIts) {
    Xs[, , i] <- rcBingUP_gibbs(Xs[, , i-1], A, B)
    
    projs[, , i] <- Xs[, 1:(P-1), i] %*% t(Conj(Xs[, 1:(P-1), i]))
}

(eigen(A)$vectors[, 1:(P-1)] %*% t(Conj(eigen(A)$vectors[, 1:(P-1)])))[1:5, 1]
apply(projs[, , (gibbsIts/2):gibbsIts], c(1, 2), mean)[1:5, 1]

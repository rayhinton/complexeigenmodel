# sampling for V matrices

# Assuming that V has positive dimensions

# library(cmvnorm)
library(rstiefel)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")

P <- 8

set.seed(9032025)
A <- rcwis(20, diag(P))
B <- runif(P, 0, 10) |> sort(decreasing = TRUE)

# generate a random initial matrix, X, which should be dist. as Complex Bing.
X <- (rcwis(20, diag(P)) |> eigen())$vector
# the matrix has orthonormal columns
t(Conj(X)) %*% X

# sample matrix using the function
newX <- rcBingUP_gibbs(X, A, B)
    
# the matrix is orthonormal
zapsmall(t(Conj(newX)) %*% newX)
zapsmall(newX %*% t(Conj(newX)))

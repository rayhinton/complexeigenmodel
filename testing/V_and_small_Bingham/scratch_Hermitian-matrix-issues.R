# Hermitian matrix issues

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

P <- 16

set.seed(18042025)
Gvecs <- matrix(rnorm(P*P) + rnorm(P*P)*1i, ncol = P) |> qr() |> qr.Q()

# try different G parameters, with large eigenvalues, implicit or explicit identity 
Gd <- Gvecs %*% (100*diag(P:1)) %*% t(Conj(Gvecs))
Gnd <- Gvecs %*% t(Conj(Gvecs))
G_I_nd <- Gvecs %*% diag(P) %*% t(Conj(Gvecs))

# show the diagonal entries - do they appear complex?
diag(Gd)
diag(Gnd)
diag(G_I_nd)

# test Hermitianity with R function
isSymmetric(Gd)
isSymmetric(Gnd)
isSymmetric(G_I_nd)

# see if the eigenvalues are complex
eigen(Gd)$values
eigen(Gnd)$values
eigen(G_I_nd)$values

# an initial unitary matrix for the sampler
X <- matrix(rnorm(P*P), ncol = P) |> qr() |> qr.Q()

# H parameter
H <- diag(P:1 + runif(P))

# returns an error
rcBingUP_gibbs(X, Gd, H, Imtol = 100*.Machine$double.eps)
# no error with a higher tolerance
rcBingUP_gibbs(X, Gd, H, Imtol = 1000*.Machine$double.eps)
rcBingUP_gibbs(X, Gnd, H)

               
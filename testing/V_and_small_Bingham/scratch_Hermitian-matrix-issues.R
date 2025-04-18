# Hermitian matrix issues

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

P <- 16

# problem with checking matrix dimensions

set.seed(18042025)
Gvecs <- matrix(rnorm(P*P) + rnorm(P*P)*1i, ncol = P) |> qr() |> qr.Q()

Gd <- Gvecs %*% (100*diag(P:1)) %*% t(Conj(Gvecs))
Gnd <- Gvecs %*% t(Conj(Gvecs))
G_I_nd <- Gvecs %*% diag(P) %*% t(Conj(Gvecs))

diag(Gd)
diag(Gnd)
diag(G_I_nd)

isSymmetric(Gd)
isSymmetric(Gnd)
isSymmetric(G_I_nd)

eigen(Gd)$values
eigen(Gnd)$values
eigen(G_I_nd)$values

# problem with checking matrix dimensions
# from the documentation on all(): 
# That all(logical(0)) is true is a useful convention: it ensures that...
dim(Gd) == dim(P:1)
all(dim(Gd) == dim(P:1))

X <- matrix(rnorm(P*P), ncol = P) |> qr() |> qr.Q()

# try drawing a sample
H <- diag(P:1 + runif(P))

# rcBingUP_gibbs(Xs[, , i-1], A, B, istatus = 0, 
#                Imtol = 10^ceiling(log10(.Machine$double.eps)))

rcBingUP_gibbs(X, Gd, H, Imtol = 100*.Machine$double.eps)
rcBingUP_gibbs(X, Gnd, H)

               
# V FCD

library(rstiefel)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")


# need
# U_k matrices, distributed as U_k ~ CGB(A, B, V), Pxd
# A diagonal matrix, PxP
# B diagonal matrix, dxd
# a true V_0 matrix, PxP

P <- 8
d <- 4

# TODO needs to be generated as a CGB, in order for the posterior to work
U_k
# it is semi-unitary
t(Conj(U_k)) %*% U_k

a <- seq(from = 1, by = 10, length.out = P)
A <- diag(a)
b <- seq(from = 1, by = 10, length.out = d)
B <- diag(b)

V_0 <- eigen(rcomplex_wishart(100, P, diag(P)))$vector
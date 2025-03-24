library(cmvnorm)
library(EigenR)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

# generate a random (semi-)unitary matrix
runitary <- function(nrows, ncols) {
    outmat <- complex(real = rnorm(nrows*ncols), 
                      imaginary = rnorm(nrows*ncols)) |> 
        matrix(nrow = nrows) |> 
        qr() |> 
        qr.Q()
    
    return(outmat)
}

new_rcmvnorm <-
    function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
              method = c("svd", "eigen", "chol"),
              tol= 100 * .Machine$double.eps){
        
        if(missing(mean) & missing(sigma)){stop("neither mean nor sigma supplied")}
        
        stopifnot(ishpd(sigma, tol)) # thus sigma known to be HPD
        
        if (length(mean) != nrow(sigma)) {
            stop("mean and sigma have non-conforming size")
        }
        
        method <- match.arg(method)
        if (method == "eigen") {
            ev <- eigen(sigma, symmetric=TRUE)
            D <- diag(sqrt(ev$values), ncol=length(ev$values))
            retval <- quad.tform(D, ev$vectors)
        } else if (method == "svd") {
            jj <- svd(sigma)
            if (!all(jj$d >= -sqrt(.Machine$double.eps) * abs(jj$d[1]))) {
                warning("sigma is numerically not positive definite")
            }
            D <- diag(sqrt(jj$d), ncol=length(jj$d))
            retval <- quad.tform(D, jj$u)
        } else {
            stop("option not recognized")
        }
        
        # old way
        out <- matrix(rcnorm(n*ncol(sigma)),ncol=n)
        # out <- Conj(crossprod(out,retval))
        # slight modification: remove Conj
        out <- crossprod(out,retval)
        
        # my new way
        # out <- matrix(rcnorm(n*ncol(sigma)),nrow=n)
        # out <- out %*% retval
        
        # original code, continued
        out <- sweep(out, 2, mean, "+")
        colnames(out) <- names(mean)
        return(out)
    }

new_rcwis <- function(n,S){
    if(length(S)==1){S <- diag(nrow=S)}
    cprod(new_rcmvnorm(n,sigma=S))
} 

n <- 100
P <- 2

its <- 2000
all_u_Ts <- array(NA, c(P, P, its))
all_u_trs <- rep(NA, its)
all_u_col1s <- matrix(NA, P, its)

all_my_Ts <- array(NA, c(P, P, its))
all_my_trs <- rep(NA, its)
all_my_col1s <- matrix(NA, P, its)

#####
# Generate a matrix parameter
#
# Cases of matrix parameter:
# - Identity: distributions seem to match
# - Diagonal: distributions seem to match
# - HPD with Identity or scaled Identity eigenvalues: dist. seem to match
# - HPD with different eigenvalues: dist. do not match
#
#####

set.seed(22032025)
S_vecs <- runitary(P, P)
S_vals <- diag(seq(50, 2, length.out = P))
Sigma_0 <- S_vecs %*% S_vals %*% t(Conj(S_vecs))

# Sigma_0 <- diag(seq(10, 1, length.out = P))

C <- Eigen_chol(solve(Sigma_0))
Eigen_chol(solve(Conj(Sigma_0)))


isSymmetric(Sigma_0)
isHermitian(Sigma_0, tol = 1e-12)
ishpd(Sigma_0)

all.equal(t(Conj(C)) %*% C, solve(Sigma_0))

set.seed(21032025)
for (i in 1:its) {
    # changing this to Conj(Sigma_0) fixes the original function 
    A_u_unt <- rcwis(n, Conj(Sigma_0))
    # A_u_unt <- new_rcwis(n, Sigma_0)
    A_u <- C %*% A_u_unt %*% t(Conj(C))
    all_u_Ts[, , i] <- t(Conj(Eigen_chol( A_u )))
    
    all_u_trs[i] <- Re(sum(diag(A_u_unt)))
    all_u_col1s[, i] <- A_u_unt[, 1]
    
    A_my_unt <- rcomplex_wishart(n, P, Sigma_0)
    A_my <- C %*% A_my_unt %*% t(Conj(C))
    all_my_Ts[, , i] <- t(Conj(Eigen_chol( A_my )))
    
    all_my_trs[i] <- Re(sum(diag(A_my_unt)))
    all_my_col1s[, i] <- A_my_unt[, 1]
}

### compare Traces
n*sum(diag(Sigma_0))
mean(all_u_trs)
mean(all_my_trs)

# plot observed trace densities
plot(density(all_u_trs))
lines(density(all_my_trs), lty = 2)

### compare means of the first column
n*Sigma_0[, 1]
rowMeans(all_u_col1s)
rowMeans(all_my_col1s)

### compare distribution of diagonal elements of Cholesky decomposition
{
ii <- 2

truemin <- min(Re(all_my_Ts[ii, ii, ])^2)
truemax <- max(Re(all_my_Ts[ii, ii, ])^2)

par(mfrow = c(1, 2))

# observed density from rcwis
plot(density( Re(all_u_Ts[ii, ii, ])^2 ), main = "rcwis")
legend("topright", legend = c("rcwis"), col = 1, lty = 1)

# observed density from rcomplex_wishart and true Gamma density
plot(density( Re(all_my_Ts[ii, ii, ])^2 ), lty = 2,
     main = "rcomplex_wishart and true Gamma")
curve(dgamma(x, shape = n-ii + 1, 
             scale = 1),
      from = 60, to = 150, add = TRUE,
      col = "red")
legend("topright", legend = c("mine", "true"), col = c(1, 2), lty = c(2, 1))

par(mfrow = c(1, 1))
}

# checking cmvnorm --------------------------------------------------------

# Real normal, first

n <- 1000

set.seed(22032025)
S_real <- matrix(rnorm(P*P), ncol = P)
Sigma_real <- S_real %*% t(S_real) + diag(P)
isSymmetric(Sigma_real)
eigen(Sigma_real)$values

# X <- mvtnorm::rmvnorm(n, mean = rep(0, P), sigma = diag(P))
X <- mvtnorm::rmvnorm(n, mean = rep(0, P), sigma = Sigma_real)
qqnorm(X[, 1])
qqline(X[, 1])

plot(X[, 1], X[, 2])
abline(lm(X[, 2] ~ X[, 1]))
abline(0, sign(Sigma_real[1, 2]) * 1/sqrt(abs(Sigma_real[1, 2])), col = "red")

abline(0, cov2cor(Sigma_real)[1, 2], col = "green")
Sigma_real[1, 2]

var(X)
cor(X)
cov2cor(Sigma_real)

# Complex normal
set.seed(21032025)
S_comp <- matrix(rnorm(P*P) + 1i*rnorm(P*P), ncol = P)
Sigma_comp <- S_comp %*% t(Conj(S_comp)) + diag(P)
isSymmetric(Sigma_comp)
eigen(Sigma_comp)$values
isHermitian(Sigma_comp)
ishpd(Sigma_comp)

mu_comp <- rep(0, P)

Z <- cmvnorm::rcmvnorm(n, mean = mu_comp, sigma = Sigma_comp)

qqnorm(Re(Z[, 1]))
qqline(Re(Z[, 1]))

qqnorm(Im(Z[, 1]))
qqline(Im(Z[, 1]))

# A statistic of a CMV should follow a certain chi-square distribution
Pmat <- Conj(solve(Conj(Sigma_comp)))
chistat <- rep(NA, n) 

for (i in 1:n) {
    Zi <- matrix(Z[i, ], ncol = 1)
    chistat[i] <- 2 * Re( t(Conj(Zi - mu_comp)) %*% Pmat %*% (Zi - mu_comp) )
}

plot(density(chistat))
curve(dchisq(x, 2*P), from = 0, to = 50, 
      add = TRUE, col = "red")


# comparing cmvnorm CN distributions --------------------------------------

# could I do something like a linear regression? compare the correlation of one column with another?

# Sigma_simp <- matrix(c(5, 3-1i, 3+1i, 4), ncol = 2, byrow = TRUE)
Sigma_simp <- Sigma_comp
t(Conj(Sigma_simp))
ishpd(Sigma_simp)
isSymmetric(Sigma_simp)
eigen(Sigma_simp)$values
eigen(Sigma_simp)$vectors

L <- t(Conj(Eigen_chol(Sigma_simp)))
L_inv <- solve(L)

Lbar <- t(Conj(Eigen_chol(Conj(Sigma_simp))))

cbind(L[, 1], Lbar[, 1])

set.seed(22032025)
Z <- rcmvnorm(1e5, sigma = Sigma_simp)

cmvnorm::var(Z)

cbind(cmvnorm::var(Z)[, 1], Sigma_simp[, 1])

# exactly what matrix do I need to multiply by here?
# t(L_inv) seems to give the covariance matrix I expect, but feels wrong somehow
# other matrices do not seem to give the right covariance matrices, though
Znorm <- Z %*% t(L_inv)

# with t(L_inv), this covariance matrix looks correct
cmvnorm::var(Znorm)

# the distribution of the Real and Im parts of cols should be N(0, sqrt(0.5))
plot(density(Im(Znorm[, 2])))
curve(dnorm(x, sd=sqrt(0.5)), -4, 4, 
      add = TRUE, col = "red")


# L and LT ----------------------------------------------------------------

Sigma_simp
all.equal(L %*% t(Conj(L)), Sigma_simp)

t(L) %*% Conj(L)

ev <- eigen(Sigma_simp, symmetric = TRUE)
jj <- svd(Sigma_simp)

ev$vectors
jj$u

D_ev <- diag(sqrt(ev$values), ncol=length(ev$values))
retval_ev <- emulator::quad.tform(D_ev, ev$vectors)

all.equal(retval_ev %*% t(Conj(retval_ev)), Sigma_simp)

ev$values
jj$d

D_svd <- diag(sqrt(jj$d), ncol=length(jj$d))
retval_svd <- quad.tform(D_svd, jj$u)

all.equal(retval_svd %*% t(Conj(retval_svd)), Sigma_simp)

all.equal(retval_svd %*% t(Conj(retval_svd)), Sigma_simp)

# quad.tform()
# claim that quad.tform(M, x) = x M x^H

# actual code:
# tcrossprod(x, tcrossprod(Conj(x), M))

# where:
# x %*% t(y) (tcrossprod)

mvtnorm::rmvnorm()
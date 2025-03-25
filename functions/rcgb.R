# sample from complex matrix Bingham distribution
# as in Abdallah et al. 2020

library(rstiefel)
# library(cmvnorm)
library(mvtnorm)

# fr <- function(b, betas) {
#     return(abs(sum(1 / (b + 2*betas)) - 1))
# }

# from Kent et al. 2016 supplementary materials
bfind <- function(lambda) {
    q <- length(lambda)
    fb <- function(b) 1 - sum(1/(b + 2*lambda))
    if (sum(lambda^2) == 0) b0 <- q
    else b0 <- uniroot(fb, interval = c(1, q))$root
    
    return(b0)
}

# rvB - sample from the real vector Bingham distribution ------------------

# input:
# A, NxN symmetric matrix parameter
# 
# output:
# u, Nx1 real-valued unit vector
rvB <- function(A, opt_upper = 1e5) {
    # based on Algorithm 4 of Abdallah 2020
    
    mat1 <- A
    N <- nrow(A)
    
    # Line 1: compute gamma
    gamm <- max(eigen(mat1)$values)
    
    # Preparation for Line 2
    mat2 <- diag(gamm, nrow = N) - mat1 
    betas <- eigen(mat2)$values
    
    # TODO can I solve this in another way?
    # TODO I should perform some sort of check that determines where I need to use a larger "upper" value
    # Line 2: solve for b
    # b <- optim(0.1, fr, method = "Brent", lower = 0, upper = opt_upper, 
    #            betas = betas)$par
    b <- bfind(betas)
    
    # Line 4: compute Omega
    Om <- diag(N) + (2/b) * mat2
    # Line 3: compute Mstar, as a function of b
    log_Mst <- .5 + gamm - (N - b)/2 +
        log((N/b)^(N/2)) + log(det(Om)^-.5)
    
    # Line 5: sample u by acceptance-rejection scheme    
    acc_prob <- 0
    u_draw <- 1
    while(u_draw >= acc_prob) {
        # Line 6: sample y from a Normal dist.
        y <- rmvnorm(1, sigma = solve(Om))
        # Line 7: normalize y
        u <- t(y) / sqrt(sum(y^2))
        # Line 8: compute the ACG density value at u
        log_f_ACG <- log(det(Om)^.5 * (t(u) %*% Om %*% u)^(-N/2))
        # Line 9: compute the vector Bingham density value at u
        log_f_vBL <- t(u) %*% A %*% u
        
        # Lines 10-11: accept u with the calculated probability
        u_draw <- runif(1)
        acc_prob <- exp(log_f_vBL - log_Mst - log_f_ACG)
    }
    
    # Line 12: return u, after the acceptance criteria is met
    return(u)
}

# rcvb - using the vBL, sample from complex VGBL --------------------------

rcvb <- function(A, opt_upper = 1e5) {
    AA <- rbind(cbind(Re(A), -Im(A)),
                cbind(Im(A), Re(A)))
    
    uu <- rvB(AA, opt_upper = opt_upper) |> matrix(ncol = 2)
    
    return(complex(real = uu[, 1], imaginary = uu[, 2]))
}


# rcvb_LN - a function to sample new column, orthogonal to others ---------

# input:
# U, an NxP matrix with orthonormal columns, for which to sample a column
# p, a scalar, the index of the column to sample
# A, NxN psd Hermitian matrix, the original matrix parameter
# B, PxP diagonal matrix, or a Px1 vector representing the diagonal entries of B

# output:
# u, Nx1 vector, to be inserted into matrix U

# Description

# This function samples a single column for U, so that the sampled column is
# orthogonal to the other columns in U. The sample column, u_p, has a vector
# Bingham distribution with parameters (A, B), restricted to the space that is
# orthogonal to the columns U[, -p]. A transformed version of the column is
# sampled, and then this result is transformed back to the originality
# dimensionality, using a basis of the left null space of U[, -p], ensuring u_p
# is orthogonal to the rest of the columns in U.
rcvb_LN <- function(U, j, A) {
    
    LN <- NullC(U[, -j])
    # this value is actually not used later
    u <- Conj(t(LN)) %*% U[, -j]
    
    Abar <- Conj(t(LN)) %*% A %*% LN
    
    u <- LN %*% rcvb(Abar)
    return(u)
}

# rcmb - sample from the random complex matrix Bingham distribution -------
# Input:

# - U, an NxP complex-valued matrix with orthonormal columns
# - A, NxN psd Hermitian matrix
# - B, PxP diagonal matrix, or a Px1 vector representing the diagonal entries of B

# Output:

# - U, an NxP complex-valued matrix with orthonormal columns, sampled from the
# complex matrix Bingham distribution

rcmb <- function(U, A, B) {
    if (nrow(U) <= ncol(U)) stop("nrow(U) must be greater than ncol(U)")
    
    if (is.matrix(B)) {
        B <- diag(B)
    }
    
    # over the columns of U, in a random order:
    js <- sample(1:ncol(U))
    for (j in js) {
        U[, j] <- rcvb_LN(U, j, B[j]*A)
    }
    
    return(U)
}


# rcBingUP ----------------------------------------------------------------

# Description: sample of a small square complex Bingham distribution on U(P)

# Input:

# - A, mxm Hermitian p.d. matrix
# - B, mxm diagonal real matrix with decreasing entries

# Output:

# - X, mxm, a sample from the complex Bingham(A, B) distribution
# - nrej, the number of rejections

# Details:

# This is a rejection-based sampler based on Hoff 2009. Currently, it works well
# for 2x2 matrices, but rejects too often for 3x3 matrices or larger.

# References:

# Hoff, P. D. (2009). Simulation of the Matrix Bingham–von Mises–Fisher Distribution, With Applications to Multivariate and Relational Data. Journal of Computational and Graphical Statistics, 18(2), 438–456. https://doi.org/10.1198/jcgs.2009.07177

rcBingUP <- function(A, B) {
    
    stopifnot("A and B must have the same dimensions" = all(dim(A) == dim(B)))
    stopifnot("A and B must be square" = dim(A)[1] == dim(A)[2])
    
    P <- nrow(A)
    
    ################################
    # original scaling of A and B from supplemental paper
    ################################
    
    # # rescale A and B
    # # TODO verify the 2 points below
    # # - shifting the diagonals does not change the distribution
    # # - shifting A and B by c, 1/c does not change the distribution
    # # shifting
    # diag(A) <- diag(A) - min(eigen(A, only.values = TRUE)$values)
    # diag(B) <- diag(B) - min(diag(B))
    # # scaling
    # # maximum eigenvalue of A is the first one
    # mA <- eigen(A, only.values = TRUE)$values[1]
    # mB <- max(B)
    # # TODO trying a different factor here
    # # gd <- (dim(A)[1] + 1) / (2*mB) + mA # this is Hoff
    # # gd <- (dim(A)[1] + 1) / (mB) + mA # this seems to lead to better acceptance rates for 2x2, but maybe worse, or at least gets stuck at one, for 3x3?
    # gd <- (dim(A)[1] + 1) / (mB) + mA
    # # TODO understand why this factor is calculated now
    # del1 <- max(eigen(A, only.values = TRUE)$values) + 0.5
    # gam <- gd/del1
    # 
    # # scale
    # A <- A / gam
    # B <- B * gam
    # 
    # As[, , i] <- A
    # Bs[, , i] <- B
    # 
    # # get and store the eigenvalues of A
    # Aevals <- eigen(A, only.values = TRUE)$values
    # # TODO calculate del
    # # del <- max(Aevals) + 0.5
    # 
    # # TODO trying a different nu
    # nu <- (dim(A)[1] + 1) # this is Hoff
    # # nu <- dim(A)[1]
    
    ################################
    # end of original scaling of A and B from supplemental paper
    ################################
    
    ################################
    # alternative scaling of A and B from rstiefel
    ################################
    
    b <- diag(B)
    bmx <- max(b)
    bmn <- min(b)
    
    A <- A*(bmx-bmn)
    b <- (b-bmn) / (bmx-bmn)
    
    Aevals <- eigen(A)$val
    diag(A) <- diag(A) - Aevals[1]
    Aevals <- eigen(A)$val
    
    nu <- max(dim(A)[1]+1,round(-Aevals[length(Aevals)]))
    del1 <- nu/2
    # TODO right now, I am also dividing this by 2 again in the sample of W
    # perhaps I should just make this, divide by 4?
    S <- solve( diag(del1,nrow=dim(A)[1] ) - A )/2
    
    ################################
    # end of alternative scaling of A and B from rstiefel
    ################################
    
    
    # S <- solve(diag(del1, nrow = P) - A)
    # TODO I am rounding this matrix so that it is Hermitian, within a tolerance that rcwis will accept.
    tol_digits <- (.Machine$double.eps * 100) |> log10() |> ceiling() |> abs()
    S <- round(S, tol_digits-1)
    
    rej <- TRUE
    nrej <- 0
    while (rej) {
        # browser()
        # W <- rcwis(nu, Conj(S)/2)
        W <- rcomplex_wishart(nu, P, S/2)
        Weigen <- eigen(W)
        
        # TODO do microbenchmark - which is faster? multiplying by a diagonal matrix of +- 1, or using vector recycling?
        X <- Weigen$vectors
        L <- Weigen$values
        
        X <- X %*% diag((-1)^rbinom(dim(A)[1], 1, .5)) 
        
        D <- sort(diag(B) - L, decreasing = TRUE)
        
        # the products must be over elements of decreasing order. Eigenvalues
        # from eigen are in decreasing order by default. D is sorted above.
        lambda <- sum(Aevals * D)
        
        # calculate the log acceptance ratio
        # part of it might be slightly complex, due to numerical issues
        # so, use Re
        # TODO but, I should double-check everything, to justify that it should be real
        # TODO double-check I am using right B-L matrix, not necessarily the sorted D
        logr <- Re(sum(diag(
            diag(B - L) %*% t(Conj(X)) %*% A %*% X))) - 
            lambda
        
        # Accept with probability r
        # draw u from Unif(0, 1)
        # TODO can I compare, say, log(u) to the logr? 
        rej <- log(runif(1)) > logr
        nrej <- nrej + rej
    }
    
    return(list(X = X, nrej = nrej))
}


# rcBingUP_gibbs ----------------------------------------------------------

# Description: sample a square complex Bingham distribution

# Input:

# X, mxm unitary matrix
# A, mxm Hermitian positive definite matrix
# B, mxm diagonal real matrix with decreasing elements

# Output:

# newX, mxm unitary matrix

# Details:

# This function is intended to be used within a Gibbs sampler, where the FCD of
# X is a complex Bingham (A, B) distribution. The columns, taken two at a time,
# can be related to a 2x2 matrix z which has a complex Bingham distribution.

# Currently, this function requires that the input X have even dimensions.

rcBingUP_gibbs <- function(X, A, B) {
    
    stopifnot("A and B must have the same dimensions" = all(dim(A) == dim(B)))
    stopifnot("A and B must be square" = dim(A)[1] == dim(A)[2])
    
    P <- nrow(A)
    
    stopifnot("A and B must have even dimensions" = P %% 2 == 0)
    
    # generate a random order of columns
    scols <- sample(1:P)
    
    for (sstep in 1:(P/2)) {
        # get the random column indices for this step, and put in increasing order. 
        ijs <- scols[c(sstep*2 - 1, sstep*2)] |> sort()
        
        # find an orthonormal basis for the left null space of X without ijs 
        N <- NullC(X[, -ijs])
        
        # transform the original parameters
        newA <- t(Conj(N)) %*% A %*% N
        newB <- diag(B[ijs])
        
        # browser()
        # sample the 2x2 columns, and transform back into X columns
        zsamp <- rcBingUP(newA, newB)
        X[, ijs] <- N %*% zsamp$X
    }
    
    return(X)
}

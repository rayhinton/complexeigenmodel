# Utility functions


# sample_gumbel -----------------------------------------------------------

# Description: draw a sample from a discrete random variable using log-scale weights that are possibly unnormalized

# Input: 

# thetas - values of the discrete random variables
# size - number of sample draws to return
# logweights - logarithm of the discrete probabilities (possibly unnormalized)

# Output:

# outsample - values from thetas sampled according to logweights

# Credit:

# Class notes (STAT 633). Dr. Anirban Bhattacharya. Fall 2024.

sample_gumbel <- function(thetas, size, logweights) { 
    stopifnot("length of thetas and logweights must be equal" = 
                  length(thetas) == length(logweights))
    
    # number of discrete values
    J <- length(thetas)
    # initialize a vector that will contain the sample
    outsample <- vector(mode(thetas), size)
    
    for (i in 1:size) {
        # generate Gumbel noise
        Gj <- -log(rexp(J))
        
        # calculate Gumbel weights = logweights + Gumbel noise
        gumweights <- logweights + Gj
        
        # perform the sampling step: return the theta corresponding to the
        # largest Gumbel weight
        outsample[i] <- thetas[which.max(gumweights)]        
    }
    
    return(outsample)
}

# rscnorm -----------------------------------------------------------------

# Description: sample from the univariate standard complex normal distribution

rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

# rFTCW -------------------------------------------------------------------

rFTCW <- function(Sigma, n, a, useEigenR = FALSE, byCholesky = FALSE) {
    W <- rcomplex_wishart(n, Sigma, useEigenR, byCholesky)
    Y <- a * W / Re(sum(diag(W)))
    return(Y)
}

# rCMACG ------------------------------------------------------------------

rCMACG <- function(nrow, ncol, Sigma) {
    # simulate X ~ CMN(0, I_P, I_d),
    # calculate Z by transforming X with Cholesky decomposition of Sigma0
    # let H_Z be orthogonal portion of polar decomposition of Z
    X <- matrix(rscnorm(nrow*ncol), nrow, ncol)
    
    # check if this is the upper or lower triang version
    # want the matrix such that Sigma = cholSigma %*% t(Conj(cholSigma))
    # Eigen_chol returns the upper triangular portion
    cholSigma <- EigenR::Eigen_chol(Sigma)
    
    Z <- t(Conj(cholSigma)) %*% X 
    sqrtZHZ <- EigenR::Eigen_sqrt( t(Conj(Z)) %*% Z )
    return(Z %*% solve(sqrtZHZ))
}

# runif_stiefel -----------------------------------------------------------

runif_stiefel <- function(P, d, betaf) {
    # real
    if (betaf == 1) {
        X1 <- matrix(rnorm(P*d), ncol = d)
        # complex
    } else if (betaf == 2) {
        X1 <- matrix(rscnorm(P*d), ncol = d)        
    }
    
    # take QR decomposition - not guaranteed to be unique due to numerical methods
    X1qr <- qr(X1)
    QX1 <- qr.Q(X1qr)
    # extract sign of the diagonal elements of R
    D <- sign(Re(diag(qr.R(X1qr))))
    # transform by Q = QS, where S = diag(D). this transformation guarantees Q is unique.
    
    # Qfinal <- t(t(QX1) * D) # faster than matrix multiplication, for larger matrices
    Qfinal <- QX1 %*% diag(D)
    
    return(Qfinal)
}

# uni.slice, univariate slice sampling ------------------------------------

uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL, ...)
{
    # Check the validity of the arguments.
    
    if (!is.numeric(x0) || length(x0)!=1
        || !is.function(g) 
        || !is.numeric(w) || length(w)!=1 || w<=0 
        || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
        || !is.numeric(lower) || length(lower)!=1 || x0<lower
        || !is.numeric(upper) || length(upper)!=1 || x0>upper
        || upper<=lower 
        || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
    { 
        stop ("Invalid slice sampling argument")
    }
    
    # Keep track of the number of calls made to this function.
    
    # uni.slice.calls <<- uni.slice.calls + 1
    
    # Find the log density at the initial point, if not already known.
    
    if (is.null(gx0)) 
    { 
        # uni.slice.evals <<- uni.slice.evals + 1
        gx0 <- g(x0, ...)
    }
    
    # Determine the slice level, in log terms.
    
    logy <- gx0 - rexp(1)
    
    # Find the initial interval to sample from.
    
    u <- runif(1,0,w)
    L <- x0 - u
    R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
    
    # Expand the interval until its ends are outside the slice, or until
    # the limit on steps is reached.
    
    if (is.infinite(m))  # no limit on number of steps
    { 
        repeat
        { if (L<=lower) break
            # uni.slice.evals <<- uni.slice.evals + 1
            if (g(L, ...)<=logy) break
            L <- L - w
        }
        
        repeat
        { if (R>=upper) break
            # uni.slice.evals <<- uni.slice.evals + 1
            if (g(R, ...)<=logy) break
            R <- R + w
        }
    }
    
    else if (m>1)  # limit on steps, bigger than one
    { 
        J <- floor(runif(1,0,m))
        K <- (m-1) - J
        
        while (J>0)
        { if (L<=lower) break
            # uni.slice.evals <<- uni.slice.evals + 1
            if (g(L, ...)<=logy) break
            L <- L - w
            J <- J - 1
        }
        
        while (K>0)
        { if (R>=upper) break
            # uni.slice.evals <<- uni.slice.evals + 1
            if (g(R, ...)<=logy) break
            R <- R + w
            K <- K - 1
        }
    }
    
    # Shrink interval to lower and upper bounds.
    
    if (L<lower) 
    { L <- lower
    }
    if (R>upper)
    { R <- upper
    }
    
    # Sample from the interval, shrinking it on each rejection.
    
    repeat
    { 
        x1 <- runif(1,L,R)
        
        # uni.slice.evals <<- uni.slice.evals + 1
        gx1 <- g(x1, ...)
        
        if (gx1>=logy) break
        
        if (x1>x0) 
        { R <- x1
        }
        else 
        { L <- x1
        }
    }
    
    # Return the point sampled, with its log density attached as an attribute.
    
    attr(x1,"log.density") <- gx1
    return (x1)
    
}

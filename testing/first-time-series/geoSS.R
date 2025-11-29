geodesic_CompSt <- function(Y, H, th) {
    n <- nrow(Y)
    p <- ncol(Y)
    
    K <- ( diag(n) - Y %*% t(Conj(Y)) ) %*% H
    qr_result <- qr(K)
    Q <- qr.Q(qr_result)
    R <- qr.R(qr_result)
    
    # Ensure R has positive diagonal entries
    signs <- sign(Re(diag(R)))
    signs[signs == 0] <- 1
    Q <- Q %*% diag(signs)
    R <- diag(signs) %*% R
    
    A <- t(Conj(Y)) %*% H
    
    # form a 2p x 2p matrix:
    toexp <- rbind(
        cbind(A, -t(Conj(R))),
        cbind(R, matrix(0, p, p))
    )
    
    # TODO change this to just selecting certain rows or columns
    MtNt <- EigenR::Eigen_exp(th * toexp) %*% diag(1, 2*p, p)
    Mt <- MtNt[1:p, ]
    Nt <- MtNt[(p+1) : (2*p), ]
    
    Yt <- Y %*% Mt + Q %*% Nt
}

step_out <- function(Y, H, logt, w, m, logd_func, ...) {
    u <- runif(1, 0, w)
    l <- -u
    r <- l + w
    iota <- sample(1:m, 1)
    
    i <- 2
    j <- 2
    
    Ytl <- geodesic_CompSt(Y, H, l)
    # logdens <- logd_CMACG(Ytl, Si)
    logdens <- logd_func(Ytl, ...)
    while (i <= iota & logdens > logt) {
        l <- l - w
        i <- i + 1
        Ytl <- geodesic_CompSt(Y, H, l)
        # logdens <- logd_CMACG(Ytl, Si) 
        logdens <- logd_func(Ytl, ...)
    }
    
    # evaluate logdens at th = r
    Ytr <- geodesic_CompSt(Y, H, r)
    # logdens <- logd_CMACG(Ytr, Si) 
    logdens <- logd_func(Ytr, ...)
    while ((j <= m + 1  - iota) & (logdens > logt)) {
        r <- r + w
        j <- j + 1
        Ytr <- geodesic_CompSt(Y, H, r)
        # logdens <- logd_CMACG(Ytr, Si) 
        logdens <- logd_func(Ytr, ...)
    }
    
    return(c(l, r))
}

shrinklr <- function(Y, H, logt, lr, logd_func, ...) {
    l <- lr[1]
    r <- lr[2]
    
    thh <- runif(1, 0, r - l)
    theta <- thh - (thh > r) * (r - l)
    thmin <- thh
    thmax <- thh
    
    Ytheta <- geodesic_CompSt(Y, H, theta)
    # logdens <- logd_CMACG(Ytheta, Si)
    logdens <- logd_func(Ytheta, ...)
    while(logdens <= logt) {
        if (thmin <= thh & thh < (r - l)) {
            thmin <- thh
        } else {
            thmax <- thh
        }
        
        thh <- runif_union(thmax, thmin, r, l)
        theta <- thh - (thh > r) * (r - l)
        
        Ytheta <- geodesic_CompSt(Y, H, theta)
        # logdens <- logd_CMACG(Ytheta, Si)
        logdens <- logd_func(Ytheta, ...)
    }
    
    return(theta)
}

drawV <- function (Y, n, p) {
    # H is generated from uniformly from the tangent space at Y
    # isomorphic to sampling on the unit sphere
    
    # simulate p*(p-1)/2 + p*(n-p) standard complex normal random variables
    Tdim <- p*(p-1)/2 + p*(n-p)
    # v <- rnorm(Tdim, 0, 1/sqrt(2)) + 1i * rnorm(Tdim, 0, 1/sqrt(2))
    # v <- v / sqrt(sum(v * Conj(v)))
    
    v <- rnorm(2*Tdim)
    v <- v / sqrt(sum(v^2))
    v <- v[1:Tdim] + 1i*v[(Tdim+1):(2*Tdim)]
    
    norm_error <- abs(1 - sqrt(sum( Mod(v)^2 )))
    
    if (norm_error > 1e-10) {
        print(sqrt(sum( Mod(v)^2 )))
        stop("generated v does not have unit norm")
    }
    
    # make a pxp skew Hermitian matrix from the first p*(p-1)/2 of them
    Pi <- matrix(0, p, p)
    idx <- 1
    for(j in 1:(p-1)) {
        for(i in (j+1):p) {
            Pi[i,j] <- v[idx]
            Pi[j,i] <- -Conj(v[idx])
            idx <- idx + 1
        }
    }
    
    # put the remaining p*(n-p) into a (n-p) x p matrix 
    Sigma <- matrix(v[(p*(p-1)/2 + 1):Tdim], nrow = n-p, ncol = p, byrow = TRUE)
    
    Yperp <- qr.Q(qr(Y), complete = TRUE)[, (p+1):n]
    
    H <- Y %*% Pi + Yperp %*% Sigma
    
    return(H)
}

geoSS <- function(Y, logd_func, w, m, ...) {
    logdens <- logd_func(Y, ...)
    # draw Unif(0, p), where p = exp(logdens)
    logt <- log(runif(1)) + logdens
    
    H <- drawV(Y, nrow(Y), ncol(Y))
    
    lr <- step_out(Y, H, logt, w, m, logd_func, ...)
    theta <- shrinklr(Y, H, logt, lr, logd_func, ...)
    
    newY <- geodesic_CompSt(Y, H, theta)
    
    return(newY)
}

# other functions ---------------------------------------------------------

# sample from a union of two intervals

runif_union <- function(theta_max, theta_min, r, l) {
    # Define intervals as a matrix (each row is [start, end])
    intervals <- rbind(c(0, theta_max),
                       c(theta_min, r - l))
    
    # Calculate interval lengths
    lengths <- intervals[,2] - intervals[,1]
    total_length <- sum(lengths)
    
    # Sample: first pick which interval (weighted by length), then uniform 
    # within it
    interval_idx <- sample(1:nrow(intervals), size = 1, 
                           prob = lengths/total_length)
    theta <- runif(1, min = intervals[interval_idx, 1], 
                   max = intervals[interval_idx, 2])
    
    return(theta)
}
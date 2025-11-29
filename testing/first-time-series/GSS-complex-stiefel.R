# geodesic function

source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/matrix_distances.R")
source("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")

logd_CGB <- function(X, A, B, V) {
    result <- Re(sum(diag(
       B %*% t(Conj(X)) %*% V %*% A %*% t(Conj(V)) %*% X 
    )))
    
    return(result)
}

# TODO write a function to evaluate a (proportional) log density
# how about starting with a simple density, like the complex Bingham?
# CMACG, since it is easy for me to actually sample from
logd_CMACG <- function(HZ, P) {
    m <- nrow(HZ)
    r <- ncol(HZ)
    result <- -m * log(EigenR::Eigen_det(t(Conj(HZ)) %*% solve(P) %*% HZ))
    return(Re(result))
} 

# Y, n x p
# H, n x p

# Q, n x p
# R, p x p

# A, p x p

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

# setup -------------------------------------------------------------------

n <- 8
p <- 2
logt <- log(0.2)

w <- 7
m <- 7

# parameter for CMACG
Si <- diag(n:1)

SiV <- eigen(Si)$vectors[, 1:p]

dataseed <- 21102025

# parameters for Complex Bingham distribution
set.seed(123)
Bb <- diag(p:1)
Aa <- diag(n:1)
Vv <- 1i * matrix(rnorm(n*n), n, n)
Vv <- qr.Q(qr(Vv), complete = TRUE)

# rcmb <- function(U, A, B) {

# generate data, Y and tangent direction H --------------------------------

set.seed(dataseed)
Y <- 1i * matrix(rnorm(n*p), n, p)
Y <- qr.Q(qr(Y), complete = TRUE)[, 1:p]
newy <- Y

H <- drawV(newy, n, p)

# norm of H should be 1
sum(diag( t(Conj(H)) %*% (diag(n) - .5* newy %*% t(Conj(newy))) %*% H ))

# check this is skew Hermitian
t(Conj(newy)) %*% H

logd_CGB(newy, Aa, Bb, Vv)
logd_CMACG(newy, Si)

# sample with Geodesic slice sampling -------------------------------------

# lr <- step_out(newy, H, logt, w, m)
# theta <- shrinklr(newy, H, logt, lr)

its <- 5000

Ys_geo <- array(NA, c(n, p, its))

Ys_proj_geo <- array(NA, c(n, n, its))
for (i in 1:its) {
    
    # logdens <- logd_CGB(newy, Aa, Bb, Vv)
    
    # logdens <- logd_CMACG(newy, Si)
    # 
    # # draw Unif(0, p), where p = exp(logdens)
    # logt <- log(runif(1)) + logdens
    # 
    # H <- drawV(newy, n, p)
    # 
    # # lr <- step_out(newy, H, logt, w, m, logd_CGB, A = Aa, B = Bb, V = Vv)
    # # theta <- shrinklr(newy, H, logt, lr, logd_CGB, A = Aa, B = Bb, V = Vv)    
    # 
    # lr <- step_out(newy, H, logt, w, m, logd_CMACG, P = Si)
    # theta <- shrinklr(newy, H, logt, lr, logd_CMACG, P = Si)
    # 
    # newy <- geodesic_CompSt(newy, H, theta)
    
    newy <- geoSS(newy, logd_CMACG, w, m, P = Si)
    Ys_geo[, , i] <- newy
    
    Ys_proj_geo[, , i] <- newy %*% t(Conj(newy))
}

# sample the same distribution with a different algorithm -----------------

Ys_r <- array(NA, c(n, p, its))
Ys_proj_r <- array(NA, c(n,n, its))

CB_A <- Vv %*% Aa %*% t(Conj(Vv))
thisy <- Y
for (i in 1:its) {
    # Ys_r[, , i] <- thisy
    # Ys_proj_r[, , i] <- thisy %*% t(Conj(thisy))
    # thisy <- rcmb(newy, CB_A, Bb)

    Ys_r[, , i] <- rCMACG(n, p, Si)
    Ys_proj_r[, , i] <- Ys_r[, , i] %*% t(Conj(Ys_r[, , i]))
}

# calculate average quantities --------------------------------------------

##### average X and average projection matrix from Geodesic slice sampling
# average X
avgYs_geo <- apply(Ys_geo[, , (its/2) : its],
                   c(1, 2), mean)
avgYs_geo <- avgYs_geo %*% 
    solve( EigenR::Eigen_sqrt( t(Conj(avgYs_geo)) %*% avgYs_geo ) )
# average projection matrix
avgYs_proj_geo <- (apply(Ys_proj_geo[, , (its/2) : its],
                         c(1, 2), mean) |> eigen())$vectors[, 1:p]
avgYs_proj_geo <- avgYs_proj_geo %*% t(Conj(avgYs_proj_geo))
# calculate some distances, for trace plotting
Ys_geo_dists <- rep(NA, its)
for (i in 1:its) {
    # Ys_geo_dists[i] <- fast_evec_Frob_stat(Vv[, 1:p], Ys_geo[, , i])
    
    Ys_geo_dists[i] <- fast_evec_Frob_stat(SiV, Ys_geo[, , i])
}

##### average X and average projection matrix from alternative sampler
# average X
avgYs_r <- apply(Ys_r,
                 c(1, 2), mean)
avgYs_r <- avgYs_r %*% 
    solve( EigenR::Eigen_sqrt( t(Conj(avgYs_r)) %*% avgYs_r ) )
# average projection matrix
avgYs_proj_r <- (apply(Ys_proj_r[, , (its/2) : its],
                       c(1, 2), mean) |> eigen())$vectors[, 1:p]
avgYs_proj_r <- avgYs_proj_r %*% t(Conj(avgYs_proj_r))
# calculate some distances, for trace plotting
Ys_r_dists <- rep(NA, its)
for (i in 1:its) {
    # Ys_r_dists[i] <- fast_evec_Frob_stat(Vv[, 1:p], Ys_r[, , i])
    
    Ys_r_dists[i] <- fast_evec_Frob_stat(SiV, Ys_r[, , i])
}

# plot summaries and display individual columns for comparison ------------

# trace plots of distances
plot(Ys_geo_dists, type = "l")
plot(Ys_r_dists, type = "l")

# summary of distances
summary(Ys_geo_dists[(its/2) : its])
summary(Ys_r_dists[(its/2) : its])

# show some columns of the average projection matrices
j <- 1
cbind(avgYs_proj_r[, j], avgYs_proj_geo[, j]) |> round(2)
# Frobenius distance of the average projection matrices
norm(avgYs_proj_r - avgYs_proj_geo, "F")


# more comparisons of the average X matrices; however, I don't think these are
# very meaningful, because of the antipodal nature of the Bingham (and CMACG?)
# distributions.

# cbind(avgYs_r[, 1], avgYs_geo[, 1]) |> round(2)
# close_avgs <- evec_Frob_stat(avgYs_r, avgYs_geo, returnMats = TRUE)
# 
# close_avgs$dist_obj
# cbind(close_avgs$Xopt[, 1], close_avgs$Yopt[, 1]) |> round(2) 
# 
# cbind(
#     (close_avgs$Xopt %*% t(Conj(close_avgs$Xopt)))[, 8],
#     (close_avgs$Yopt %*% t(Conj(close_avgs$Yopt)))[, 8]) |> round(2)
# 
# rcmb_and_Vv <- evec_Frob_stat(avgYs_r, Vv[, 1:p], returnMats = TRUE)
# rcmb_and_Vv$dist_obj
# 
# cbind(rcmb_and_Vv$Xopt[, 1], Vv[, 1]) |> round(2)
# 
# geo_and_Vv <- evec_Frob_stat(avgYs_geo, Vv[, 1:p], returnMats = TRUE)
# geo_and_Vv$dist_obj
# 
# cbind(geo_and_Vv$Xopt[, 1], geo_and_Vv$Yopt[, 1]) |> round(2)

# check real matrix Bingham distribution ----------------------------------

set.seed(24102025)
Ar <- rWishart(1, 5, diag(4))[, , 1]
Vr <- eigen(Ar)$vectors
Ar_Evals <- diag(eigen(Ar)$values)
Br <- diag(4:1)

# X <- rstiefel::rustiefel(4, 4)
# rstiefel::rbing.matrix.gibbs(X, Ar, Br)

Bingits <- 10000
Bings <- array(NA, c(4, 4, Bingits))

Ar_Bings <- array(NA, c(4, 4, Bingits))

for (i in 1:Bingits) {
    if (i %% 500 == 0) {print(i)}
    thisB <- rstiefel::rbing.Op(Ar, Br)
    Bings[, , i] <- thisB
    Ar_Bings[, , i] <- thisB %*% Ar_Evals %*% t(thisB)
}

avgBings <- apply(Bings,
                 c(1, 2), mean)
avgBings <- avgBings %*% 
    solve( EigenR::Eigen_sqrt( t(Conj(avgBings)) %*% avgBings ) )

cbind(Vr[, 1], avgBings[, 1])

avgBings |> round(2)
Vr |> round(2)

evec_Frob_stat(avgBings, Vr, returnMats = TRUE)

avg_ArBings <- apply(Ar_Bings, c(1, 2), mean)
cbind(Ar[, 1], avg_ArBings[, 1])
round(Ar, 2)
round(avg_ArBings, 2)
norm(Ar - avg_ArBings, "F")

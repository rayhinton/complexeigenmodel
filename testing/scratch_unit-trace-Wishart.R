# Unit trace Wishart


# functions ---------------------------------------------------------------

mvGamma <- function(c, m, betaf, logscale = FALSE) {
    term1 <- .25*m*(m-1)*betaf*log(pi)
    term2 <- sum(lgamma(c - ((1:m)-1)/2 * betaf))
    
    if (logscale) {
        return(term1 + term2)
    } else {
        return(exp(term1 + term2))
    }
}

# generate a Uniform unit trace sym. p.d. matrix

runif_unittr_spd <- function(num_samples, n) {
    
    ls <- 1:(0.5*n*(n+1) - 1)
    is <- 1:(n-1)
    
    pl <- rep(0, length(ls))
    kappal <- rep(0, length(ls))
    lambdal <- rep(0, length(ls))
    ql <- rep(0, length(ls))
    cl <- rep(0, length(ls))
    Intl <- rep(pi, length(ls))
    mul <- rep(pi/2, length(ls))
    dl <- rep(1, length(ls))
    
    for (l in ls) {
        if (any( (0.5 * is * (is + 1)) == l )) {
            il <- is[(0.5 * is * (is + 1)) == l]
            # pl is verified correct by Table II
            pl[l] <- n + 1 - il
        }
        
        for (i in 1:(n-1)) {
            for (m in 0:i) {
                if (l == .5*i*(i + 1) + m) {
                    # print(c(l, .5*i*(i + 1) + m, i, m))
                    # kappal is verified correct by Table II
                    kappal[l] <- i - 1
                    # lambdal is verified correct by Table II
                    lambdal[l] <- i + 1 + m
                    # ql is verified correct by Table II
                    ql[l] <- n^2 - kappal[l]*n - lambdal[l]
                }
            }
        }
        
        if (any( (0.5 * is * (is + 1)) == l )) {
            # cl is verified correct by Table II
            cl[l] <- 2 * beta(.5*pl[l] + .5, .5*ql[l] + .5)^-1
            
            # TODO cannot verify mul or dl against a table
            mul[l] <- atan(sqrt(ql[l] / pl[l]))
            dl[l] <- (1 + ql[l]/pl[l])^(pl[l]/2) / 
                (1 + pl[l]/ql[l])^(-ql[l]/2)
        } else {
            # cl is verified correct by Table II
            cl[l] <- pi^-.5 * gamma(.5*ql[l] + 1) / gamma(.5*ql[l] + .5)
        }
    }
    
    As <- array(NA, c(n, n, num_samples))
    
    for(s in 1:num_samples) {
    
        yl <- rbeta(length(ql), (ql + 1)/2, (pl + 1)/2)
        bl <- rep(0, length(yl))
        (bl[pl == 0] <- rbinom(sum(pl == 0), 1, 0.5))
        
        phil <- (1 - bl)*yl + bl*(pi - yl)
        
        ks <- 1:(.5*n*(n+1))
        xk <- rep(0, length(ks))
        
        # Continue Algorithm 1
        # compute x from (11)
        
        for (k in ks) {
            if (k == 1) {
                xk[k] <- 1    
            } else {
                xk[k] <- prod(sin(phil[1:(k-1)]))
            }
        }
        xk <- xk * c(cos(phil), 1)
        
        # With vector x generate matrix U from (9).
        U <- matrix(0, n, n)
        
        Uind <- 1
        for (j in 1:n) {
            for (i in 1:j) {
                # print(c(i, j, Uind))
                # U[i, j] <- Uind
                U[i, j] <- xk[Uind]
                Uind <- Uind + 1
            }
        }

        As[, , s] <- t(U) %*% U
    } # end of for loop to generate Ais
    
    return(As)
}

# distribution of elements ------------------------------------------------

set.seed(30052025)
# V <- rWishart(1, 20, diag(4))[, , 1]
V <- diag(4)

nu <- 50

Ws <- rWishart(100000, nu, V)
cholWs <- array(NA, dim(Ws))

for (s in 1:dim(Ws)[3]) {
    Ws[, , s] <- Ws[, , s] / sum(diag(Ws[, , s]))
    cholWs[, , s] <- chol(Ws[, , s])
}

apply(Ws, c(1, 2), mean)
V / (sum(diag(V)))

summary(cholWs[1, 1, ])

plot(density(cholWs[1, 1, ]))
# curve(dbeta(x, nu-1 + 1, 1), from = 0, to = 1, add = TRUE, col = "red")

modew <- .50
kappaw <- 125
alphaw <- modew * (kappaw - 2) + 1
betaw <- (1 - modew) * (kappaw - 2) + 1

abline(v = modew, col = "red")
curve(dbeta(x, alphaw, betaw), from = 0, to = 1, n = 1001, add = TRUE, col = "red")


# Uniform matrices on constant trace manifold -----------------------------

n <- 4

# Algorithm 3: Generate a sample of phil
# - compute pl, ql with equations 29-31

ls <- 1:(0.5*n*(n+1) - 1)
is <- 1:(n-1)

pl <- rep(0, length(ls))
kappal <- rep(0, length(ls))
lambdal <- rep(0, length(ls))
ql <- rep(0, length(ls))
cl <- rep(0, length(ls))
Intl <- rep(pi, length(ls))
mul <- rep(pi/2, length(ls))
dl <- rep(1, length(ls))

for (l in ls) {
    if (any( (0.5 * is * (is + 1)) == l )) {
        il <- is[(0.5 * is * (is + 1)) == l]
        # print(c(l, il))
        # pl is verified correct by Table II
        pl[l] <- n + 1 - il
        
        Intl[l] <- pi/2 
    }
    
    for (i in 1:(n-1)) {
        for (m in 0:i) {
            if (l == .5*i*(i + 1) + m) {
                # print(c(l, .5*i*(i + 1) + m, i, m))
                # kappal is verified correct by Table II
                kappal[l] <- i - 1
                # lambdal is verified correct by Table II
                lambdal[l] <- i + 1 + m
                # ql is verified correct by Table II
                ql[l] <- n^2 - kappal[l]*n - lambdal[l]
            }
        }
    }
    
    if (any( (0.5 * is * (is + 1)) == l )) {
        # cl is verified correct by Table II
        cl[l] <- 2 * beta(.5*pl[l] + .5, .5*ql[l] + .5)^-1
        
        # TODO cannot verify mul or dl against a table
        mul[l] <- atan(sqrt(ql[l] / pl[l]))
        dl[l] <- (1 + ql[l]/pl[l])^(pl[l]/2) / 
            (1 + pl[l]/ql[l])^(-ql[l]/2)
    } else {
        # cl is verified correct by Table II
        cl[l] <- pi^-.5 * gamma(.5*ql[l] + 1) / gamma(.5*ql[l] + .5)
    }
}

cbind(pl, kappal, lambdal, ql, cl)
# c(131072/143/pi, 3003/2048, 60, 
#   128/35/pi, 35/32, 32/pi, 
#   3/4, 2/pi, 1/2)

yl <- rbeta(length(ql), (ql + 1)/2, (pl + 1)/2)
bl <- rep(0, length(yl))
(bl[pl == 0] <- rbinom(sum(pl == 0), 1, 0.5))

phil <- (1 - bl)*yl + bl*(pi - yl)
# TODO based on fix suggested by Gemini

ks <- 1:(.5*n*(n+1))
xk <- rep(0, length(ks))

# Continue Algorithm 1
# compute x from (11)

for (k in ks) {
    if (k == 1) {
        xk[k] <- 1    
    } else {
        xk[k] <- prod(sin(phil[1:(k-1)]))
    }
}
xk <- xk * c(cos(phil), 1)

# constraint, by equation (10)
.5*n*(n+1)
sum(xk^2)

# With vector x generate matrix U from (9).
U <- matrix(0, n, n)

Uind <- 1
for (j in 1:n) {
    for (i in 1:j) {
        # print(c(i, j, Uind))
        # U[i, j] <- Uind
        U[i, j] <- xk[Uind]
        Uind <- Uind + 1
    }
}
# U
# xk

Ai <- t(U) %*% U

sum(diag(Ai))


# sanity check of x calculations ------------------------------------------

# y could all be 0.5
# b could all be 0
phi_test <- c(0.5, 0.5)

x1 <- c(cos(.5)*sin(.5), cos(.5)*sin(.5), sin(.5)*sin(.5))
x2 <- c(cos(.5), cos(.5)*sin(.5), sin(.5)*sin(.5))

sum(x1^2)
sum(x2^2)

# alternative method ------------------------------------------------------

# calculate cl and Il as in (29) - (32) and (16)

# compute mul, sigmal2 and nul from (36)-(39)

sigmal2 <- 1 / (sqrt(pl) + sqrt(ql))^2
nul <- sqrt(2 * pi * sigmal2) * cl / dl

# fxil, fphil as in 28, 35
# fxil
cl * cos(phil)^pl * sin(phil)^ql * 
    (phil > 0 & phil < Intl)

# fphil
dnorm(phil, mul, sqrt(sigmal2))

phil <- rep(0, length(ls))
for (l in ls) {
    phi_rej <- TRUE
    
    while(phi_rej) {
        u <- runif(1)
        newphi <- rnorm(1, mul[l], sqrt(sigmal2[l]))
        
        u * nul[l]
        fxil <- cl[l] * cos(newphi)^pl[l] * sin(newphi)^ql[l] *
            (0 < newphi & newphi < Intl[l])
        fphil <- dnorm(newphi, mul[l], sqrt(sigmal2[l]))
        
        if (u * nul[l] * fphil < fxil) {
            phi_rej <- FALSE
        }
    }
    
    phil[l] <- newphi
}

ks <- 1:(.5*n*(n+1))
xk <- rep(0, length(ks))

# Continue Algorithm 1
# compute x from (11)

for (k in ks) {
    if (k == 1) {
        xk[k] <- 1    
    } else {
        xk[k] <- prod(sin(phil[1:(k-1)]))
    }
}
xk <- xk * c(cos(phil), 1)

# constraint, by equation (10)
.5*n*(n+1)
sum(xk^2)

# With vector x generate matrix U from (9).
U <- matrix(0, n, n)

Uind <- 1
for (j in 1:n) {
    for (i in 1:j) {
        # print(c(i, j, Uind))
        # U[i, j] <- Uind
        U[i, j] <- xk[Uind]
        Uind <- Uind + 1
    }
}
# U

Ai <- t(U) %*% U

sum(diag(Ai))


# Density evaluation ------------------------------------------------------

nu <- 6
p <- 4

# V <- rWishart(1, p+1, diag(p))[, , 1]
V <- diag(p)
V <- V / sum(diag(V))

invV <- solve(V)

U <- runif_unittr_spd(1, p)[ , , 1]

MC_its <- 1e6
Us <- runif_unittr_spd(MC_its, p)

HypergeoMat::mvgamma(nu/2, p)^(-1) *
    2^(-nu*p/2) *
    # ignore V for now - assume it is the identity
    det(U)^((nu-p-1)/2) *
    gamma(p*nu/2) *
    2^(p*nu/2)

intgd <- rep(NA, MC_its)
for (s in 1:MC_its) {
    Uu <- Us[, , s]
    intgd[s] <- HypergeoMat::mvgamma(nu/2, p)^(-1) *
        2^(-nu*p/2) *
        # ignore V for now - assume it is the identity
        det(V)^(-nu/2) *
        det(Uu)^((nu-p-1)/2) *
        gamma(p*nu/2) *
        2^(p*nu/2) *
        sum(diag(invV %*% Uu))^(-p*nu/2)
}

intgd[1:5]
summary(intgd)
plot(density(intgd))

summary(log(intgd))
plot(density(log(intgd)))

mean(intgd)

# volume?
Rvol <- pi^(.25*p*(p-1)) *
    prod(gamma( ((2:p) + 1)/2 )) / gamma(.5*p*(p+1))

Rvol * mean(intgd)

1/Rvol

c(mean(intgd) + 1.96*sd(intgd) / sqrt(MC_its), 
  mean(intgd) - 1.96*sd(intgd) / sqrt(MC_its))

# MH proposal calculations

U1 <- runif_unittr_spd(1, p)[, , 1]
U2 <- runif_unittr_spd(1, p)[, , 1]

det(U1)^(-nu/2) * det(U2)^((nu-p-1)/2) * sum(.5 * diag(solve(U1) %*% U2))^(-p*nu/2)
det(U2)^(-nu/2) * det(U1)^((nu-p-1)/2) * sum(.5 * diag(solve(U2) %*% U1))^(-p*nu/2)



# Cayley for rectangular matrices -----------------------------------------

m <- 4

S <- matrix(0, m, m)
# S[upper.tri(S)] <- c(1, rep(0, choose(m, 2)-1))
S[upper.tri(S)] <- c(1.5, -0.5, rep(0, choose(m, 2)-2))
# S[upper.tri(S)] <- rep(-0.1, choose(m, 2))

for(j in 1:(m-1)) {
    for (i in (j+1):m) {
        S[i, j] <- -S[j, i]
    }
}
S
eigen(diag(m) + S)
solve(diag(m) + S)

U <- (diag(m) - S) %*% solve(diag(m) + S)
Ur <- U %*% diag(1, nrow = p, ncol = 2)

Ur %*% t(Ur)


# degrees of freedom of unit trace complex Wishart ------------------------

p <- 4

set.seed(4062025)
V <- rcomplex_wishart(p+1, p, diag(p))
V <- V / Re(sum(diag(V)))

nu <- 10

Us <- array(NA, c(p, p, 1e5))

for(i in 1:1e5) {
    Uu <- rcomplex_wishart(nu, p, V)
    Us[, , i] <- Uu / Re(sum(diag(Uu)))
}

dCTCW <- function(X, df, V, a = 1, logscale = FALSE) {
    p <- nrow(X)
    
    logdens <- -mvGamma(df, p, beta = 2, logscale = TRUE) +
        -df * log(Re(EigenR::Eigen_det(V))) +
        (df - p) * log(Re(EigenR::Eigen_det(X))) +
        lgamma(p*df) +
        (p^2 - 1)*log(a) + 
        -p*df * log(Re(sum(diag( solve(V) %*% X ))))
    
    if (logscale) {
        return(logdens)
    } else {
        return(exp(logdens))
    }
}

apply(Us, c(1, 2), mean)

nu_test <- nu
U <- V
U <- U / Re(sum(diag(U)))

dCTCW(U, nu_test, V)
dCTCW(U, nu_test, V)

dCTCW(Us[, , 1], nu_test, V)
dCTCW(V, nu_test, Us[, , 1])

# dCTCW
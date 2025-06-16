# compare calculations of 1F1 from Butler and Wood 2005

# HGF 0F0 by MC integration functions

### functions to simulate from complex Stiefel manifold
# return n independent samples from the standard complex normal distribution
rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

# draw matrix distributed uniformly on V_d^P(F_betaf), where 
# - betaf = 1 --> Real
# - betaf = 2 --> Complex
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
    # transform by Q = QS, where S = diag(D)
    
    # Qfinal <- t(t(QX1) * D) # faster than matrix multiplication, for larger matrices
    Qfinal <- QX1 %*% diag(D)
    
    return(Qfinal)
}

# 1F1 one parameter function ----------------------------------------------

hgf_1F1_1p <- function(a, b, X) {
    # Ensure X is a diagonal matrix or a vector of diagonal elements
    if (is.matrix(X)) {
        x <- diag(X)
    } else {
        x <- X
    }
    
    m <- length(x)
    
    ys <- 2*a / (b - x + sqrt( (x - b)^2 + 4*a*x ))
    
    R11 <- 1
    for (i in 1:m) {
        for (j in i:m) {
            R11 <- R11 * (ys[i]*ys[j]/a + (1 - ys[i])*(1 - ys[j])/(b - a))
        }
    }
    
    output <- b^(b*m - m*(m+1)/4) *
        R11^(-1/2) *
        prod((ys/a)^a * ((1-ys)/(b-a))^(b-a) * exp(x*ys))
    
    return(output)
}

# a function needed for calculating derivatives of log(1F1(a, b, X))

# - x: the one value in X to be differentiated against and evaluated at
# - i: the index of x
# - X: a vector containing all other values
# - a, b: the usual arguments

hgf_1F1_1p_x <- function(x, i, X, a, b) {
    X[i] <- x
    return(log(hgf_1F1_1p(a, b, X)))
}


# example parameters ------------------------------------------------------

# Butler and Wood example matrices and parameters
# m <- 3
# n1 <- 3
# n2 <- 5
# Theta <- c(2, 1, 1/2)
# B1 <- c(3/4, 1/2, 1/4)
# # B1 <- c(1/2, 3/10, 1/10)
# a <- (n1 + n2)/2
# b <- n1 / 2

m <- 4
n1 <- 10
n2 <- 60
Theta <- (4:1)/10
B1 <- c(8, 6, 4, 2)/10
# B1 <- c(7, 5, 3, 1)/10
a <- (n1 + n2)/2
b <- n1 / 2

# m <- 8
# n1 <- m
# n2 <- m*2-1
# Theta <- seq(.5, by = .5, length.out = m) |> rev()
# B1 <- seq(.25, by = .25, length.out = m) |> rev()
# # B1 <- c(1/2, 3/10, 1/10)
# a <- (n1 + n2)/2
# b <- n1 / 2

# 1F1 approximation by MC integral ----------------------------------------

# Monte Carlo integral approximation

# Qi: generate a uniform O(m) matrix
# - Wi: generate Wishart_m(2m, I_m)
# - define Qi as eigenvectors of Wi

# calculate single-parameter 1F1

# - try HypergeoMat function
# - try to implement equation (6) in Butler and Wood (somewhat complicated, but I have all the information I need)

MCits <- 1e5
intgd <- rep(NA, MCits)

for (i in 1:MCits) {

    if (i %% 1e4 == 0) {
        print(paste0("i = ", i))
    }
    
    # Wi <- rWishart(1, 2*m, diag(m))[, , 1]
    # Qi <- eigen(Wi)$vectors
    Qi <- runif_stiefel(m, m, 1)
    
    # argument to the 1-parameter HGF
    Xi <- eigen(n2/2 * diag(Theta) %*% Qi %*% diag(B1) %*% t(Qi))$values
    
    # 1F1(X) approximation
    intgd[i] <- hgf_1F1_1p(a, b, Xi)
    # intgd[i] <- hypergeometric_1F1(a, b, Xi)
    # intgd[i] <- HypergeoMat::hypergeomPFQ(15, a, b, Xi)
}

# first "Exact" value in Butler and Wood 2005 is 8569.0 +- 9.8

mean(intgd)
2*sd(intgd)/(sqrt(MCits))

# 1F1 Laplace approximation -----------------------------------------------

# - pi stuff
# - Omega: volume, dimensions of matrices, multivariate Gamma, etc.
# - J: product of as, bs, and Psi, which also depends on fi
# - fi: partial derivative of one-parameter HGF
# - 1F1(X): one-parameter HGF

# fi, partial derivatives of 1F1(a, b, X)

A <- .5*n2*Theta
B <- B1
X <- A * B
ss <- m*(m-1)/2

fs <- rep(NA, m)
for (i in 1:m) {
    fs[i] <- pracma::fderiv(hgf_1F1_1p_x, X[i], i = i, X = X, a = a, b = b)
}

ijs <- combn(m, 2)
ii <- ijs[1, ]
jj <- ijs[2, ]

# J(A, B)
Psi_ij <- 2 * (fs[ii]*A[ii]*B[ii] - fs[jj]*A[jj]*B[jj]) / 
    (A[ii]*B[ii] - A[jj]*B[jj])
cij <- (A[ii] - A[jj]) * (B[ii] - B[jj])
JAB <- prod(cij * Psi_ij)

# Omega
omegami <- 2^1 * pi^(1^2 / 2) / HypergeoMat::mvgamma(1/2, 1)
omegam <- 2^m * pi^(m^2 / 2) / HypergeoMat::mvgamma(m/2, m)
Omegams <- (omegami^m) / omegam

(2*pi)^(ss/2) *
    Omegams *
    JAB^(-1/2) *
    hgf_1F1_1p(a, b, X)

# try 0F0(A, B) approximation ---------------------------------------------

MCits <- 1e5
intgd_0F0 <- rep(NA, MCits)

for (i in 1:MCits) {
    
    if (i %% 1e4 == 0) {
        print(paste0("i = ", i))
    }
    
    # Wi <- rWishart(1, 2*m, diag(m))[, , 1]
    # Qi <- eigen(Wi)$vectors
    Qi <- runif_stiefel(m, m, 1)
    
    # argument to the 1-parameter HGF
    Xi <- eigen(n2/2 * diag(Theta) %*% Qi %*% diag(B1) %*% t(Conj(Qi)))$values
    
    # 1F1(X) approximation
    # intgd[i] <- hgf_1F1_1p(a, b, Xi)
    # 0F0(X)
    intgd_0F0[i] <- exp(sum(Xi))
}

mean(intgd_0F0)
2*sd(intgd_0F0)/sqrt(MCits)

# Laplace approximation

# J(A, B); Psi is just 2 for 0F0
Psi_0F0 <- 2
cij <- (A[ii] - A[jj]) * (B[ii] - B[jj])
JAB_0F0 <- prod(cij * Psi_0F0)

(2*pi)^(ss/2) *
    Omegams *
    JAB_0F0^(-1/2) *
    exp(sum(X))

# understanding MACG and CMACG


# functions ---------------------------------------------------------------

dMACG <- function(X, Sigma, k, logscale = FALSE) {
    n <- nrow(X)
    logdens <- -(k/2) * log(det(Sigma)) +
        (-n/2) * log(det( t(X) %*% solve(Sigma) %*% X ))
    if (logscale) {
        return(logdens)
    } else {
        return(exp(logdens))
    }
}

frob_dist <- function(A, B, returnDists = FALSE) {
    diffAB <- A - B
    sqDists <- Re(diag( t(Conj(diffAB)) %*% diffAB ))
    Fdist <- Re(sum(sqDists))
    if (returnDists) {
        return(list(Fdist = sqrt(Fdist),
                    sqDists = sqDists))
    } else {
        return(sqrt(Fdist))
    }
}

evec_Frob_stat <- function(X, Y, returnDists = FALSE) {
    Rmat <- diag(complex(modulus = 1, 
                         argument = -Arg(diag(t(Conj(X)) %*% Y))))    
    return(frob_dist(X, Y %*% Rmat, returnDists = returnDists))
}


# test of uniqueness of eigenvectors and rotations ------------------------

P <- 4

set.seed(7062025)
V0 <- pracma::randortho(4)
L <- diag(P:1)

Sigma <- V0 %*% L %*% t(V0)

eigen(Sigma)

Lr <- diag(c(4, 4, 2, 1))

Vr <- eigen(V0 %*% Lr %*% t(V0))$vectors

R <- t(V0[, 1:2]) %*% Vr[, 1:2]

t(R) %*% R
R %*% t(R)

V0[, 1:2] %*% t(V0[, 1:2])
Vr[, 1:2] %*% t(Vr[, 1:2])

# k and k + 1 eigenvalues are not unique
Lnu <- diag(c(4, 3, 3, 1))
Vnu <- eigen(V0 %*% Lnu %*% t(V0))$vectors

Vnu_swapped <- Vnu[, c(1, 3, 2, 4)]

V0 %*% Lnu %*% t(V0)
Vnu %*% Lnu %*% t(Vnu)
Vnu_swapped %*% Lnu %*% t(Vnu_swapped)


Rnu <- t(V0[, 2:3]) %*% Vnu[, 2:3]

t(Rnu) %*% Rnu
Rnu %*% t(Rnu)


# comparing modes and concentration ---------------------------------------

set.seed(7062025)
V <- pracma::randortho(4)
# a bit larger
A1 <- diag(c(1.5, 1.5, .5, .5))
# smaller
A2 <- diag(c(1.25, 1.25, .75, .75))
# much larger
A3 <- diag(c(1.95, 1.95, .05, .05))
# unique first 2
A4 <- diag(c(2.1, 1.8, .05, .05))


Sigma1 <- V %*% A1 %*% t(V)
Sigma2 <- V %*% A2 %*% t(V)
Sigma3 <- V %*% A3 %*% t(V)
Sigma4 <- V %*% A4 %*% t(V)

V1 <- eigen(Sigma1)$vectors
V2 <- eigen(Sigma2)$vectors
V4 <- eigen(Sigma4)$vectors


dMACG(V[, 1:2], Sigma1, 2, TRUE)
dMACG(V1[, 1:2], Sigma1, 2, TRUE)

dMACG(V[, 1:2], Sigma2, 2, TRUE)
dMACG(V2[, 1:2], Sigma2, 2, TRUE)

dMACG(V[, 1:2], Sigma3, 2, TRUE)

# what if the initial eigenvalues are not the same?
dMACG(V[, 1:2], Sigma4, 2, TRUE)


# upper bounds to rotation distances --------------------------------------

RL <- pracma::randortho(4, "unitary")

U <- diag(1, nrow = 4, ncol = 2)

RL %*% U

t(U) %*% RL
t(U) %*% U

U %*% t(U) %*% RL


# average of unit modulus complex numbers ---------------------------------

unitmods <- complex(modulus = 1, argument = runif(1e4, -pi, pi))

mean(unitmods)


# is Bingham distribution invariant to Right rotations? -------------------

A <- diag(c(4, 3, 2, 1))
B <- diag(c(4, 2.5, 1.5, 1))
V <- pracma::randortho(4)

U <- pracma::randortho(4)
R <- pracma::randortho(4)

# proportional log density value
sum(diag(B %*% t(U) %*% V%*%A%*%t(V) %*% U ))
sum(diag(B %*% t(U%*%R) %*% V%*%A%*%t(V) %*% U%*%R))

# Wishart with constant eigenvalues ---------------------------------------

set.seed(8062025)
Ww <- rWishart(1, 5, diag(4))[, , 1]
Wevals <- eigen(Ww)$values
Wevecs <- eigen(Ww)$vectors

Wctrue <- Wevecs %*% diag(4:1) %*% t(Wevecs)

newevecs <- 1/Wevals * (4:1)

diag(sqrt(newevecs)) %*% Ww %*% diag(sqrt(newevecs))

Ww |> det()
Ww |> qr() |> qr.Q()
Ww |> qr() |> qr.R()
Ww |> qr() |> qr.R() |> solve()

Ww |> solve() |> qr() |> qr.Q()
Ww |> solve() |> qr() |> qr.R()

# checking for positive elements ------------------------------------------

set.seed(9062025)
P <- rWishart(1, 7, diag(6))[, , 1]

X <- matrix(rnorm(6*4), ncol = 4)
qrX <- qr(X)
Tt <- qrX |> qr.R()
H <- (qrX |> qr.Q()) %*% diag(sign(diag(Tt)))
Tt <- diag(sign(diag(Tt))) %*% Tt

all.equal(H %*% Tt, X)

t(H) %*% solve(P) %*% H

t(H) %*% H

# randomly 0 Cayley transformations ---------------------------------------

# choose only some of the triangular elements to be non-zer0

m <- 4

S <- matrix(0, m, m)
# S[upper.tri(S)] <- c(1, rep(0, choose(m, 2)-1))
# S[upper.tri(S)] <- c(1.5, -0.5, rep(0, choose(m, 2)-2))
# S[upper.tri(S)] <- rep(-0.1, choose(m, 2))
S[upper.tri(S)] <- runif(choose(m, 2), -1, 1) * rbinom(choose(m, 2), 1, 0.5)

# S[upper.tri(S)] <- runif(choose(m, 2), -1, 1) + runif(choose(m, 2), -1, 1) * 1i

for(j in 1:(m-1)) {
    for (i in (j+1):m) {
        S[i, j] <- -Conj(S[j, i])
    }
}
S

U <- (diag(m) - S) %*% solve(diag(m) + S)

SU <- (diag(m) - U) %*% solve(diag(m) + U)
SU[upper.tri(SU)]

SUt <- (diag(m) - t(U)) %*% solve(diag(m) + t(U))
SUt[upper.tri(SUt)]
S[upper.tri(S)]


# complex numbers and vectors ---------------------------------------------

z1 <- complex(modulus = 1, argument = -0.25)
plot(z1,
     xlim = c(-1, 1),
     ylim = c(-1, 1),
     asp = 1)
DescTools::DrawCircle(0, 0)
points(c(z1 * complex(modulus = 1, argument = 2),
         z1 * complex(modulus = 1, argument = 1)),
       col = 2)


# estimate of true eigenvectors from data ---------------------------------

# E (Y / n) = s2 * (U L U^H + I_P)
Uk <- pracma::randortho(4, "orthonormal")[, 1:2]
LL <- diag(seq(20, 1, length.out = 2))
s2 <- 5

Gammak <- s2 * (Uk %*% LL %*% t(Uk) + diag(4))
eigen(Gammak)


# distance of flipped eigenvectors ----------------------------------------

a <- pracma::randortho(4, "u")[, 1:2]
b <- -a

sum(diag(t(Conj(a - b)) %*% (a - b)))

# optimal eigenvector matrix distances ------------------------------------

set.seed(11062025)
X <- pracma::randortho(8, "unitary")[, 1:4]
Y <- pracma::randortho(8, "unitary")[, 1:4]
# # Y <- X
Z <- pracma::randortho(8, "unitary")[, 1:4]

# A <- pracma::randortho(8, "unitary")
# X <- A[, 1:4]
# Y <- pracma::randortho(8, "unitary")[, 1:4]
# Y <- X %*% pracma::randortho(4, "unitary")
# Z <- A[, 5:8]

t(Conj(X)) %*% Z

diag(t(Conj(X)) %*% Y)

Arg(diag(t(Conj(X)) %*% Y))
(svd(t(Conj(X)) %*% Y)$d) |> atan()

Mod(diag(t(Conj(X)) %*% Y))

Rmat <- diag(complex(modulus = 1, 
                     argument = -Arg(diag(t(Conj(X)) %*% Y))))

Re(sum(diag( t(Conj(X)) %*% Y %*% Rmat )))
sum(Mod(diag( t(Conj(X)) %*% Y )))

sum(Mod(diag( t(Conj(Y)) %*% Z )))

sum(Mod(diag( t(Conj(X)) %*% Z )))

sqrt(2 * 4 - 2 * sum(Mod(diag(t(Conj(X)) %*% Y))))
evec_Frob_stat(X, Y, returnDists = TRUE)

# top should be less than or equal to bottom
(-sum(Mod(diag( t(Conj(X)) %*% Z ))))
(-sum(Mod(diag( t(Conj(X)) %*% Y ))) + -sum(Mod(diag( t(Conj(Y)) %*% Z ))) + 4)

# top should be greater than or equal to bottom
(sum(Mod(diag( t(Conj(X)) %*% Z ))))
(sum(Mod(diag( t(Conj(X)) %*% Y ))) + sum(Mod(diag( t(Conj(Y)) %*% Z ))) - 4)

t(Conj(X[, 1])) %*% Y[, 1]
Mod( t(Conj(X[, 1])) %*% Y[, 1] )

triTrue <- rep(NA, 1e4)
# triangle inequality?
for (i in 1:1e4) {
    X <- pracma::randortho(8, "unitary")[, 1:4]
    Y <- pracma::randortho(8, "unitary")[, 1:4]
    Z <- pracma::randortho(8, "unitary")[, 1:4]
    
    triTrue[i] <- (sum(Mod(diag( t(Conj(X)) %*% Z )))) >=
    (sum(Mod(diag( t(Conj(X)) %*% Y ))) + sum(Mod(diag( t(Conj(Y)) %*% Z ))) - 4)
    
}
mean(triTrue)


# distance between equivalence classes ------------------------------------

set.seed(11062025)
X <- pracma::randortho(8, "unitary")[, 1:4]
Y <- pracma::randortho(8, "unitary")[, 1:4]
# # Y <- X
Z <- pracma::randortho(8, "unitary")[, 1:4]

R <- complex(modulus = 1, argument = runif(4, -pi, pi)) |> diag()
Ro <- complex(modulus = 1, argument = runif(4, -pi, pi)) |> diag()

evec_Frob_stat(X, Y)
evec_Frob_stat(Y, X)
evec_Frob_stat(Y%*%R, X)
evec_Frob_stat(X%*%R, Y)
evec_Frob_stat(X%*%R, Y%*%Ro)

evec_Frob_stat(X, X)
evec_Frob_stat(X%*%R, X)
evec_Frob_stat(X%*%R, X%*%Ro)


# examples of axis estimates ----------------------------------------------

trueX <- matrix(sample(20, 6), ncol = 2)
trueX <- qr(trueX) |> qr.Q()
smallR <- diag( complex(modulus = 1, argument = runif(2, -0.1, 0.1)) )
bigR <- diag( complex(modulus = 1, argument = runif(2, -pi, pi)) )

trueX %*% diag(c(-1, -1))
frob_dist(trueX, trueX %*% diag(c(-1, -1)), returnDists = TRUE)

trueX %*% smallR
frob_dist(trueX, trueX %*% smallR, returnDists = TRUE)

trueX %*% bigR
frob_dist(trueX, trueX %*% bigR, returnDists = TRUE)

evec_Frob_stat(trueX, trueX %*% smallR)
evec_Frob_stat(trueX, trueX %*% bigR)

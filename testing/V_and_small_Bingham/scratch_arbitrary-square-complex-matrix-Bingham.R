# testing generating Complex matrix Bingham matrices

library(rstiefel)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")

# an example of orthogonal 2D complex vectors
x <- matrix(c(1, 1i), ncol = 1)
y <- matrix(c(-1, 1i), ncol = 1)

t(Conj(x)) %*% y
t(Conj(y)) %*% x


# parameters, A and B -----------------------------------------------------

istatus <- 0
gibbsIts <- 2000

P <- 4

set.seed(9032025)
Av <- runitary(P, P)
Aevals <- (P:1)*20
# A <- rcomplex_wishart(P, P, diag(P))
A0 <- Av %*% diag(Aevals) %*% t(Conj(Av))
# B <- runif(P, 0, 10) |> sort(decreasing = TRUE)
B0 <- (P:1)

eigen(A0)$values
eigen(A0)
B0


# a random starting matrix ------------------------------------------------

X <- runitary(P, P)
Xs <- array(NA, c(P, P, gibbsIts))

Xs[, , 1] <- X 

# function for generating the Bingham matrix ------------------------------


for (s in 2:gibbsIts) {
    
    if (s %% 50 == 0) print(paste("Gibbs", s))
    
    X <- Xs[, , s-1]
    
    stopifnot("A and B must have the same dimensions" = all(dim(A0) == dim(B0)))
    stopifnot("A and B must be square" = dim(A0)[1] == dim(A0)[2])
    
    P <- nrow(A0)
    
    stopifnot("A and B must have even dimensions" = P %% 2 == 0)
    
    # generate a random order of columns
    # perhaps a neater way: matrix(sample(1:4), ncol = 2, byrow = TRUE)
    scols <- sample(1:P)
    
    for (sstep in 1:(P/2)) {
        # get the random column indices for this step, and put in increasing order. 
        ijs <- scols[c(sstep*2 - 1, sstep*2)] |> sort()
        
        # find an orthonormal basis for the left null space of X without ijs 
        N <- NullC(X[, -ijs])
        
        # transform the original parameters
        newA <- t(Conj(N)) %*% A0 %*% N
        newB <- diag(B0[ijs])
        
        # sample the 2x2 columns, and transform back into X columns
        # zsamp <- rcBingUP(newA, newB)
        
        zsamp <- my.rCbing.Op(newA, newB, istatus = istatus)
        X[, ijs] <- N %*% zsamp$X
    }

    Xs[, , s] <- X
}

isSymmetric(A0)
isSymmetric(newA)
eigen(newA)

newA
newB

A <- newA
B <- newB

Xs[, , gibbsIts]

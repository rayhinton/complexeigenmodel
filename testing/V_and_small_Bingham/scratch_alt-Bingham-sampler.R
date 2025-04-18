# alternative Bingham sampler

library(rstiefel)

Gvecs <- 1/sqrt(2)*matrix(c(1, -1, 1, 1), ncol = 2, byrow = TRUE)
Gvals <- c(4, 1)
G <- Gvecs %*% diag(Gvals) %*% t(Gvecs)

Hvals <- c(100, 1)
H <- diag(Hvals)
Hri <- diag(1/(sqrt(Hvals)))

nreps <- 100

plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, col = 1:2)

set.seed(2)
Us <- array(NA, c(2, 2, nreps))
Ys <- array(NA, c(2, 2, nreps))

S <- solve(diag(7, nrow = 2) - G)/2
eigen(S)

for (i in 1:nreps) {
    Us[, , i] <- rbing.Op(G, H)
    
    points(t(Us[, , i]), 
           xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, col = 1:2)
    
    W <- rWishart(1, 4, S)[, , 1]
    Weig <- eigen(W)
    L <- Weig$values
    Lroot <- diag(sqrt(L))
    
    Ys[, , i] <- Weig$vectors %*% Lroot %*% Hri
    
    points(t(Ys[, , i]), 
           xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, col = 3:4)
}

crossprod(Ys[, , 1])

#

set.seed(7042025)
X <- matrix(rnorm(12, 0, 2), ncol = 3)
sds <- apply(X, 2, function(x) (3/4)*sd(x))
Xs <- t(t(X) / sds)
Xs <- scale(Xs, center = TRUE, scale = FALSE)

eigen(t(Xs) %*% Xs)

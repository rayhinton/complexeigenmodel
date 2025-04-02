# investigate 2x2 Bingham matrices from rejection samplers

library(rstiefel)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

###
# multiplying matrices column wise
dd <- 64
ZZ <- matrix(1:dd^2, ncol = dd)

diagpart <- c(-1, rep(1, dd-1))

# by matrix multiplication
# ZZ %*% diag(diagpart)
# by transposing
# t(t(ZZ) * diagpart)

microbenchmark::microbenchmark(ZZ %*% diag(diagpart), 
                               t(t(ZZ) * diagpart))

###

?rbing.O2
?rbing.Op

A <- rWishart(1, 3, diag(3))[, , 1]

(Av <- 1/sqrt(2) * matrix(c(1, 1, -1, 1), ncol = 2))
(A <- Av %*% diag(c(2, 1)) %*% t(Av))
eigen(A)

B <- diag(c(20, 1))

# reference vector
ur <- matrix(c(1, 0), ncol = 1)

# Acv
Acv <- 1/sqrt(2) * matrix(c(1, 1i, -1i, 1), ncol = 2)
t(Conj(Acv)) %*% Acv

plot2colvec <- function(U, add = FALSE, col = c("black", "red")) {
    if (add == FALSE) {
        plot(x = NA, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
    }
    arrows(x0 = 0, y0 = 0, x1 = U[1, ], y1 = U[2, ], col = col,
           xlim = c(-1, 1))
}

plot2colvec(Av, add = FALSE, col = c("blue", "green"))
abline(a = 0, b = 1, col = "blue")
abline(a = 0, b = -1, col = "green")
legend(x = "topright",
       legend = c("true 1st", "true 2nd", "rand. 1st", "rand. 2nd"),
       col = c("blue", "green", "black", "red"),
       lwd = 2)

{
U <- rbing.Op(A, B)
plot2colvec(U, add = TRUE)
}

colmod <- sign(t(ur) %*% U)
t(t(U) * c(colmod))

# if we change the sign of 1 column, the density is the same
# if we changes the sign of all columns, obviously the density would be the same (since that would be equivalent to multiplying the matrix by -1, and thus multiplying by 1 in the density)
sum(diag(B %*% t(U) %*% A %*% U))
U1 <- U %*% diag(c(-1, 1))
sum(diag(B %*% t(U1) %*% A %*% U1))

eigen(A[, 1] %*% t(A[, 1]))

its <- 2000
Us <- array(NA, c(2, 2, its))
U_fs <- array(NA, dim = c(2, 2, its))
U_flipped <- array(FALSE, dim = c(2, its))
Ps <- array(NA, c(2, 2, 2, its))

set.seed(1042025)
for (i in 1:its) {
    
    # generate a random matrix
    U <- rbing.Op(A, B)
    # store the raw sample
    Us[, , i] <- U
    
    # flip the column signs, relative to reference vector
    # if (Re(t(U[, 1]) %*% ur) < 0) {
    #     U[, 1] <- -U[, 1]
    #     U_flipped[1, i] <- TRUE
    # }
    # if (Re(t(U[, 2]) %*% ur) < 0) {
    #     U[, 2] <- -U[, 2]
    #     U_flipped[2, i] <- TRUE
    # }
    
    colmod <- sign(t(ur) %*% U)
    U <- t(t(U) * c(colmod))
    
    U_flipped[, i] <- colmod == -1
    
    # store the modified sample
    U_fs[, , i] <- U
}

(avgUs <- apply(Us, c(1, 2), mean))
(avgU_fs <- apply(U_fs, c(1, 2), mean))
Av

rowMeans(U_flipped)

Us[, , 1] |> tcrossprod()
U_fs[, , 1]

plot2colvec(Av, add = FALSE, col = c("blue", "green"))
abline(a = 0, b = 1, col = "blue")
abline(a = 0, b = -1, col = "green")
legend(x = "topright",
       legend = c("true 1st", "true 2nd", "rand. 1st", "rand. 2nd"),
       col = c("blue", "green", "black", "red"),
       lwd = 2)

plot2colvec(avgU_fs, add = TRUE)


# complex vectors ---------------------------------------------------------

u <- c(1 +0i, 0+0i)

Ucs <- array(NA, dim = c(2, 2, its))
Uc_fs <- array(NA, dim = c(2, 2, its))

# reference vector
ur <- matrix(c(1 +0i, 0+0i), ncol = 1)
for (i in 1:its) {
    # generate a random matrix
    U <- my.rCbing.Op(A, B)$X
    # store the raw sample
    Ucs[, , i] <- U
    
    # flip the column signs, relative to reference vector
    if (Re(t(U[, 1]) %*% ur) < 0) {
        U[, 1] <- -U[, 1]
    }
    if (Re(t(U[, 2]) %*% ur) < 0) {
        U[, 2] <- -U[, 2]
    }
    
    # store the modified sample
    Uc_fs[, , i] <- U
}

apply(Uc_fs, c(1, 2), mean)
Av

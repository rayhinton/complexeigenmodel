# compare true SDFs to estimated SDFs

sdf_ar2 <- function(omega, phis, sigma2 = 1) {
    phi1 <- phis[1]
    phi2 <- phis[2]
    
    denom <- 1 + phi1^2 + phi2^2 + 2*phi1*(phi2 - 1) * cos(2*pi*omega) - 
        2*phi2*cos(4*pi*omega)
    
    return(sigma2 / denom)
}

P <- 4
d <- 2
Tt <- 1024

# phis <- c(1, -.25, -.75, -.5,
#           .75, -.65, -1, -.5) |>
#     array(c(2, d, K))

sigmak2 <- 5

set.seed(10172025)
U <- pracma::randortho(P)[, 1:d]

Xt <- ts(unname(cbind(arima.sim(list(ar = c(1, -.25)), n = Tt),
               arima.sim(list(ar = c(-.75, -.5)), n = Tt))))

Yt <- ts( sqrt(sigmak2) * t(U %*% t(Xt) + rnorm(Tt * P)) )

Y_sdfest <- astsa::mvspec(Yt, kernel("modified.daniell", 
                                     c(32,32,32)), taper = .5)

# calculate true SDMs -----------------------------------------------------

omegas <- seq(1/Tt, 1/2, by = 1/Tt)
trueSDMs <- array(NA, c(P, P, Tt/2))

for (l in 1:(Tt/2)) {
    Lambdaw <- diag(c(
        sdf_ar2(omegas[l], c(1, -.25)),
        sdf_ar2(omegas[l], c(-.75, -.5)))
    )
    
    trueSDMs[, , l] <- sigmak2 * (U %*% Lambdaw %*% t(Conj(U)) + diag(P))
}

plot(trueSDMs[1, 1, ], type = "l", lty = 2, ylim = c(0, 50))
for (j in 2:P) {
    lines(trueSDMs[j, j, ], lty = 2, col = j)
}

for (j in 1:P) {
    lines(pi*Re(Y_sdfest$fxx[j, j, ]), lty = 1, col = j)
    # lines(Re(Y_sdfest$spec[, j]), lty = 1, col = j)
}

Re(Y_sdfest$fxx[1, 1, 1:5]) / Y_sdfest$spec[1:5, 1]

Y_sdfest$fxx[, , 1:5]

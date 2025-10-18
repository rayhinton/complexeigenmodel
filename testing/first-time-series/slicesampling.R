# slice sampling

w <- 0.1
m <- 10
its <- 100000

xs <- rep(NA, its)

x0 <- 0

for (i in 1:its) {
    
    y <- runif(1, 0, dnorm(x0))
    
    U <- runif(1)
    L <- x0 - w*U
    R <- L + w
    V <- runif(1)
    J <- floor(m*V)
    K <- (m-1) - J
    
    while (J > 0 & y < dnorm(L)) {
        L <- L - w
        J <- J - 1
    }
    
    while (K > 0 & y < dnorm(R)) {
        R <- R + w
        K <- K - 1
    }
    
    x0 <- runif(1, L, R)
    xs[i]<- x0
}

plot(density(xs[seq(from = its/2, to = its, by = 10)]))
curve(dnorm(x), from = -4, to = 4, n = 501, col = "red", add = TRUE)

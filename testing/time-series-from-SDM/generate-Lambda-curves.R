# generate N points

N <- 510
P <- 4
d <- 2

tau1 <- .1
tau2 <- .1

Lambdas <- array(NA, c(d, N))

om1 <- runif(d) |> sort(decreasing = TRUE)

Lambdas[, 1] <- om1 / (1 - om1)

for (j in 1:d) { 
    # Lambdas[j, 2] <- rnorm(1, Lambdas[j, 1], sqrt(tau1))
    Lambdas[j, 2] <- truncdist::rtrunc(1, "norm", a = 0, b = Inf, 
                                       mean = Lambdas[j, 1], sd = sqrt(tau1))
}

for (n in 3:N) {
    for (j in 1:d) {
        Lambdas[j, n] <- 
            truncdist::rtrunc(1, "norm", a = 0, b = Inf, 
                              mean = 2 * Lambdas[j, n-1] - Lambdas[j, n-2], 
                              sd = sqrt(tau2))
    }
}

plot(Lambdas[1, ], type = "l", )
lines(Lambdas[2, ])

fxs <- rep(NA, N)
fxs[1] <- 1
fxs[N] <- 10

fxs <- seq(-1, 2, length.out = N)           
fxs[2:(N - 1)] <- rnorm(N-2, fxs[2:(N - 1)], sd = sqrt(tau2))

fxs <- exp(fxs)
plot(fxs, type = "l")

# Generate a smooth positive random function using splines
library(splines)

generate_smooth_positive_function <- function(N = 100, n_knots = 7, 
                                              log_mean = 0, log_sd = 1) {
    # Define evaluation points
    x <- 1:N
    
    # Sample knot locations (including boundaries)
    knot_locs <- sort(c(1, sample(2:(N-1), n_knots - 2, replace = FALSE), N))
    
    # Sample log-values at knots
    log_y_knots <- rnorm(n_knots, mean = log_mean, sd = log_sd)
    
    # Fit smooth spline through knots
    spline_fit <- splinefun(knot_locs, log_y_knots, method = "natural")
    
    # Evaluate at all points and exponentiate
    log_f <- spline_fit(x)
    f <- exp(log_f)
    
    return(list(x = x, f = f, knot_locs = knot_locs, log_y_knots = log_y_knots))
}

# Example usage
set.seed(123)
# result <- generate_smooth_positive_function(N = 100, n_knots = 5, 
                                            log_mean = 2, log_sd = 0.25)

# Plot
plot(result$x, result$f, type = 'l', xlab = 'x', ylab = 'f(x)',
     main = 'Smooth Positive Random Function')
points(result$knot_locs, exp(result$log_y_knots), col = 'red', pch = 19)

set.seed(30102025)
for (j in 1:d) {
    Lambdas[j, ] <- 
        generate_smooth_positive_function(N, n_knots = 6, 
                                          log_mean = 1, log_sd = 0.25)$f
}

plot(Lambdas[1, ], type = "l", ylim = c(min(Lambdas), max(Lambdas)))
lines(Lambdas[2,], col = "red")

order(Lambdas[, 1], decreasing = TRUE)

Lambdas <- Lambdas[order(Lambdas[, 1], decreasing = TRUE), ]
plot(Lambdas[1, ], type = "l", ylim = c(min(Lambdas), max(Lambdas)))
lines(Lambdas[2,], col = "red")

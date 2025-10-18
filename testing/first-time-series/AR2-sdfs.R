# need spectral density function for AR(2)
sdf_ar2 <- function(omega, phis, sigma2 = 1) {
    phi1 <- phis[1]
    phi2 <- phis[2]
    
    denom <- 1 + phi1^2 + phi2^2 + 2*phi1*(phi2 - 1) * cos(2*pi*omega) - 
        2*phi2*cos(4*pi*omega)
    
    return(sigma2 / denom)
}

# validate conditions for AR(2) causality:

validate_ar2 <- function(phis) {
    pars_valid <- (sum(phis) < 1) & (phis[2] - phis[1] < 1) & 
        (abs(phis[2]) < 1)
    return(pars_valid) 
}

# randomly generate AR(2) parameters

random_AR2_pars <- function() {
    if (sample(0:1, 1)) {
        phi1 <- runif(1, -2, 2)
        phi2 <- runif(1, -1, -abs(phi1) + 1)
    } else {
        phi2 <- runif(1, -1, 1)
        phi1 <- runif(1, -abs(1-phi2), abs(1-phi2))
    }
    
    return(c(phi1, phi2))
}

# choose phi1, phi2 randomly

# phis1 <- random_AR2_pars()
# phis2 <- random_AR2_pars()

# phis1 <- c(1, -.25)
# phis2 <- c(-.75, -0.5)
phis1 <- c(.75, -0.65)
phis2 <- c(-1, -0.5)

validate_ar2(phis1)
validate_ar2(phis2)

curve(sdf_ar2(x, phis1), from = 0, to = 1/2, n = 501)
curve(sdf_ar2(x, phis2), from = 0, to = 1/2, n = 501, add = TRUE, col = "red")

optimize(sdf_ar2, phis = phis1, c(0, 1/2), maximum = TRUE)
optimize(sdf_ar2, phis = phis2, c(0, 1/2), maximum = TRUE)

sdf_ar2(2/2024, phis1)
sdf_ar2(2/2024, phis2)


# compare mvspec estimates to true ----------------------------------------

Xtest <- arima.sim(list(ar = c(1, -.9)), n = 1024)
plot(Xtest)

sdfXest <- astsa::mvspec(Xtest, kernel("daniell", 4))

Re(sdfXest$fxx[1, 1, 1:5])
sdfXest$spec[1:5]

curve(abs(1 - exp(-2i*pi*x) + .9 * exp(-4i*pi*x))^-2, from = 0, to = 1/2,
      n = 1001,
      add = TRUE, col = "red")

# same as the previous calculation
curve(sdf_ar2(x, c(1, -.9)), from = 0, to = 1/2,
      n = 1001,
      add = TRUE, col = "green")

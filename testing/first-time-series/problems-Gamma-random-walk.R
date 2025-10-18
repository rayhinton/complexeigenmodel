# proportional FCD using 2nd order random walk prior

# define using runuran

#distr <- Runuran::udgamma(shape = g_shape, scale = g_scale, 0, 1)
#gen <- Runuran::pinvd.new(distr)
#gen <- Runuran::arsd.new(distr)
#xi_jk <- Runuran::ur(gen, 1)

# functions ---------------------------------------------------------------

# unnormalized PDF
unnorm_logPDF <- function(x, tau2, ajk, l1, l2, N, logscale = FALSE) {
    logf1 <- -N * log(1 + x)
    # logf2 <- a * x / (1+x) - (.5/tau2) * (x^2 + 2*x*(l1 - 2*l2))
    logf2 <- -ajk / (1+x) - (.5/tau2) * (x - (2*l2 - l1))^2
    
    result <- ifelse(x <= 0, -Inf, logf1 + logf2)
    
    if (logscale) {
        return(result)
    } else {
        return(exp(result))
    }
}

d_unnorm_logPDF <- function(x, tau2, ajk, l1, l2, N, logscale = TRUE) {
    f1 <- -N/(1 + x) 
    f2 <- ajk*(1 / (1+x) -x/(1 + x)^2)
    f3 <- -(.5/tau2) * (2*x + 2*(l1 - 2*l2))
    
    return(f1 + f2 + f3)
}

# problems with current Random Walk prior ---------------------------------

tau2 <- 1

# prior, in terms of actual original variable

xi1 <- 0.1
lp1 <- 1/xi1 - 1
alpha <- xi1^2/tau2
beta <- xi1/tau2

xs <- seq(0, 50, by = 0.1)
fxs <- 1/(xs+1)^(alpha + 1) * exp(-beta/(1 + xs))

plot(xs, fxs, type = "l")
abline(v = lp1)

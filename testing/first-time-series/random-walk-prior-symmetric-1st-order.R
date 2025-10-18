# proportional FCD using 2nd order random walk prior

# define using runuran

#distr <- Runuran::udgamma(shape = g_shape, scale = g_scale, 0, 1)
#gen <- Runuran::pinvd.new(distr)
#gen <- Runuran::arsd.new(distr)
#xi_jk <- Runuran::ur(gen, 1)

# functions ---------------------------------------------------------------

# unnormalized PDF
unnorm_logPDF <- function(x, tau2, ajk, l1, N, logscale = FALSE) {
    logf1 <- -N * log(1 + x)
    logf2 <- -ajk / (1+x) - (.5/tau2) * (x - l1)^2
    
    result <- ifelse(x <= 0, -Inf, logf1 + logf2)
    
    if (logscale) {
        return(result)
    } else {
        return(exp(result))
    }
}

# try "normalizing" the PDF, so it is not huge ----------------------------

# problem with these parameters:
# > tau2_Lambda
# [1] 10
# > ajk
# [,1]
# [1,] 546.339
# > l1
# [1] 41.03722
# > l2
# [1] 47.54683
# > LL
# [1] 32

tau2 <- 0.5
ajk <- 50
l1 <- 1
N <- 32

# tau2 <- 10
# ajk <- 546.339
# l1 <- 41.03722
# l2 <- 47.54683
# N <- 32

xs <- seq(-.99, 100, by = .001)

log_fxs <- unnorm_logPDF(xs, tau2, ajk, l1, N, logscale = TRUE)
log_fxs <- log_fxs - max(log_fxs)
# log_fxs <- log_fxs - log(normval)

fxs <- exp(log_fxs)

which.max(fxs)
fxs[which.max(fxs)]
xs[which.max(fxs)]

plot(xs, fxs, type = "l", xlim = c(0, 5))
abline(v = xs[which.max(fxs)], col = "red")

# other cubic solutions ---------------------------------------------------

mu <- l1

a <- 1
b <- 2 - mu
c <- N*tau2 + 1 - 2*mu
d <- N*tau2 - ajk*tau2 - mu

xroots <- polyroot(c(d, c, b, a))

microbenchmark::microbenchmark(
    xroots <- polyroot(c(d, c, b, a)),
    times = 1000
)

microbenchmark::microbenchmark( {
    D0 <- b^2 - 3*a*c
    D1 <- 2*b^3 - 9*a*b*c + 27*a^2*d
    C <- ( (D1 + sqrt(D1^2 - 4*D0^3)) / 2 )^(1/3)
    xmode <- -1/3/a * (b + C + D0/C)
}, times = 1000)


which(zapsmall(Im(xroots)) == 0)

# discriminant < 0 means 1 real root; > 0 means 3 real roots

DISC <- 18*a*b*c*d - 4*b^3*d + b^2*c^2 - 4*a*c^3 - 27*a^2*d^2

if (DISC < 0) {
    # xmode <- max(0, Re(xroots[zapsmall(Im(xroots)) == 0]))
    xmode <- max(1e-12, Re(xroots[zapsmall(Im(xroots)) == 0]))
    # D0 <- b^2 - 3*a*c
    # D1 <- 2*b^3 - 9*a*b*c + 27*a^2*d
    # C <- ( (D1 + sqrt(D1^2 - 4*D0^3)) / 2 )^(1/3)
    # xmode <- -1/3/a * (b + C + D0/C)
} else {
    # TODO just report an error for now
    # but, there could be ways to handle if there are "multiple roots;
    # - good: it could that they are 3 repeated real roots
    # - maybe good: 1 and 2 mult. roots; maybe one is <= 0, outside domain?
    stop(paste0("There are multiple real roots. DISC = ", DISC))
}

print(xmode)

# draw samples ------------------------------------------------------------

microbenchmark::microbenchmark(
    gen <- 
        Runuran::tabl.new(unnorm_logPDF, 
                          lb = 0, ub = Inf, mode = xmode,
                          islog = TRUE,
                          tau2 = tau2, ajk = ajk, l1 = l1, N = N, 
                          logscale = TRUE),
    times = 1000
)

# create the generator object
gen <- 
    Runuran::pinv.new(unnorm_logPDF, lb = 0, ub = Inf,
                      center = xmode, islog = TRUE,
                      uresolution = 1.0e-10, # 1.0e-10 is default
                      smooth = FALSE, # FALSE is default
                      tau2 = tau2, ajk = ajk, l1 = l1, N = N, 
                      logscale = TRUE)

gen <- 
    Runuran::tabl.new(unnorm_logPDF, 
                      lb = 0, ub = Inf, mode = xmode,
                      islog = TRUE,
                      tau2 = tau2, ajk = ajk, l1 = l1, N = N, 
                      logscale = TRUE)

microbenchmark::microbenchmark(Runuran::ur(gen, 1),
                               times = 1000)

# draw a random sample
random_sample <- Runuran::ur(gen, 100000)
# plot the sample density
plot(density(random_sample, from = 0), xlim = c(0, 1))
lines(xs, fxs*1.4, type = "l", col = "red")

# comparing time to set up Runuran distributions --------------------------

distr <- Runuran::udgamma(shape = 1, scale = 1, 0, 1)
#gen <- Runuran::pinvd.new(distr)

microbenchmark::microbenchmark(gen <- Runuran::arsd.new(distr),
                               times = 1000)

# tau2 FCD ----------------------------------------------------------------

Nw <- 384
K <- 2
d <- 2

alpha <- 1
beta <- 1

lambda_jkl <- array(rgamma(Nw*K*d, 1, 1), c(d, K, Nw))

sumall <- sum( (lambda_jkl[, , 2:Nw] - lambda_jkl[, , 1:(Nw-1)])^2 )

alpha + .5*Nw*K*d
beta + .5*sumall

dinvgamma <- function(x, alpha, beta, logscale = FALSE) {
    logf1 <- alpha*log(beta) - lgamma(alpha)
    logf2 <- (-alpha-1) * log(x)
    logf3 <- -beta/x
    
    result <- ifelse(x <= 0, -Inf, logf1 + logf2 + logf3)
    
    if (logscale) {
        return(result)
    } else {
        return(exp(result))
    }
}

curve(dinvgamma(x, 1, 1), from = 0, to = 3, n = 501, ylim = c(0, 5), col = "red")
curve(dinvgamma(x, 3, 0.5), from = 0, to = 3, n = 501, add = TRUE, col = "lightblue")
curve(dinvgamma(x, 2, 1), from = 0, to = 3, n = 501, add = TRUE, col = "green")
curve(dinvgamma(x, 3, 1), from = 0, to = 3, n = 501, add = TRUE, col = "blue")

curve(dinvgamma(x, 1, 1), from = 0, to = 3, n = 501, col = "red", ylim = c(0, 5))
curve(dinvgamma(x, 501, 1086), from = 0, to = 3, n = 1001, 
      add = TRUE, col = "black")

IGsamples <- 1/rgamma(10000, 501, rate = 1086)
lines(density(IGsamples), lty = 2)


# rejection sampling, trying to estimate variance -------------------------

lognormval <- log(sqrt(2*pi/N /
                           (1 / (1+xmode)^2) - ajk/N/(1+xmode)^3 - 1/(N*tau2))) +
    unnorm_logPDF(xmode, tau2, ajk, l1, N, logscale = TRUE)

exp(lognormval)


estvar <- sqrt(2*pi/N /
         abs(-(2/N)/xmode^2 + 1/(1+xmode)^2 - ajk/N/(1 + xmode)^3 - 1/tau2/N)) /
    exp(lognormval) * xmode^2 *
    unnorm_logPDF(xmode, tau2, ajk, l1, N, logscale = FALSE) - 
    xmode^2    

curve(dnorm(x, xmode, sd = sqrt(estvar)), from = 0, to = 5, add = TRUE, 
      lty = 2)

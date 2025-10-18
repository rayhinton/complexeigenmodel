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

tau2 <- 1
ajk <- 5
l1 <- 5
l2 <- 5.1
N <- 32

# tau2 <- 10
# ajk <- 546.339
# l1 <- 41.03722
# l2 <- 47.54683
# N <- 32

xs <- seq(-.99, 100, by = .001)

log_fxs <- unnorm_logPDF(xs, tau2, ajk, l1, l2, N, logscale = TRUE)
log_fxs <- log_fxs - max(log_fxs)
# log_fxs <- log_fxs - log(normval)

fxs <- exp(log_fxs)

which.max(fxs)
fxs[which.max(fxs)]
xs[which.max(fxs)]

plot(xs, fxs, type = "l")
abline(v = xs[which.max(fxs)], col = "red")

# other cubic solutions ---------------------------------------------------

mu <- 2*l2 - l1

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
                          tau2 = tau2, ajk = ajk, l1 = l1, l2 = l2, N = N, 
                          logscale = TRUE),
    times = 1000
)

# create the generator object
gen <- 
    Runuran::pinv.new(unnorm_logPDF, lb = 0, ub = Inf,
                      center = xmode, islog = TRUE,
                      uresolution = 1.0e-10, # 1.0e-10 is default
                      smooth = FALSE, # FALSE is default
                      tau2 = tau2, ajk = ajk, l1 = l1, l2 = l2, N = N, 
                      logscale = TRUE)

gen <- 
    Runuran::tabl.new(unnorm_logPDF, 
                      lb = 0, ub = Inf, mode = xmode,
                      islog = TRUE,
                      tau2 = tau2, ajk = ajk, l1 = l1, l2 = l2, N = N, 
                      logscale = TRUE)

# draw a random sample
random_sample <- Runuran::ur(gen, 10000000)
# plot the sample density
plot(density(random_sample, from = 0), xlim = c(0, 1))
lines(xs, fxs*18.5, type = "l", col = "red")

# rejection sampling, instead? --------------------------------------------

fxsprop <- dnorm(xs, mean = xmode, sd = sqrt(80))

{
plot(xs, fxs, col = "red", type = "l")
lines(xs, fxsprop)
}

# normval <- integrate(unnorm_logPDF, 0, Inf, 
#                      tau2 = tau2, a = ajk, l1 = l1, l2 = l2, N = N, logscale = FALSE)$value

pracma::quadgk(unnorm_logPDF, a = 0, b = Inf,
               tau2 = tau2, ajk = ajk, l1 = l1, l2 = l2, N = N, logscale = FALSE)

meanfun <- function(x) {
    x * unnorm_logPDF(x, tau2, ajk, l1, l2, N, logscale = FALSE)/normval
}

integrate(meanfun, 0.0001, Inf)

1/(1 / (1+xmode)^2) - a/N/(1+xmode)^3 - 1/(N*tau2)

lognormval <- log(sqrt(2*pi/N /
                           (1 / (1+xmode)^2) - a/N/(1+xmode)^3 - 1/(N*tau2))) +
    unnorm_logPDF(xmode, tau2, ajk, l1, l2, N, logscale = TRUE)

exp(lognormval)


sqrt(2*pi/N /
         abs(-(2/N)/xmode^2 + 1/(1+xmode)^2 - ajk/N/(1 + xmode)^3 - 1/tau2/N)) /
    exp(lognormval) * xmode^2 *
    unnorm_logPDF(xmode, tau2, ajk, l1, l2, N, logscale = FALSE) - 
    xmode^2    

# tilted exponential sampling ---------------------------------------------

if (ajk > N) {
    M <- (N/a)^N * exp(-N)
} else {
    M <- exp(-ajk)
}

Y <- truncdist::rtrunc(1, "norm", 0, Inf, mean = mu, sd = sqrt(tau2))
U <- runif(1)
(logR <- -N * log(1 + Y) - ajk/(1+Y) - log(M))
if (log(U) <= logR) {
    print(Y)
}

# beta <- (N-2)/xmode
-xmode * (N/(1 + xmode)^2 - 2*ajk/(1 + xmode)^3 - 1/tau2)

beta <- 2
alpha <- beta*xmode + 1

am <- -1/tau2
bm <- -2/tau2 + beta
cm <- -N - 2/tau2 + mu/tau2 - (alpha-1) + 2*beta
dm <- -N + ajk + 2*mu/tau2 - 2/tau2 - 2*(alpha-1) + beta
em <- mu/tau2 - 1/tau2 - (alpha-1)

polyroot(c(em, dm, cm, bm, am))

{
    plot(xs, fxs, type = "l", xlim = c(37, 62))
    abline(v = xs[which.max(fxs)], col = "red")
    curve(dgamma(x, alpha, beta)*13, 
          from = 37, to = 62, n = 1001, add = TRUE, col = "blue")
}

Mfunc <- function(x) {
    -N*log(1+x) - ajk/(1+x) - .5/tau2 * (x - mu)^2 - (alpha-1)*log(x) + beta*x
}

Mfunc(0.6789101+ 0.6362002i)
curve(Mfunc, 0, 101, n = 20001)


# problems with current Random Walk prior ---------------------------------

# prior, in terms of actual original variable

xi1 <- 0.1
lp1 <- 1/xi1 - 1
alpha <- xi1^2/tau2
beta <- xi1/tau2

xs <- seq(0, 50, by = 0.1)
fxs <- 1/(xs+1)^(alpha + 1) * exp(-beta/(1 + xs))

plot(xs, fxs, type = "l")
abline(v = lp1)

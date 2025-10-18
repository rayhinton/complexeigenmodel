# rejection sampling for Lambda FCDs


# functions ---------------------------------------------------------------

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

dtruncnorm <- function(x, a, b, mean, sd, log = FALSE) {
    logZ <- log(pnorm((b - mean)/sd) - pnorm((a - mean)/sd))
    ds <- dnorm((x - mean)/sd, log = TRUE) - log(sd) - logZ
    result <- ifelse(a <= x & x <= b, ds, -Inf)
    
    if (log) {
        return(result)
    } else {
        return(exp(result))
    }
}

FCD_rejsampler <- function(n, mode, sigmap2, tau2, ajk, l1, L, 
                           print_ntries = FALSE) {
    ap <- 1/sigmap2 - 1/tau2
    bp <- (2 - xmode)/sigmap2 - (2 - l1)/tau2
    cp <- (1 - 2*xmode)/sigmap2 - (1 - 2*l1)/tau2 - L
    dp <- -L + ajk + l1/tau2 - xmode/sigmap2
    
    proots <- zapsmall(polyroot(c(dp, cp, bp, ap)))
    
    (xM <- max(1e-12, max(Re(proots[Im(proots) == 0]))) )
    # xM <- 1e-12
    # xM <- 5.24962
    
    logM <- unnorm_logPDF(xM, tau2, ajk, l1, L, logscale = TRUE) -
        dtruncnorm(xM, 0, Inf, xmode, sqrt(sigmap2), log = TRUE)
    
    FCDsamples <- rep(NA, n)
    ntries <- rep(NA, n)
    logrvals <- c()
    
    for (i in 1:n) {
        logr <- -Inf
        logu <- 1
        ntry <- 0
        while (logu >= logr) {
            (propx <- truncdist::rtrunc(1, "norm", 0, Inf,
                                        mean = xmode, sd = sqrt(sigmap2)))
            
            loggpropx <- dtruncnorm(propx, 0, Inf, xmode, sqrt(sigmap2), log = TRUE)
            
            logu <- log(runif(1))
            logr <- unnorm_logPDF(propx, tau2, ajk, l1, L, TRUE) - loggpropx - logM
            ntry <- ntry + 1
            
            logrvals <- c(logrvals, logr)
        }
        
        ntries[i] <- ntry
        FCDsamples[i] <- propx
    }
    
    if (print_ntries) {
        print(mean(ntries))
        print("wow")
    }
    
    return(FCDsamples)
}

# compare plots -----------------------------------------------------------

# a <- -10
# b <- 10
# 
# mu <- -8
# sigma <- 2
# 
# Z <- pnorm((b - mu)/sigma) - pnorm((a - mu)/sigma)
# NC <- sqrt(sigma^2*2*pi)*Z
# 
# xs <- seq(a, b, by = 0.01)
# 
# plot(xs,
#      truncdist::dtrunc(xs, "norm", a, b, mean = mu, sd = sigma),
#      type = "l", ylab = "density")
# 
# manual_den <- 1/NC * exp(-.5/sigma^2 * (xs - mu)^2)
# lines(xs, manual_den, col = "red")
# lines(xs, dtruncnorm(xs, a, b, mu, sigma), col = "blue")

# setup -------------------------------------------------------------------

tau2 <- 1
sigmap2 <- tau2*1.5
l1 <- 10
ajk <- 25
L <- 32

# find mode ---------------------------------------------------------------

mu <- l1

a <- 1
b <- 2 - mu
c <- L*tau2 + 1 - 2*mu
d <- L*tau2 - ajk*tau2 - mu

xroots <- polyroot(c(d, c, b, a))

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

# proposal ----------------------------------------------------------------

ap <- 1/sigmap2 - 1/tau2
bp <- (2 - xmode)/sigmap2 - (2 - l1)/tau2
cp <- (1 - 2*xmode)/sigmap2 - (1 - 2*l1)/tau2 - L
dp <- -L + ajk + l1/tau2 - xmode/sigmap2

DISCp <- 18*ap*bp*cp*dp - 4*bp^3*dp + bp^2*cp^2 - 4*ap*cp^3 - 27*ap^2*dp^2

# TODO not the best way to do this, in general
proots <- zapsmall(polyroot(c(dp, cp, bp, ap)))
print(proots)

(xM <- max(1e-12, max(Re(proots[Im(proots) == 0]))) )
# xM <- 1e-12
# xM <- 5.24962

logM <- unnorm_logPDF(xM, tau2, ajk, l1, L, logscale = TRUE) -
    dtruncnorm(xM, 0, Inf, xmode, sqrt(sigmap2), log = TRUE)

xmaxer <- xM

-L/(1+xmaxer) + ajk/(1+xmaxer)^2 - (xmaxer - l1)/tau2 + (xmaxer - xmode)/sigmap2

sum(xmaxer^(3:0) * c(ap, bp, cp, dp))

# plot to compare FCD and proposal ----------------------------------------

xs <- seq(1e-12, xmode+400, by = 0.01)
logfxs <- unnorm_logPDF(xs, tau2, ajk, l1, L, logscale = TRUE)

logprops <- dtruncnorm(xs, 0, Inf, xmode, sqrt(sigmap2), log = TRUE)
fprops <- exp(logprops)

fxs <- exp(logfxs)

# is the target density ever larger than the scaled proposal density?
sum((logfxs - logprops) > logM)
which((logfxs - logprops) > logM)[1:5]

xs[(logfxs - logprops) > logM][1:5]

max(abs((logfxs - logprops)))

sort(abs(logfxs - logprops), decreasing = TRUE)[1:5]

max(abs(logfxs))
max(abs(logprops))
min(logprops)
min(logfxs)

plot(xs, fxs, type = "l", xlim = c(1e-12, xmode+4))
abline(v = xmode, lty = 2)
lines(xs, fprops*exp(logM), col = 2)
legend("topright", legend = c("FCD", "proposal", "mode"), col = c(1, 2, 1), 
       lwd = c(2, 2, 1), lty = c(1, 1, 2))

unnorm_logPDF(xmode, tau2, ajk, l1, L, logscale = TRUE)
logM + dtruncnorm(xmode, 0, Inf, xmode, sqrt(sigmap2), log = TRUE)

# M <- 1.43e-55
# do rejection sampling ---------------------------------------------------

nsamps <- 25000

FCDsamples <- FCD_rejsampler(nsamps, xmode, sigmap2, tau2, ajk, l1, L, 
                             print_ntries = TRUE)

# FCDsamples <- rep(NA, nsamps)
# ntries <- rep(NA, nsamps)
# logrvals <- c()
# 
# for (i in 1:nsamps) {
#     logr <- -Inf
#     logu <- 1
#     ntry <- 0
#     while (logu >= logr) {
#         (propx <- truncdist::rtrunc(1, "norm", 0, Inf,
#                                    mean = xmode, sd = sqrt(sigmap2)))
#         
#         loggpropx <- dtruncnorm(propx, 0, Inf, xmode, sqrt(sigmap2), log = TRUE)
#         
#         logu <- log(runif(1))
#         logr <- unnorm_logPDF(propx, tau2, ajk, l1, L, TRUE) - loggpropx - logM
#         ntry <- ntry + 1
#         
#         logrvals <- c(logrvals, logr)
#     }
#     
#     ntries[i] <- ntry
#     FCDsamples[i] <- propx
# }

mean(ntries)
quantile(ntries, probs = c(0.5, 0.75, 0.95, 0.99, 0.999))
max(ntries)

max(logrvals)
sort(logrvals, decreasing = TRUE)[1:10]

plot(density(FCDsamples, bw = "SJ", n = 4096), col = "black")
lines(xs, fxs/exp(logM)*1.2, col = 2)
legend("topright", legend = c("sample", "true"), col = 1:2, lwd = 2)

mean(FCDsamples)
var(FCDsamples)

xmode
tau2


# compare time of rejection sampler and Runuran ---------------------------

microbenchmark::microbenchmark(
    {
        gen <- 
            Runuran::tabl.new(unnorm_logPDF, 
                              lb = 0, ub = Inf, mode = xmode, 
                              islog = TRUE,
                              tau2 = tau2, 
                              ajk = ajk, l1 = l1, N = L, 
                              logscale = TRUE)
        # draw a random sample
        Runuran::ur(gen, 1)
    },
    times = 5000
)

microbenchmark::microbenchmark(
    FCD_rejsampler(1, xmode, sigmap2, tau2, ajk, l1, L),
    times = 5000
)

FCDsamples <- FCD_rejsampler(5000, xmode, sigmap2, tau2, ajk, l1, L)

plot(density(FCDsamples))

# try using different Runuran functions -----------------------------------

gen <- 
    Runuran::tabl.new(unnorm_logPDF, 
                      lb = 0, ub = Inf, mode = xmode, 
                      islog = TRUE,
                      tau2 = tau2, 
                      ajk = ajk, l1 = l1, N = L, 
                      logscale = TRUE)
# draw a random sample
ru_samples <- Runuran::ur(gen, 120000)

plot(density(ru_samples))

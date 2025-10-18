# proportional FCD using 2nd order random walk prior and previous, next neighbors

# define using runuran

#distr <- Runuran::udgamma(shape = g_shape, scale = g_scale, 0, 1)
#gen <- Runuran::pinvd.new(distr)
#gen <- Runuran::arsd.new(distr)
#xi_jk <- Runuran::ur(gen, 1)

# functions ---------------------------------------------------------------

# slice sampling, by Neal
# https://github.com/radfordneal/unislice/blob/master/slice.r

uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL, ...)
{
    # Check the validity of the arguments.
    
    if (!is.numeric(x0) || length(x0)!=1
        || !is.function(g) 
        || !is.numeric(w) || length(w)!=1 || w<=0 
        || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
        || !is.numeric(lower) || length(lower)!=1 || x0<lower
        || !is.numeric(upper) || length(upper)!=1 || x0>upper
        || upper<=lower 
        || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
    { 
        stop ("Invalid slice sampling argument")
    }
    
    # Keep track of the number of calls made to this function.
    
    uni.slice.calls <<- uni.slice.calls + 1
    
    # Find the log density at the initial point, if not already known.
    
    if (is.null(gx0)) 
    { uni.slice.evals <<- uni.slice.evals + 1
    gx0 <- g(x0, ...)
    }
    
    # Determine the slice level, in log terms.
    
    logy <- gx0 - rexp(1)
    
    # Find the initial interval to sample from.
    
    u <- runif(1,0,w)
    L <- x0 - u
    R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
    
    # Expand the interval until its ends are outside the slice, or until
    # the limit on steps is reached.
    
    if (is.infinite(m))  # no limit on number of steps
    { 
        repeat
        { if (L<=lower) break
            uni.slice.evals <<- uni.slice.evals + 1
            if (g(L, ...)<=logy) break
            L <- L - w
        }
        
        repeat
        { if (R>=upper) break
            uni.slice.evals <<- uni.slice.evals + 1
            if (g(R, ...)<=logy) break
            R <- R + w
        }
    }
    
    else if (m>1)  # limit on steps, bigger than one
    { 
        J <- floor(runif(1,0,m))
        K <- (m-1) - J
        
        while (J>0)
        { if (L<=lower) break
            uni.slice.evals <<- uni.slice.evals + 1
            if (g(L, ...)<=logy) break
            L <- L - w
            J <- J - 1
        }
        
        while (K>0)
        { if (R>=upper) break
            uni.slice.evals <<- uni.slice.evals + 1
            if (g(R, ...)<=logy) break
            R <- R + w
            K <- K - 1
        }
    }
    
    # Shrink interval to lower and upper bounds.
    
    if (L<lower) 
    { L <- lower
    }
    if (R>upper)
    { R <- upper
    }
    
    # Sample from the interval, shrinking it on each rejection.
    
    repeat
    { 
        x1 <- runif(1,L,R)
        
        uni.slice.evals <<- uni.slice.evals + 1
        gx1 <- g(x1, ...)
        
        if (gx1>=logy) break
        
        if (x1>x0) 
        { R <- x1
        }
        else 
        { L <- x1
        }
    }
    
    # Return the point sampled, with its log density attached as an attribute.
    
    attr(x1,"log.density") <- gx1
    return (x1)
    
}

# unnormalized PDF
unnorm_logPDF <- function(x, tau2, ajk, mu, N, logscale = FALSE) {
    logf1 <- -N * log(1 + x)
    logf2 <- -ajk / (1+x) - (.5/tau2) * (x - mu)^2
    
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

tau2 <- 1
ajk <- 25
lp <- 10
ln <- 10.2
mu <- (lp + ln)/2
N <- 32

# tau2 <- 10
# ajk <- 546.339
# l1 <- 41.03722
# l2 <- 47.54683
# N <- 32

xs <- seq(-.99, 100, by = .001)

log_fxs <- unnorm_logPDF(xs, tau2, ajk, mu, N, logscale = TRUE)
log_fxs <- log_fxs - max(log_fxs)
# log_fxs <- log_fxs - log(normval)

fxs <- exp(log_fxs)

which.max(fxs)
fxs[which.max(fxs)]
xs[which.max(fxs)]

plot(xs, fxs, type = "l", xlim = c(0, 10))
abline(v = xs[which.max(fxs)], col = "red")

# other cubic solutions ---------------------------------------------------

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
    {gen <- 
        Runuran::tabl.new(unnorm_logPDF, 
                          lb = 0, ub = Inf, mode = xmode,
                          islog = TRUE,
                          tau2 = tau2, ajk = ajk, mu = mu, N = N, 
                          logscale = TRUE)
    rr <- Runuran::ur(gen, 1)},
    times = 10000
)

# create the generator object
# gen <- 
#     Runuran::pinv.new(unnorm_logPDF, lb = 0, ub = Inf,
#                       center = xmode, islog = TRUE,
#                       uresolution = 1.0e-10, # 1.0e-10 is default
#                       smooth = FALSE, # FALSE is default
#                       tau2 = tau2, ajk = ajk, l1 = l1, N = N, 
#                       logscale = TRUE)

gen <- 
    Runuran::tabl.new(unnorm_logPDF, 
                      lb = 0, ub = Inf, mode = xmode,
                      islog = TRUE,
                      tau2 = tau2, ajk = ajk, mu = mu, N = N, 
                      logscale = TRUE)

microbenchmark::microbenchmark(Runuran::ur(gen, 1),
                               times = 1000)

# draw a random sample
random_sample <- Runuran::ur(gen, 25000)
# plot the sample density
plot(density(random_sample, from = 0), xlim = c(0, 10))
lines(xs, fxs*0.25, type = "l", col = "red")

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


# try slice sampling ------------------------------------------------------

# uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL, ...)
# unnorm_logPDF(xs, tau2, ajk, mu, N, logscale = TRUE)

uni.slice.calls <- 0	# Number of calls of the slice sampling function
uni.slice.evals <- 0	# Number of density evaluations done in these calls

uni.slice(3, unnorm_logPDF, w = 1, m = Inf, lower = 0, upper = Inf, 
          tau2 = tau2, ajk = ajk, mu = mu, N = N, logscale = TRUE)

print(c(uni.slice.calls, uni.slice.evals))

x0 <- 3

slicesamples <- rep(NA, 50000)

for (i in 1:50000) {
    x0 <- uni.slice(x0, unnorm_logPDF, w = 1, m = Inf, lower = 0, upper = Inf, 
                    tau2 = tau2, ajk = ajk, mu = mu, N = N, logscale = TRUE)
    slicesamples[i] <- x0
}

plot(density(slicesamples, from = 0), col = 2)
lines(xs, fxs*0.25)
abline(v = xs[which.max(fxs)], col = "red")

microbenchmark::microbenchmark(
    {
        x0 <- 3
        x0 <- uni.slice(x0, unnorm_logPDF, w = 1, m = Inf, lower = 0, upper = Inf, 
                        tau2 = tau2, ajk = ajk, mu = mu, N = N, logscale = TRUE)
    },
    times = 10000
)

# test Lambda sampler and truncated distribution issues

logDiffExp <- function(a, b) {
    return(a + log(1 - exp(b - a)))
}

a <- .233
b <- 1
shape <- 508
rate <- 4035

x <- pgamma(b, shape, rate, log.p = TRUE)
y <- pgamma(a, shape, rate, log.p = TRUE)


# test logDiffExp ---------------------------------------------------------

x <- 3
y <- 2

log(exp(x) - exp(y))
x + log(1 - exp(y - x))

lb <- 4.354526e-03
ub <- 9.147465e-01

lp <- pgamma(lb, shape = 496, 2.475665e+05)
up <- pgamma(ub, shape = 496, 2.475665e+05)
mlog_lp <- -pgamma(lb, shape = 496, 2.475665e+05, log.p = TRUE)
mlog_up <- -pgamma(ub, shape = 496, 2.475665e+05, log.p = TRUE)

logu <- -truncdist::rtrunc(1, "exp", mlog_up, mlog_lp)
xi_jk <- qgamma(logu, shape = 496, 2.475665e+05, log.p = TRUE)

1/xi_jk - 1

lb <- 0
ub <- .6405
n1 <- 494
tjk <- 434.7

mlog_lp <- -pgamma(ub, shape = n1, rate = tjk, log.p = TRUE)
mlog_up <- -pgamma(lb, shape = n1, rate = tjk, log.p = TRUE)

# logu <- -truncdist::rtrunc(1, "exp", mlog_lp, mlog_up)
logu <- -LaplacesDemon::rtrunc(1, "exp", mlog_lp, mlog_up)
qgamma(logu, shape = n1, rate = tjk, log.p = TRUE)


theta_t <- -tjk + (n1 - 1)/mean(c(lb, ub))
rgamma(1, n1, tjk-theta_t)

lpdf <- function(x, shape, rate) {(shape-1)*log(x) - rate*x}
dlpdf <- function(x, shape, rate) {(shape-1)/x - rate}
Runuran::arou.new(lpdf, dpdf = dlpdf, lb = 0, ub = .6405, islog = TRUE, 
                  shape = n1, rate = tjk)
gen <- Runuran::ars.new(lpdf, dlogpdf = dlpdf, lb = 0, ub = .6405, 
                 shape = n1, rate = tjk)
Runuran::ur(gen, 1)

microbenchmark::microbenchmark(
    {gen <- Runuran::ars.new(lpdf, dlogpdf = dlpdf, lb = 0, ub = .6405, 
                             shape = n1, rate = tjk)
    Runuran::ur(gen, 1)},
    rgamma(1, 1, 1), times = 1000)


distr <- Runuran::udgamma(shape = n1, scale = 1/tjk, lb, ub)
# gen <- Runuran::pinvd.new(distr)
gen <- Runuran::arsd.new(distr)
Runuran::ur(gen, 1)

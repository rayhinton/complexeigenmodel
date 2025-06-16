intgd1 <- function(x, nu, p) {
    return(sin(x)^(nu-p))
}

nu <- 6
p <- 2

# demonstrate known integral identity
integrate(intgd1, 0, pi/2, nu = nu, p = p)
.5 * beta((nu-p + 1)/2, 1/2)

curve(intgd1(x, 3.75, p), 0, pi)
abline(v = pi/2)

# proposed integral identity
integrate(intgd1, 0, pi, nu = nu, p = p)
beta((nu-p + 1)/2, 1/2)




sqrt(pi) * gamma(5 - 8*nu + nu^2) * gamma(nu) / gamma(13/2 - 8*nu + nu^2)

plot(1:1000, (1:1000 - p + 1)^2 / 4 / (p - 1:1000))


xs <- seq(0, 1, length.out = 201)
plot(xs, xs*sqrt(1 - xs^(2/(nu - p + 1))))

lims <- c(0, pi/2)
sin(lims)
cos(lims)
Arg(cos(lims))
Arg(sin(lims))


curve((sin(x)*cos(x))^(nu-p+1), from = 0, to = pi/2, n = 201)
integrate(
    function(x) {(sin(x)*cos(x))^(nu-p+1)},
    0, pi/2)

curve((sin(x+pi/4)*cos(x+pi/4))^(nu-p+1), from = -pi/4, to = pi/4, n = 201)
integrate(
    function(x) {(sin(x+pi/4)*cos(x+pi/4))^(nu-p+1)},
    -pi/4, pi/4)


1 / (2*(nu - p + 1)) 

nu-p
2*(nu - p + 1)

integrate(
    function(x) (x^(nu - p + 1) * (1 - x^2)^( (nu - p)/2 )),
    0, 1)

integrate(
    function(x) {0.5 * x^((nu - p)/2) * (1 - x)^((nu - p)/2) },
    0, 1)

beta((nu-p)/2 + 1, (nu-p)/2 + 1) * 0.5


# non-Identity V ----------------------------------------------------------

nu <- 5
p <- 2

V <- matrix(c(.75, .1, .1, .25), ncol = 2)
eigen(V)
sum(diag(V))

x <- .1

# sin(x)*cos(x)^(nu - p + 1) *
#     (cos(x)^2 * (1 - 2*V[1,1]) -
#          )^(-p*nu/2)

xs <- seq(0, 1, length.out = 201)
plot(xs , xs * (1 - xs))

AA <- matrix(c(1.25, -1, 
               -1, -.25), 
             ncol = 2, byrow = TRUE)

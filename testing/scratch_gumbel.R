# Gumbel trick sampling

?sample

# number of discrete possibilities
J <- 101
# outcomes
thetas <- seq(0.01, 10, length.out = J)

# have a set of weights - these are not normalized to sum to 1
set.seed(8032025)
pweights <- dgamma(thetas, 3, 1)
summary(pweights)

points(thetas, pweights)

# weights on the log scale
logweights <- log(pweights) - 1000
summary(logweights)

# Sampling function, sample_gumbel
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")

gsample <- sample_gumbel(thetas, 1000, logweights)
plot(density(gsample))
points(thetas, pweights)
lines(thetas, dgamma(thetas, 3, 1), col = "red")

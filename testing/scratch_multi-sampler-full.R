# test multiple samplers, but with A, B constant

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Lambdak.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_sigmak2.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_V.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_Uk.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/FCD_AB.R")

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/matrix_distances.R")

# generate true parameters ------------------------------------------------

P <- 8
d <- 4
K <- 3
n_k <- rep(500, K)

Uk_warmup_its <- 5000
Uk_thin_its <- 100

gibbsIts <- 100
gibbsPrint <- 10

set.seed(27042025)

# A, B - constant
A0s <- seq(2*P, P, length.out = P)
B0s <- seq(2*P, P, length.out = d)

# V0
V0s <- runitary(P, P)
G0 <- V0s %*% diag(A0s) %*% t(Conj(V0s))

# Uk ~ CGB(A, B, V)
Ukinit <- runitary(P, d)

for (i in 1:Uk_warmup_its) {
    if (i %% 500 == 0) print(i)
    Ukinit <- rcmb(Ukinit, G0, B0s)
}
Uk0s <- array(NA, c(P, d, K))
for (i in 1:(Uk_thin_its * K)) {
    Ukinit <- rcmb(Ukinit, G0, B0s)
    if (i %% Uk_thin_its == 0) {
        print(i/Uk_thin_its)
        Uk0s[, , i/Uk_thin_its] <- Ukinit
    }
}

# Lambdak
Lambdak0s <- array(NA, c(d, K))
for (k in 1:K) {
    Lambdak0s[, k] <- rgamma(d, 1, 1) |> sort(decreasing = TRUE)
}

# sigmak2
sigmak20s <- rgamma(K, 1, 1)

# generate data according to the true parameters --------------------------

# Data: Yk ~ CW(nk, Gammak)
# Gammak = sigmak2 * (Uk * Lk * UkH + I_P)

data_list <- list()
for (k in 1:K) {
    Gammak0 <- sigmak20s[k] * (Uk0s[, , k] %*% diag(Lambdak0s[, k]) %*% 
                                   t(Conj(Uk0s[, , k])) + diag(P))
    data_list[[k]] <- rcomplex_wishart(n_k[k], Gammak0)
}

# initialize arrays to store posterior samples ----------------------------

U_ks_init <- array(NA, c(P, d, K))
Lambda_ks_init <- array(NA, c(d, K))

for (k in 1:K) {
    U_ks_init[, , k] <- runitary(P, d)
    Lambda_ks_init[, k] <- rgamma(d, 1, 1) |> sort(decreasing = TRUE)
}

U_ks <- array(NA, c(P, d, K, gibbsIts))
Lambda_ks <- array(NA, c(d, K, gibbsIts))
Vs <- array(NA, c(P, P, gibbsIts))
sigma_k2s <- matrix(NA, K, gibbsIts)
As <- matrix(NA, P, gibbsIts)
Bs <- matrix(NA, d, gibbsIts)

U_ks[, , , 1] <- U_ks_init
Lambda_ks[, , 1] <- Lambda_ks_init
Vs[, , 1] <- runitary(P, P)
sigma_k2s[, 1] <- rgamma(K, 1, 1)

# TODO think of a better way to initialize this
As[, 1] <- P*c(2, sort(runif(P-2, 1, 2), decreasing = TRUE), 1)
Bs[, 1] <- P*c(2, sort(runif(d-2, 1, 2), decreasing = TRUE), 1)

# put data, parameters into lists -----------------------------------------

param_list <- list(
    P = P,
    d = d,
    K = K,
    n_k = n_k,
    U_ks = U_ks[, , , 1],
    Lambda_ks = Lambda_ks[, , 1],
    sigma_k2s = sigma_k2s[, 1],
    Vs = Vs[, , 1],
    As = As[, 1],
    Bs = Bs[, 1]
)

# Gibbs sampling ----------------------------------------------------------
set.seed(28042025)
print(Sys.time())
for (s in 2:gibbsIts) {
    
    if (s %% gibbsPrint == 0) {
        print(paste0("s = ", s, ": ", Sys.time()))
    }
    
    param_list$U_ks <- U_ks[, , , s] <- 
        Uk_gibbs_densCovar(data_list, param_list)
    
    param_list$Lambda_ks <- Lambda_ks[, , s] <-
        Lambdak_gibbs_densCovar(data_list, param_list)
    
    param_list$sigma_k2s <- sigma_k2s[, s] <-
        sigmak2_gibbs_densCovar(data_list, param_list)
    
    param_list$Vs <- Vs[, , s] <-
        V_gibbs_densCovar(data_list, param_list, Imtol = 1e5 * .Machine$double.eps)
    
    ABs <- AB_gibbs_densCovar(data_list, param_list)
    param_list$As <- As[, s] <- ABs$As
    param_list$Bs <- Bs[, s] <- ABs$Bs
}

gibbsKeep <- seq(gibbsIts/2, gibbsIts, by = 1)

# assess sigmak samples ---------------------------------------------------

sigmak20s
summary(sigma_k2s[2, gibbsKeep])
quantile(sigma_k2s[2, gibbsKeep], c(.025, .975))

# assess V samples --------------------------------------------------------

avgVcovs <- matrix(0 + 0i, P, P)
for (s in 1:length(gibbsKeep)) {
    avgVcovs <- avgVcovs + Vs[, , gibbsKeep[s]] %*% diag(A0s) %*% t(Conj(Vs[, , gibbsKeep[s]]))
}
avgVcovs <- avgVcovs / length(gibbsKeep)
avgVvecs <- eigen(avgVcovs)$vectors

coli <- 2
cbind(avgVvecs[, coli], V0s[, coli])
cbind(avgVvecs[, coli], eigen(V0s %*% diag(A0s) %*% t(Conj(V0s)))$vectors[, coli])

eigen(V0s %*% diag(A0s) %*% t(Conj(V0s)))$vectors[, coli]

frame_distance(avgVvecs, V0s)

# assess Uk samples -------------------------------------------------------

avgUkcovs <- array(0 + 0i, c(P, P, K))
avgUkevecs <- array(0 + 0i, c(P, d, K))

for (k in 1:K) {
    for (s in 1:length(gibbsKeep)) {
        avgUkcovs[, , k] <- avgUkcovs[, , k] + U_ks[, , k, gibbsKeep[s]] %*% 
            diag(Lambda_ks[, k, gibbsKeep[s]]) %*% t(Conj(U_ks[, , k, gibbsKeep[s]]))
    }
    avgUkcovs[, , k] <- avgUkcovs[, , k] / length(gibbsKeep)
    avgUkevecs[, , k] <- eigen(avgUkcovs[, , k])$vectors[, 1:d]
}

coli <- 1
k <- 3

cbind(avgUkevecs[, coli, k], Uk0s[, coli, k])
frame_distance(avgUkevecs[, , k], Uk0s[, , k])

# save.image("/home/rayhinton/Documents/PhD_research/RA_time-series/code-experiments/multi-sampler-20280425.Rdata")


# assess A samples --------------------------------------------------------

# TODO troubleshooting: just plot or extract all the non-NA values - don't worry
# about burn-in

nonNAs <- !is.na(As[1, ])

plot(As[6, nonNAs])
summary(As[2, nonNAs])
summary(As[3, nonNAs])
summary(As[4, nonNAs])
summary(As[7, nonNAs])

apply(As[2:(P-1), gibbsKeep-1], 1, 
      quantile, probs = c(0.025, 0.975)) |> 
    t() |>
    cbind("true" = A0s[2:(P-1)])

# assess B samples --------------------------------------------------------

apply(Bs[2:(d-1), gibbsKeep-1], 1, 
      quantile, probs = c(0.025, 0.975)) |> 
    t() |>
    cbind("true" = B0s[2:(d-1)])

# assess Lambdak samples --------------------------------------------------

ki <- 1
Lambda_ks[, ki, gibbsKeep-1] |> 
    apply(1, quantile, probs = c(.025, .975)) |> 
    t() |> 
    cbind("true" = Lambdak0s[, ki])

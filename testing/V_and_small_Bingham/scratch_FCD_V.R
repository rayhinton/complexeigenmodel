# V FCD

# library(rstiefel)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/matrix_distances.R")


V_gibbs_densCovar <- function(data_list, param_list, 
                              Imtol = 2*.Machine$double.eps) {
    
    Vpar <- matrix(0 + 0i, param_list$P, param_list$P)
    for (k in 1:param_list$K) {
        Vpar <- Vpar + param_list$U_ks[, , k] %*% param_list$Bs %*% 
            t(Conj(param_list$U_ks[, , k]))
    }
    
    outV <- rcBingUP_gibbs(param_list$Vs, Vpar, param_list$As,
                           Imtol = Imtol)
    
    return(outV)
}

# set up parameters ------------------------------------------------------

# need
# U_k matrices, distributed as U_k ~ CGB(A, B, V), Pxd
# A diagonal matrix, PxP
# B diagonal matrix, dxd
# a true V_0 matrix, PxP

P <- 8
d <- 4
gibbsIts <- 1000

a <- seq(from = P*2, to = 1, length.out = P)
As <- diag(a)
b <- seq(from = P*2, to = 1, length.out = d)
Bs <- diag(b)

set.seed(19042025)
V_0 <- eigen(rcomplex_wishart(100, diag(P)))$vector
G_0 <- V_0 %*% As %*% t(Conj(V_0))

# multiple Uk matrices ----------------------------------------------------

# needs to be generated as a CGB, in order for the posterior to work

# TODO to generate multiple Uk from the same distribution, I could just run the same Gibbs sampler
# it would probably be best to take them, say, every 100th sample, to make sure they are not too correlate

K <- 3
initIts <- 5000

U_ks <- array(NA, c(P, d, K))

set.seed(19042025)
# burn in the CGB sampler
Uinit <- runitary(P, d)
Uinit <- rcmb(Uinit, G_0, Bs)
for (i in 1:initIts) {
    if (i %% 500 == 0) print(i)
    Uinit <- rcmb(Uinit, G_0, Bs)
}

# sample the U_k matrices every 100 iterations, to avoid correlation
for (i in 1:(K*100)) {
    Uinit <- rcmb(Uinit, G_0, Bs)
    if (i %% 100 == 0) {
        print(i)
        U_ks[, , i/100] <- Uinit
    }
}

data_list <- list()
param_list <- list(P = P,
                   d = d,
                   K = K,
                   U_ks = U_ks,
                   As = As,
                   Bs = Bs)

# compare the mode of FCD to the V_0 parameter - how "biased" is the FCD?
Vpar <- matrix(0 + 0i, P, P)
for (k in 1:K) {
    Vpar <- Vpar + U_ks[, , k] %*% Bs %*% t(Conj(U_ks[, , k]))
}

Vparevecs <- eigen(Vpar)$vectors
eigen(Vpar)$values

coli <- 2
cbind(Vparevecs[, coli], V_0[, coli])

frame_distance(Vparevecs, V_0)
procrustes_distance(Vparevecs, V_0)


# how large are the eigenvalues for the 2x2 parameters? -------------------

# N <- MASS::Null(X[, -ijs])
# 
# # transform the original parameters
# newA <- t(Conj(N)) %*% A %*% N
# newB <- B[ijs, ijs]

ijs <- c(7,8)
N <- MASS::Null(Vpar[, -ijs])
newG <- t(Conj(N)) %*% Vpar %*% N
newH <- As[ijs, ijs]

eigen(newG)
newH

my.rCbing.Op(newG, newH, Imtol = 100*.Machine$double.eps)

# initialize V samples ----------------------------------------------------

set.seed(20042025)
Vs <- array(NA, c(P, P, gibbsIts))
Vs[, , 1] <- eigen(rcomplex_wishart(P+1, diag(P)))$vector
Vcovs <- array(NA, c(P, P, gibbsIts))

param_list$Vs <- Vs[, , 1]

# run Gibbs sampler function ----------------------------------------------

print(Sys.time())
for (s in 2:gibbsIts) {
    if (s %% 100 == 0) print(paste0(s, ", ", Sys.time()))
    
    # Vs[, , s] <- rcBingUP_gibbs(Vs[, , s-1], Vpar, A,
                                # Imtol = 10000*.Machine$double.eps)
    
    Vs[, , s] <- V_gibbs_densCovar(data_list, param_list, 
                                   Imtol = 100000*.Machine$double.eps)
    param_list$Vs <- Vs[, , s]
    
    Vcovs[, , s] <- Vs[, , s] %*% param_list$As %*% t(Conj(Vs[, , s]))
}

# summarize samples -------------------------------------------------------

avgVcovs <- apply(Vcovs[, , (gibbsIts/2) : gibbsIts], c(1, 2), mean)
avgVs <- eigen(avgVcovs)$vectors

coli <- 4
cbind(avgVs[, coli], V_0[, coli])
cbind(avgVs[, coli], Vparevecs[, coli])

frame_distance(avgVs, V_0)
procrustes_distance(avgVs, V_0)

frame_distance(avgVs, Vparevecs)
procrustes_distance(avgVs, Vparevecs)

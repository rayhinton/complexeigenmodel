# scratch for rejection sampler for small complex square matrix Bingham dist.

# first, do it with no scaling

# output: a PxP matrix

# parameters:
# - A, PxP Hermitian matrix
# - B, PxP real diagonal matrix 

??wishart
# TODO perhaps write a separate complex Wishart sampler
library(cmvnorm)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

P <- 2

set.seed(26022025)
# A <- rcwis(P, diag(P))
A <- rcomplex_wishart(P, P, diag(P))
# B is a real diagonal matrix, with descending values
B <- diag(sort(rgamma(P, 1, 1), decreasing = TRUE))


# exp(logr)

rcBingUP(A, B)

reps <- 2000
nrejs <- rep(NA, reps)
tracestat <- rep(NA, reps)
As <- array(NA, c(P, P, reps))
Bs <- array(NA, c(P, P, reps))
# set.seed(26022025)
set.seed(8032025)

for (i in 1:reps) {
    if (i%% 20 == 0) print(i)

    # A <- rcwis(P, diag(P))
    A <- rcomplex_wishart(P, P, diag(P))
    # B is a real diagonal matrix, with descending values
    B <- diag(sort(rgamma(P, 1, 1), decreasing = TRUE))    
    
    ################################
    ### function
    ################################
    
    i_sample <- rcBingUP(A, B)
    

    ### end of function
    # i_sample <- list(X = X, nrej = nrej)
        
    nrejs[i] <- i_sample$nrej
    
    tracestat[i] <- Re(sum(diag(B %*% t(Conj(i_sample$X)) %*% A %*% i_sample$X)))
}

# investigate S, if there are issues with the most recently sampled value after an error
# are the eigenvalues real and positive?
eigen(S)$values
# is S Hermitian?
S - t(Conj(S))
isSymmetric(S)

S_round <- round(S, digits = tol_digits-1)
eigen(S_round)$values
S_round - t(Conj(S_round))
isSymmetric(S_round)

eigen(A)$values
isSymmetric(A)


# investigate last A, B parameters for issues -----------------------------

As[, , which.min(nrejs)] |> eigen()
As[, , which.max(nrejs)] |> eigen()

Bs[, , which.min(nrejs)] |> diag()
Bs[, , which.max(nrejs)] |> diag()

plot(ifelse(nrejs == 0, 0, log(nrejs)) + runif(reps, -.1, .1), log(Bs[1, 1, ] / Bs[2, 2, ]))

# look at summaries -------------------------------------------------------

summary(nrejs)
quantile(nrejs, c(0.95, 0.99))
plot(density(Re(tracestat)))
plot(log(nrejs), log(tracestat))


eigen(i_sample$X)
t(Conj(i_sample$X[1:4])) %*% i_sample$X[1:4]
# Hoff's function from paper
### 
rwish<-function(nu,M,cholM=chol(M)) 
{ 
    Z<-matrix(rnorm(nu*dim(M)[1]),nrow=nu,ncol=dim(M)[1]) 
    Y<-Z%*%cholM 
    t(Y)%*%Y 
} 
### 
### 
rZ.AB<-function(A,B) { 
    #rescale A and B, and select gamma, to improve acceptance rate 
    diag(A)<-diag(A)-min(eigen(A)$val) 
    diag(B)<-diag(B)-min(diag(B)) 
    mA<-max(eigen(A)$val) ; mB<-max(B) 
    gd<-(dim(A)[1]+1)/(2*mB) + mA 
    del<-max(eigen(A)$val) +.5 ; gam<-gd/del 
    A<-A/gam ; B<-B*gam ; 
    vA<-diag(eigen(A)$val) 
    
    S<- solve( diag(del,nrow=dim(A)[1] ) - A )/2 
    nu<-dim(A)[1]+1 

    rej<-TRUE ; nrej<-0 
    while(rej) { 
        tmp<-eigen(rwish(nu,S)) 
        Z<-tmp$vec ; L<-diag(tmp$val) ; Z<-Z%*%diag((-1)^rbinom(dim(A)[1],1,.5)) 
        D<-B-L 
        lrr<- sum(diag(( D%*%t(Z)%*%A%*%Z)) ) - sum( -sort(diag(-D))*diag(vA)) 
        rej<- ( log(runif(1))> lrr ) 
        nrej<-nrej+rej 
    } 
    
    list(Z=Z,nrej=nrej) 
    
} 

P <- 4

set.seed(26022025)
A <- rWishart(1, P, diag(P))[ , , 1]
# A is Hermitian
A <- A %*% t(Conj(A))
# B is a real diagonal matrix, with descending values
B <- diag(sort(rgamma(P, 1, 1), decreasing = TRUE))

rZ.AB(A, B)

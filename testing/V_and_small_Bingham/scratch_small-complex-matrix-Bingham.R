# scratch for rejection sampler for small complex square matrix Bingham dist.

# output: a PxP matrix

# parameters:
# - A, PxP Hermitian matrix
# - B, PxP real diagonal matrix 

??wishart
# library(cmvnorm)
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

P <- 2

### this A matrix has eigenvalues with a large difference
Aevals <- c(8000, 20)
(Av <- 1/sqrt(2) * matrix(c(1, 1, -1, 1), ncol = 2))
(A <- Av %*% diag(Aevals) %*% t(Av))
eigen(A)

B <- diag(c(8000, 20))

# a matrix that will be eigenvalues of sampled "covariance" matrix
C <- diag(c(2, 1))

# set.seed(26022025)
# # A <- rcwis(P, diag(P))
# A <- rcomplex_wishart(P, P, diag(P))
# # B is a real diagonal matrix, with descending values
# B <- diag(sort(rgamma(P, 1, 1), decreasing = TRUE))

my.rCbing.Op(A, B)

reps <- 1000
nrejs <- rep(NA, reps)
tracestat <- rep(NA, reps)
As <- array(NA, c(P, P, reps))
Bs <- array(NA, c(P, P, reps))
Us <- array(NA, c(P, P, reps))
Covs <- array(NA, c(P, P, reps))
# set.seed(26022025)
set.seed(8032025)

for (i in 1:reps) {
    if (i%% 50 == 0) print(i)

    # A <- rcomplex_wishart(P, P, diag(P))
    # B is a real diagonal matrix, with descending values
    # B <- diag(sort(rgamma(P, 1, 1), decreasing = TRUE))
    
    As[, , i] <- A
    Bs[, , i] <- B
    
    ################################
    ### function
    ################################
    
    # i_sample <- rcBingUP(A, B)
    i_sample <- my.rCbing.Op(A, B)
    
    ### end of function
    # i_sample <- list(X = X, nrej = nrej)
        
    nrejs[i] <- i_sample$nrej
    
    tracestat[i] <- Re(sum(diag(B %*% t(Conj(i_sample$X)) %*% A %*% i_sample$X)))
    
    Us[, , i] <- i_sample$X
    Covs[, , i] <- i_sample$X %*% C %*% t(Conj(i_sample$X))
}

# look at summaries -------------------------------------------------------

summary(nrejs)
quantile(nrejs, c(0.95, 0.99))
quantile(nrejs, c(0.95, 0.99), na.rm = TRUE)

eigen(apply(Covs, c(1, 2), mean))
eigen(A)

eigen(i_sample$X)
t(Conj(i_sample$X[1:4])) %*% i_sample$X[1:4]

# Hoff's function from paper --------------------------------------------
### 
rwish<-function(nu,M,cholM=chol(M)) 
{ 
    Z<-matrix(rnorm(nu*dim(M)[1]),nrow=nu,ncol=dim(M)[1]) 
    Y<-Z%*%cholM 
    t(Y)%*%Y 
} 
### 

###
my.rbing.Op <- function(A,B) {
        #simulate from the bingham distribution on O(p) 
        #having density proportional to etr(B t(U)%*%A%*%U ) 
        #using the rejection sampler described in Hoff(2009)
        #this only works for small matrices, otherwise the sampler
        #will reject too frequently
        
        ### assumes B is a diagonal matrix with *decreasing* entries 
        
        b<-diag(B) ; bmx<-max(b) ; bmn<-min(b)  
        if(bmx>bmn)
        { 
            A<-A*(bmx-bmn) ; b<-(b-bmn)/(bmx -bmn)
            vlA<-eigen(A)$val  
            diag(A)<-diag(A)-vlA[1]
            vlA<-eigen(A)$val  
            
            nu<- max(dim(A)[1]+1,round(-vlA[length(vlA)]))
            del<- nu/2
            M<- solve( diag(del,nrow=dim(A)[1] ) - A )/2
            
            rej<-TRUE
            cholM<-chol(M)
            nrej<-0
            while(rej)
            {
                Z<-matrix(rnorm(nu*dim(M)[1]),nrow=nu,ncol=dim(M)[1])
                Y<-Z%*%cholM ; tmp<-eigen(t(Y)%*%Y)
                U<-tmp$vec%*%diag((-1)^rbinom(dim(A)[1],1,.5)) ; L<-diag(tmp$val)
                D<-diag(b)-L
                lrr<- sum(diag(( D%*%t(U)%*%A%*%U)) ) - sum( -sort(diag(-D))*vlA)
                rej<- ( log(runif(1))> lrr )
                # nrej<-nrej+1
                nrej <- nrej + rej
            }
        }
        if(bmx==bmn) { U<-rustiefel(dim(A)[1],dim(A)[1]) } 
        return(list(X=U, nrej = nrej))
    }

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

P <- 2

set.seed(26022025)
# A <- rWishart(1, P, diag(P))[ , , 1]
# A is Hermitian
# A <- A %*% t(Conj(A))
# B is a real diagonal matrix, with descending values
# B <- diag(sort(rgamma(P, 1, 1), decreasing = TRUE))

rZ.AB(A, B)
my.rbing.Op(A, B)

# repeated testing
reps <- 1000
nrejs <- rep(NA, reps)
tracestat <- rep(NA, reps)
As <- array(NA, c(P, P, reps))
Bs <- array(NA, c(P, P, reps))
Us <- array(NA, c(P, P, reps))
Covs <- array(NA, c(P, P, reps))
# set.seed(26022025)
set.seed(8032025)

for (i in 1:reps) {
    if (i%% 100 == 0) print(i)
    
    # A <- rcwis(P, diag(P))
    # A <- rcomplex_wishart(P, P, diag(P))
    # A <- rWishart(1, P+1, diag(P))[, , 1]
    # B is a real diagonal matrix, with descending values
    # B <- diag(sort(rgamma(P, 1, 1), decreasing = TRUE))    
    
    ################################
    ### function
    ################################
    
    # i_sample <- rZ.AB(A, B)
    i_sample <- my.rbing.Op(A, B)
    # As[, , i] <- A
    # Bs[, , i] <- B
    
    ### end of function
    # i_sample <- list(X = X, nrej = nrej)
    
    nrejs[i] <- i_sample$nrej
    
    Us[, , i] <- i_sample$X
    Covs[, , i] <- i_sample$X %*% C %*% t(i_sample$X)
}

summary(nrejs)
quantile(nrejs, c(0.95, 0.99))

eigen(apply(Covs, c(1, 2), mean))
eigen(A)

# show the 10 worst numbers of rejections
nrejs[order(nrejs, decreasing = TRUE)[1:10]]

# show the A and B matrices with the worst rejections
whichworst <- 1
nrejs[order(nrejs, decreasing = TRUE)[whichworst]]
As[, , order(nrejs, decreasing = TRUE)[whichworst]]
Bs[, , order(nrejs, decreasing = TRUE)[whichworst]]

eigen(As[, , order(nrejs, decreasing = TRUE)[whichworst]])$values

B12ratios <- rep(NA, reps)
for (i in 1:reps) {
    B12ratios[i] <- Bs[1, 1, i] / Bs[2, 2, i]
}
B12ratios[order(nrejs, decreasing = TRUE)[1:20]]

# 3x3 test ----------------------------------------------------------------

##### define A and B parameters

A3v <- matrix(c(1/sqrt(2), -1/sqrt(2), 0,
                1/sqrt(2), 1/sqrt(2), 0,
                0, 0, 1), ncol = 3, byrow = TRUE)
# A3evals <- c(140, 80, 20)
A3evals <- c(80, 50, 20)

A3 <- A3v %*% diag(A3evals) %*% t(Conj(A3v))
B3 <- diag(c(7, 4, 1))

my.rCbing.Op(A3, B3, istatus = 100)

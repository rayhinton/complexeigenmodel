# scratch for rejection sampler for small complex square matrix Bingham dist.

# first, do it with no scaling

# output: a PxP matrix

# parameters:
# - A, PxP Hermitian matrix
# - B, PxP real diagonal matrix 

??wishart
# TODO perhaps write a separate complex Wishart sampler
library(cmvnorm)

P <- 2

set.seed(26022025)
A <- rcwis(P, diag(P))
# A is Hermitian
A <- A %*% t(Conj(A))
# B is a real diagonal matrix, with descending values
B <- diag(sort(rgamma(P, 1, 1), decreasing = TRUE))

rcBing <- function(A, B) {

    # rescale A and B
    # TODO verify the 2 points below
    # - shifting the diagonals does not change the distribution
    # - shifting A and B by c, 1/c does not change the distribution
    # shifting
    diag(A) <- diag(A) - min(eigen(A, only.values = TRUE)$values)
    diag(B) <- diag(B) - min(diag(B))
    # scaling
    # maximum eigenvalue of A is the first one
    mA <- eigen(A, only.values = TRUE)$values[1]
    mB <- max(B)
    # TODO trying a different factor here
    gd <- (dim(A)[1] + 1) / (2*mB) + mA # this is Hoff
    # gd <- (dim(A)[1] + 1) / (2*mB) + mA/2
    # TODO understand why this factor is calculated now
    del1 <- max(eigen(A, only.values = TRUE)$values) + 0.5
    gam <- gd/del1
    
    # scale
    A <- A / gam
    B <- B * gam
    
    # get and store the eigenvalues of A
    Aevals <- eigen(A, only.values = TRUE)$values
    # TODO calculate del
    # del <- max(Aevals) + 0.5
    
    # TODO trying a different nu
    nu <- (dim(A)[1] + 1) # this is Hoff
    # nu <- dim(A)[1]
    
    S <- solve(diag(del1, nrow = P) - A)
    
    rej <- TRUE
    nrej <- 0
    while (rej) {
        W <- rcwis(nu, S)
        Weigen <- eigen(W)
        
        
        # TODO need to randomly multiply columns by +- 1
        # TODO do microbenchmark - which is faster? multiplying by a diagonal matrix of +- 1, or using vector recycling?
        X <- Weigen$vectors
        L <- Weigen$values
        
        X <- X %*% diag((-1)^rbinom(dim(A)[1], 1, .5)) 
        
        D <- sort(diag(B) - L, decreasing = TRUE)
        
        # the products must be over elements of decreasing order. Eigenvalues from eigen are in decreasing order by default. D is sorted above.
        lambda <- sum(Aevals * D)
        
        # calculate the log acceptance ratio
        # part of it might be slightly complex, due to numerical issues
        # so, use Re
        # TODO but, I should double-check everything, to justify that it should be real
        # TODO double-check I am using right B-L matrix, not necessarily the sorted D
        logr <- Re(sum(diag(
            diag(B - L) %*% t(Conj(X)) %*% A %*% X))) - 
            lambda
        
        # Accept with probability r
        # draw u from Unif(0, 1)
        # TODO can I compare, say, log(u) to the logr? 
        rej <- log(runif(1)) > logr
        nrej <- nrej + rej
    }
    
    return(list(X = X, nrej = nrej))
}
# exp(logr)

rcBing(A, B)

reps <- 1000
nrejs <- rep(NA, reps)
tracestat <- rep(NA, reps)
set.seed(26022025)

for (i in 1:reps) {
    if (i%% 20 == 0) print(i)

    A <- rcwis(P, diag(P))
    # A is Hermitian
    A <- A %*% t(Conj(A))
    # B is a real diagonal matrix, with descending values
    B <- diag(sort(rgamma(P, 1, 1), decreasing = TRUE))    
    
    i_sample <- rcBing(A, B)
    
    nrejs[i] <- i_sample$nrej
    
    tracestat[i] <- Re(sum(diag(B %*% t(Conj(i_sample$X)) %*% A %*% i_sample$X)))
}

summary(nrejs)
plot(density(Re(tracestat)))

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

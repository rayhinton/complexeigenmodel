# test Hoff's code for A, B estimation


# Hoff functions ----------------------------------------------------------

hadj<-function(alpha,beta,w) {
    tmp<-outer(alpha,alpha,"-")*outer(beta,beta,"-")
    sm1<-sum(1/ tmp[upper.tri(tmp)] )/4
    sm2<-sum(1/ tmp[upper.tri(tmp)]^2 )/4 + sm1^2/2
    1/( 1 + sm1/w + sm2/w^2 )
}
#####


#####
rabw.mh<-function(a,b,w,UA,V,gs=200,wa=1/p,wb=1/p^2 ) {
    K<-dim(UA)[3] ; p<-dim(UA)[1] 
    M<-diag(rep(K,p))
    for(k in 1:K) { M<- M - (t(V)%*%UA[,,k])^2 }
    
    # Original Hoff code
    wp<-rgamma(1,choose(p,2)*K/2+wa, t(a)%*%M%*%b + wb )
    # RJH adjustment: change wb to wa*wb in the 2nd parameter
    # wp<-rgamma(1,choose(p,2)*K/2+ wa, t(a)%*%M%*%b + wa*wb )
    if(log(runif(1))< K*( log(hadj(a,b,wp)) - log(hadj(a,b,w)) ) ) {w<-wp}
    # RJH the sampler seems to always be rejecting new samples; try without MH for w
    # w <- wp
    
    
    for(j in sample(2:(p-1))) {
        ap<-a
        ub<-min( a[ (1:p)<j ] )
        lb<-max( a[ (1:p)>j ] )
        as<-seq(lb,ub,length=gs)
        lp<-  -w*as*(M[j,]%*%b)+.5*K*apply(log(abs(outer(as,a[-j],"-"))),1,sum)
        ap[j]<-sample(as,1,prob=exp(lp-max(lp)))
        if(log(runif(1))< K*( log(hadj(ap,b,w)) - log(hadj(a,b,w)) ) ){a<-ap }
        
    }
    
    for(j in sample(2:(p-1))) {
        bp<-b
        ub<-min( b[ (1:p)<j ] )
        lb<-max( b[ (1:p)>j ] )
        bs<-seq(lb,ub,length=gs)
        lp<-  -w*bs*(M[,j]%*%a)+.5*K*apply(log(abs(outer(bs,b[-j],"-"))),1,sum)
        bp[j]<-sample(bs,1,prob=exp(lp-max(lp)))
        if(log(runif(1))< K*( log(hadj(a,bp,w)) - log(hadj(a,b,w)) ) ) {b<-bp}
    }
    
    list(a=a, b=b,w=w )
}

# fix Hoff bing -----------------------------------------------------------

my.rbing.matrix.gibbs <- function(A,B,X) {
        #simulate from the matrix bmf distribution as described in Hoff(2009) 
        #this is one Gibbs step, and must be used iteratively
        
        ### assumes B is a diagonal matrix with *decreasing* entries 
        
        m<-dim(X)[1] ;  R<-dim(X)[2]
        if(m>R)
        {
            for(r in sample( seq(1,R,length=R)))
            {
                N<- rstiefel::NullC(X[,-r])
                An<-B[r,r]*t(N)%*%(A)%*%N 
                X[,r]<-N%*% rstiefel::rbing.vector.gibbs(An,t(N)%*%X[,r])
            }
        }
        
        #If m=R then the fc of one vector given all the others is 
        #just +-1 times the vector in the null space. In this case, 
        #the matrix needs to be updated at least two columns at a 
        #time. 
        if(m==R)
        {
            rs <- sample(1:R) |> matrix(ncol = 2)
            # for(s in seq(1,R,length=R))
            for (s in 1:nrow(rs))
            {
                # r<-sort( sample(seq(1,R,length=R),2) )
                r <- sort(rs[s, ])
                N<- rstiefel::NullC( X[,-r]  )
                An<- t(N)%*%A%*%N
                #X[,r]<-N%*%rbing.O2(An,B[r,r]) 
                X[,r]<-N%*% rstiefel::rbing.Op(An,B[r,r]) 
            }
        }
        X
    }

# setup -------------------------------------------------------------------

P <- 4
K <- 100

eta0 <- 2
tau02 <- 1/1000

bing_burn <- 5000
bing_status <- 500

gibbsIts <- 10000

set.seed(3052025)
# alpha0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)
# beta0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)

alpha0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)
beta0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)

w0 <- rgamma(1, 1, 1)

A0 <- diag(sqrt(w0) * alpha0)
B0 <- diag(sqrt(w0) * beta0)

V0 <- rstiefel::rustiefel(P, P)

G0 <- V0 %*% A0 %*% t(V0) 

Uks <- array(NA, c(P, P, K))

Ukinit <- rstiefel::rustiefel(P, P)

for (i in 1:bing_burn) {
    if (i %% bing_status == 0) print(i)
    Ukinit <- rstiefel::rbing.matrix.gibbs(G0, B0, Ukinit)
    # Ukinit <- my.rbing.matrix.gibbs(G0, B0, Ukinit)
}

for (k in 1:(K*100)) {
    Ukinit <- rstiefel::rbing.matrix.gibbs(G0, B0, Ukinit)
    # Ukinit <- my.rbing.matrix.gibbs(G0, B0, Ukinit)
    if (k %% 100 == 0) {
        print(k/100)
        Uks[, , k/100] <- Ukinit
    }
}

avgCovs <- matrix(0, P, P)

for (k in 1:K) {
    # Uks[, , k] <- rstiefel::rbing.Op(G0, B0)
    avgCovs <- avgCovs + Uks[, , k] %*% diag(P:1) %*% t(Uks[, , k])
}

avgCovs <- avgCovs/K
avgVecs <- eigen(avgCovs)$vectors

G0vecs <- eigen(G0)$vectors

G0vecs
V0

coli <- 4
cbind(avgVecs[, coli], G0vecs[, coli])

Uks[, , 1]
Uks[, , K]

# posterior mean of w, based on A0, B0 being fixed ------------------------

M <- matrix(0, P, P)

for (k in 1:K) {
    M <- M + (t(V0) %*% Uks[, , k])^2
}

aw <- eta0/2 + K/2 * choose(P, 2)
bw <- eta0*tau02/2 + t(matrix(alpha0, ncol=1)) %*% (diag(K, nrow = P) - M) %*% matrix(beta0, ncol=1)

# "posterior" mean of w, given A0 and B0
aw / bw
w0
# sampling ----------------------------------------------------------------

alphas <- array(NA, c(P, gibbsIts))
betas <- array(NA, c(P, gibbsIts))
ws <- rep(NA, gibbsIts)

# alphas[, 1] <- seq(1, 0, length.out = P)
# betas[, 1] <- seq(1, 0, length.out = P)

# alphas[, 1] <- c(1, sort(runif(P-2), decreasing = TRUE), 0)
# betas[, 1] <- c(1, sort(runif(P-2), decreasing = TRUE), 0)

alphas[, 1] <- alpha0
betas[, 1] <- beta0

ws[1] <- rgamma(1, 1, 1)

for (s in 2:gibbsIts) {
    if (s %% 500 == 0) print(s)
    
    abwsample <- rabw.mh(alphas[, s-1], 
                         betas[, s-1], 
                         ws[s-1], 
                         Uks, V0, gs = 200, 
                         # wa = 1/P, wb = 1/P^2)
                         # wa = 1, wb = 1/1000)
                         wa = eta0/2, wb = tau02)
    
    # alphas[, s] <- abwsample$a
    # betas[, s] <- abwsample$b
    
    alphas[, s] <- alpha0
    betas[, s] <- beta0
    
    ws[s] <- abwsample$w
}

noburn <- (gibbsIts/2):gibbsIts

plot(ws, type = "l")
abline(h = w0, col = "red")

apply(alphas[2:(P-1), noburn], 1,
      quantile, probs = c(.025, .975)) |> 
    t() |> 
    cbind("true" = alpha0[2:(P-1)])

alphas[, 1]

apply(betas[2:(P-1), noburn], 1,
      quantile, probs = c(.025, .975)) |> 
    t() |> 
    cbind("true" = beta0[2:(P-1)])

betas[, 1]

summary(ws[noburn])
w0

curve(dgamma(x, 1, 1/1), from = 0, to = 1000, n = 1001)

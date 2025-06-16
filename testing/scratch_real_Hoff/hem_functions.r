#####
tr<-function(X) { sum(diag(X)) }
#####


#####
Null<-function(M){
    tmp <- qr(M)
    set <- if (tmp$rank == 0)
            1:nrow(M)  #modified from package MASS
    else -(1:tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
                   }
#####


#####
rmvnorm<-function(n,mu,Sigma) {
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
                               }
#####


#####
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
 
  wp<-rgamma(1,choose(p,2)*K/2+wa, t(a)%*%M%*%b + wb )
  
  if (doMH) {
      if(log(runif(1))< K*( log(hadj(a,b,wp)) - log(hadj(a,b,w)) ) ) {w<-wp}
  } else {
      w <- wp
  }


    
  for(j in sample(2:(p-1))) {
    ap<-a
    ub<-min( a[ (1:p)<j ] )
    lb<-max( a[ (1:p)>j ] )
    as<-seq(lb,ub,length=gs)
    lp<-  -w*as*(M[j,]%*%b)+.5*K*apply(log(abs(outer(as,a[-j],"-"))),1,sum)

    # RJH troubleshooting
    aNAs <- sum(is.na(exp(lp-max(lp))))
    if (aNAs > 0) {
        print( paste0("s = ", s, "; a NAs:", aNAs) )
    }
    
    ap[j]<-sample(as,1,prob=exp(lp-max(lp)))
    if (doMH) { 
     if(log(runif(1))< K*( log(hadj(ap,b,w)) - log(hadj(a,b,w)) ) ){a<-ap }
    } else {
        a <- ap
    }

                             }

  for(j in sample(2:(p-1))) {
    bp<-b
    ub<-min( b[ (1:p)<j ] )
    lb<-max( b[ (1:p)>j ] )
    bs<-seq(lb,ub,length=gs)
    lp<-  -w*bs*(M[,j]%*%a)+.5*K*apply(log(abs(outer(bs,b[-j],"-"))),1,sum)
    
    # RJH troubleshooting
    bNAs <- sum(is.na(exp(lp-max(lp))))
    if (bNAs > 0) {
        print( paste0("s = ", s, "; b NAs:", bNAs) )
    }
    
    bp[j]<-sample(bs,1,prob=exp(lp-max(lp)))
    if (doMH) {
        if(log(runif(1))< K*( log(hadj(a,bp,w)) - log(hadj(a,b,w)) ) ) {b<-bp}
    } else {
        b <- bp
    }

                             }

  list(a=a, b=b,w=w )
                              }
#####



#####
rUL.gibbs<-function(YYA,A,B,V,UA,LA,upool=TRUE,nu0=2,s20=1)
{
 M<-V%*%A%*%t(V)*upool 
  for(k in 1:K) 
  {
    YY<-YYA[,,k] ; U<-UA[,,k] ; L<-LA[,,k]
    a<- (N[k]-1)/2 + nu0/2 ; b<-diag(t(U)%*%YY%*%U)/2 + nu0*s20/2
    for(j in sample(1:p))
    {
      lb<-suppressWarnings( max(c(0, 1/diag(L)[ (1:p)<j ]) ) )
      ub<-suppressWarnings( min( 1/diag(L)[ (1:p)>j ]) )
      L[j,j]<-1/qgamma(runif(1,pgamma(lb,a,b[j]),pgamma(ub,a,b[j])),a,b[j])
    }
    for(ss in 1:p) 
    {
      r<-sample(1:p,2)
      Nr<-Null(U[,-r])
      M1<-t(Nr)%*%( B[r[1],r[1]]*M - .5*YY/L[r[1],r[1]] )%*%Nr
      M2<-t(Nr)%*%( B[r[2],r[2]]*M - .5*YY/L[r[2],r[2]] )%*%Nr
      U[,r]<-Nr%*%rU1U2.M1M2.th(M1,M2)
    }
    UA[,,k]<-U ; LA[,,k]<-L
  }
  list(UA=UA, LA=LA)
}
#####


#####
rV.fc<-function(S,A,V){

    for(ss in 1:p) {
      r<-sample(1:p,2)
      Nr<-Null(V[,-r])
      M1<-t(Nr)%*%( A[r[1],r[1]]*S  )%*%Nr
      M2<-t(Nr)%*%( A[r[2],r[2]]*S  )%*%Nr
      V[,r]<-Nr%*%rU1U2.M1M2.th(M1,M2)
                     }
      V                 } 
#####


#####
THETA<-seq(0,2*pi,length=200)
C2T<- cos(THETA)^2
S2T<- sin(THETA)^2
CST<- cos(THETA)*sin(THETA)

rU1U2.M1M2.th<-function(M1,M2) {
  lpz<-C2T*(M1[1,1]+M2[2,2])+S2T*(M1[2,2]+M2[1,1])+
       CST*(M1[1,2]+M1[2,1]-M2[1,2]-M2[2,1] )

    # RJH troubleshooting
    thetaNAs <- sum(is.na(exp(lpz-max(lpz))))
    if (thetaNAs > 0) {
        print( paste0("s = ", s, "; theta NAs: ", thetaNAs) )
    }
    
  theta<-sample(THETA,1,prob=exp(lpz-max(lpz)) )
  x<-(-1)^rbinom(1,1,.5)
  matrix( c(cos(theta),sin(theta),x*sin(theta),-x*cos(theta)),2,2)
                                }
#####

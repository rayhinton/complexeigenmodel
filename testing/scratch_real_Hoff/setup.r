##### starting values
A<-B<-diag( (p:1-1),nrow=p) 
UA<-array(dim=c(p,p,K)) ; LA<-array(diag(p),dim=c(p,p,K)) ; V<-matrix(0,p,p)
for(k in 1:K){
  tmp<-eigen(YYA[,,k]/(N[k]-1));UA[,,k]<-tmp$vec 
  lv<-tmp$val
  ns<-sum(lv<1/p) 
  lv[ (lv<1/p)]<-seq(1/p,1/p^2,length=ns) 
  LA[,,k]<-diag(lv)
  V<- V + YYA[,,k]/(N[k]-1)
              }
V<-eigen(V)$vec
a<-diag(A)/max(A*1)
b<-diag(B)/max(B*1)
w<-max(A)

UPE<-V
UAE<-UA
#####


##### storage
ABPS<-array(dim=c(p,p,(NSCAN+NBURN)/ODENS))
VPS<-array(0,dim=c(p,p,(NSCAN+NBURN)/ODENS)) 
MSD<-matrix(0,0,2)
PRES<-matrix(nrow=0,ncol=p)
X2PS<-matrix(0,p,p)
X<-diag(p)
SDZ<-matrix(nrow=0,ncol=p)
W<-TSTAT<-NULL

LPS<- array(dim=c(K,p,(NSCAN+NBURN)/ODENS))
#####


##### functions for output

tstat<-function(VUA) {

  VU2<-VUA$V*0
  for(k in 1:K) { VU2<-VU2 + ( t(VUA$V)%*%VUA$UA[,,k] )^2 }
  range(diag( VU2/K ))
                        }
##
gof.stat<-function(V,UA ) {
  sm<-matrix(0,p,p)
  for(k in 1:K) { sm<-sm + (t(V)%*%UA[,,k])^2 }
  diag( sm/K )
                           }



sUV<-function(U,V) {
  G<- (t(U)%*%V)^2
  mxs<-NULL
  for(l in 1:(dim(G)[1]-1)) {
    rm<-apply(G,1,max)
    cm<-apply(G,2,max)
    mx<-max(rm)
    G<-G[ (rm!=mx),(cm!=mx) ]
    mxs<-c(mxs,mx)
                              }
  c(mxs,G)
                       }




YYA.sim<-function(N,UA,LA ) {
 K<-length(N) ; p<-dim(UA)[1]
 YYA<-array(dim=c(p,p,K)) 
 
    for(k in 1:K){
    Y<-rmvnorm(N[k],rep(0,p),UA[,,k]%*%LA[,,k]%*%t(UA[,,k]))
    Y<-t( t(Y)-apply(Y,2,mean))
    YYA[,,k]<-t(Y)%*%Y
                  }
 YYA           }


VUA.emp<-function(YYA,N) {
  K<-length(N) ; p<-dim(YYA)[1]
  V<-matrix(0,p,p) ; UA<-array(dim=c(p,p,K))

   for(k in 1:K){
    V<-V + YYA[,,k]/(N[k]-1) 
    UA[,,k]<- eigen( YYA[,,k] )$vec
                }
    V<-eigen(V)$vec    
  list(V=V,UA=UA ) 
                            }

VVPS<-matrix(0,p,p)

 VUA.obs<-VUA.emp(YYA,N)
 ostat<-gof.stat(VUA.obs$V,VUA.obs$UA)


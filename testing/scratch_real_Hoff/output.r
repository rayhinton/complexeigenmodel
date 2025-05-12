###
SO<- matrix(0,p,p) ;SO[upper.tri(SO,diag=T)]<-1 ; DO<-solve(SO)
dab<- c(  (DO%*%diag(A))[-p] , (DO%*%diag(B))[-p] )
if( all(dab!=0) ) {MSD<-rbind(MSD, c(mean(log(dab)),sd(log(dab)))  ) }
if( any(dab==0)) {MSD<-rbind(MSD, c(mean(dab),sd(dab))  ) }
###


###
ABPS[,,s/ODENS]<- outer(diag(A),diag(B))
VPS[,,s/ODENS]<- V
LPS[,,s/ODENS]<- t(apply(LA,3,diag))

print(apply(LPS,c(1,2),mean,na.rm=T))
###

###
X<- rV.fc(A,B,X)
X2PS<-X2PS + X^2
###

###

VUA.sim<-VUA.emp( YYA.sim(N,UA,LA), N) 
sstat<- gof.stat(VUA.sim$V,VUA.sim$UA)
PRES<-rbind(PRES,sstat)
TSTAT<-rbind(TSTAT,range(sstat))
###


###  text
cat(s," ",round(apply(PRES,2,median),2)," ", 
          round(diag(X2PS)/dim(MSD)[1],2)," ", 
          round(apply(MSD,2,median),2),"\n" )


cat(w, round(diag(A),2),"  ",round(diag(B),2),"\n")
W<-c(W,w) 



### plots
if(plt==TRUE) {
par(mfrow=c(3,2))

plot(MSD[,1],type="l",col="green",ylim=range(MSD));lines(MSD[,2],col="blue")
abline(h=apply(MSD,2,median),col=c("light green","light blue"))
  
image( 1:p,1:p, t(X2PS[p:1,1:p]),col=gray( (0:20)/20 ))

lpemp<-log( ostat/(1-ostat))

LPRES<-log(PRES/(1-PRES))
if(dim(LPRES)[1]>1)  {
  plot(density(LPRES[,1]),xlim=range(c(LPRES,lpemp))) 
  abline(v=lpemp[1])
  for(j in 2:p) { 
    lines(density(LPRES[,j]),col=j)  ; abline(v=lpemp[j],col=j)
                 }
                       }
###

LTSTAT<-log( TSTAT/(1-TSTAT))
ltobs<-log(range(ostat)/(1-range(ostat)))

plot(c(1,length(TSTAT)/2), range(c(LTSTAT),ltobs),type="n")
#points(LTSTAT[,1],pch="-")
#points(LTSTAT[,2],pch="+")
points(LTSTAT[,1],pch=16,col="blue")
points(LTSTAT[,2],pch=16,col="red")
abline(h=ltobs)

plot(log(W)) ; abline(h=median(log(W)))
hist(W,prob=T) ; lines(1:max(W),dgamma(1:max(W),wa,wb) )

cat(mean(LTSTAT[,1]>ltobs[1]), mean(LTSTAT[,2]<ltobs[2]),   "\n")


print( round( X2PS/(s/ODENS),2)) 
      }

if(s%%(ODENS*10) == 0 ) {
  ss<-s/ODENS
  dput(list(PRES=PRES,W=W,VPS=VPS[,,1:ss],ABPS=ABPS[,,1:ss],TSTAT=TSTAT,
            LPS=LPS[,,1:ss]), fname)
                         }



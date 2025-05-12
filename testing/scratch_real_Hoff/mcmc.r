#####
seed<-1
set.seed(seed)
#####


#####
source("hem_functions.r")
#####


##### data
#dat<-dget("data_vole")

# Hoff defaults
wa <- 1
wb <- 1/1000 
nu0 <- 1
s20 <- 1

P <- 8
K <- 4
n_k <- c(82, 70, 58, 54)

eta0 <- 2*wa
tau02 <- wb/wa

nu0 <- 1
sigma02 <- s20

bing_burn <- 5000
bing_status <- 500

# gibbsIts <- 10000

set.seed(3052025)
# alpha0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)
# beta0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)

alpha0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)
beta0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)

# w0 <- rgamma(1, 1, 1)
w0 <- rgamma(1, eta0/2, eta0 * tau02 / 2)

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

Lambdak0 <- array(NA, c(P, K))

for (k in 1:K) {
	Lambdak0[, k] <- 1 / sort(rgamma(P, nu0/2, sigma02/2))
}

YYA <- array(NA, c(P, P, K))

for (k in 1:K) {
	YYA[, , k] <- rWishart(1, n_k[k]-1, 
		Uks[, , k] %*% diag(Lambdak0[, k]) %*% t(Uks[, , k]))[, , 1]
}

N <- n_k
# YYA <- dat$YYA
p <- P
# K <- dim(YYA)[3]

#####


##### MCMC output
NSCAN<-10000
NBURN<-5000
ODENS<-1000
source("setup.r")
UPOOL<-TRUE
LPOOL<-FALSE
plt<-FALSE

# RJH I copied this earlier in the code, so that the data-generating parameters 
# are based on these defaults.
# wa<-1 ; wb<-1/1000 ; nu0=1 ; s20=1

fname<-paste("out_seed",seed,"upool",1*UPOOL,"_",Sys.Date(),sep="")
#####

ws <- rep(0, NBURN+NSCAN)
as <- array(NA, c(P, NBURN+NSCAN))
bs <- array(NA, c(P, NBURN+NSCAN))

##### MCMC
for(s in 1:(NBURN+NSCAN)) {
    if (s %% 1000 == 0) print(s)

  if( UPOOL==TRUE )
  {
    ### update V
    S<-matrix(0,p,p) ; for(k in 1:K) { S<-S+UA[,,k]%*%B%*%t(UA[,,k]) }
    V<-rV.fc(S,A,V)
    ###
  
    ### update A and B 
    tmp<-rabw.mh(a,b,w,UA,V,gs=200,wa=wa,wb=wb)
    a<-tmp$a ; b<-tmp$b ; w<-tmp$w
    ws[s] <- w
    as[, s] <- a
    bs[, s] <- b
    A<-diag(a*sqrt(w)) ; B<-diag(b*sqrt(w))
    ###  
  }

  ### update U and L
  tmp<-rUL.gibbs(YYA,A,B,V,UA,LA,upool=UPOOL,nu0=nu0,s20=s20) 
  UA<-tmp$UA ; LA<-tmp$LA
  ###  


  ### output 
  if(s%%ODENS==0) {  source("output.r") }
  ###

}


save.image(paste("image_",fname,sep=""))

gibbsKeep <- (NBURN+1):(NBURN+NSCAN)

# assess ws ----------------------------------------------------------

w0

mean(ws[gibbsKeep])
quantile(ws[gibbsKeep], probs = c(.025, .975))


# assess as ---------------------------------------------------------------

apply(as[2:(P-1), gibbsKeep], 1, 
      quantile, probs = c(0.025, 0.975)) |> 
    t() |>
    cbind("true" = alpha0[2:(P-1)])

# assess bs ---------------------------------------------------------------
apply(bs[2:(P-1), gibbsKeep], 1, 
      quantile, probs = c(0.025, 0.975)) |> 
    t() |>
    cbind("true" = beta0[2:(P-1)])

# assess Lambdas ----------------------------------------------------------

ki <- 2

Lambdak0[, ki]
diag(LA[, , ki])

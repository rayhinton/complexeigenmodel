
# functions from rstiefel -------------------------------------------------

rustiefel <-
    function(m,R) {
        #simulate uniformly from V_{R,m}  
        #see Chikuse 
        #note R is given second, the result is an m*R matrix
        X<-matrix(rnorm(m*R),m,R)
        tmp<-eigen(t(X)%*%X)
        
        return(X%*%( tmp$vec%*%sqrt(diag(1/tmp$val,nrow=R))%*%t(tmp$vec) ))
    }

NullC <-
    function(M) {
        #modified from package "MASS" 
        #MASS version : Null(matrix(0,4,2))  returns a 4*2 matrix
        #this version : NullC(matrix(0,4,2)) returns diag(4)
        
        tmp <- qr(M)
        set <- if (tmp$rank == 0L)
            1L:nrow(M)
        else -(1L:tmp$rank)
        
        return(qr.Q(tmp, complete = TRUE)[, set, drop = FALSE])
    }

rbing.vector.gibbs <-
    function(A,x) {
        #simulate from the vector bmf distribution as described in Hoff(2009) 
        #this is one Gibbs step, and must be used iteratively
        evdA<-eigen(A,symmetric=TRUE)
        E<-evdA$vec
        l<-evdA$val
        
        y<-t(E)%*%x
        x<-E%*%ry_bing(y,l)
        
        return(x/sqrt(sum(x^2)))
        #One improvement might be a rejection sampler 
        #based on a mixture of vector mf distributions. 
        #The difficulty is finding the max of the ratio.
    }

rbing.matrix.gibbs <-
    function(A,B,X) {
        #simulate from the matrix bmf distribution as described in Hoff(2009) 
        #this is one Gibbs step, and must be used iteratively
        
        ### assumes B is a diagonal matrix with *decreasing* entries 
        
        m<-dim(X)[1] ;  R<-dim(X)[2]
        if(m>R)
        {
            for(r in sample( seq(1,R,length=R)))
            {
                N<-NullC(X[,-r])
                An<-B[r,r]*t(N)%*%(A)%*%N 
                X[,r]<-N%*%rbing.vector.gibbs(An,t(N)%*%X[,r])
            }
        }
        
        #If m=R then the fc of one vector given all the others is 
        #just +-1 times the vector in the null space. In this case, 
        #the matrix needs to be updated at least two columns at a 
        #time. 
        if(m==R)
        {
            for(s in seq(1,R,length=R))
            {
                r<-sort( sample(seq(1,R,length=R),2) )
                N<-NullC( X[,-r]  )
                An<- t(N)%*%A%*%N
                #X[,r]<-N%*%rbing.O2(An,B[r,r]) 
                X[,r]<-N%*%rbing.Op(An,B[r,r]) 
            }
        }
        
        return(X)
    }

rbing.Op <-
    function(A,B) {
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
                nrej<-nrej+1
            }
        }
        if(bmx==bmn) { U<-rustiefel(dim(A)[1],dim(A)[1]) } 
        
        return(U)
    }

#####
seed<-2
set.seed(seed)
#####


#####
source("hem_functions.r")
#####

NSCAN<-50000
NBURN<-5000
ODENS<-5000
UPOOL<-TRUE
LPOOL<-FALSE
plt<-FALSE

constantLk <- TRUE
constantUk <- FALSE
constantV <- TRUE
constantw <- TRUE
constantalpha <- TRUE
constantbeta <- TRUE

##### data
#dat<-dget("data_vole")

doMH <- TRUE

# Hoff defaults
wa <- 1
wb <- 1/1000 
nu0 <- 1
s20 <- 1

P <- 8
K <- 2
# n_k <- c(82, 70, 58, 54)
n_k <- sample(c(82, 70, 58, 54), K, replace = TRUE) + 
    sample(-5:5, K, replace = TRUE)

eta0 <- 2*wa
tau02 <- wb/wa

nu0 <- 1
sigma02 <- s20

bing_burn <- 5000
bing_status <- 500

# gibbsIts <- 10000

set.seed(14052025)
# alpha0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)
# beta0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)

alpha0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)
beta0 <- c(1, sort(runif(P-2), decreasing = TRUE), 0)

# w0 <- rgamma(1, 1, 1)
w0 <- rgamma(1, eta0/2, eta0 * tau02 / 2)

A0 <- diag(sqrt(w0) * alpha0)
B0 <- diag(sqrt(w0) * beta0)

# V0 <- rstiefel::rustiefel(P, P)
V0 <- rustiefel(P, P)

G0 <- V0 %*% A0 %*% t(V0) 

Uks <- array(NA, c(P, P, K))

# Ukinit <- rstiefel::rustiefel(P, P)
Ukinit <- rustiefel(P, P)

for (i in 1:bing_burn) {
    if (i %% bing_status == 0) print(i)
    Ukinit <- rbing.matrix.gibbs(G0, B0, Ukinit)
    # Ukinit <- rstiefel::rbing.matrix.gibbs(G0, B0, Ukinit)
    # Ukinit <- my.rbing.matrix.gibbs(G0, B0, Ukinit)
}

for (k in 1:(K*100)) {
    Ukinit <- rbing.matrix.gibbs(G0, B0, Ukinit)
    # Ukinit <- rstiefel::rbing.matrix.gibbs(G0, B0, Ukinit)
    # Ukinit <- my.rbing.matrix.gibbs(G0, B0, Ukinit)
    if (k %% 100 == 0) {
        print(k/100)
        Uks[, , k/100] <- Ukinit
    }
}

Lambdak0 <- array(NA, c(P, K))
Lk0 <- array(NA, c(P, P, K))

for (k in 1:K) {
	# Lambdak0[, k] <- 1 / sort(rgamma(P, nu0/2, sigma02/2))
	Lambdak0[, k] <- 1 / sort(rexp(P))
	Lk0[, , k] <- diag(Lambdak0[, k])
}

YYA <- array(NA, c(P, P, K))

for (k in 1:K) {
	YYA[, , k] <- rWishart(1, n_k[k]-1, 
		Uks[, , k] %*% diag(Lambdak0[, k]) %*% t(Uks[, , k]))[, , 1]
}

YYA[, , 1] |> eigen()

N <- n_k
# YYA <- dat$YYA
p <- P
# K <- dim(YYA)[3]

#####

source("setup.r")

##### MCMC output

# RJH I copied this earlier in the code, so that the data-generating parameters 
# are based on these defaults.
# wa<-1 ; wb<-1/1000 ; nu0=1 ; s20=1

fname<-paste("out_seed",seed,"upool",1*UPOOL,"_",Sys.Date(),sep="")
#####

ws <- rep(0, NBURN+NSCAN)
as <- array(NA, c(P, NBURN+NSCAN))
bs <- array(NA, c(P, NBURN+NSCAN))
UA_S <- array(NA, c(P, P, K, NBURN+NSCAN))

##### MCMC
for(s in 1:(NBURN+NSCAN)) {
    
    if( UPOOL==TRUE )
    {
        ### update V
        if (constantV) {
            V <- V0
        } else { 
            S<-matrix(0,p,p) ; for(k in 1:K) { S<-S+UA[,,k]%*%B%*%t(UA[,,k]) }
            V<-rV.fc(S,A,V)
        }
        ###
        
        ### update A and B 
        tmp<-rabw.mh(a,b,w,UA,V,gs=200,wa=wa,wb=wb)
        if (constantalpha) { a <- alpha0 } else {a<-tmp$a}
        if (constantbeta) {b <- beta0} else {b<-tmp$b}
        if (constantw) {w <- w0} else {w<-tmp$w}
        # a<-tmp$a ; b<-tmp$b ; w<-tmp$w
        ws[s] <- w
        as[, s] <- a
        bs[, s] <- b
        A<-diag(a*sqrt(w)) ; B<-diag(b*sqrt(w))
        ###  
    }
    
    ### update U and L
    tmp<-rUL.gibbs(YYA,A,B,V,UA,LA,upool=UPOOL,nu0=nu0,s20=s20) 
    if (constantLk) {LA <- Lk0} else {LA<-tmp$LA}
    if (constantUk) {UA <- Uks} else {UA<-tmp$UA}
    # UA<-tmp$UA ; LA<-tmp$LA
    UA_S[, , , s] <- UA
    ###  
    
    ### output 
    if(s%%ODENS==0) {  source("output.r") }
    ###
    
}


# save.image(paste("image_",fname,sep=""))

load("image_out_seed2upool1_2025-05-18_2005h")

gibbsKeep <- (NBURN+1):(NBURN+NSCAN)
gibbsKeep <- seq(NBURN+1, NBURN+NSCAN, by = 25)

# assess Uks --------------------------------------------------------------

# manysamples <- list(list(Uk_S = Uks[, , whichk, ], accCount = 1))
# new_param_list <- list(P = P, d = d, Lambda_ks = diag(Lambdaks[, , whichk]), 
                       # Uk0 = Uk0s[, , whichk])
whichk <- 1

manysamples <- list(list(Uk_S = UA_S[, , whichk, ], accCount = whichk))
new_param_list <- list(P = P, d = P, Lambda_ks = Lambdak0[, whichk],
                       Uk0 = Uks[, , whichk])

# summarize_Uk <- function(sample_list, param_list, burnin) {
summarystuff <- summarize_Uk(manysamples, new_param_list, NBURN)

# avg_tracePlots <- function(sample_list, avgs_list, param_list, tracePlotEvery = 10) {
# dist_vals <- avg_tracePlots(manysamples, summarystuff, new_param_list, 
                            # tracePlotEvery = 5)

dev.new()
dist_vals <- avg_tracePlots(manysamples, summarystuff, new_param_list,
                            tracePlotEvery = 1)
somemcmc <- coda::mcmc(dist_vals$d_to_avgUk[, 1], thin = 1)
coda::effectiveSize(somemcmc)

dev.off()

# assess ws ----------------------------------------------------------

w0

mean(ws[gibbsKeep])
quantile(ws[gibbsKeep], probs = c(.025, .5, .975))

plot(ws, type = "l")

ws_mcmc <- coda::mcmc(ws, start = 1, thin = 1)
coda::effectiveSize(ws_mcmc)
# assess as ---------------------------------------------------------------

as[, 1]

apply(as[2:(P-1), gibbsKeep], 1, 
      quantile, probs = c(0.025, .5, 0.975)) |> 
    t() |>
    cbind("true" = alpha0[2:(P-1)])

plot(as[2, gibbsKeep], type = "l", main = "alpha2")
plot(as[3, gibbsKeep], type = "l", main = "alpha3")
plot(as[6, gibbsKeep], type = "l", main = "alpha6")
plot(as[7, gibbsKeep], type = "l", main = "alpha7")

# assess bs ---------------------------------------------------------------

bs[, 1]

apply(bs[2:(P-1), gibbsKeep], 1, 
      quantile, probs = c(0.025, .5, 0.975)) |> 
    t() |>
    cbind("true" = beta0[2:(P-1)])

plot(bs[3, gibbsKeep], type = "l")

# assess Lambdas ----------------------------------------------------------

ki <- 2

Lambdak0[, ki]
diag(LA[, , ki])

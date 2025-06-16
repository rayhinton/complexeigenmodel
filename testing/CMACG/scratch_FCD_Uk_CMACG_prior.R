# scratch CMACG prior for Uk

library(EigenR)
library(MASS)
library(foreach)

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_rcomplex_wishart.R")

# functions ---------------------------------------------------------------

mvGamma <- function(c, m, betaf, logscale = FALSE) {
    term1 <- .25*m*(m-1)*betaf*log(pi)
    term2 <- sum(lgamma(c - ((1:m)-1)/2 * betaf))
    
    if (logscale) {
        return(term1 + term2)
    } else {
        return(exp(term1 + term2))
    }
}

rscnorm <- function(n) {
    return(rnorm(n, 0, 1/sqrt(2)) + 
               1i * rnorm(n, 0, 1/sqrt(2)))
}

runif_stiefel <- function(P, d, betaf) {
    # real
    if (betaf == 1) {
        X1 <- matrix(rnorm(P*d), ncol = d)
        # complex
    } else if (betaf == 2) {
        X1 <- matrix(rscnorm(P*d), ncol = d)        
    }
    
    # take QR decomposition - not guaranteed to be unique due to numerical methods
    X1qr <- qr(X1)
    QX1 <- qr.Q(X1qr)
    # extract sign of the diagonal elements of R
    D <- sign(Re(diag(qr.R(X1qr))))
    # transform by Q = QS, where S = diag(D). this transformation guarantees Q is unique.
    
    # Qfinal <- t(t(QX1) * D) # faster than matrix multiplication, for larger matrices
    Qfinal <- QX1 %*% diag(D)
    
    return(Qfinal)
}

rCMACG <- function(nrow, ncol, Sigma) {
    # simulate X ~ CMN(0, I_P, I_d),
    # calculate Z by transforming X with Cholesky decomposition of Sigma0
    # let H_Z be orthogonal portion of polar decomposition of Z
    X <- matrix(rscnorm(nrow*ncol), nrow, ncol)
    
    # check if this is the upper or lower triang version
    # want the matrix such that Sigma = cholSigma %*% t(Conj(cholSigma))
    # Eigen_chol returns the upper triangular portion
    cholSigma <- EigenR::Eigen_chol(Sigma)
    
    Z <- t(Conj(cholSigma)) %*% X 
    sqrtZHZ <- EigenR::Eigen_sqrt( t(Conj(Z)) %*% Z )
    return(Z %*% solve(sqrtZHZ))
}

dCMACG <- function(X, Sigma, logscale = TRUE) {
    m <- nrow(X)
    r <- ncol(X)
    
    invSigma <- solve(Sigma)
    
    term1 <- -r * log( EigenR::Eigen_det(Sigma) )
    term2 <- -m * log( EigenR::Eigen_det( t(Conj(X)) %*% invSigma %*% X ) )
    
    log_density <- Re(term1 + term2)
    
    if (logscale) {
        return(log_density)
    } else {
        return(exp(log_density))
    }
}

# generate a proposal via a small rotation parameterized with Cayley transform
Uk_MH_Cayley <- function(Us, tau_U, doCayleyZeros, CayleyZeroProb) {
    P <- nrow(Us)
    
    Sp <- matrix(0, P, P)
    Sp[upper.tri(Sp)] <- tau_U*rscnorm(choose(P, 2))

    if (doCayleyZeros) {
        Sp[upper.tri(Sp)] <- 
            Sp[upper.tri(Sp)] * rbinom(choose(P, 2), 1, 1 - CayleyZeroProb)
    }

    for(j in 1:(P-1)) {
        for (i in (j+1):P) {
            Sp[i, j] <- -Conj(Sp[j, i])
        }
    }
    Rp <- (diag(P) - Sp) %*% solve(diag(P) + Sp)
    Up <- Rp %*% Us
    
    return(Up)
}

frob_dist <- function(A, B, returnDists = FALSE) {
    diffAB <- A - B

    if (returnDists) {
        sqDists <- Re(diag( t(Conj(diffAB)) %*% diffAB ))
        Fdist <- Re(sum(sqDists))
        return(list(Fdist = sqrt(Fdist),
                    sqDists = sqDists))
    } else {
        return(norm(diffAB, "F"))
    }
}

fast_evec_Frob_stat <- function(X, Y) {
    k <- ncol(X)
    return(
        sqrt(2*k - 2*sum(Mod( diag(t(Conj(X)) %*% Y)) ))
    )
}

evec_Frob_stat <- function(X, Y, returnDists = FALSE, returnMats = FALSE) {
    Rmat <- diag(complex(modulus = 1, 
                         argument = -Arg(diag(t(Conj(X)) %*% Y))))  
    dist_list <- frob_dist(X, Y %*% Rmat, returnDists = returnDists)
    
    if (returnMats) {
        return(list(dist_obj = dist_list,
               Xopt = X,
               Yopt = Y %*% Rmat))
    } else {
        return(dist_list)
    }
}

# data generation ---------------------------------------------------------

genData_Uk_CMACG <- function(P, d, n_k, tracenorm, customEvals, sigmaevals, 
                             sortLambda, 
                             parameterseed = NULL, dataseed = NULL) {
    
    if (is.null(parameterseed)) {
        parameterseed <- (Sys.time() |> as.integer()) - 17e8 
    }
    if (is.null(dataseed)) {
        dataseed <- (Sys.time() |> as.integer()) - 17e8 - sample(10, 1)
    }
    

    set.seed(parameterseed)
    # Sigma0 ~ diffuse Complex Wishart, PxP
    Sigma0 <- rcomplex_wishart(P+1, P, diag(P))
    if (tracenorm) {
        Sigma0 <- P * Sigma0 / Re(sum(diag(Sigma0)))
    }
    if (customEvals) {
        SigmaV <- eigen(Sigma0)$vectors
        Sigma0 <- SigmaV %*% diag(sigmaevals) %*% t(Conj(SigmaV))
    }
    
    invSigma0 <- solve(Sigma0)
    
    ###
    # Uk0 ~ CMACG(Sigma0)
    ###
    Uk0 <- rCMACG(P, d, Sigma0)
    
    # Lambdak0: diagonal matrix of decreasing positive entries
    # Lambdak0 <- diag(d:1)
    # Omegak0 <- solve( solve(Lambdak0) + diag(d) )
    
    Omegak0 <- runif(d, 0, 1)
    Lambdak0 <- Omegak0 / (1 - Omegak0)
    if (sortLambda) {
        Lambdak0 <- diag( sort(Lambdak0, decreasing = TRUE) )
    } else {
        Lambdak0 <- diag(Lambdak0)
    }
    Omegak0 <- solve( solve(Lambdak0) + diag(d) )
    
    # sigmak02: random positive scalar
    # Jeffreys prior would use lower values of shape and rate, to approach 1/x
    sigmak02 <- rgamma(1, shape = 1, rate = 1)
    
    # Yk ~ CW(n_k, sigmak02 * (Uk0 %*% Lambdak0 %*% t(Conj(Uk0))) )
    Gamma0 <- sigmak02 * (Uk0 %*% Lambdak0 %*% t(Conj(Uk0)) + diag(P))
    set.seed(dataseed)
    Yk <- rcomplex_wishart(n_k, P, Gamma0)
    
    # create output
    data_list <- list(Yk)
    param_list <- list(
        Sigmas = Sigma0,
        invSigmas = invSigma0,
        Lambda_ks = diag(Lambdak0),
        Omega_ks = diag(Omegak0),
        sigma_k2s = sigmak02,
        P = P,
        d = d,
        n_k = n_k,
        Uk0 = Uk0
    )
    
    return(list(data_list = data_list, param_list = param_list,
                parameterseed = parameterseed, dataseed = dataseed))
}

# MH proposal: small rotation with Cayley transform and CMACG prior -------

sample_Uk_MH_Cayley <- 
    function(data_list, param_list, S_its, tau_U, initialUkseed = NULL, 
             initialUk = NULL, doCayleyZeros = FALSE, CayleyZeroProb = 0.5,
             doStatus = FALSE, statusNum = 100) {
    
    # extract stuff from param_list
    P <- param_list$P
    d <- param_list$d
    n_k <- param_list$n_k
    Sigma0 <- param_list$Sigmas
    # invSigma0 <- param_list$invSigmas
    invSigma0 <- solve(Sigma0)
    # Omegak0 <- diag(param_list$Omega_ks)
    Lambdak0 <- diag(param_list$Lambda_ks)
    Omegak0 <- solve( solve(Lambdak0) + diag(d) )
    sigmak02 <- param_list$sigma_k2s
    # extract data
    Yk <- data_list[[1]]
    
    # initialize array to hold Uk samples
    Uk_S <- array(NA, c(P, d, S_its))
    # set seed, if provided
    if (!is.null(initialUkseed)) { 
        set.seed(initialUkseed) 
    } else {set.seed((Sys.time() |> as.integer()) - 17e8)}
    
    if (is.null(initialUk)) {
        Uk_S[, , 1] <- runif_stiefel(P, d, betaf = 2)
    } else {
        Uk_S[, , 1] <- initialUk
    }
    Us <- Uk_S[, , 1]
    
    accCount <- rep(TRUE, S_its)
    for (s in 2:S_its) {
        if (doStatus & (s %% (S_its/statusNum) == 0)) print(paste0("s = ", s))
        
        # propose U by small rotation of Us
        Up <- Uk_MH_Cayley(Us, tau_U, doCayleyZeros, CayleyZeroProb)
        
        # calculate terms based on priors
        tracep <- Re(sum(diag( Up %*% Omegak0 %*% t(Conj(Up)) %*% Yk )))
        logdetp <- log(Re(EigenR::Eigen_det( t(Conj(Up)) %*% invSigma0 %*% Up )))
        traces <- Re(sum(diag( Us %*% Omegak0 %*% t(Conj(Us)) %*% Yk )))
        logdets <- log(Re(EigenR::Eigen_det( t(Conj(Us)) %*% invSigma0 %*% Us )))
        
        # calculate acceptance ratio
        (logr <- tracep/sigmak02 - P*logdetp - traces/sigmak02 + P*logdets)
        
        # accept or reject
        if (log(runif(1)) <= logr) {
            Uk_S[, , s] <- Up
            Us <- Up
        } else {
            Uk_S[, , s] <- Us
            accCount[s] <- FALSE
        }
    } # end of s sampling loop
    
    return(list(Uk_S = Uk_S, accCount = accCount))
}

# calculate summaries of samples ------------------------------------------

summarize_Uk <- function(sample_list, param_list, burnin) {

    S_its <- dim(sample_list[[1]]$Uk_S)[3]    
    P <- param_list$P
    d <- param_list$d
    
    if (burnin > 1) {
        itsKeep <- (burnin+1) : S_its 
    } else {
        itsKeep <- (floor(S_its/2) + 1):S_its
    }
    
    avgs_df <- data.frame(matrix(NA, length(sample_list), 5))
    names(avgs_df) <- c("accRate", "d_avgCov", "d_avgCovEvec", "d_avgProj", 
                        "d_avgUk")
    
    avgCovs <- list()
    avgCovEvecs <- list()
    avgProjs <- list()
    avgUks <- list()
    
    for (j in 1:length(sample_list)) {
        avgs_df$accRate[j] <- mean(sample_list[[j]]$accCount)
        
        Lambdak0 <- diag(param_list$Lambda_ks)
        Uk0 <- param_list$Uk0
        
        projSum <- matrix(0, P, P)
        covSum <- matrix(0, P, P)
        for (s in itsKeep) {
            N <- length(itsKeep)
            projSum <- projSum + 
                (1/N) * sample_list[[j]]$Uk_S[, , s] %*% 
                t(Conj(sample_list[[j]]$Uk_S[, , s]))
            covSum <- covSum +
                (1/N) * sample_list[[j]]$Uk_S[, , s] %*% Lambdak0 %*% 
                t(Conj(sample_list[[j]]$Uk_S[, , s]))
        }
        
        # eigenvectors of average projection and covariance matrices
        evec_projSum <- eigen(projSum)$vectors[, 1:d]
        evec_covSum <- eigen(covSum)$vectors[, 1:d]
        
        avgs_df$d_avgCov[j] <- frob_dist(covSum, 
                                         Uk0 %*% Lambdak0 %*% t(Conj(Uk0)))
        avgs_df$d_avgCovEvec[j] <- evec_Frob_stat(evec_covSum, Uk0)
        
        # average projection matrix
        avgProj <- evec_projSum %*% t(Conj(evec_projSum))
        avgs_df$d_avgProj[j] <- frob_dist(avgProj, Uk0 %*% t(Conj(Uk0)))
        
        ### natural mean of eigenvectors, from Chikuse 2012, p. 237
        S <- apply(sample_list[[j]]$Uk_S[, , itsKeep], c(1, 2), mean)
        HS <- S %*% EigenR::Eigen_sqrt(solve(t(Conj(S)) %*% S))
        avgs_df$d_avgUk[j] <- evec_Frob_stat(HS, Uk0)
        
        avgCovs[[j]] <- covSum
        avgCovEvecs[[j]] <- evec_covSum
        avgProjs[[j]] <- projSum
        avgUks[[j]] <- HS
    } # end of sample list average loop
    
    return(list(avgs_df = avgs_df, avgCovs = avgCovs, avgCovEvecs = avgCovEvecs, 
                avgProjs = avgProjs, avgUks = avgUks))
}

# separate plotting function

avg_tracePlots <- function(sample_list, avgs_list, param_list, burnin,
                           tracePlotEvery = 10) {
    
    numchains <- length(sample_list)
    
    S_its <- dim(sample_list[[1]]$Uk_S)[3] 
    if (burnin < 1) {
        burnin <- floor(S_its * burnin)
    }
    
    Lambdak0 <- diag(param_list$Lambda_ks)
    
    S_its <- dim(sample_list[[1]]$Uk_S)[3]
    # metric by metric
    par(mfrow = c(2, 2), cex.axis = 1.5, cex.lab = 1.5)
    
    toPlot <- seq(1, S_its, by = tracePlotEvery)
    
    d_to_avgCov <- d_to_avgCovEvec <- d_to_avgProj <- d_to_avgUk <- 
        matrix(NA, nrow = length(toPlot), ncol = numchains)
    for (j in 1:numchains) {
        for (sind in 1:length(toPlot)) {
            s <- toPlot[sind]
            thisUk <- sample_list[[j]]$Uk_S[, , s]
            
            # Covariances
            thisCov <- thisUk %*% Lambdak0 %*% t(Conj(thisUk))
            d_to_avgCov[sind, j] <- norm(avgs_list$avgCovs[[j]] - thisCov, "F")
            
            # Covariance eigenvectors
            d_to_avgCovEvec[sind, j] <- 
                fast_evec_Frob_stat(avgs_list$avgCovEvecs[[j]],
                               eigen(thisCov)$vectors[, 1:(param_list$d)])
            
            # Projection matrices
            d_to_avgProj[sind, j] <-
                norm(avgs_list$avgProjs[[j]] - thisUk %*% t(Conj(thisUk)), "F")
            
            # Uk matrices
            d_to_avgUk[sind, j] <- 
                fast_evec_Frob_stat(avgs_list$avgUks[[j]], thisUk)
        }
    }
    
    # print quantiles
    # want to do it with respect to some burn-in choice
    # the index in each d_to matrix that corresponds to a raw sample index greater than the burn-in is
    burninLine <- (which(toPlot > burnin)[1])
    itsKeep <- burninLine:nrow(d_to_avgCov)
    q_probs <- c(.025, .25, .5, .75, .975)
    
    print(paste0("dim of d_to_avgCov: ", dim(d_to_avgCov)))
    
    chain_quantiles <- rbind(
        apply(d_to_avgCov[itsKeep, , drop = FALSE], 2, quantile, 
              probs = q_probs) |> t(),
        apply(d_to_avgCovEvec[itsKeep, , drop = FALSE], 2, quantile, 
              probs = q_probs) |> t(),
        apply(d_to_avgProj[itsKeep, , drop = FALSE], 2, quantile, 
              probs = q_probs) |> t(),
        apply(d_to_avgUk[itsKeep, , drop = FALSE], 2, quantile, 
              probs = q_probs) |> t())
    
    rownames(chain_quantiles) <- 
        rep(c("d_to_avgCov", "d_to_avgCovEvec", "d_to_avgProj", "d_to_avgUk"),
            each = numchains)
    
    # plot Covariance distances
    plot(NA, xlab = "", ylab = "Frob. dist.",
         xlim = c(1, length(toPlot)),
         ylim = c(min(d_to_avgCov), max(d_to_avgCov)),
         main = "dists. to Avg. Cov.")
    for(j in 1:numchains) {
        points(1, d_to_avgCov[1, j], col = j, cex = 2)
        lines(x = 1:length(toPlot), y = d_to_avgCov[, j], col = j)
    }
    abline(v = burninLine, lty = 2)
    legend("topright", legend = c(1:numchains,"burn-in"), 
           col = c(1:numchains, 1), lwd = 2, lty = c(rep(1, numchains), 2),
           cex = 1.25)
    # plot Covariance eigenvector distances
    plot(NA, xlab = "", ylab = "Axis Frob. dist.",
         xlim = c(1, length(toPlot)),
         ylim = c(min(d_to_avgCovEvec), max(d_to_avgCovEvec)),
         main = "dists. to Avg. Cov. evec.")
    for(j in 1:numchains) {
        points(1, d_to_avgCovEvec[1, j], col = j, cex = 2)
        lines(x = 1:length(toPlot), y = d_to_avgCovEvec[, j], col = j)
    }
    abline(v = burninLine, lty = 2)
    # plot Projection matrix distances
    plot(NA, xlab = "", ylab = "Frob. dist.",
         xlim = c(1, length(toPlot)),
         ylim = c(min(d_to_avgProj), max(d_to_avgProj)),
         main = "dists. to Avg. Proj. matrix")
    for(j in 1:numchains) {
        points(1, d_to_avgProj[1, j], col = j, cex = 2)
        lines(x = 1:length(toPlot), y = d_to_avgProj[, j], col = j)
    }
    abline(v = burninLine, lty = 2)
    # plot Uk matrix distances
    plot(NA, xlab = "", ylab = "Axis Frob. dist.",
         xlim = c(1, length(toPlot)),
         ylim = c(min(d_to_avgUk), max(d_to_avgUk)),
         main = "dists. to Avg. Uk matrix")
    for(j in 1:numchains) {
        points(1, d_to_avgUk[1, j], col = j, cex = 2)
        lines(x = 1:length(toPlot), y = d_to_avgUk[, j], col = j)
    }
    abline(v = burninLine, lty = 2)
    
    par(mfrow = c(1, 1), cex.axis = 1, cex.lab = 1)
    
    return(list(d_to_avgCov = d_to_avgCov, d_to_avgCovEvec = d_to_avgCovEvec,
                d_to_avgProj = d_to_avgProj, d_to_avgUk = d_to_avgUk,
                chain_quantiles = chain_quantiles))
}

# make MCMC lists for coda
make_mcmclist <- function(dist_vals, S_its, burnin, thinBy) {
    if (burnin < 1) {
        startAt <- floor(S_its * burnin) + 1
    } else {
        startAt <- burnin + 1
    }
    
    numchains <- dim(dist_vals$d_to_avgCov)[2]
    
    toPlot <- seq(1, S_its, by = thinBy)
    itsKeep <- (which(toPlot >= startAt)[1]) : nrow(dist_vals$d_to_avgCov)
    print(length(itsKeep))
    
    out_mcmcs <- list()
    
    for (j in 1:numchains) {
        out_mcmcs[[j]] <- 
            cbind(d_to_avgCov     = dist_vals$d_to_avgCov[itsKeep, j],
                  d_to_avgCovEvec = dist_vals$d_to_avgCovEvec[itsKeep, j],
                  d_to_avgProj    = dist_vals$d_to_avgProj[itsKeep, j],
                  d_to_avgUk      = dist_vals$d_to_avgUk[itsKeep, j]) |> 
            coda::mcmc(start = startAt, thin = thinBy)
    }
    
    return(coda::mcmc.list(out_mcmcs))
}

# set up ------------------------------------------------------------------

# data parameters
P <- 8
d <- 4
n_k <- 1000

tracenorm <- TRUE
customEvals <- FALSE
sigmaevals <- c(rep(1.5, d),
                rep(.5, d))

parameterseed <- 27052025
dataseed <- 10062025

sortLambda <- TRUE

# sampling parameters
numchains <- 1
S_its <- 1e5
burnin <- 33000
thinBy <- 5
tau_U <- .01 # SD of the entries of the Cayley transformation

doCayleyZeros <- FALSE
CayleyZeroProb <- 0.75

# make data and do sampling -----------------------------------------------

datastuff <- genData_Uk_CMACG(P, d, n_k, tracenorm, customEvals, sigmaevals, 
                              sortLambda, parameterseed, dataseed)

Ukest <- eigen(datastuff$data_list[[1]])$vectors[, 1:d]

# attempts to generate "overdispersed" starting points
initialUks <- list(
    NULL, 
    Uk_MH_Cayley(Ukest, .1, FALSE, NULL),
    Uk_MH_Cayley(Ukest, .2, FALSE, NULL),
    Ukest, MASS::Null(Ukest)[, 1:d]
)

cluster <- parallel::makeCluster(numchains)
doParallel::registerDoParallel(cluster)

manysamples <- list()
{
    print(Sys.time())
    manysamples <- foreach::foreach(i=1:numchains) %do% {
        manysamples[[i]] <- 
            sample_Uk_MH_Cayley(datastuff$data_list, datastuff$param_list,
                                S_its = S_its, tau_U = tau_U, 
                                initialUk = initialUks[[i]],
                                doCayleyZeros = doCayleyZeros,
                                CayleyZeroProb = CayleyZeroProb,
                                doStatus = FALSE, statusNum = 10)
    }
    print(Sys.time())
}

parallel::stopCluster(cluster)

summarizeStuff <- summarize_Uk(manysamples, datastuff$param_list, burnin)

summarizeStuff$avgs_df

dev.new()
pdf("Uk_CMACG_sampler_trace.pdf", width = 9, height = 7.5)
dist_vals <- avg_tracePlots(manysamples, summarizeStuff, datastuff$param_list, 
                            burnin = burnin, tracePlotEvery = thinBy)
dev.off()

dist_vals$chain_quantiles
with(dist_vals,
     chain_quantiles[rownames(chain_quantiles) == "d_to_avgUk", ])

# max distances -----------------------------------------------------------

Lambda0 <- datastuff$param_list$Lambda_ks
# maximum distance of the covariance matrices
# actually, not quite this: 2 * sqrt( sum(Lambda0^2) )
# maximum distance of projection and eigenvector matrices
sqrt(2*d)

# mcmc lists --------------------------------------------------------------

bigoutmcmcs <- make_mcmclist(dist_vals, S_its, burnin, thinBy)

coda::gelman.diag(bigoutmcmcs)

coda::effectiveSize(bigoutmcmcs)
lapply(bigoutmcmcs, coda::effectiveSize)

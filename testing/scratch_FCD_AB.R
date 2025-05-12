# FCD sampler function for A, B matrices

source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/generatedata.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/rcgb.R")
# source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/testing/scratch_data-covar-dense.R")
source("~/Documents/PhD_research/RA_time-series/code-experiments/complexeigenmodel/functions/utility.R")

AB_gibbs_densCovar <- function(data_list, param_list, gs = 201, al_be_lower = 1,
                               al_be_upper = 2,
							   eta0 = 2, tau02 = 1/1000) {
	P <- param_list$P
	d <- param_list$d
	K <- param_list$K
	Vs <- param_list$Vs
	U_ks <- param_list$U_ks
	As <- param_list$As
	Bs <- param_list$Bs
	
	# all i < j
	firstpairs <- combn(1:P, 2) |> t()
	# such that i < d + 1, equivalent to i <= d
	# TODO the choice of pairs should depend on the smallest value in A, B matrices
	# - if alpha, beta range from 2 to 1: use i <= d
	# - if alpha, beta range from 1 to 0: use i <= d - 1
	Ppairs <- firstpairs[firstpairs[, 1] <= d, ]
	
	# calculate matrix M, sum of function of V and U_k matrices
	M <- matrix(0, nrow = P, ncol = d)
	for (k in 1:K) {
	    M <- M + Conj( t(Conj(Vs))%*%U_ks[, , k] ) * t(Conj(Vs))%*%U_ks[, , k]
	}
	# M should be exactly real, but it comes out complex with 0 imaginary values
	M <- Re(M)
	
	ws <- (param_list$As[1] / al_be_upper)^2
	
	alphas <- As / sqrt(ws)
	betas <- Bs / sqrt(ws)
	
	##########
	# sample alpha
	##########
	# this vector stays constant, when sampling alpha=
	Mbeta <- M %*% betas
	
	# TODO sample a_ind randomly, so we do it in a different order each time
	for (a_ind in sample(2:(P-1))) {

		# for the a index: select all the rows where the 1st or 2nd column is a_ind
		# a_ind <- 3
		# these are the combined columns
		# Ppairs[Ppairs[, 1] == a_ind | Ppairs[, 2] == a_ind, ]
		# these have the indices I actually need
		others <- Ppairs[Ppairs[, 1] == a_ind, 2]
		some <- Ppairs[Ppairs[, 2] == a_ind, 1]

		# for each value in grid xs = seq(0, 1, .1) (or some other sequence),
		# calculate the product of (a[some] - x) and (x - a[others]). the
		# sampled alpha value should be between the surrounding alpha values.
		# since the alphas are sorted decreasing, that means it should be
		# between the alpha values with a_ind+1 and a_ind-1, minimum and maximum
		# respectively. However, since the prior is that the alpha values are
		# order statistics from a Uniform(1, 2) distribution, and thus they are
		# equal with probability 0, and further, since the approximation we are
		# using requires that the alpha values are not equal, we should not give
		# positive probability to the values possibly being equal in the
		# posterior. Therefore, we should not use a sequence starting and ending
		# exactly with the respective min and max alpha values. Instead, we will
		# begin and end the sequences on a small offset from the exact
		# surrounding alpha values.

		# xs <- seq(0, 1, .1)
		# create the min and max values for the sequence. Give them an offset, so that we do not possibly sample equal alpha values.
		min_al <- alphas[a_ind + 1]
		max_al <- alphas[a_ind - 1]
		gint <- (max_al - min_al)/(gs-1)
		min_al <- min_al + gint/2
		max_al <- max_al - gint/2

		xs <- seq(min_al, max_al, length.out = gs)

		# TODO check that this sum of logs is actually calculating the right log density
		aprod1 <- sapply(xs, function(x) {
		    sum(log((alphas[some] - x))) + sum(log(x - alphas[others]))
		    })

		# aprod2
		aprod2 <- ws * xs * Mbeta[a_ind]
		# calculate the proportional density, with a possible additional factor
		proplogdens <- K * aprod1 + aprod2
		if (2 <= a_ind & a_ind <= d) {
		    # optional: aprod3 (only needed if 2 <= a_ind <= d)
		    aprod3 <- -K * ws * alphas[a_ind] * betas[a_ind]
		    proplogdens <- proplogdens + aprod3
		}

		# draw one value from xs (i.e. an alpha value), using the logweights
		alphas[a_ind] <- sample_gumbel(xs, 1, proplogdens)
	} # end of alphas sample for loop
	
	
	##########
	# sample beta
	##########
	
	alphaM <- t(alphas) %*% M

	for (b_ind in sample(2:(d-1))) {
        
        # select all the rows where the 1st or 2nd column is b_ind
        # indexes that are after b_ind (i.e. 2nd)
        b2s <- Ppairs[Ppairs[, 1] == b_ind & Ppairs[, 2] <= d, 2]
        # indexes that are before b_ind (i.e. 1st)
        b1s <- Ppairs[Ppairs[, 2] == b_ind, 1]
        
        # create a grid of beta values, depending on the next highest and lowest
        # betas values. first, determine the upper and lower limits of the
        # grid.
        min_bl <- betas[b_ind + 1]
        max_bl <- betas[b_ind - 1]
        gint <- (max_bl - min_bl)/(gs-1)
        min_bl <- min_bl + gint/2
        max_bl <- max_bl - gint/2
        # make the grid
        xs <- seq(min_bl, max_bl, length.out = gs)
        
        # calculate terms in the log density
        bprod1 <- sapply(xs, function(x) {
            sum(log((betas[b1s] - x))) + sum(log(x - betas[b2s]))
            })
        bprod2 <- (P-d)*log(betas[b_ind])
        bprod3 <- xs * (-K * ws * alphas[b_ind] + alphaM[b_ind])
        
        # calculate the log density
        proplogdens_b <- K*bprod1 + K*bprod2 + bprod3
        
        # draw one value from xs (i.e. beta value), using the logweights
        betas[b_ind] <- sample_gumbel(xs, 1, proplogdens_b)
    }
	
	##########
	# sample w
	##########
	
	fstar <- choose(P, 2) - choose(P-d, 2)
	betatilde <- c(betas, rep(0, P-d))

	w_shape <- K*fstar + eta0/2
	w_rate <- K * (t(alphas) %*% betatilde) - t(alphas) %*% M %*% betas + tau02

	ws1 <- rgamma(1, shape = w_shape, rate = w_rate)
	
	return(list(As = sqrt(ws1) * alphas, 
	            Bs = sqrt(ws1) * betas))
}

# test the sampler function ----------------------------------------------

# P <- param_list$P
# d <- param_list$d
# K <- param_list$K
# Vs <- param_list$Vs
# U_ks <- param_list$U_ks
# As <- param_list$As
# Bs <- param_list$Bs

P <- 8
d <- 4
K <- 3
Vs <- runitary(P, P)

# TODO properly generate Uk - this is just a placeholder to make sure the function runs
U_ks <- array(NA, c(P, d, K))
for (k in 1:K) {
	U_ks[, , k] <- runitary(P, d)
}

As <- seq(2, 1, length.out = P)
Bs <- seq(2, 1, length.out = d)

# TODO placeholder data - the AB FCD does not even require this data
data_list <- list()

param_list <- list(
	P = P,
	d = d,
	K = K,
	Vs = Vs,
	U_ks = U_ks,
	As = As,
	Bs = Bs
)

# test running the function
AB_gibbs_densCovar(data_list, param_list)

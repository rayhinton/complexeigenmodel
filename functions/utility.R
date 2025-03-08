# Utility functions


# sample_gumbel -----------------------------------------------------------

# Description: draw a sample from a discrete random variable using log-scale weights that are possibly unnormalized

# Input: 

# thetas - values of the discrete random variables
# size - number of sample draws to return
# logweights - logarithm of the discrete probabilities (possibly unnormalized)

# Output:

# outsample - values from thetas sampled according to logweights

# Credit:

# Class notes (STAT 633). Dr. Anirban Bhattacharya. Fall 2024.

sample_gumbel <- function(thetas, size, logweights) { 
    stopifnot("length of thetas and logweights must be equal" = 
                  length(thetas) == length(logweights))
    
    # number of discrete values
    J <- length(thetas)
    # initialize a vector that will contain the sample
    outsample <- vector(mode(thetas), size)
    
    for (i in 1:size) {
        # generate Gumbel noise
        Gj <- -log(rexp(J))
        
        # calculate Gumbel weights = logweights + Gumbel noise
        gumweights <- logweights + Gj
        
        # perform the sampling step: return the theta corresponding to the
        # largest Gumbel weight
        outsample[i] <- thetas[which.max(gumweights)]        
    }
    
    return(outsample)
}
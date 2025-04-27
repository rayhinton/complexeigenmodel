# matrix distances and statistics

frob_norm <- function(X) {
    outdiag <- Re(diag(crossprod(Conj(X), X)))
    outlist <- list(fnorm = sum(outdiag), diags = outdiag)
    
    return(outlist)
}

frame_distance <- function(A, B) {
    # return(norm(A - B, "F"))
    return(frob_norm(A - B))
}

procrustes_distance <- function(A, B) {
    # Get orthogonal bases
    Q_A <- qr.Q(qr(A))
    Q_B <- qr.Q(qr(B))
    
    # Find optimal rotation between the frames
    svd_result <- svd(t(Conj(Q_B)) %*% Q_A)
    R_opt <- svd_result$v %*% t(Conj(svd_result$u))
    
    # Compute distance after optimal alignment
    # return(norm(Q_A - Q_B %*% R_opt, "F"))
    return(frob_norm(Q_A - Q_B %*% R_opt))
}

# Description: calculate Grassmann Distance between subspaces

# Inputs:
# A, p x k matrix, semi-unitary
# B, p x k matrix, semi-unitary
# r, the number of columns of A and B to use for calculating subspaces
# s_tol, tolerance to use for checking if singular values are greater than 1

# Output:
# sub_dist, the Grassmann distance between the spaces spanned by A[, 1:d] and B[, 1:d]

grass_dist <- function (A, B, r = ncol(A), s_tol = 2*.Machine$double.eps) {
    stopifnot("A and B must have the same dimensions" = all(dim(A) == dim(B)))
    
    Sigma <- t(Conj(A)) %*% B
    svd_Sigma <- svd(Sigma)
    
    # check for numerical issues with SVDs
    max_diff <- max(svd_Sigma$d[1:r] - 1)
    if (any(svd_Sigma$d[1:r] > (1 + s_tol))) {
        # warning()
        stop(paste0("Singular values exceed 1 by more than s_tol = ", s_tol, 
                    ". Maximum excess: ", max_diff))
    } else {
        sigmas <- pmin(svd_Sigma$d[1:r], 1)
    }
    
    thetas <- acos(sigmas)
    sq_thetas <- thetas^2
    sub_dist <- sqrt(sum(sq_thetas))
    
    outlist <- list(sq_thetas = sq_thetas, sub_dist = sub_dist)
    
    return(outlist)
}

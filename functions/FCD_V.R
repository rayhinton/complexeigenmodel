V_gibbs_densCovar <- function(data_list, param_list, 
                              Imtol = 2*.Machine$double.eps) {
    
    Vpar <- matrix(0 + 0i, param_list$P, param_list$P)
    for (k in 1:param_list$K) {
        Vpar <- Vpar + param_list$U_ks[, , k] %*% diag(param_list$Bs) %*% 
            t(Conj(param_list$U_ks[, , k]))
    }
    
    outV <- rcBingUP_gibbs(param_list$Vs, Vpar, diag(param_list$As),
                           Imtol = Imtol)
    
    return(outV)
}

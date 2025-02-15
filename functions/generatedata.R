# generate sample data

# generate a random unitary matrix
# V matrix, PxP
runitary <- function(nrows, ncols) {
    outmat <- complex(real = rnorm(nrows*ncols), 
                      imaginary = rnorm(nrows*ncols)) |> 
        matrix(nrow = nrows) |> 
        qr() |> 
        qr.Q()
    
    return(outmat)
}

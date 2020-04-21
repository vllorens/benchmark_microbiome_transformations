#' Impute zeros
#' 
#' Function to impute zeros prior to CLR transformation
#' 
#' @param matrix.test matrix with rows as taxa and columns as samples
#' 
#' @return matrix with zeros imputed

estimate0.min <- function(matrix.test) { 
    print(paste("You have",ncol(matrix.test),"samples. If this is not correct, transpose matrix!"))
    matrix.test.p <- sweep(matrix.test,MARGIN=2,FUN="/",STATS=colSums(matrix.test)) # divide taxa in each sample by the total counts
    matrix.test.p[is.nan(matrix.test.p)] <- 0
    samplesums <- colSums(matrix.test)
    
    matrix.f.n0 <- matrix.test
    for (i in 1:nrow(matrix.f.n0)) {
        min <- min(matrix.test.p[i,][matrix.test.p[i,] > 0])
        for (j in 1:ncol(matrix.f.n0 )) {
            if (matrix.f.n0 [i,j] == 0)
                matrix.f.n0 [i,j] <- min*samplesums[j]
        }
    }
    return(matrix.f.n0)
}

#' Calculates taxon-taxon correlations 
#' 
#' Function to calculate Spearman correlation coefficients and p-values from the provided (transformed) matrix
#' 
#' @param taxon_matrix matrix with rows as taxa and columns as samples, can be already transformed
#' 
#' @return list with two elements: the matrix of proportionality coefficients and the matrix of p-values. Since this matrix is symmetric, only the upper triangle is calculated

## cORRECT NEGATIVE CORRELATIONS


propr_core <- function(ind1, ind2, taxon_matrix){
    cor <- 1-var(taxon_matrix[ind1,]-taxon_matrix[ind2,])/(var(taxon_matrix[ind1,])+var(taxon_matrix[ind2,]))
    permutations   <- 1000
    permutedY <- replicate(permutations, sample(c(unlist(taxon_matrix[ind2,]))))
    cor.perm <- apply(permutedY, 2, function(Y){
        return(1-var(taxon_matrix[ind1,]-Y)/(var(taxon_matrix[ind1,])+var(Y)))
    })
    p.value <- length(abs(cor.perm)[abs(cor.perm) >= abs(cor)])/permutations
    return(as.list(set_names(c(cor,p.value), c("cor", "pvalue"))))
}

propr_custom <- function(taxon_matrix){
    if(is.data.frame(taxon_matrix)){
        taxon_matrix <- as.matrix(taxon_matrix)
    }
    df_rho = matrix(nrow = nrow(taxon_matrix), ncol = nrow(taxon_matrix), NA)
    df_pvalue = matrix(nrow = nrow(taxon_matrix), ncol = nrow(taxon_matrix), NA)
    rownames(df_rho) <- rownames(taxon_matrix)
    colnames(df_rho) <- rownames(taxon_matrix)
    rownames(df_pvalue) <- rownames(taxon_matrix)
    colnames(df_pvalue) <- rownames(taxon_matrix)
    
    ttnames <- t(combn(rownames(taxon_matrix),2, simplify = T)) %>% 
        as.data.frame(., stringsAsFactors=F)
    colnames(ttnames) <- c("ind1", "ind2")
    
    ttnames <- ttnames %>% 
        bind_cols(future_pmap_dfr(., propr_core, taxon_matrix, .progress = TRUE))
    
    df_rho[lower.tri(df_rho, diag=F)] <- ttnames$cor
    df_rho <- t(df_rho)
    df_pvalue[lower.tri(df_pvalue, diag=F)] <- ttnames$pvalue
    df_pvalue <- t(df_pvalue)
    
    return(list(df_rho, df_pvalue))
}
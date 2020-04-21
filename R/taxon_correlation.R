#' Calculates taxon-taxon correlations 
#' 
#' Function to calculate Spearman correlation coefficients and p-values from the provided (transformed) matrix
#' 
#' @param taxon_matrix matrix with rows as taxa and columns as samples, can be already transformed
#' 
#' @return list with two elements: the matrix of correlation coefficients and the matrix of p-values. Since this matrix is symmetric, only the upper triangle is calculated

taxon_correlation <- function(taxon_matrix){
    df_rho = matrix(nrow = nrow(taxon_matrix), ncol = nrow(taxon_matrix), 0)
    df_rho[lower.tri(df_rho, diag = TRUE)] <- NA
    df_pvalue = matrix(nrow = nrow(taxon_matrix), ncol = nrow(taxon_matrix), 0)
    df_pvalue[lower.tri(df_pvalue, diag = TRUE)] <- NA
    for (l in 1:nrow(taxon_matrix)){
        for (k in 1:nrow(taxon_matrix)){
            if (is.na(df_rho[l,k]) == FALSE){
                cor = cor.test(taxon_matrix[l,], taxon_matrix[k,], alternative = "two.sided", method = "spearman", exact = F)
                df_rho[l,k] = cor$estimate
                df_pvalue[l,k] = cor$p.value
            }
        }
    }
    rownames(df_rho) <- rownames(taxon_matrix)
    colnames(df_rho) <- rownames(taxon_matrix)
    rownames(df_pvalue) <- rownames(taxon_matrix)
    colnames(df_pvalue) <- rownames(taxon_matrix)
    return(list(df_rho, df_pvalue))
}

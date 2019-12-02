#' Calculates taxon-metadata correlations 
#' 
#' Function to calculate Spearman correlation coefficients and p-values from the provided (transformed) taxonomy and metadata matrices
#' 
#' @param taxon_matrix matrix with rows as taxa and columns as samples, can be already transformed
#' @param metadata_matrix matrix with metadata, rows are samples and columns are metadata categories
#' 
#' @return list with two elements: the matrix of correlation coefficients and the matrix of p-values. Since this matrix is symmetric, only the upper triangle is calculated

taxon_metadata_correlation <- function(taxon_matrix, metadata_matrix){
    df_rhometa = matrix(nrow = nrow(taxon_matrix), ncol = ncol(metadata_matrix), 0)
    df_pvaluemeta = matrix(nrow = nrow(taxon_matrix), ncol = ncol(metadata_matrix), 0)
    for(l in 1:ncol(metadata_matrix)){
        if(is.numeric(metadata_matrix[,l])){
            for(k in 1:nrow(taxon_matrix)){
                cor <- cor.test(taxon_matrix[k,], metadata_matrix[,l], alternative="t", method="sp", exact = F)
                df_rhometa[k,l] <- cor$estimate
                df_pvaluemeta[k,l] <- cor$p.value
            }
        } else if(is.factor(metadata_matrix[,l]) | is.character(metadata_matrix[,l])){
            for(k in 1:nrow(taxon_matrix)){
                cor <- kruskal.test(taxon_matrix[k,] ~ as.factor(metadata_matrix[,l]))
                df_rhometa[k,l] <- NA
                df_pvaluemeta[k,l] <- cor$p.value
            }
        }
    }
    rownames(df_rhometa) <- rownames(taxon_matrix)
    colnames(df_rhometa) <- colnames(metadata_matrix)
    rownames(df_pvaluemeta) <- rownames(taxon_matrix)
    colnames(df_pvaluemeta) <- colnames(metadata_matrix)
    return(list(df_rhometa, df_pvaluemeta))
}

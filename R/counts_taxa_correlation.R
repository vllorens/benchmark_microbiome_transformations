#' Counts-taxa correlation
#' Function to calculate Spearman correlation coefficients and p-values from the provided (transformed) taxonomy and metadata matrices
#' 
#' @param taxon_matrix matrix with rows as taxa and columns as samples, can be already transformed
#' @param total_counts vector with total counts per sample
#' @return list with two elements: the matrix of correlation coefficients and the matrix of p-values. Since this matrix is symmetric, only the upper triangle is calculated

counts_taxa_correlation <- function(taxon_matrix, total_counts){
    tot_rhometa <- c()
    tot_pvaluemeta <- c()
    for(l in 1:nrow(taxon_matrix)){
        cor <- cor.test(total_counts, taxon_matrix[l,], alternative="t", method="sp", exact = F)
        tot_rhometa[l] <- cor$estimate
        tot_pvaluemeta[l] <- cor$p.value
    }
    names(tot_pvaluemeta) <- rownames(taxon_matrix)
    names(tot_rhometa) <- rownames(taxon_matrix)
    return(list(tot_rhometa, tot_pvaluemeta))
}

# rarefaction to even sampling depth #
# author: Doris Vandeputte           #
######################################
# this script doesn't include copy number correction, a function for copy number correction is included in RDP classifier 2.12 
# this script uses function rarefy_even_depth from phyloseq 1.20.0, it needs package phyloseq to be installed and loaded in order to work.
# with cnv_corrected_abundance_table: taxa are rows, samples are columns
# with cell_counts_table: a table with sample-identifiers as rows, cell counts as columns 
library(phyloseq)
rarefy_even_sampling_depth <- function(cnv_corrected_abundance_table, cell_counts_table) {
    cnv_corrected_abundance_table <- t(cnv_corrected_abundance_table)
    try(if(all(row.names(cnv_corrected_abundance_table) == row.names(cell_counts_table))==FALSE) stop("Cnv_corrected_abundance_table and cell_counts_table do not have the same sample-names, Please check!"))
    cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations
    cell_counts_table = t(cell_counts_table[row.names(cnv_corrected_abundance_table),]) # make sure the order of the samples is the same in both files  
    sample_sizes = rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
    sampling_depths = sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
    minimum_sampling_depth = min(sampling_depths) # minimum of all sampling depths
    rarefy_to = (cell_counts_table * minimum_sampling_depth) - 1 # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
    cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
    rarefied_matrix=matrix(nrow = nrow(cnv_corrected_abundance_table_phyloseq), ncol = ncol(cnv_corrected_abundance_table_phyloseq), dimnames = list(rownames(cnv_corrected_abundance_table_phyloseq), colnames(cnv_corrected_abundance_table_phyloseq)))
    for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq))
    {
        x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i], rngseed = 711, replace = FALSE, trimOTUs = F, verbose = FALSE)
        rarefied_matrix[i,] = x
    }
    normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
    QMP = t(normalised_rarefied_matrix*cell_counts_table[1,])
    return(QMP)
}
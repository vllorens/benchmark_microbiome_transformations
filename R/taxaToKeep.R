#' Select taxa to keep
#' 
#' Function to impute zeros prior to CLR transformation
#' 
#' @param filename matrix with rows as taxa and columns as samples
#' 
#' @return indices of the taxa matrix that should be kept

taxaToKeep <- function(filename){
    rmp_name <- gsub(filename, pattern="[[:print:]]*taxonomy", replacement="seqOut_taxonomy")
    rmp_matrix <- read.table(paste0("data/seq_matrices/", rmp_name), header=T, stringsAsFactors = F, sep="\t")
    taxa_to_keep = as.vector(which(apply(rmp_matrix == 0 ,1, sum) <= 0.50*ncol(rmp_matrix)))
    return(taxa_to_keep)
}

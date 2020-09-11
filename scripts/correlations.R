# Taxon-taxon and taxon-metadata correlations
# Mon Dec  2 15:01:23 2019 ------------------------------

# We will use different data and transformation methods to calculate all the correlations. 
# For comparisons, we will only use Spearman correlations for taxon-taxon and either Spearman correlation or Kruskal-wallis test for taxon-metadata associations 
# (Spearman for numeric and Kruskal-Wallis for factors in the metadata)

# The following methods will be tested:
# * REAL: Original data (only used as reference). Not transformed.

# * Seq: Sequencing data. Not transformed.
# * Rel: Relative abundances based on the sequencing data. 
# * AST: Relative abundances based on the sequencing data + arcsine squared-root transformation. 
# * RMP: Relative microbiome profiling (rarefying data)
# * CLR: Centered Log Ratio transformation (CoDaSeq implementation)
# * CSS: Data nrmalized using cumulative sum scaling (MetagenomeSeq implementation).
# * GMPR: Size factors for normalization are calculated using the Geometric Mean of Pairwise Ratios (GMPR package).
# * TMM: Normalization using trimmed mean of M-values (EdgeR implementation).
# * RLE: Normalization using relative log-expression (EdgeR implementation).
# * UQ: Normalization using upper quartile (EdgeR implementation).
# * VST: Using variance stabilizing transformation from DESeq2.
# * QMP: QMP data. Normalized by rarefying to even sampling depth and scaling to total counts. 
# * ACS: Absolute count scaling. Normalized by scaling to total counts. 
# 
# In all of these methods, we will only evaluate those taxa with prevalence > 50%. We will use the sequencing matrix to determine this threshold on the zeros. 
# 
# The correlations that we will calculate, for each type of transformation, are:
# * Taxon-counts correlations (i.e. with total microbial loads)
# * Taxon-taxon correlations
# * Taxon-metadata correlations


#### Configure the environment ####

# load required packages
suppressPackageStartupMessages(library(CoDaSeq))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(metagenomeSeq))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(rstatix))

# load functions
source("R/estimateZeros.R")
source("R/taxaToKeep.R")
source("R/taxon_correlation.R")
source("R/taxon_metadata_correlation.R")
source("R/counts_taxa_correlation.R")
source("R/GMPR.R")
source("R/plot_comparisons.R")

# If it has not been run yet, it is necessary to run the script to generate the data (scripts/generate_data.R)
# source("scripts/generate_data.R")

# create folders, if not existing
system("mkdir -p data/correlations_taxontaxon")
system("mkdir -p data/correlations_taxonmetadata")
system("mkdir -p data/correlations_taxoncounts")
system("mkdir -p data/all_matrices")
system("mkdir -p output/taxontaxon")
system("mkdir -p output/taxonmetadata")
system("mkdir -p output/taxoncounts")
system("mkdir -p output/raw_transformed_correlation")

# clean these folders and remove all previous results (if any) 
system("rm -r data/correlations_taxontaxon/*")
system("rm -r data/all_matrices/*")
system("rm -r data/correlations_taxonmetadata/*")
system("rm -r data/correlations_taxoncounts/*")
system("rm -r output/taxontaxon/*")
system("rm -r output/taxonmetadata/*")
system("rm -r output/taxoncounts/*")
system("rm -r output/raw_transformed_correlation/*")

# create additional folders for the reference data (on the real matrices) and the special taxa
system("mkdir -p data/correlations_taxontaxon/reference")
system("mkdir -p data/all_matrices/reference")
system("mkdir -p data/correlations_taxonmetadata/reference")
system("mkdir -p data/correlations_taxoncounts/reference")
system("mkdir -p output/taxonmetadata/specialtaxon")
system("mkdir -p output/taxontaxon/specialtaxon")
system("mkdir -p output/taxoncounts/specialtaxon")
set.seed(777)

#### Calculate correlations for all methods ####

# REAL
for(file in list.files("data/tax_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # read metadata file
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename))
    # get counts and select only those taxa to keep
    counts_data <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_rho_REAL")
    outputname_pval <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_pvalue_REAL")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/reference/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/reference/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="taxonomy", replacement="taxonmetadata_rho_REAL")
    outputname_pval <- gsub(filename, pattern="taxonomy", replacement="taxonmetadata_pvalue_REAL")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/reference/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/reference/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="taxonomy", replacement="taxoncounts_rho_REAL")
    outputname_pval <- gsub(filename, pattern="taxonomy", replacement="taxoncounts_pvalue_REAL")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/reference/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/reference/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/reference/REAL_", filename))
}
set.seed(777)

# Seq
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix)
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_Seq")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_Seq")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_Seq")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_Seq")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_Seq")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_Seq")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="Seq")))
}
set.seed(777)

# Rel
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    total_counts <- apply(tax_matrix,2,sum)
    tax_matrix <- sweep(tax_matrix, MARGIN = 2, total_counts, '/')
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix)
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_Rel")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_Rel")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_Rel")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_Rel")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_Rel")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_Rel")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="Rel")))
}
set.seed(777)

# CLR
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    tax_matrix <- t(zCompositions::cmultRepl(X = t(tax_matrix), output="p-counts")) # estimates zeros
    tax_matrix <- t(codaSeq.clr(tax_matrix, samples.by.row=F))
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix)
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_CLR")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_CLR")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_CLR")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_CLR")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_CLR")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_CLR")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="CLR")))
    
}
set.seed(777)


# GMPR
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    gmpr.size.factor <- GMPR(comm=tax_matrix)
    tax_matrix <- sweep(tax_matrix, MARGIN = 2, gmpr.size.factor, '/')
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_GMPR")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_GMPR")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_GMPR")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_GMPR")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_GMPR")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_GMPR")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="GMPR")))

}

set.seed(777)

# CSS-metagenomeSeq
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    tax_matrix <- newMRexperiment(tax_matrix)
    p <- cumNormStatFast(tax_matrix)
    tax_matrix <- cumNorm(tax_matrix, p = p)
    tax_matrix <- MRcounts(tax_matrix, norm=T)
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix)
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix)
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data)
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_CSS")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_CSS")
    write.table(taxon_correlation_object[[1]],
                paste0("data/correlations_taxontaxon/", outputname_correlation),
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]],
                paste0("data/correlations_taxontaxon/", outputname_pval),
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_CSS")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_CSS")
    write.table(metadata_correlation_object[[1]],
                paste0("data/correlations_taxonmetadata/", outputname_correlation),
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]],
                paste0("data/correlations_taxonmetadata/", outputname_pval),
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_CSS")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_CSS")
    write.table(counts_correlation_object[[1]],
                paste0("data/correlations_taxoncounts/", outputname_correlation),
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]],
                paste0("data/correlations_taxoncounts/", outputname_pval),
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="CSS")))
}
set.seed(777)


# TMM-edgeR
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    tax_matrix <- DGEList(counts=tax_matrix)
    tax_matrix <- calcNormFactors(tax_matrix, method = "TMM")
    tax_matrix <- cpm(tax_matrix)
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_TMM")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_TMM")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_TMM")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_TMM")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_TMM")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_TMM")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="TMM")))
}

set.seed(777)

# UQ-edgeR
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    tax_matrix <- DGEList(counts=tax_matrix)
    tax_matrix <- calcNormFactors(tax_matrix, method = "upperquartile")
    tax_matrix <- cpm(tax_matrix)
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_UQ")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_UQ")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_UQ")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_UQ")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_UQ")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_UQ")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="UQ")))
}
set.seed(777)


# RLE-edgeR
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    tax_matrix <- DGEList(counts=tax_matrix)
    tax_matrix <- calcNormFactors(tax_matrix, method = "RLE")
    tax_matrix <- cpm(tax_matrix)
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_RLE")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_RLE")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_RLE")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_RLE")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_RLE")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_RLE")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="RLE")))
    
}
set.seed(777)

# TSS-AST
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    tax_matrix = apply(tax_matrix, 2, function(x){x/sum(x)}) # TSS
    tax_matrix <- asin(sqrt(tax_matrix)) # AST
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_AST")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_AST")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_AST")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_AST")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_AST")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_AST")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="AST")))
    
}
set.seed(777)


# VST-DESeq2
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    tax_matrix <- varianceStabilizingTransformation(tax_matrix)
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_rho_VST")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxontaxon_pvalue_VST")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_rho_VST")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxonmetadata_pvalue_VST")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_rho_VST")
    outputname_pval <- gsub(filename, pattern="seqOut_taxonomy", replacement="taxoncounts_pvalue_VST")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", gsub(filename, pattern="seqOut", replacement="VST")))
    
}
set.seed(777)


# RMP
for(file in list.files("data/rmp_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # read metadata file
    filename_metadata <- gsub(filename, pattern="RMP_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="RMP_taxonomy", replacement="taxontaxon_rho_RMP")
    outputname_pval <- gsub(filename, pattern="RMP_taxonomy", replacement="taxontaxon_pvalue_RMP")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="RMP_taxonomy", replacement="taxonmetadata_rho_RMP")
    outputname_pval <- gsub(filename, pattern="RMP_taxonomy", replacement="taxonmetadata_pvalue_RMP")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="RMP_taxonomy", replacement="taxoncounts_rho_RMP")
    outputname_pval <- gsub(filename, pattern="RMP_taxonomy", replacement="taxoncounts_pvalue_RMP")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", filename))
    
}
set.seed(777)


# QMP
for(file in list.files("data/qmp_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # read metadata file
    filename_metadata <- gsub(filename, pattern="QMP_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="QMP_taxonomy", replacement="taxontaxon_rho_QMP")
    outputname_pval <- gsub(filename, pattern="QMP_taxonomy", replacement="taxontaxon_pvalue_QMP")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="QMP_taxonomy", replacement="taxonmetadata_rho_QMP")
    outputname_pval <- gsub(filename, pattern="QMP_taxonomy", replacement="taxonmetadata_pvalue_QMP")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="QMP_taxonomy", replacement="taxoncounts_rho_QMP")
    outputname_pval <- gsub(filename, pattern="QMP_taxonomy", replacement="taxoncounts_pvalue_QMP")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", filename))
    
}
set.seed(777)

# ACS
for(file in list.files("data/acs_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # read metadata file
    filename_metadata <- gsub(filename, pattern="ACS_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # get counts and select only those taxa to keep
    origfile <- paste0("data/tax_matrices/", filename_metadata)
    counts_data <- read.table(origfile, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(., 2, sum)
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # calculate correlations and pvalues (taxon-counts)
    counts_correlation_object <- counts_taxa_correlation(tax_matrix, counts_data) 
    counts_correlation_object[[2]] <- p.adjust(counts_correlation_object[[2]], method="BH")
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="ACS_taxonomy", replacement="taxontaxon_rho_ACS")
    outputname_pval <- gsub(filename, pattern="ACS_taxonomy", replacement="taxontaxon_pvalue_ACS")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="ACS_taxonomy", replacement="taxonmetadata_rho_ACS")
    outputname_pval <- gsub(filename, pattern="ACS_taxonomy", replacement="taxonmetadata_pvalue_ACS")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-counts)
    outputname_correlation <- gsub(filename, pattern="ACS_taxonomy", replacement="taxoncounts_rho_ACS")
    outputname_pval <- gsub(filename, pattern="ACS_taxonomy", replacement="taxoncounts_pvalue_ACS")
    write.table(counts_correlation_object[[1]], 
                paste0("data/correlations_taxoncounts/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(counts_correlation_object[[2]], 
                paste0("data/correlations_taxoncounts/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write matrix
    write.table(tax_matrix, paste0("data/all_matrices/", filename))
}





#### Evaluate taxon-counts correlations ####

# define parameters for significance and initialize vectors for the evaluation data
significance_level <- 0.05

true_positive <- c()
false_positive <- c()
true_negative <- c()
false_negative <- c()
discordant <- c()
method <- c()
spread <- c()
matrixnum <- c()
datatable <- c()
scen <- c()
discordant_tb <- tibble()

for(file in list.files("data/correlations_taxoncounts/", recursive = F, pattern = "pvalue", full.names = T)){
    # define file name, method used and reference to be compared against
    filename <- basename(file)
    methodname <- strsplit(filename, split = "_")[[1]][3]
    matrixname <- strsplit(filename, split = "_")[[1]][4]
    spreadname <- strsplit(filename, split = "_")[[1]][5]
    scenarioname <- strsplit(filename, split="_")[[1]][7] %>% gsub(., pattern="\\.tsv", replacement="")
    referencename <- gsub(filename, pattern=methodname, replacement="REAL")
    rhoname <- gsub(filename, pattern="pvalue", replacement="rho")
    referencerhoname <- gsub(referencename, pattern="pvalue", replacement="rho")
    # read files
    test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
    reference_file <- read.table(paste0("data/correlations_taxoncounts/reference/", referencename),
                                 header=T, stringsAsFactors=F, sep="\t")
    rho_file <- read.table(paste0("data/correlations_taxoncounts/", rhoname),
                           header=T, stringsAsFactors = F, sep="\t")
    referencerho_file <- read.table(paste0("data/correlations_taxoncounts/reference/", referencerhoname),
                                    header=T, stringsAsFactors = F, sep = "\t")
    
    # all
    test_significant <- c(unlist(test_file < significance_level))
    reference_significant <- c(unlist(reference_file < significance_level))
    test_sign <- sign(c(unlist(rho_file)))
    reference_sign <- sign(c(unlist(referencerho_file)))
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "all")
    true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
    false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
    true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
    false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
    discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
    discordant_names <- tibble(taxa=rownames(test_file)[which(test_significant & reference_significant & test_sign!=reference_sign)], matrix=matrixname, method=methodname)
    discordant_tb <- bind_rows(discordant_tb, discordant_names)
    
    # only positive correlations
    reference_sign_pos <- reference_sign
    reference_sign_pos[reference_sign<0 | reference_significant==F] <- NA
    test_sign_pos <- test_sign
    test_sign_pos[reference_sign<0 | reference_significant==F] <- NA
    reference_significant_pos <- reference_significant
    reference_significant_pos[reference_sign<0 | reference_significant==F] <- NA
    test_significant_pos <- test_significant
    test_significant_pos[reference_sign<0 | reference_significant==F] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "positive")
    true_positive <- c(true_positive, length(which(test_significant_pos & reference_significant_pos & test_sign_pos==reference_sign_pos)))
    false_positive <- c(false_positive, length(which(test_significant_pos & !reference_significant_pos)))
    true_negative <- c(true_negative, length(which(!test_significant_pos & !reference_significant_pos)))
    false_negative <- c(false_negative, length(which(!test_significant_pos & reference_significant_pos)))
    discordant <- c(discordant, length(which(test_significant_pos & reference_significant_pos & test_sign_pos!=reference_sign_pos)))
    
    # only negative correlations
    reference_sign_neg <- reference_sign
    reference_sign_neg[reference_sign>0 | reference_significant==F] <- NA
    test_sign_neg <- test_sign
    test_sign_neg[reference_sign>0 | reference_significant==F] <- NA
    reference_significant_neg <- reference_significant
    reference_significant_neg[reference_sign>0 | reference_significant==F] <- NA
    test_significant_neg <- test_significant
    test_significant_neg[reference_sign>0 | reference_significant==F] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "negative")
    true_positive <- c(true_positive, length(which(test_significant_neg & reference_significant_neg & test_sign_neg==reference_sign_neg)))
    false_positive <- c(false_positive, length(which(test_significant_neg & !reference_significant_neg)))
    true_negative <- c(true_negative, length(which(!test_significant_neg & !reference_significant_neg)))
    false_negative <- c(false_negative, length(which(!test_significant_neg & reference_significant_neg)))
    discordant <- c(discordant, length(which(test_significant_neg & reference_significant_neg & test_sign_neg!=reference_sign_neg)))
    
    # only not correlated taxa
    reference_sign_neg <- reference_sign
    reference_sign_neg[reference_significant==T] <- NA
    test_sign_neg <- test_sign
    test_sign_neg[reference_significant==T] <- NA
    reference_significant_neg <- reference_significant
    reference_significant_neg[reference_significant==T] <- NA
    test_significant_neg <- test_significant
    test_significant_neg[reference_significant==T] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "non_correlated")
    true_positive <- c(true_positive, length(which(test_significant_neg & reference_significant_neg & test_sign_neg==reference_sign_neg)))
    false_positive <- c(false_positive, length(which(test_significant_neg & !reference_significant_neg)))
    true_negative <- c(true_negative, length(which(!test_significant_neg & !reference_significant_neg)))
    false_negative <- c(false_negative, length(which(!test_significant_neg & reference_significant_neg)))
    discordant <- c(discordant, length(which(test_significant_neg & reference_significant_neg & test_sign_neg!=reference_sign_neg)))
    
}

# write results table
# discordant associations (DA, significant association but with opposite sign) have been calculated separately from FP and FN, 
# but in the calculations for FPR, precision and recall, they are considered both as FP and FN
results <- tibble(method, spread, datatable, matrixnum, scen, true_positive, false_positive, true_negative, false_negative, discordant)
results <- results %>%
    mutate(FDR=100*(false_positive+discordant)/(false_positive+discordant+true_positive)) %>% 
    mutate(Recall=100*true_positive/(true_positive+false_negative+discordant)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive+discordant)) %>% 
    mutate(Accuracy=100*(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative+discordant))%>%
    mutate(TP=100*true_positive/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FP=100*false_positive/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(TN=100*true_negative/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FN=100*false_negative/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(DA=100*discordant/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(false_positive_rate=100*(false_positive+discordant)/(true_negative+false_positive+discordant))

write_tsv(discordant_tb, path="output/taxoncounts/discordant_taxa.tsv", col_names=T)

write_tsv(results, "output/taxoncounts/statistics_taxoncounts_correlation.tsv", col_names = T)
method_type <- tibble(method=c("Seq", "RMP", "Rel", "AST", "CLR", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "VST", "QMP", "ACS"),
                      method_type=c("Sequencing - non-transformed", 
                                    rep("Traditional transformations",  times=2),
                                    rep("Compositional transformations", times=8),
                                    rep("Transformations incorporating microbial loads", times=2)))


results <- results %>% 
    left_join(method_type, by="method")

results$method_type <- factor(results$method_type, 
                              levels=c("Sequencing - non-transformed", "Traditional transformations", 
                                       "Compositional transformations",
                                       "Transformations incorporating microbial loads"))
results$method <- factor(results$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

results$spread <- factor(results$spread, levels=c("low", "medium", "high"))
results$scen <- factor(results$scen, levels=c("Healthy", "Blooming","Dysbiosis"))
results_all <- results %>% dplyr::filter(datatable=="all")

results_all$F1 <- 2*((results_all$Precision * results_all$Recall)/(results_all$Precision + results_all$Recall))

# plot results
toplot <- results_all %>% 
    dplyr::select(method, scen, TP, TN, FP, FN, DA, method_type) %>% 
    gather(key="metric", value="value", -method, -scen, -method_type) %>% 
    group_by(method, scen, method_type, metric) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot$metric <- factor(toplot$metric, levels=c("TP", "FP", "TN", "FN", "DA"))

p1 <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by="scen", palette="Dark2",
          legend.title="Metric", title="Associations of transformed taxon abundances with original microbial loads",
          xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p1, filename = "output/taxoncounts/plot_performance_taxoncounts.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)

# now only scaled to the number of positives in the samples
results_all <- results_all %>% 
    mutate(actual_positives=true_positive+false_negative+discordant) 

toplot2 <- results_all %>% 
    mutate(TP_scaled=100*true_positive/actual_positives) %>% 
    mutate(FP_scaled=100*false_positive/actual_positives) %>% 
    mutate(DA_scaled=100*discordant/actual_positives) %>% 
    dplyr::select(method, scen, method_type, TP_scaled, FP_scaled, DA_scaled) %>% 
    gather(key="metric", value="value", -method, -scen, -method_type) %>% 
    group_by(method, scen, method_type, metric) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot2$metric <- factor(toplot2$metric, levels=c("FP_scaled", "DA_scaled", "TP_scaled"))

p1b <- ggbarplot(toplot2, x="method", y="value", fill="metric", facet.by="scen", palette="Dark2", 
                 subtitle="True Positives, False Positives and Discordant Associations scaled to the total actual positives (100%)",
                legend.title="Metric", title="Associations of transformed taxon abundances with original microbial loads",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold")) +
    geom_hline(yintercept = 100, col="black", lwd=1)

ggsave(p1b, filename = "output/taxoncounts/plot_performance_taxoncounts_scaledpositives.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)


# calculate sensitivity, precision and FPR 
sensitivity_kruskal <- results_all %>% 
     group_by(scen) %>% 
     kruskal_test(Recall ~ method)

sensitivity_dunn <- results_all %>%
     group_by(scen) %>%
     dunn_test(Recall ~ method, p.adjust.method="BH")

precision_kruskal <- results_all %>% 
    group_by(scen) %>% 
    kruskal_test(Precision ~ method)

precision_dunn <- results_all %>%
    group_by(scen) %>%
    dunn_test(Precision ~ method, p.adjust.method="BH")

fpr_kruskal <- results_all %>% 
    group_by(scen) %>% 
    kruskal_test(false_positive_rate ~ method)

fpr_dunn <- results_all %>%
    group_by(scen) %>%
    dunn_test(false_positive_rate ~ method, p.adjust.method="BH")

# plot
p1 <- plot_comparisons(dunn=sensitivity_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=sensitivity_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=sensitivity_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

psensitivity <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")

p1 <- plot_comparisons(dunn=precision_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=precision_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=precision_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

pprecision <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")

p1 <- plot_comparisons(dunn=fpr_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=fpr_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=fpr_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

pfpr <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")

all_res <- ggarrange(psensitivity, pprecision, pfpr, nrow=3, ncol=1)
ggsave(all_res, filename = "output/taxoncounts/dunn_tests_plots.pdf", device="pdf", height=11, width=9)

# save statistics 
kruskal_table <- bind_rows(sensitivity_kruskal, precision_kruskal, fpr_kruskal)
dunn_table <- bind_rows(sensitivity_dunn, precision_dunn, fpr_dunn)
write_tsv(kruskal_table, path="output/taxoncounts/kruskal_wallis_method_comparison.txt", col_names=T)
write_tsv(dunn_table, path="output/taxoncounts/dunn_test_method_comparison.txt", col_names=T)



# associations classified by type (positive, negative or non significant)
results_class <- results %>% dplyr::filter(datatable!="all")

toplot <- results_class %>% 
    dplyr::select(method, scen, TP, TN, FP, FN, DA, method_type, datatable) %>% 
    gather(key="metric", value="value", -method, -scen, -method_type, -datatable) %>% 
    drop_na() %>% 
    group_by(method, scen, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup 

toplot$metric <- factor(toplot$metric, levels=c("TP", "FP", "TN", "FN", "DA"))
toplot$datatable <- factor(toplot$datatable, levels=c("positive", "negative", "non_correlated"))

p2 <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by=c("datatable", "scen"), palette="Dark2",
                legend.title="Metric", title="Associations of transformed taxon abundances with original microbial loads (classified by association type)",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p2, filename = "output/taxoncounts/plot_performance_taxoncounts_bytype.pdf", device = "pdf", width = 8, height=7, useDingbats=FALSE)

 
# now only scaled to the number of positives in the samples
actual_positives <- results_all %>% 
    dplyr::select(method, matrixnum, actual_positives)

results_class <- results_class %>% 
    left_join(actual_positives, by=c("method", "matrixnum"))

toplot2 <- results_class %>% 
    mutate(TP_scaled=100*true_positive/actual_positives) %>% 
    mutate(FP_scaled=100*false_positive/actual_positives) %>% 
    mutate(DA_scaled=100*discordant/actual_positives) %>% 
    dplyr::select(method, scen, method_type, TP_scaled, FP_scaled, DA_scaled, datatable) %>%
    drop_na() %>% 
    gather(key="metric", value="value", -method, -scen, -method_type, -datatable) %>% 
    group_by(method, scen, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot2$metric <- factor(toplot2$metric, levels=c("FP_scaled", "DA_scaled", "TP_scaled"))
toplot2$datatable <- factor(toplot2$datatable, levels=c("positive", "negative", "non_correlated"))

p2b <- ggbarplot(toplot2, x="method", y="value", fill="metric", facet.by=c("datatable", "scen"), palette="Dark2", 
                 subtitle="True Positives, False Positives and Discordant Associations scaled to the total actual positives (100%)",
                 legend.title="Metric", title="Associations of transformed taxon abundances with original microbial loads",
                 xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold")) +
    geom_hline(yintercept = 100, col="black", lwd=1)

ggsave(p2b, filename = "output/taxoncounts/plot_performance_taxoncounts_bytype_scaledpositives.pdf", device = "pdf", width = 8, height=7, useDingbats=FALSE)


#### Evaluate taxon-counts correlations only in special taxa ####
# define parameters for significance and initialize vectors for the evaluation data
significance_level <- 0.05

true_positive <- c()
false_positive <- c()
true_negative <- c()
false_negative <- c()
method <- c()
spread <- c()
matrixnum <- c()
datatable <- c()
discordant <- c()
scen <- c()
specialtaxon <- c()
matrix_stats <- read_tsv("data/raw/matrix_stats_3scenarios.tsv")

for(file in list.files("data/correlations_taxoncounts/", recursive = F, pattern = "pvalue", full.names = T)){
    # define file name, method used and reference to be compared against
    filename <- basename(file)
    methodname <- strsplit(filename, split = "_")[[1]][3]
    matrixname <- strsplit(filename, split = "_")[[1]][4]
    spreadname <- strsplit(filename, split = "_")[[1]][5]
    scenarioname <- strsplit(filename, split="_")[[1]][7] %>% gsub(., pattern="\\.tsv", replacement="")
    referencename <- gsub(filename, pattern=methodname, replacement="REAL")
    rhoname <- gsub(filename, pattern="pvalue", replacement="rho")
    referencerhoname <- gsub(referencename, pattern="pvalue", replacement="rho")
    if(scenarioname=="Blooming"){
        # read files
        test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
        reference_file <- read.table(paste0("data/correlations_taxoncounts/reference/", referencename),
                                     header=T, stringsAsFactors=F, sep="\t")
        rho_file <- read.table(paste0("data/correlations_taxoncounts/", rhoname),
                               header=T, stringsAsFactors = F, sep="\t")
        referencerho_file <- read.table(paste0("data/correlations_taxoncounts/reference/", referencerhoname),
                                        header=T, stringsAsFactors = F, sep = "\t")
        
        # special taxon (bloomer)
        sptax <- matrix_stats %>% dplyr::filter(matrix==matrixname) %>% pull(special_taxon)
        test_file <- test_file[sptax, ,drop=F]
        reference_file <- reference_file[sptax, ,drop=F]
        rho_file <- rho_file[sptax, ,drop=F]
        referencerho_file <- referencerho_file[sptax, ,drop=F]
        
        # all
        test_significant <- c(unlist(test_file < significance_level))
        reference_significant <- c(unlist(reference_file < significance_level))
        test_sign <- sign(c(unlist(rho_file)))
        reference_sign <- sign(c(unlist(referencerho_file)))
        
        # populate evaluation vectors
        method <- c(method, methodname)
        spread <- c(spread, spreadname)
        matrixnum <- c(matrixnum, matrixname)
        scen <- c(scen, scenarioname)
        datatable <- c(datatable, "all")
        true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
        false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
        true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
        false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
        discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
        specialtaxon <- c(specialtaxon, "Bloomer")
        
    }
    if(scenarioname=="Dysbiosis"){
        # read files
        test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
        reference_file <- read.table(paste0("data/correlations_taxoncounts/reference/", referencename),
                                     header=T, stringsAsFactors=F, sep="\t")
        rho_file <- read.table(paste0("data/correlations_taxoncounts/", rhoname),
                               header=T, stringsAsFactors = F, sep="\t")
        referencerho_file <- read.table(paste0("data/correlations_taxoncounts/reference/", referencerhoname),
                                        header=T, stringsAsFactors = F, sep = "\t")
        
        # special taxon (opportunist)
        sptax <- matrix_stats %>% dplyr::filter(matrix==matrixname) %>% pull(special_taxon)
        test_file <- test_file[sptax,, drop=F]
        reference_file <- reference_file[sptax, ,drop=F]
        rho_file <- rho_file[sptax,, drop=F]
        referencerho_file <- referencerho_file[sptax,, drop=F]
        
        # all
        test_significant <- c(unlist(test_file < significance_level))
        reference_significant <- c(unlist(reference_file < significance_level))
        test_sign <- sign(c(unlist(rho_file)))
        reference_sign <- sign(c(unlist(referencerho_file)))
        
        # populate evaluation vectors
        method <- c(method, methodname)
        spread <- c(spread, spreadname)
        matrixnum <- c(matrixnum, matrixname)
        scen <- c(scen, scenarioname)
        datatable <- c(datatable, "all")
        true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
        false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
        true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
        false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
        discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
        specialtaxon <- c(specialtaxon, "Opportunist")
        
        
        # read files
        test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
        reference_file <- read.table(paste0("data/correlations_taxoncounts/reference/", referencename),
                                     header=T, stringsAsFactors=F, sep="\t")
        rho_file <- read.table(paste0("data/correlations_taxoncounts/", rhoname),
                               header=T, stringsAsFactors = F, sep="\t")
        referencerho_file <- read.table(paste0("data/correlations_taxoncounts/reference/", referencerhoname),
                                        header=T, stringsAsFactors = F, sep = "\t")
        
        # special taxon (unresponsive)
        sptax <- matrix_stats %>% dplyr::filter(matrix==matrixname) %>% pull(flat_taxon_dysbiosis)
        sptax <- strsplit(sptax, split=",")[[1]]
        sptax <- intersect(sptax, rownames(test_file))
        test_file <- na.omit(c(unlist(test_file[sptax,, drop=F])))
        reference_file <- na.omit(c(unlist(reference_file[sptax,, drop=F])))
        rho_file <- na.omit(c(unlist(rho_file[sptax,, drop=F])))
        referencerho_file <- na.omit(c(unlist(referencerho_file[sptax,, drop=F])))
        
        # all
        test_significant <- c(unlist(test_file < significance_level))
        reference_significant <- c(unlist(reference_file < significance_level))
        test_sign <- sign(c(unlist(rho_file)))
        reference_sign <- sign(c(unlist(referencerho_file)))
        
        # populate evaluation vectors
        method <- c(method, methodname)
        spread <- c(spread, spreadname)
        matrixnum <- c(matrixnum, matrixname)
        scen <- c(scen, scenarioname)
        datatable <- c(datatable, "all")
        true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
        false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
        true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
        false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
        discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
        specialtaxon <- c(specialtaxon, "Unresponsive")
    }
}

# write results table
results <- tibble(method, spread, datatable, matrixnum, scen, true_positive, false_positive, true_negative, false_negative, discordant, specialtaxon)
results <- results %>%
    mutate(FDR=100*(false_positive+discordant)/(false_positive+discordant+true_positive)) %>% 
    mutate(Recall=100*true_positive/(true_positive+false_negative+discordant)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive+discordant)) %>% 
    mutate(Accuracy=100*(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative+discordant))%>%
    mutate(TP=100*true_positive/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FP=100*false_positive/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(TN=100*true_negative/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FN=100*false_negative/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(DA=100*discordant/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(false_positive_rate=100*(false_positive+discordant)/(true_negative+false_positive+discordant)) 


write_tsv(results, "output/taxoncounts/specialtaxon/statistics_taxoncounts_correlation_onlyspecialtaxa.tsv", col_names = T)
method_type <- tibble(method=c("Seq", "RMP", "Rel", "AST", "CLR", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "VST", "QMP", "ACS"),
                      method_type=c("Sequencing - non-transformed", 
                                    rep("Relative transformations",  times=2),
                                    rep("Compositional transformations", times=8),
                                    rep("Quantitative transformations", times=2)))


results <- results %>% 
    left_join(method_type, by="method")

results$method_type <- factor(results$method_type, 
                              levels=c("Sequencing - non-transformed", "Relative transformations", 
                                       "Compositional transformations",
                                       "Quantitative transformations"))
results$method <- factor(results$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

results$spread <- factor(results$spread, levels=c("low", "medium", "high"))
results$scen <- factor(results$scen, levels=c("Healthy", "Blooming","Dysbiosis"))
results_all <- results %>% dplyr::filter(datatable=="all")

results_all$F1 <- 2*((results_all$Precision * results_all$Recall)/(results_all$Precision + results_all$Recall))

# plot results
toplot <- results_all %>% 
    dplyr::select(method, specialtaxon, TP, TN, FP, FN, DA, method_type, datatable) %>% 
    gather(key="metric", value="value", -method, -specialtaxon, -method_type, -datatable) %>% 
    drop_na() %>% 
    group_by(method, specialtaxon, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup 

toplot$metric <- factor(toplot$metric, levels=c("TP", "FP", "TN", "FN", "DA"))
toplot$specialtaxon <- factor(toplot$specialtaxon, levels=c("Bloomer", "Opportunist", "Unresponsive"))

p3 <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by=c("specialtaxon"), palette="Dark2",
                legend.title="Metric", title="Associations of transformed taxon abundances with original microbial loads (special taxa)",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p3, filename = "output/taxoncounts/specialtaxon/plot_performance_taxoncounts_special.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)



# now only scaled to the number of positives in the samples
stats_matrix <- read_tsv("data/raw/matrix_stats_3scenarios.tsv", col_names=T)
actual_positives <- read_tsv("output/taxoncounts/statistics_taxoncounts_correlation.tsv") %>% 
    dplyr::filter(datatable=="all") %>% 
    mutate(actual_positives=true_positive+false_negative+discordant) %>% 
    dplyr::select(method, matrixnum, actual_positives) %>% 
    distinct()

results_all <- results_all %>% 
    left_join(actual_positives, by=c("method", "matrixnum"))

toplot <- results_all %>% 
    left_join(stats_matrix, by=c("matrixnum" = "matrix", "scen" = "scenario", "specialtaxon" = "special_taxon")) %>% 
    mutate(TP_scaled=100*true_positive/(actual_positives)) %>% 
    mutate(FP_scaled=100*false_positive/(actual_positives)) %>% 
    mutate(DA_scaled=100*discordant/(actual_positives)) %>% 
    dplyr::select(method, specialtaxon, method_type, TP_scaled, FP_scaled, DA_scaled, datatable) %>%
    drop_na() %>% 
    gather(key="metric", value="value", -method, -specialtaxon, -method_type, -datatable) %>% 
    group_by(method, specialtaxon, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot$metric <- factor(toplot$metric, levels=c("FP_scaled", "DA_scaled", "TP_scaled"))
toplot$specialtaxon <- factor(toplot$specialtaxon, levels=c("Bloomer", "Opportunist", "Unresponsive"))
toplot$method <- factor(toplot$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

p3b <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by=c("specialtaxon"), palette="Dark2",
                subtitle="True Positives, False Positives and Discordant Associations, scaled to the number of true associations",
                legend.title="Metric", title="Associations of transformed taxon abundances with original microbial loads (special taxa)",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p3b, filename = "output/taxoncounts/specialtaxon/plot_performance_taxoncounts_special_scaledpositives.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)




#### Evaluate correlations between raw and transformed data ####

results <- tibble()
for(file in list.files("data/all_matrices/", recursive = F, pattern = "taxonomy", full.names = T)){
    # define file name, method used and reference to be compared against
    filename <- basename(file)
    methodname <- strsplit(filename, split = "_")[[1]][1]
    matrixname <- strsplit(filename, split = "_")[[1]][3]
    spreadname <- strsplit(filename, split = "_")[[1]][4]
    scenarioname <- strsplit(filename, split="_")[[1]][6] %>% gsub(., pattern="\\.tsv", replacement="")
    referencename <- gsub(filename, pattern=methodname, replacement="REAL")
    
    # read files
    test_file <- read.table(file, header=T, stringsAsFactors = F)
    reference_file <- read.table(paste0("data/all_matrices/reference/", referencename),
                                 header=T, stringsAsFactors=F)
    
    # calculate correlation between raw and transformed data
    correlation_taxa <- sapply(1:nrow(test_file), function(i) cor(c(unlist(reference_file[i,])), c(unlist(test_file[i,]))))
    correlation_taxa_p <- sapply(1:nrow(test_file), function(i) cor.test(c(unlist(reference_file[i,])), c(unlist(test_file[i,])))[[3]])
    
    # results table
    results_tmp <- tibble(method=methodname, matrixname, spread=spreadname, scen=scenarioname, taxon=rownames(test_file), correlation=correlation_taxa, pvalue=correlation_taxa_p)
    results <- bind_rows(results, results_tmp)
}

results$pvalue <- p.adjust(results$pvalue, method="BH")

method_type <- tibble(method=c("Seq", "RMP", "Rel", "AST", "CLR", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "VST", "QMP", "ACS"),
                      method_type=c("Sequencing - non-transformed", 
                                    rep("Relative transformations",  times=2),
                                    rep("Compositional transformations", times=8),
                                    rep("Quantitative transformations", times=2)))


results <- results %>% 
    left_join(method_type, by="method")

results$method_type <- factor(results$method_type, 
                              levels=c("Sequencing - non-transformed", "Relative transformations", 
                                       "Compositional transformations",
                                       "Quantitative transformations"))
results$method <- factor(results$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

results$spread <- factor(results$spread, levels=c("low", "medium", "high"))
results$scen <- factor(results$scen, levels=c("Healthy", "Blooming","Dysbiosis"))

p1 <- ggboxplot(results, x="method", y="correlation", fill="method_type", 
          facet.by="scen", palette="Dark2", ylab="R", xlab="Method", 
          title="Correlation between raw and transformed taxa abundances", legend.title="Method type") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p1, filename = "output/raw_transformed_correlation/plot_raw_transformed_correlation_values.pdf", device = "pdf", width = 9, height=3.5, useDingbats=FALSE)



# compare methods according to the correlation between raw and transformed taxa
correlation_kruskal <- results %>% 
    group_by(scen) %>% 
    kruskal_test(correlation ~ method)

correlation_dunn <- results %>%
    group_by(scen) %>%
    dunn_test(correlation ~ method, p.adjust.method="BH")

# plot
p1 <- plot_comparisons(dunn=correlation_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=correlation_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=correlation_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

pcorrelation <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")
ggsave(pcorrelation, filename = "output/raw_transformed_correlation/dunn_test_plots.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)

write_tsv(correlation_kruskal, path = "output/raw_transformed_correlation/kruskal_test_correlations.txt", col_names=T)
write_tsv(correlation_dunn, path="output/raw_transformed_correlation/dunn_test_correlations.txt", col_names=T)


# plot by categories
results <- results %>% 
    mutate(category="High correlation (R>0.8)")

results[results$correlation<0.8 & results$pvalue<0.05,"category"] <- "Moderate correlation (R>0.5)"
results[results$correlation<0.5 & results$pvalue<0.05,"category"] <- "Mild correlation (R<0.5)"
results[results$pvalue>=0.05,"category"] <- "Non-significant correlation"
results[results$correlation<0 & results$pvalue<0.05,"category"] <- "Negative correlation (R<0)"

results_table <- results %>% 
    dplyr::select(method, scen, category) %>% 
    table %>% 
    as_tibble() %>% 
    group_by(method, scen) %>% 
    mutate(n=100* n/sum(n))

results_table$category <- factor(results_table$category,  levels=c("High correlation (R>0.8)",
                                                                   "Moderate correlation (R>0.5)",
                                                                   "Mild correlation (R<0.5)",
                                                                   "Non-significant correlation",
                                                                   "Negative correlation (R<0)"))
results_table$scen <- factor(results_table$scen, levels=c("Healthy", "Blooming","Dysbiosis"))

p2 <- ggbarplot(results_table, x="method", y="n", fill="category", facet.by="scen",
          palette="Dark2", title="Classification of transformed taxa according to their correlation with the raw values",
          legend.title="Category", ylab="Percentage of taxa (average per method)", xlab="Method") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p2, filename = "output/raw_transformed_correlation/plot_raw_transformed_correlation_values_category.pdf", device = "pdf", width = 9, height=3.5, useDingbats=FALSE)


# focus on the special taxa only
matrix_stats <- read_tsv("data/raw/matrix_stats_3scenarios.tsv")
bloomer_taxa <- matrix_stats %>% 
    dplyr::filter(scenario=="Blooming") %>% 
    dplyr::select(matrix, special_taxon) %>% 
    mutate(type="Bloomer")
opportunist_taxa <- matrix_stats %>% 
    dplyr::filter(scenario=="Dysbiosis") %>% 
    dplyr::select(matrix, special_taxon) %>% 
    mutate(type="Opportunist")
flat_taxa <- matrix_stats %>% 
    dplyr::filter(scenario=="Dysbiosis") %>% 
    dplyr::select(matrix, flat_taxon_dysbiosis) %>% 
    separate(col="flat_taxon_dysbiosis", into=paste0("tax", 1:11), sep = ",") %>% 
    gather(key="t", value="special_taxon", -matrix) %>% 
    dplyr::select(-t) %>% 
    mutate(type="Flat")

special_tax <- bind_rows(bloomer_taxa,
                         opportunist_taxa,
                         flat_taxa) %>% 
    drop_na() %>% 
    as.data.frame()

results <- results %>% mutate(special_taxa="none")
for(i in 1:nrow(special_tax)){
    results[results$matrixname==special_tax$matrix[i] & results$taxon==special_tax$special_taxon[i], "special_taxa"] <- special_tax$type[i]
}

results_special <- results %>% 
    dplyr::filter(special_taxa!="none")

p3 <- ggboxplot(results_special, x="method", y="correlation", fill="method_type", 
                facet.by="special_taxa", palette="Dark2", ylab="R", xlab="Method", 
                title="Correlation between raw and transformed taxa abundances", legend.title="Method type") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p3, filename = "output/raw_transformed_correlation/plot_raw_transformed_correlation_values_specialtaxa.pdf", device = "pdf", width = 10, height=5, useDingbats=FALSE)

write_tsv(results, path = "output/raw_transformed_correlation/correlation_raw_transformed_data.tsv", col_names = T)






#### Evaluate taxon-taxon correlations ####

# define parameters for significance and initialize vectors for the evaluation data
significance_level <- 0.05

true_positive <- c()
true_positive_discordant <- c()
false_positive <- c()
true_negative <- c()
false_negative <- c()
discordant <- c()
method <- c()
spread <- c()
matrixnum <- c()
datatable <- c()
scen <- c()

discordant_taxa <- read_tsv("output/taxoncounts/discordant_taxa.tsv")

for(file in list.files("data/correlations_taxontaxon", recursive = F, pattern = "pvalue", full.names = T)){
    # define file name, method used and reference to be compared against
    filename <- basename(file)
    methodname <- strsplit(filename, split = "_")[[1]][3]
    matrixname <- strsplit(filename, split = "_")[[1]][4]
    spreadname <- strsplit(filename, split = "_")[[1]][5]
    scenarioname <- strsplit(filename, split="_")[[1]][7] %>% gsub(., pattern="\\.tsv", replacement="")
    referencename <- gsub(filename, pattern=methodname, replacement="REAL")
    rhoname <- gsub(filename, pattern="pvalue", replacement="rho")
    referencerhoname <- gsub(referencename, pattern="pvalue", replacement="rho")
    # read files
    test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
    reference_file <- read.table(paste0("data/correlations_taxontaxon/reference/", referencename),
                                 header=T, stringsAsFactors=F, sep="\t")
    rho_file <- read.table(paste0("data/correlations_taxontaxon/", rhoname),
                           header=T, stringsAsFactors = F, sep="\t")
    referencerho_file <- read.table(paste0("data/correlations_taxontaxon/reference/", referencerhoname),
                                    header=T, stringsAsFactors = F, sep = "\t")
    
    # all
    test_significant <- na.omit(c(unlist(t(test_file < significance_level))))
    reference_significant <- na.omit(c(unlist(t(reference_file < significance_level))))
    test_sign <- na.omit(sign(c(unlist(t(rho_file)))))
    reference_sign <- na.omit(sign(c(unlist(t(referencerho_file)))))
    
    #combn_names
    pairs_taxa <- combn(rownames(test_file), 2, simplify = T) %>% t()
    
    # discordant_taxa in this setting
    discordant_names <- discordant_taxa %>% 
        filter(matrix==matrixname, method==methodname) %>% 
        pull(taxa)
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "all")
    # take into account those TPs that are TPs because both taxa are discordant
    tptax <- pairs_taxa[which(test_significant & reference_significant & test_sign==reference_sign & test_sign==1),,drop=F ] 
    tptax2 <- pairs_taxa[which(test_significant & reference_significant & test_sign==reference_sign & test_sign==-1),,drop=F ] 
    tptax_discordant <- tptax[which(tptax[,1] %in% discordant_names & tptax[,2] %in% discordant_names),,drop=F ]
    tptax_discordant2 <- tptax2[which(tptax2[,1] %in% discordant_names | tptax2[,2] %in% discordant_names),, drop=F]
    true_positive <- c(true_positive, nrow(tptax)+nrow(tptax2)-nrow(tptax_discordant)-nrow(tptax_discordant2))
    true_positive_discordant <- c(true_positive_discordant, nrow(tptax_discordant)+nrow(tptax_discordant2))
    
    false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
    true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
    false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
    discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
    
    # only positive correlations
    reference_sign_pos <- reference_sign
    reference_sign_pos[reference_sign<0 | reference_significant==F] <- NA
    test_sign_pos <- test_sign
    test_sign_pos[reference_sign<0 | reference_significant==F] <- NA
    reference_significant_pos <- reference_significant
    reference_significant_pos[reference_sign<0 | reference_significant==F] <- NA
    test_significant_pos <- test_significant
    test_significant_pos[reference_sign<0 | reference_significant==F] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "positive")
    # take into account those TPs that are TPs because both taxa are discordant
    tptax <- pairs_taxa[which(test_significant_pos & reference_sign_pos & test_sign_pos==reference_sign_pos),,drop=F ] 
    tptax_discordant <- tptax[which(tptax[,1] %in% discordant_names & tptax[,2] %in% discordant_names),,drop=F ]
    true_positive <- c(true_positive, nrow(tptax)-nrow(tptax_discordant))
    true_positive_discordant <- c(true_positive_discordant, nrow(tptax_discordant))

    false_positive <- c(false_positive, length(which(test_significant_pos & !reference_significant_pos)))
    true_negative <- c(true_negative, length(which(!test_significant_pos & !reference_significant_pos)))
    false_negative <- c(false_negative, length(which(!test_significant_pos & reference_significant_pos)))
    discordant <- c(discordant, length(which(test_significant_pos & reference_significant_pos & test_sign_pos!=reference_sign_pos)))
    
    # only negative correlations
    reference_sign_neg <- reference_sign
    reference_sign_neg[reference_sign>0 | reference_significant==F] <- NA
    test_sign_neg <- test_sign
    test_sign_neg[reference_sign>0 | reference_significant==F] <- NA
    reference_significant_neg <- reference_significant
    reference_significant_neg[reference_sign>0 | reference_significant==F] <- NA
    test_significant_neg <- test_significant
    test_significant_neg[reference_sign>0 | reference_significant==F] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "negative")
    # take into account those TPs that are TPs because both taxa are discordant
    tptax <- pairs_taxa[which(test_significant_neg & reference_significant_neg & test_sign_neg==reference_sign_neg),,drop=F ] 
    tptax_discordant <- tptax[which(tptax[,1] %in% discordant_names | tptax[,2] %in% discordant_names),,drop=F ]
    true_positive <- c(true_positive, nrow(tptax)-nrow(tptax_discordant))
    true_positive_discordant <- c(true_positive_discordant, nrow(tptax_discordant))
    
    false_positive <- c(false_positive, length(which(test_significant_neg & !reference_significant_neg)))
    true_negative <- c(true_negative, length(which(!test_significant_neg & !reference_significant_neg)))
    false_negative <- c(false_negative, length(which(!test_significant_neg & reference_significant_neg)))
    discordant <- c(discordant, length(which(test_significant_neg & reference_significant_neg & test_sign_neg!=reference_sign_neg)))
    
    # only not correlated taxa
    reference_sign_neg <- reference_sign
    reference_sign_neg[reference_significant==T] <- NA
    test_sign_neg <- test_sign
    test_sign_neg[reference_significant==T] <- NA
    reference_significant_neg <- reference_significant
    reference_significant_neg[reference_significant==T] <- NA
    test_significant_neg <- test_significant
    test_significant_neg[reference_significant==T] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "non_correlated")
    # take into account those TPs that are TPs because both taxa are discordant
    tptax <- pairs_taxa[which(test_significant_neg & reference_significant_neg & test_sign_neg==reference_sign_neg),,drop=F ] 
    tptax_discordant <- tptax[which(tptax[,1] %in% discordant_names & tptax[,2] %in% discordant_names),,drop=F ]
    true_positive <- c(true_positive, nrow(tptax)-nrow(tptax_discordant))
    true_positive_discordant <- c(true_positive_discordant, nrow(tptax_discordant))
    
    false_positive <- c(false_positive, length(which(test_significant_neg & !reference_significant_neg)))
    true_negative <- c(true_negative, length(which(!test_significant_neg & !reference_significant_neg)))
    false_negative <- c(false_negative, length(which(!test_significant_neg & reference_significant_neg)))
    discordant <- c(discordant, length(which(test_significant_neg & reference_significant_neg & test_sign_neg!=reference_sign_neg)))
    
}

# write results table
results <- tibble(method, spread, datatable, matrixnum, scen, true_positive, true_positive_discordant,false_positive, true_negative, false_negative, discordant)
results <- results %>%
    mutate(FDR=100*(false_positive+discordant)/(false_positive+discordant+true_positive+true_positive_discordant)) %>% 
    mutate(Recall=100*(true_positive+true_positive_discordant)/(true_positive+true_positive_discordant+false_negative+discordant)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive+discordant)) %>% 
    mutate(Accuracy=100*(true_positive+true_positive_discordant+true_negative)/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant))%>%
    mutate(TP=100*true_positive/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(TPD=100*true_positive_discordant/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FP=100*false_positive/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(TN=100*true_negative/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FN=100*false_negative/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(DA=100*discordant/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(false_positive_rate=100*(false_positive+discordant)/(true_negative+false_positive+discordant)) 


write_tsv(results, "output/taxontaxon/statistics_taxontaxon_correlation.tsv", col_names = T)
method_type <- tibble(method=c("Seq", "RMP", "Rel", "AST", "CLR", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "VST", "QMP", "ACS"),
                      method_type=c("Sequencing - non-transformed", 
                                    rep("Traditional transformations",  times=2),
                                    rep("Compositional transformations", times=8),
                                    rep("Transformations incorporating microbial loads", times=2)))


results <- results %>% 
    left_join(method_type, by="method")

results$method_type <- factor(results$method_type, 
                                  levels=c("Sequencing - non-transformed", "Traditional transformations", 
                                           "Compositional transformations",
                                           "Transformations incorporating microbial loads"))
results$method <- factor(results$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

results$spread <- factor(results$spread, levels=c("low", "medium", "high"))
results$scen <- factor(results$scen, levels=c("Healthy", "Blooming","Dysbiosis"))
results_all <- results %>% dplyr::filter(datatable=="all")

results_all$F1 <- 2*((results_all$Precision * results_all$Recall)/(results_all$Precision + results_all$Recall))

# plot results
toplot <- results_all %>% 
    dplyr::select(method, scen, TP, TPD, TN, FP, FN, DA, method_type) %>% 
    gather(key="metric", value="value", -method, -scen, -method_type) %>% 
    group_by(method, scen, method_type, metric) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot$metric <- factor(toplot$metric, levels=c("TP", "TPD", "FP", "TN", "FN", "DA"))

p1 <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by="scen", palette="Dark2",
                legend.title="Metric", title="Taxon-taxon associations",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p1, filename = "output/taxontaxon/plot_performance_taxontaxon.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)

# now only scaled to the number of positives in the samples
results_all <- results_all %>% 
    mutate(actual_positives=true_positive+false_negative+discordant+true_positive_discordant) 

toplot2 <- results_all %>% 
    mutate(TP_scaled=100*true_positive/actual_positives) %>% 
    mutate(TPD_scaled=100*true_positive_discordant/actual_positives) %>% 
    mutate(FP_scaled=100*false_positive/actual_positives) %>% 
    mutate(DA_scaled=100*discordant/actual_positives) %>% 
    dplyr::select(method, scen, method_type, TP_scaled,TPD_scaled, FP_scaled, DA_scaled) %>% 
    gather(key="metric", value="value", -method, -scen, -method_type) %>% 
    group_by(method, scen, method_type, metric) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot2$metric <- factor(toplot2$metric, levels=c("FP_scaled", "DA_scaled", "TPD_scaled", "TP_scaled"))

palette_custom <-  get_palette("Dark2",4)[c(1,2,4,3)]
p1b <- ggbarplot(toplot2, x="method", y="value", fill="metric", facet.by="scen", palette=palette_custom, 
                 subtitle="True Positives, False Positives and Discordant Associations scaled to the total positives (100%)",
                 legend.title="Metric", title="Taxon-taxon associations",
                 xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold")) +
    geom_hline(yintercept = 100, col="black", lwd=1)

ggsave(p1b, filename = "output/taxontaxon/plot_performance_taxontaxon_scaledpositives.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)


# associations classified by type (positive, negative or non significant)
results_class <- results %>% dplyr::filter(datatable!="all")

toplot <- results_class %>% 
    dplyr::select(method, scen, TP, TPD, TN, FP, FN, DA, method_type, datatable) %>% 
    gather(key="metric", value="value", -method, -scen, -method_type, -datatable) %>% 
    drop_na() %>% 
    group_by(method, scen, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup 

toplot$metric <- factor(toplot$metric, levels=c("TP", "TPD", "FP", "TN", "FN", "DA"))
toplot$datatable <- factor(toplot$datatable, levels=c("positive", "negative", "non_correlated"))

p2 <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by=c("datatable", "scen"), palette="Dark2",
                legend.title="Metric", title="Taxon-taxon associations (classified by true association type)",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p2, filename = "output/taxontaxon/plot_performance_taxontaxon_bytype.pdf", device = "pdf", width =8, height=7, useDingbats=FALSE)


# now only scaled to the number of positives in the samples
actual_positives <- results_all %>% 
    dplyr::select(method, matrixnum, actual_positives)

results_class <- results_class %>% 
    left_join(actual_positives, by=c("method", "matrixnum"))

toplot2 <- results_class %>% 
    mutate(TP_scaled=100*true_positive/actual_positives) %>% 
    mutate(TPD_scaled=100*true_positive_discordant/actual_positives) %>% 
    mutate(FP_scaled=100*false_positive/actual_positives) %>% 
    mutate(DA_scaled=100*discordant/actual_positives) %>% 
    dplyr::select(method, scen, method_type, TP_scaled, TPD_scaled, FP_scaled, DA_scaled, datatable) %>%
    drop_na() %>% 
    gather(key="metric", value="value", -method, -scen, -method_type, -datatable) %>% 
    group_by(method, scen, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot2$metric <- factor(toplot2$metric, levels=c("FP_scaled", "DA_scaled", "TPD_scaled", "TP_scaled"))
toplot2$datatable <- factor(toplot2$datatable, levels=c("positive", "negative", "non_correlated"))

p2b <- ggbarplot(toplot2, x="method", y="value", fill="metric", facet.by=c("datatable", "scen"), palette=palette_custom, 
                 subtitle="True Positives, False Positives and Discordant Associations scaled to the total actual positives (100%)",
                 legend.title="Metric", title="Taxon-taxon associations (classified by true association type)",
                 xlab="Method", ylab="Percentage", scales="free_y") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p2b, filename = "output/taxontaxon/plot_performance_taxontaxon_bytype_scaledpositives.pdf", device = "pdf", width = 8, height=7, useDingbats=FALSE)

# plot comparison of positive and negative associations
comp <- results %>% 
    dplyr::filter(datatable %in% c("positive", "negative")) %>% 
    dplyr::select(method, datatable, scen, Recall, Precision, false_positive_rate)
comp$method <- factor(comp$method, levels=c("Seq", "RMP", "Rel", "AST", "CLR", "CSS", "GMPR",
                                      "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))
comp$scen <- factor(comp$scen, levels=c("Healthy", "Blooming", "Dysbiosis"))


p1 <- ggboxplot(comp, x="datatable", y="Recall", fill="datatable", palette="Dark2",
          facet.by=c("scen", "method"), ylim=c(0,120), ylab="Sensitivity") + 
    theme_bw() + 
    rotate_x_text(45) + 
    stat_compare_means(comparisons = list(c("negative", "positive")), 
                       label = "p.signif")
p2 <- ggboxplot(comp, x="datatable", y="Precision", fill="datatable", palette="Dark2",
                facet.by=c("scen", "method"), ylim=c(0,120), ylab="Precision") + 
    theme_bw() + 
    rotate_x_text(45) + 
    stat_compare_means(comparisons = list(c("negative", "positive")), 
                       label = "p.signif")
pfin <- ggarrange(p1, p2, nrow=2)
ggsave(pfin, filename = "output/taxontaxon/comparison_positive_negative_interactions.pdf", device="pdf", height=8, width=9,useDingbats=FALSE)

# calculate sensitivity, precision and FPR 
sensitivity_kruskal <- results_all %>% 
    group_by(scen) %>% 
    kruskal_test(Recall ~ method)

sensitivity_dunn <- results_all %>%
    group_by(scen) %>%
    dunn_test(Recall ~ method, p.adjust.method="BH")

precision_kruskal <- results_all %>% 
    group_by(scen) %>% 
    kruskal_test(Precision ~ method)

precision_dunn <- results_all %>%
    group_by(scen) %>%
    dunn_test(Precision ~ method, p.adjust.method="BH")

fpr_kruskal <- results_all %>% 
    group_by(scen) %>% 
    kruskal_test(false_positive_rate ~ method)

fpr_dunn <- results_all %>%
    group_by(scen) %>%
    dunn_test(false_positive_rate ~ method, p.adjust.method="BH")

# plot
p1 <- plot_comparisons(dunn=sensitivity_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=sensitivity_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=sensitivity_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

psensitivity <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")

p1 <- plot_comparisons(dunn=precision_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=precision_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=precision_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

pprecision <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")

p1 <- plot_comparisons(dunn=fpr_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=fpr_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=fpr_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

pfpr <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")

all_res <- ggarrange(psensitivity, pprecision, pfpr, nrow=3, ncol=1)
ggsave(all_res, filename = "output/taxontaxon/dunn_tests_plots.pdf", device="pdf", height=11, width=9)
# save statistics 
kruskal_table <- bind_rows(sensitivity_kruskal, precision_kruskal, fpr_kruskal)
dunn_table <- bind_rows(sensitivity_dunn, precision_dunn, fpr_dunn)

dunn_table <- dunn_table %>% 
    dplyr::filter(!(.y.=="false_positive_rate" & scen=="Healthy")) # remove non-significant data from kruskal - only kept for plotting 
write_tsv(kruskal_table, path="output/taxontaxon/kruskal_wallis_method_comparison.txt", col_names=T)
write_tsv(dunn_table, path="output/taxontaxon/dunn_test_method_comparison.txt", col_names=T)




#### Evaluate taxon-taxon correlations only in special taxa ####
# define parameters for significance and initialize vectors for the evaluation data
significance_level <- 0.05

true_positive <- c()
true_positive_discordant <- c()
false_positive <- c()
true_negative <- c()
false_negative <- c()
discordant <- c()
method <- c()
spread <- c()
matrixnum <- c()
datatable <- c()
scen <- c()
specialtaxon <- c()
matrix_stats <- read_tsv("data/raw/matrix_stats_3scenarios.tsv")

for(file in list.files("data/correlations_taxontaxon/", recursive = F, pattern = "pvalue", full.names = T)){
    # define file name, method used and reference to be compared against
    filename <- basename(file)
    methodname <- strsplit(filename, split = "_")[[1]][3]
    matrixname <- strsplit(filename, split = "_")[[1]][4]
    spreadname <- strsplit(filename, split = "_")[[1]][5]
    scenarioname <- strsplit(filename, split="_")[[1]][7] %>% gsub(., pattern="\\.tsv", replacement="")
    referencename <- gsub(filename, pattern=methodname, replacement="REAL")
    rhoname <- gsub(filename, pattern="pvalue", replacement="rho")
    referencerhoname <- gsub(referencename, pattern="pvalue", replacement="rho")
    if(scenarioname=="Blooming"){
        # read files
        test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
        reference_file <- read.table(paste0("data/correlations_taxontaxon/reference/", referencename),
                                     header=T, stringsAsFactors=F, sep="\t")
        rho_file <- read.table(paste0("data/correlations_taxontaxon/", rhoname),
                               header=T, stringsAsFactors = F, sep="\t")
        referencerho_file <- read.table(paste0("data/correlations_taxontaxon/reference/", referencerhoname),
                                        header=T, stringsAsFactors = F, sep = "\t")
        
        # special taxon (bloomer)
        nm <- rep(colnames(test_file),times=2) # keep the names for the generated vectors
        sptax <- matrix_stats %>% dplyr::filter(matrix==matrixname) %>% pull(special_taxon)
        test_file <-c(unlist(test_file[sptax,, drop=F]), unlist(test_file[,sptax, drop=F]))
        names(test_file) <- nm
        test_file <- na.omit(test_file)
        reference_file <- na.omit(c(unlist(reference_file[sptax,, drop=F]), unlist(reference_file[,sptax, drop=F])))
        rho_file <- na.omit(c(unlist(rho_file[sptax,, drop=F]), unlist(rho_file[,sptax, drop=F])))
        referencerho_file <- na.omit(c(unlist(referencerho_file[sptax,, drop=F]), unlist(referencerho_file[,sptax, drop=F])))
        names(reference_file) <- names(test_file)
        names(rho_file) <- names(test_file)
        names(referencerho_file) <- names(test_file)
        
        # all
        test_significant <- c(unlist(test_file < significance_level))
        reference_significant <- c(unlist(reference_file < significance_level))
        test_sign <- sign(c(unlist(rho_file)))
        reference_sign <- sign(c(unlist(referencerho_file)))

        
        # discordant_taxa in this setting
        discordant_names <- discordant_taxa %>% 
            filter(matrix==matrixname, method==methodname) %>% 
            pull(taxa)
        
        # populate evaluation vectors
        method <- c(method, methodname)
        spread <- c(spread, spreadname)
        matrixnum <- c(matrixnum, matrixname)
        scen <- c(scen, scenarioname)
        datatable <- c(datatable, "all")
        
        # take into account those TPs that are TPs because both taxa are discordant
        tptax <- names(which(test_significant & reference_significant & test_sign==reference_sign & test_sign==1))
        tptax2 <- names(which(test_significant & reference_significant & test_sign==reference_sign & test_sign==-1))
        tptax_discordant <- tptax[which(tptax %in% discordant_names & sptax %in% discordant_names)]
        tptax_discordant2 <- tptax2[which(xor(tptax2 %in% discordant_names, sptax %in% discordant_names))]
        true_positive <- c(true_positive, length(tptax)+length(tptax2)-length(tptax_discordant)-length(tptax_discordant2))
        true_positive_discordant <- c(true_positive_discordant, length(tptax_discordant)+length(tptax_discordant2))
        false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
        true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
        false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
        discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
        specialtaxon <- c(specialtaxon, "Bloomer")

    }
    if(scenarioname=="Dysbiosis"){
        # read files
        test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
        reference_file <- read.table(paste0("data/correlations_taxontaxon/reference/", referencename),
                                     header=T, stringsAsFactors=F, sep="\t")
        rho_file <- read.table(paste0("data/correlations_taxontaxon/", rhoname),
                               header=T, stringsAsFactors = F, sep="\t")
        referencerho_file <- read.table(paste0("data/correlations_taxontaxon/reference/", referencerhoname),
                                        header=T, stringsAsFactors = F, sep = "\t")
        
        # special taxon (opportunist)
        nm <- rep(colnames(test_file),times=2) # keep the names for the generated vectors
        sptax <- matrix_stats %>% dplyr::filter(matrix==matrixname) %>% pull(special_taxon)
        test_file <-c(unlist(test_file[sptax,, drop=F]), unlist(test_file[,sptax, drop=F]))
        names(test_file) <- nm
        test_file <- na.omit(test_file)
        reference_file <- na.omit(c(unlist(reference_file[sptax,, drop=F]), unlist(reference_file[,sptax, drop=F])))
        rho_file <- na.omit(c(unlist(rho_file[sptax,, drop=F]), unlist(rho_file[,sptax, drop=F])))
        referencerho_file <- na.omit(c(unlist(referencerho_file[sptax,, drop=F]), unlist(referencerho_file[,sptax, drop=F])))
        names(reference_file) <- names(test_file)
        names(rho_file) <- names(test_file)
        names(referencerho_file) <- names(test_file)
        
        # all
        test_significant <- c(unlist(test_file < significance_level))
        reference_significant <- c(unlist(reference_file < significance_level))
        test_sign <- sign(c(unlist(rho_file)))
        reference_sign <- sign(c(unlist(referencerho_file)))
        
        # discordant_taxa in this setting
        discordant_names <- discordant_taxa %>% 
            filter(matrix==matrixname, method==methodname) %>% 
            pull(taxa)
        
        # populate evaluation vectors
        method <- c(method, methodname)
        spread <- c(spread, spreadname)
        matrixnum <- c(matrixnum, matrixname)
        scen <- c(scen, scenarioname)
        datatable <- c(datatable, "all")
        
        # take into account those TPs that are TPs because both taxa are discordant
        tptax <- names(which(test_significant & reference_significant & test_sign==reference_sign & test_sign==1))
        tptax2 <- names(which(test_significant & reference_significant & test_sign==reference_sign & test_sign==-1))
        tptax_discordant <- tptax[which(tptax %in% discordant_names & sptax %in% discordant_names)]
        tptax_discordant2 <- tptax2[which(xor(tptax2 %in% discordant_names, sptax %in% discordant_names))]
        true_positive <- c(true_positive, length(tptax)+length(tptax2)-length(tptax_discordant)-length(tptax_discordant2))
        true_positive_discordant <- c(true_positive_discordant, length(tptax_discordant)+length(tptax_discordant2))
        false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
        true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
        false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
        discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
        specialtaxon <- c(specialtaxon, "Opportunist")
        
        
        # read files
        test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
        reference_file <- read.table(paste0("data/correlations_taxontaxon/reference/", referencename),
                                     header=T, stringsAsFactors=F, sep="\t")
        rho_file <- read.table(paste0("data/correlations_taxontaxon/", rhoname),
                               header=T, stringsAsFactors = F, sep="\t")
        referencerho_file <- read.table(paste0("data/correlations_taxontaxon/reference/", referencerhoname),
                                        header=T, stringsAsFactors = F, sep = "\t")
        
        # special taxon (unresponsive)
        sptax <- matrix_stats %>% dplyr::filter(matrix==matrixname) %>% pull(flat_taxon_dysbiosis)
        sptax <- strsplit(sptax, split=",")[[1]]
        sptax <- intersect(sptax, rownames(test_file))
        test_file <- na.omit(c(unlist(test_file[sptax,, drop=F])))
        reference_file <- na.omit(c(unlist(reference_file[sptax,, drop=F])))
        rho_file <- na.omit(c(unlist(rho_file[sptax,, drop=F])))
        referencerho_file <- na.omit(c(unlist(referencerho_file[sptax,, drop=F])))
        
        # all
        test_significant <- c(unlist(test_file < significance_level))
        reference_significant <- c(unlist(reference_file < significance_level))
        test_sign <- sign(c(unlist(rho_file)))
        reference_sign <- sign(c(unlist(referencerho_file)))
        
        # populate evaluation vectors
        method <- c(method, methodname)
        spread <- c(spread, spreadname)
        matrixnum <- c(matrixnum, matrixname)
        scen <- c(scen, scenarioname)
        datatable <- c(datatable, "all")
        true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
        true_positive_discordant <- c(true_positive_discordant, 0) # because this does not correlate with cell counts, there cannot be true associations 
        false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
        true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
        false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
        discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
        specialtaxon <- c(specialtaxon, "Unresponsive")

    }
}

# write results table
results <- tibble(method, spread, datatable, matrixnum, scen, true_positive, true_positive_discordant,false_positive, true_negative, false_negative, discordant, specialtaxon)
results <- results %>%
    mutate(FDR=100*(false_positive+discordant)/(false_positive+discordant+true_positive+true_positive_discordant)) %>% 
    mutate(Recall=100*(true_positive+true_positive_discordant)/(true_positive+true_positive_discordant+false_negative+discordant)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive+discordant)) %>% 
    mutate(Accuracy=100*(true_positive+true_positive_discordant+true_negative)/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant))%>%
    mutate(TP=100*true_positive/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(TPD=100*true_positive_discordant/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FP=100*false_positive/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(TN=100*true_negative/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FN=100*false_negative/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(DA=100*discordant/(true_positive+true_positive_discordant+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(false_positive_rate=100*(false_positive+discordant)/(true_negative+false_positive+discordant)) 


write_tsv(results, "output/taxontaxon/specialtaxon/statistics_taxontaxon_correlation_onlyspecialtaxa.tsv", col_names = T)
method_type <- tibble(method=c("Seq", "RMP", "Rel", "AST", "CLR", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "VST", "QMP", "ACS"),
                      method_type=c("Sequencing - non-transformed", 
                                    rep("Relative transformations",  times=2),
                                    rep("Compositional transformations", times=8),
                                    rep("Quantitative transformations", times=2)))


results <- results %>% 
    left_join(method_type, by="method")

results$method_type <- factor(results$method_type, 
                              levels=c("Sequencing - non-transformed", "Relative transformations", 
                                       "Compositional transformations",
                                       "Quantitative transformations"))
results$method <- factor(results$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

results$spread <- factor(results$spread, levels=c("low", "medium", "high"))
results$scen <- factor(results$scen, levels=c("Healthy", "Blooming","Dysbiosis"))
results_all <- results %>% dplyr::filter(datatable=="all")

results_all$F1 <- 2*((results_all$Precision * results_all$Recall)/(results_all$Precision + results_all$Recall))

# plot results
toplot <- results_all %>% 
    dplyr::select(method, specialtaxon, TP, TPD, TN, FP, FN, DA, method_type, datatable) %>% 
    gather(key="metric", value="value", -method, -specialtaxon, -method_type, -datatable) %>% 
    drop_na() %>% 
    group_by(method, specialtaxon, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup 

toplot$metric <- factor(toplot$metric, levels=c("TP", "TPD", "FP", "TN", "FN", "DA"))
toplot$specialtaxon <- factor(toplot$specialtaxon, levels=c("Bloomer", "Opportunist", "Unresponsive"))

p3 <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by=c("specialtaxon"), palette="Dark2",
                legend.title="Metric", title="Taxon-taxon associations (special taxa)",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p3, filename = "output/taxontaxon/specialtaxon/plot_performance_taxontaxon_special.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)


# now only scaled to the number of positives in the samples

stats_matrix <- read_tsv("data/raw/matrix_stats_3scenarios.tsv", col_names=T)
actual_positives <- read_tsv("output/taxontaxon/statistics_taxontaxon_correlation.tsv") %>% 
    dplyr::filter(datatable=="all") %>% 
    mutate(actual_positives=true_positive+true_positive_discordant+false_negative+discordant) %>% 
    dplyr::select(method, matrixnum, actual_positives) %>% 
    distinct()

results_all <- results_all %>% 
    left_join(actual_positives, by=c("method", "matrixnum"))


toplot <- results_all %>% 
    mutate(TP_scaled=100*true_positive/actual_positives) %>% 
    mutate(TPD_scaled=100*true_positive_discordant/actual_positives) %>% 
    mutate(FP_scaled=100*false_positive/actual_positives) %>% 
    mutate(DA_scaled=100*discordant/actual_positives) %>% 
    dplyr::select(method, specialtaxon, method_type, TP_scaled, TPD_scaled, FP_scaled, DA_scaled, datatable) %>%
    drop_na() %>% 
    gather(key="metric", value="value", -method, -specialtaxon, -method_type, -datatable) %>% 
    group_by(method, specialtaxon, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()


toplot$metric <- factor(toplot$metric, levels=c("FP_scaled", "DA_scaled","TPD_scaled", "TP_scaled"))
toplot$specialtaxon <- factor(toplot$specialtaxon, levels=c("Bloomer", "Opportunist", "Unresponsive"))
toplot$method <- factor(toplot$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

p3b <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by=c("specialtaxon"), palette=palette_custom,
                 subtitle="True Positives, False Positives and Discordant Associations scaled to the total actual positives (100%)",
                 legend.title="Metric", title="Taxon-taxon associations (special taxa)",
                 xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p3b, filename = "output/taxontaxon/specialtaxon/plot_performance_taxontaxon_special_scaledpositives.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)





#### Evaluate taxon-metadata correlations ####

# define parameters for significance and initialize vectors for the evaluation data
significance_level <- 0.05

true_positive <- c()
false_positive <- c()
true_negative <- c()
false_negative <- c()
method <- c()
spread <- c()
matrixnum <- c()
datatable <- c()
scen <- c()
discordant <- c()

for(file in list.files("data/correlations_taxonmetadata/", recursive = F, pattern = "pvalue", full.names = T)){
    # define file name, method used and reference to be compared against
    filename <- basename(file)
    methodname <- strsplit(filename, split = "_")[[1]][3]
    matrixname <- strsplit(filename, split = "_")[[1]][4]
    spreadname <- strsplit(filename, split = "_")[[1]][5]
    scenarioname <- strsplit(filename, split="_")[[1]][7] %>% gsub(., pattern="\\.tsv", replacement="")
    referencename <- gsub(filename, pattern=methodname, replacement="REAL")
    rhoname <- gsub(filename, pattern="pvalue", replacement="rho")
    referencerhoname <- gsub(referencename, pattern="pvalue", replacement="rho")
    # read files
    test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
    reference_file <- read.table(paste0("data/correlations_taxonmetadata/reference/", referencename),
                                 header=T, stringsAsFactors=F, sep="\t")
    rho_file <- read.table(paste0("data/correlations_taxonmetadata/", rhoname),
                           header=T, stringsAsFactors = F, sep="\t")
    referencerho_file <- read.table(paste0("data/correlations_taxonmetadata/reference/", referencerhoname),
                                    header=T, stringsAsFactors = F, sep = "\t")
    
    
    # all
    test_significant <- c(unlist(test_file < significance_level))
    reference_significant <- c(unlist(reference_file < significance_level))
    test_sign <- sign(c(unlist(rho_file)))
    reference_sign <- sign(c(unlist(referencerho_file)))
    test_sign[is.na(test_sign)] <- 1  # for the joint analyses we consider the "categorical" interactions (with NA sign), to be positive, otherwise they won't count
    reference_sign[is.na(reference_sign)] <- 1
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "all")
    true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
    false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
    true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
    false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
    discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))     
    
    # only positive correlations
    test_significant <- c(unlist(test_file[,1:50] < significance_level))
    reference_significant <- c(unlist(reference_file[,1:50] < significance_level))
    test_sign <- sign(c(unlist(rho_file[,1:50])))
    reference_sign <- sign(c(unlist(referencerho_file[,1:50])))
    
    reference_sign_pos <- reference_sign
    reference_sign_pos[reference_sign<0 | reference_significant==F] <- NA
    test_sign_pos <- test_sign
    test_sign_pos[reference_sign<0 | reference_significant==F] <- NA
    reference_significant_pos <- reference_significant
    reference_significant_pos[reference_sign<0 | reference_significant==F] <- NA
    test_significant_pos <- test_significant
    test_significant_pos[reference_sign<0 | reference_significant==F] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "positive")
    true_positive <- c(true_positive, length(which(test_significant_pos & reference_significant_pos & test_sign_pos==reference_sign_pos)))
    false_positive <- c(false_positive, length(which(test_significant_pos & !reference_significant_pos)))
    true_negative <- c(true_negative, length(which(!test_significant_pos & !reference_significant_pos)))
    false_negative <- c(false_negative, length(which(!test_significant_pos & reference_significant_pos)))
    discordant <- c(discordant, length(which(test_significant_pos & reference_significant_pos & test_sign_pos!=reference_sign_pos)))     
    
    # only negative correlations
    reference_sign_neg <- reference_sign
    reference_sign_neg[reference_sign>0 | reference_significant==F] <- NA
    test_sign_neg <- test_sign
    test_sign_neg[reference_sign>0 | reference_significant==F] <- NA
    reference_significant_neg <- reference_significant
    reference_significant_neg[reference_sign>0 | reference_significant==F] <- NA
    test_significant_neg <- test_significant
    test_significant_neg[reference_sign>0 | reference_significant==F] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "negative")
    true_positive <- c(true_positive, length(which(test_significant_neg & reference_significant_neg & test_sign_neg==reference_sign_neg)))
    false_positive <- c(false_positive, length(which(test_significant_neg & !reference_significant_neg)))
    true_negative <- c(true_negative, length(which(!test_significant_neg & !reference_significant_neg)))
    false_negative <- c(false_negative, length(which(!test_significant_neg & reference_significant_neg)))
    discordant <- c(discordant, length(which(test_significant_neg & reference_significant_neg & test_sign_neg!=reference_sign_neg)))     
    
    # only categorical correlations
    test_significant <- c(unlist(test_file[,51:100] < significance_level))
    reference_significant <- c(unlist(reference_file[,51:100] < significance_level))
    reference_significant_cat <- reference_significant
    reference_significant_cat[reference_significant_cat==F] <- NA
    test_significant_cat <- test_significant
    test_significant_cat[test_significant_cat==F] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "categorical")
    true_positive <- c(true_positive, length(which(test_significant_cat & reference_significant_cat)))
    false_positive <- c(false_positive, length(which(test_significant_cat & !reference_significant_cat)))
    true_negative <- c(true_negative, length(which(!test_significant_cat & !reference_significant_cat)))
    false_negative <- c(false_negative, length(which(!test_significant_cat & reference_significant_cat)))
    discordant <- c(discordant, 0) # no discordant associations with categorical metadata
    
    # only not correlated taxa-metadata pairs
    test_significant <- c(unlist(test_file < significance_level))
    reference_significant <- c(unlist(reference_file < significance_level))
    test_sign <- sign(c(unlist(rho_file)))
    reference_sign <- sign(c(unlist(referencerho_file)))
    test_sign[is.na(test_sign)] <- 1  # for the joint analyses we consider the "categorical" interactions (with NA sign), to be positive, otherwise they won't count
    reference_sign[is.na(reference_sign)] <- 1
    
    reference_sign_no <- reference_sign
    reference_sign_no[reference_significant==T] <- NA
    test_sign_no <- test_sign
    test_sign_no[reference_significant==T] <- NA
    reference_significant_no <- reference_significant
    reference_significant_no[reference_significant==T] <- NA
    test_significant_no <- test_significant
    test_significant_no[reference_significant==T] <- NA
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "non_correlated")
    true_positive <- c(true_positive, length(which(test_significant_no & reference_significant_no & test_sign_no==reference_sign_no)))
    false_positive <- c(false_positive, length(which(test_significant_no & !reference_significant_no)))
    true_negative <- c(true_negative, length(which(!test_significant_no & !reference_significant_no)))
    false_negative <- c(false_negative, length(which(!test_significant_no & reference_significant_no)))
    discordant <- c(discordant, length(which(test_significant_no & reference_significant_no & test_sign_no!=reference_sign_no)))
    
}


# write results table
results <- tibble(method, spread, datatable, matrixnum, scen, true_positive, false_positive, true_negative, false_negative, discordant)
results <- results %>%
    mutate(FDR=100*(false_positive+discordant)/(false_positive+discordant+true_positive)) %>% 
    mutate(Recall=100*true_positive/(true_positive+false_negative+discordant)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive+discordant)) %>% 
    mutate(Accuracy=100*(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative+discordant))%>%
    mutate(TP=100*true_positive/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FP=100*false_positive/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(TN=100*true_negative/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FN=100*false_negative/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(DA=100*discordant/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(false_positive_rate=100*(false_positive+discordant)/(true_negative+false_positive+discordant)) 


write_tsv(results, "output/taxonmetadata/statistics_taxonmetadata_correlation.tsv", col_names = T)
method_type <- tibble(method=c("Seq", "RMP", "Rel", "AST", "CLR", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "VST", "QMP", "ACS"),
                      method_type=c("Sequencing - non-transformed", 
                                    rep("Traditional transformations",  times=2),
                                    rep("Compositional transformations", times=8),
                                    rep("Transformations incorporating microbial loads", times=2)))


results <- results %>% 
    left_join(method_type, by="method")

results$method_type <- factor(results$method_type, 
                              levels=c("Sequencing - non-transformed", "Traditional transformations", 
                                       "Compositional transformations",
                                       "Transformations incorporating microbial loads"))
results$method <- factor(results$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

results$spread <- factor(results$spread, levels=c("low", "high"))
results$scen <- factor(results$scen, levels=c("Healthy", "Blooming","Dysbiosis"))
results_all <- results %>% dplyr::filter(datatable=="all")

results_all$F1 <- 2*((results_all$Precision * results_all$Recall)/(results_all$Precision + results_all$Recall))

# plot results
toplot <- results_all %>% 
    dplyr::select(method, scen, TP, TN, FP, FN, DA, method_type) %>% 
    gather(key="metric", value="value", -method, -scen, -method_type) %>% 
    group_by(method, scen, method_type, metric) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot$metric <- factor(toplot$metric, levels=c("TP", "FP", "TN", "FN", "DA"))

p1 <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by="scen", palette="Dark2",
                legend.title="Metric", title="Taxon-taxon associations",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p1, filename = "output/taxonmetadata/plot_performance_taxonmetadata.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)

# now only scaled to the number of positives in the samples
results_all <- results_all %>% 
    mutate(actual_positives=true_positive+false_negative+discordant) 

toplot2 <- results_all %>% 
    mutate(TP_scaled=100*true_positive/actual_positives) %>% 
    mutate(FP_scaled=100*false_positive/actual_positives) %>% 
    mutate(DA_scaled=100*discordant/actual_positives) %>% 
    dplyr::select(method, scen, method_type, TP_scaled, FP_scaled, DA_scaled) %>% 
    gather(key="metric", value="value", -method, -scen, -method_type) %>% 
    group_by(method, scen, method_type, metric) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot2$metric <- factor(toplot2$metric, levels=c("FP_scaled", "DA_scaled", "TP_scaled"))

p1b <- ggbarplot(toplot2, x="method", y="value", fill="metric", facet.by="scen", palette="Dark2", 
                 subtitle="True Positives, False Positives and Discordant Associations scaled to the total positives (100%)",
                 legend.title="Metric", title="Taxon-metadata associations",
                 xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold")) +
    geom_hline(yintercept = 100, col="black", lwd=1)

ggsave(p1b, filename = "output/taxonmetadata/plot_performance_taxonmetadata_scaledpositives.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)


# associations classified by type (positive, negative or non significant)
results_class <- results %>% dplyr::filter(datatable!="all")

toplot <- results_class %>% 
    dplyr::select(method, scen, TP, TN, FP, FN, DA, method_type, datatable) %>% 
    gather(key="metric", value="value", -method, -scen, -method_type, -datatable) %>% 
    drop_na() %>% 
    group_by(method, scen, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup 

toplot$metric <- factor(toplot$metric, levels=c("TP", "FP", "TN", "FN", "DA"))
toplot$datatable <- factor(toplot$datatable, levels=c("positive", "negative", "categorical", "non_correlated"))

p2 <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by=c("datatable", "scen"), palette="Dark2",
                legend.title="Metric", title="Taxon-metadata associations (classified by true association type)",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p2, filename = "output/taxonmetadata/plot_performance_taxonmetadata_bytype.pdf", device = "pdf", width = 8, height=7, useDingbats=FALSE)


# now only scaled to the number of positives in the samples
actual_positives <- results_all %>% 
    dplyr::select(method, matrixnum, actual_positives)

results_class <- results_class %>% 
    left_join(actual_positives, by=c("method", "matrixnum"))

toplot2 <- results_class %>% 
    mutate(TP_scaled=100*true_positive/actual_positives) %>% 
    mutate(FP_scaled=100*false_positive/actual_positives) %>% 
    mutate(DA_scaled=100*discordant/actual_positives) %>% 
    dplyr::select(method, scen, method_type, TP_scaled, FP_scaled, DA_scaled, datatable) %>%
    drop_na() %>% 
    gather(key="metric", value="value", -method, -scen, -method_type, -datatable) %>% 
    group_by(method, scen, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot2$metric <- factor(toplot2$metric, levels=c("FP_scaled", "DA_scaled", "TP_scaled"))
toplot2$datatable <- factor(toplot2$datatable, levels=c("positive", "negative", "categorical", "non_correlated"))

p2b <- ggbarplot(toplot2, x="method", y="value", fill="metric", facet.by=c("datatable", "scen"), scales="free_y", palette="Dark2", 
                 subtitle="True Positives, False Positives and Discordant Associations scaled to the total actual positives (100%)",
                 legend.title="Metric", title="Taxon-metadata associations (classified by true association type)",
                 xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold")) +
    geom_hline(yintercept = 100, col="black", lwd=1)

ggsave(p2b, filename = "output/taxonmetadata/plot_performance_taxonmetadata_bytype_scaledpositives.pdf", device = "pdf", width = 8, height=7, useDingbats=FALSE)


# calculate sensitivity, precision and FPR 
sensitivity_kruskal <- results_all %>% 
    group_by(scen) %>% 
    kruskal_test(Recall ~ method)

sensitivity_dunn <- results_all %>%
    group_by(scen) %>%
    dunn_test(Recall ~ method, p.adjust.method="BH")

precision_kruskal <- results_all %>% 
    group_by(scen) %>% 
    kruskal_test(Precision ~ method)

precision_dunn <- results_all %>%
    group_by(scen) %>%
    dunn_test(Precision ~ method, p.adjust.method="BH")

fpr_kruskal <- results_all %>% 
    group_by(scen) %>% 
    kruskal_test(false_positive_rate ~ method)

fpr_dunn <- results_all %>%
    group_by(scen) %>%
    dunn_test(false_positive_rate ~ method, p.adjust.method="BH")

# plot
p1 <- plot_comparisons(dunn=sensitivity_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=sensitivity_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=sensitivity_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

psensitivity <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")

p1 <- plot_comparisons(dunn=precision_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=precision_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=precision_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

pprecision <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")

p1 <- plot_comparisons(dunn=fpr_dunn %>% filter(scen=="Healthy"), title="Healthy succession")
p2 <- plot_comparisons(dunn=fpr_dunn %>% filter(scen=="Blooming"), title="Blooming")
p3 <- plot_comparisons(dunn=fpr_dunn %>% filter(scen=="Dysbiosis"), title="Dysbiosis")

pfpr <- ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="top")

all_res <- ggarrange(psensitivity, pprecision, pfpr, nrow=3, ncol=1)
ggsave(all_res, filename = "output/taxonmetadata/dunn_tests_plots.pdf", device="pdf", height=11, width=9)
# save statistics 
kruskal_table <- bind_rows(sensitivity_kruskal, precision_kruskal, fpr_kruskal)
dunn_table <- bind_rows(sensitivity_dunn, precision_dunn, fpr_dunn)

dunn_table <- dunn_table %>% 
    dplyr::filter(!(.y.=="false_positive_rate" & scen=="Healthy")) # remove non-significant data from kruskal - only kept for plotting 
write_tsv(kruskal_table, path="output/taxonmetadata/kruskal_wallis_method_comparison.txt", col_names=T)
write_tsv(dunn_table, path="output/taxonmetadata/dunn_test_method_comparison.txt", col_names=T)








#### Evaluate taxon-metadata correlations only in special taxa ####

# define parameters for significance and initialize vectors for the evaluation data
significance_level <- 0.05

true_positive <- c()
false_positive <- c()
true_negative <- c()
discordant <- c()
false_negative <- c()
method <- c()
spread <- c()
matrixnum <- c()
datatable <- c()
scen <- c()
specialtaxon <- c()
matrix_stats <- read_tsv("data/raw/matrix_stats_3scenarios.tsv")

for(file in list.files("data/correlations_taxonmetadata/", recursive = F, pattern = "pvalue", full.names = T)){
    # define file name, method used and reference to be compared against
    filename <- basename(file)
    methodname <- strsplit(filename, split = "_")[[1]][3]
    matrixname <- strsplit(filename, split = "_")[[1]][4]
    spreadname <- strsplit(filename, split = "_")[[1]][5]
    scenarioname <- strsplit(filename, split="_")[[1]][7] %>% gsub(., pattern="\\.tsv", replacement="")
    referencename <- gsub(filename, pattern=methodname, replacement="REAL")
    rhoname <- gsub(filename, pattern="pvalue", replacement="rho")
    referencerhoname <- gsub(referencename, pattern="pvalue", replacement="rho")
    if(scenarioname=="Blooming"){
        # read files
        test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
        reference_file <- read.table(paste0("data/correlations_taxonmetadata/reference/", referencename),
                                     header=T, stringsAsFactors=F, sep="\t")
        rho_file <- read.table(paste0("data/correlations_taxonmetadata/", rhoname),
                               header=T, stringsAsFactors = F, sep="\t")
        referencerho_file <- read.table(paste0("data/correlations_taxonmetadata/reference/", referencerhoname),
                                        header=T, stringsAsFactors = F, sep = "\t")
        
        # special taxon (bloomer)
        sptax <- matrix_stats %>% dplyr::filter(matrix==matrixname) %>% pull(special_taxon)
        test_file <- test_file[sptax,, drop=F]
        reference_file <- reference_file[sptax,, drop=F]
        rho_file <- rho_file[sptax,, drop=F]
        referencerho_file <- referencerho_file[sptax,, drop=F]
        
        # all
        test_significant <- c(unlist(test_file < significance_level))
        reference_significant <- c(unlist(reference_file < significance_level))
        test_sign <- sign(c(unlist(rho_file)))
        reference_sign <- sign(c(unlist(referencerho_file)))
        test_sign[is.na(test_sign)] <- 1  # for the joint analyses we consider the "categorical" interactions (with NA sign), to be positive, otherwise they won't count
        reference_sign[is.na(reference_sign)] <- 1
        
        # populate evaluation vectors
        method <- c(method, methodname)
        spread <- c(spread, spreadname)
        matrixnum <- c(matrixnum, matrixname)
        scen <- c(scen, scenarioname)
        datatable <- c(datatable, "all")
        true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
        false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
        true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
        false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
        discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
        specialtaxon <- c(specialtaxon, "Bloomer")
        
    }
    if(scenarioname=="Dysbiosis"){
        # read files
        test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
        reference_file <- read.table(paste0("data/correlations_taxonmetadata/reference/", referencename),
                                     header=T, stringsAsFactors=F, sep="\t")
        rho_file <- read.table(paste0("data/correlations_taxonmetadata/", rhoname),
                               header=T, stringsAsFactors = F, sep="\t")
        referencerho_file <- read.table(paste0("data/correlations_taxonmetadata/reference/", referencerhoname),
                                        header=T, stringsAsFactors = F, sep = "\t")
        
        # special taxon (opportunist)
        sptax <- matrix_stats %>% dplyr::filter(matrix==matrixname) %>% pull(special_taxon)
        test_file <- test_file[sptax,, drop=F]
        reference_file <- reference_file[sptax,, drop=F]
        rho_file <- rho_file[sptax,, drop=F]
        referencerho_file <- referencerho_file[sptax,, drop=F]
        
        # all
        test_significant <- c(unlist(test_file < significance_level))
        reference_significant <- c(unlist(reference_file < significance_level))
        test_sign <- sign(c(unlist(rho_file)))
        reference_sign <- sign(c(unlist(referencerho_file)))
        test_sign[is.na(test_sign)] <- 1  # for the joint analyses we consider the "categorical" interactions (with NA sign), to be positive, otherwise they won't count
        reference_sign[is.na(reference_sign)] <- 1
        
        # populate evaluation vectors
        method <- c(method, methodname)
        spread <- c(spread, spreadname)
        matrixnum <- c(matrixnum, matrixname)
        scen <- c(scen, scenarioname)
        datatable <- c(datatable, "all")
        true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
        false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
        true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
        false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
        discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
        specialtaxon <- c(specialtaxon, "Opportunist")
        
        # read files
        test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
        reference_file <- read.table(paste0("data/correlations_taxonmetadata/reference/", referencename),
                                     header=T, stringsAsFactors=F, sep="\t")
        rho_file <- read.table(paste0("data/correlations_taxonmetadata/", rhoname),
                               header=T, stringsAsFactors = F, sep="\t")
        referencerho_file <- read.table(paste0("data/correlations_taxonmetadata/reference/", referencerhoname),
                                        header=T, stringsAsFactors = F, sep = "\t")
        
        # special taxon (unresponsive)
        sptax <- matrix_stats %>% dplyr::filter(matrix==matrixname) %>% pull(flat_taxon_dysbiosis)
        sptax <- strsplit(sptax, split=",")[[1]]
        sptax <- intersect(sptax, rownames(test_file))
        test_file <- test_file[sptax,, drop=F]
        reference_file <- reference_file[sptax,, drop=F]
        rho_file <- rho_file[sptax,, drop=F]
        referencerho_file <- referencerho_file[sptax,, drop=F]
        
        # all
        test_significant <- c(unlist(test_file < significance_level))
        reference_significant <- c(unlist(reference_file < significance_level))
        test_sign <- sign(c(unlist(rho_file)))
        reference_sign <- sign(c(unlist(referencerho_file)))
        test_sign[is.na(test_sign)] <- 1  # for the joint analyses we consider the "categorical" interactions (with NA sign), to be positive, otherwise they won't count
        reference_sign[is.na(reference_sign)] <- 1
        
        # populate evaluation vectors
        method <- c(method, methodname)
        spread <- c(spread, spreadname)
        matrixnum <- c(matrixnum, matrixname)
        scen <- c(scen, scenarioname)
        datatable <- c(datatable, "all")
        true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
        false_positive <- c(false_positive, length(which(test_significant & !reference_significant)))
        true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
        false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
        discordant <- c(discordant, length(which(test_significant & reference_significant & test_sign!=reference_sign)))
        specialtaxon <- c(specialtaxon, "Unresponsive")
    }
}


# write results table
results <- tibble(method, spread, datatable, matrixnum, scen, true_positive, false_positive, true_negative, false_negative, discordant, specialtaxon)
results <- results %>%
    mutate(FDR=100*(false_positive+discordant)/(false_positive+discordant+true_positive)) %>% 
    mutate(Recall=100*true_positive/(true_positive+false_negative+discordant)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive+discordant)) %>% 
    mutate(Accuracy=100*(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative+discordant))%>%
    mutate(TP=100*true_positive/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FP=100*false_positive/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(TN=100*true_negative/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(FN=100*false_negative/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(DA=100*discordant/(true_positive+true_negative+false_positive+false_negative+discordant)) %>% 
    mutate(false_positive_rate=100*(false_positive+discordant)/(true_negative+false_positive+discordant)) 


write_tsv(results, "output/taxonmetadata/specialtaxon/statistics_taxonmetadata_correlation_onlyspecialtaxa.tsv", col_names = T)
method_type <- tibble(method=c("Seq", "RMP", "Rel", "AST", "CLR", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "VST", "QMP", "ACS"),
                      method_type=c("Sequencing - non-transformed", 
                                    rep("Relative transformations",  times=2),
                                    rep("Compositional transformations", times=8),
                                    rep("Quantitative transformations", times=2)))


results <- results %>% 
    left_join(method_type, by="method")

results$method_type <- factor(results$method_type, 
                              levels=c("Sequencing - non-transformed", "Relative transformations", 
                                       "Compositional transformations",
                                       "Quantitative transformations"))
results$method <- factor(results$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

results$spread <- factor(results$spread, levels=c("low", "medium", "high"))
results$scen <- factor(results$scen, levels=c("Healthy", "Blooming","Dysbiosis"))
results_all <- results %>% dplyr::filter(datatable=="all")

results_all$F1 <- 2*((results_all$Precision * results_all$Recall)/(results_all$Precision + results_all$Recall))

# plot results
toplot <- results_all %>% 
    dplyr::select(method, specialtaxon, TP, TN, FP, FN, DA, method_type, datatable) %>% 
    gather(key="metric", value="value", -method, -specialtaxon, -method_type, -datatable) %>% 
    drop_na() %>% 
    group_by(method, specialtaxon, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup 

toplot$metric <- factor(toplot$metric, levels=c("TP", "FP", "TN", "FN", "DA"))
toplot$specialtaxon <- factor(toplot$specialtaxon, levels=c("Bloomer", "Opportunist", "Unresponsive"))

p3 <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by=c("specialtaxon"), palette="Dark2",
                legend.title="Metric", title="Taxon-metadata associations (special taxa)",
                xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p3, filename = "output/taxonmetadata/specialtaxon/plot_performance_taxonmetadata_special.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)



# now only scaled to the number of positives in the samples
stats_matrix <- read_tsv("data/raw/matrix_stats_3scenarios.tsv", col_names=T)
actual_positives <- read_tsv("output/taxonmetadata/statistics_taxonmetadata_correlation.tsv") %>% 
    dplyr::filter(datatable=="all") %>% 
    mutate(actual_positives=true_positive+false_negative+discordant) %>% 
    dplyr::select(method, matrixnum, actual_positives) %>% 
    distinct()

results_all <- results_all %>% 
    left_join(actual_positives, by=c("method", "matrixnum"))

toplot <- results_all %>% 
    left_join(stats_matrix, by=c("matrixnum" = "matrix", "scen" = "scenario", "specialtaxon" = "special_taxon")) %>% 
    mutate(TP_scaled=100*true_positive/(actual_positives)) %>% 
    mutate(FP_scaled=100*false_positive/(actual_positives)) %>% 
    mutate(DA_scaled=100*discordant/(actual_positives)) %>% 
    dplyr::select(method, specialtaxon, method_type, TP_scaled, FP_scaled, DA_scaled, datatable) %>%
    drop_na() %>% 
    gather(key="metric", value="value", -method, -specialtaxon, -method_type, -datatable) %>% 
    group_by(method, specialtaxon, method_type, metric, datatable) %>% 
    summarize(value=mean(value)) %>% 
    ungroup %>% 
    drop_na()

toplot$metric <- factor(toplot$metric, levels=c("FP_scaled", "DA_scaled", "TP_scaled"))
toplot$specialtaxon <- factor(toplot$specialtaxon, levels=c("Bloomer", "Opportunist", "Unresponsive"))
toplot$method <- factor(toplot$method, levels=c("Seq", "RMP", "Rel","AST", "CLR", "CSS", "GMPR",
                                                "UQ", "RLE", "TMM", "VST", "QMP", "ACS"))

p3b <- ggbarplot(toplot, x="method", y="value", fill="metric", facet.by=c("specialtaxon"), palette="Dark2",
                 subtitle="True Positives, False Positives and Discordant Associations, scaled to the number of true associations",
                 legend.title="Metric", title="Taxon-metadata associations (special taxa)",
                 xlab="Method", ylab="Percentage") +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(p3b, filename = "output/taxonmetadata/specialtaxon/plot_performance_taxonmetadata_special_scaledpositives.pdf", device = "pdf", width = 8, height=3.5, useDingbats=FALSE)

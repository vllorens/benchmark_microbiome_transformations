# Taxon-taxon and taxon-metadata correlations
# Mon Dec  2 15:01:23 2019 ------------------------------

# We will use different data and transformation methods to calculate all the correlations. 
# For comparisons, we will only use Spearman correlations for taxon-taxon and either Spearman correlation or Kruskal-wallis test for taxon-metadata associations 
# (Spearman for numeric and Kruskal-Wallis for factors in the metadata)

# The following methods will be tested:
# * REAL: Original data (as reference). Not transformed.
# * CLR: Sequencing data. CLR transformed.
# * CSS: Sequencing data. Normalized using cumulative sum scaling (MetagenomeSeq implemetnation).
# * AST: Sequencing data. Normalized using total sum scaling + arcsine squared-root transformation. 
# * GMPR: Sequencing data. Size factors for normalization are calculated using the Geometric Mean of Pairwise Ratios.
# * TMM: Sequencing data. Normalization using trimmed mean of M-values (EdgeR implementation).
# * RLE: Sequencing data. Normalization using relative log-expression (EdgeR implementation).
# * UQ: Sequencing data. Normalization using upper quartile (EdgeR implementation).
# * VST: Sequencing data. Using variance stabilizing transformation from DESeq2.
# * RMP: RMP data. Rarefied (already done above).
# * QMP: QMP data. Normalized by rarefying to even sampling depth and scaling to total counts (already done above). 
# * QMP-NR: QMP non-rarefied data. Normalized by scaling to total counts (already done above). 
# 
# In all of these methods, we will only evaluate those taxa with prevalence > 50%. We will use the sequencing matrix to determine this threshold on the zeros. 
# 
# The correlations that we will calculate, for each type of transformation, are:
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
suppressPackageStartupMessages(library(dabestr))
suppressPackageStartupMessages(library(gdata))

# load functions
source("R/estimateZeros.R")
source("R/taxaToKeep.R")
source("R/taxon_correlation.R")
source("R/taxon_metadata_correlation.R")
source("R/GMPR.R")

# If it has not been run yet, it is necessary to run the script to generate the data (scripts/generate_data.R)
# source("scripts/generate_data.R")

# create folders, if not existing
system("mkdir -p data/correlations_taxontaxon")
system("mkdir -p data/correlations_taxonmetadata")
system("mkdir -p output/taxontaxon")
system("mkdir -p output/taxonmetadata")

# clean these folders and remove all previous results (if any) 
system("rm -r data/correlations_taxontaxon/*")
system("rm -r data/correlations_taxonmetadata/*")
system("rm -r output/taxontaxon/*")
system("rm -r output/taxonmetadata/*")

# create additional folders for the reference data (on the real matrices)
system("mkdir -p data/correlations_taxontaxon/reference")
system("mkdir -p data/correlations_taxonmetadata/reference")

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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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

}

# CLR
for(file in list.files("data/seq_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # transform data
    tax_matrix <- t(zCompositions::cmultRepl(X = t(tax_matrix), output="counts")) # estimates zeros
    tax_matrix <- t(codaSeq.clr(tax_matrix, samples.by.row=F))
    # read metadata file
    filename_metadata <- gsub(filename, pattern="seqOut_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix)
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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
}


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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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

}


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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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
}


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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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
}


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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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
}


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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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
}

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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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
}


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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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
}


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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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
}


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
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
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
}

# QMP2
for(file in list.files("data/qmp2_matrices", full.names = T)){
    # read file and select only those taxa to keep
    filename <- basename(file)
    taxa_to_keep <- taxaToKeep(filename)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    tax_matrix <- tax_matrix[taxa_to_keep,]
    # read metadata file
    filename_metadata <- gsub(filename, pattern="QMP2_", replacement="")
    metadata_matrix <- read.table(paste0("data/metadata_matrices/metadata_", filename_metadata))
    # calculate correlations and pvalues (taxon-taxon)
    taxon_correlation_object <- taxon_correlation(tax_matrix) 
    taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])] <- p.adjust(taxon_correlation_object[[2]][upper.tri(taxon_correlation_object[[2]])], method="BH")
    # calculate correlations and pvalues (taxon-metadata)
    metadata_correlation_object <- taxon_metadata_correlation(tax_matrix, metadata_matrix) 
    metadata_correlation_object[[2]] <- matrix(p.adjust(as.vector(as.matrix(metadata_correlation_object[[2]])), method='fdr'),ncol=ncol(metadata_correlation_object[[2]]))
    colnames(metadata_correlation_object[[2]]) <- colnames(metadata_correlation_object[[1]])
    rownames(metadata_correlation_object[[2]]) <- rownames(metadata_correlation_object[[1]])
    # write output (taxon-taxon)
    outputname_correlation <- gsub(filename, pattern="QMP2_taxonomy", replacement="taxontaxon_rho_QMP-NR")
    outputname_pval <- gsub(filename, pattern="QMP2_taxonomy", replacement="taxontaxon_pvalue_QMP-NR")
    write.table(taxon_correlation_object[[1]], 
                paste0("data/correlations_taxontaxon/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(taxon_correlation_object[[2]], 
                paste0("data/correlations_taxontaxon/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output (taxon-metadata)
    outputname_correlation <- gsub(filename, pattern="QMP2_taxonomy", replacement="taxonmetadata_rho_QMP-NR")
    outputname_pval <- gsub(filename, pattern="QMP2_taxonomy", replacement="taxonmetadata_pvalue_QMP-NR")
    write.table(metadata_correlation_object[[1]], 
                paste0("data/correlations_taxonmetadata/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(metadata_correlation_object[[2]], 
                paste0("data/correlations_taxonmetadata/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
}


#### Evaluate taxon-taxon correlations (Figure 2a) ####

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
    
    lowerTriangle(test_file) <- upperTriangle(test_file, byrow = T)
    lowerTriangle(reference_file) <- upperTriangle(reference_file, byrow = T)
    lowerTriangle(rho_file) <- upperTriangle(rho_file, byrow = T)
    lowerTriangle(referencerho_file) <- upperTriangle(referencerho_file, byrow = T)
    
    
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
    false_positive <- c(false_positive, length(which(test_significant & !reference_significant))+length(which(test_significant & reference_significant & test_sign!=reference_sign)))
    true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
    false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
}

# write results table
results <- tibble(method, spread, datatable, matrixnum, scen, true_positive, false_positive, true_negative, false_negative)
results <- results %>%
    mutate(FDR=100*false_positive/(false_positive+true_positive)) %>% 
    mutate(Recall=100*true_positive/(true_positive+false_negative)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive)) %>% 
    mutate(Accuracy=100*(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative))%>%
    mutate(true_positive_percent=100*true_positive/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_positive_percent=100*false_positive/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(true_negative_percent=100*true_negative/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_negative_percent=100*false_negative/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_positive_rate=100*false_positive/(true_negative+false_positive)) %>% 
    dplyr::select(-c(true_positive,true_negative,false_positive,false_negative)) %>% 
    drop_na()


write_tsv(results, "output/taxontaxon/statistics_taxontaxon_correlation.tsv", col_names = T)
method_type <- tibble(method=c("RMP", "AST", "CLR", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "VST", "QMP", "QMP-NR"),
                      method_type=c("Traditional transformations", 
                                    rep("Compositional transformations", times=8),
                                    rep("Quantitative transformations", times=2)))

results <- results %>% 
    left_join(method_type, by="method")

results$method_type <- factor(results$method_type, 
                                  levels=c("Sequencing", "Traditional transformations", 
                                           "Compositional transformations",
                                           "Quantitative transformations"))
results$method <- factor(results$method, levels=c("RMP", "AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "QMP-NR"))

results$spread <- factor(results$spread, levels=c("low", "high"))
results$scen <- factor(results$scen, levels=c("Healthy", "Dysbiosis", "Blooming"))
results_all <- results %>% dplyr::filter(datatable=="all")


# plot results
custom_palette=get_palette("Spectral",11)[c(2,5,10)]
ggboxplot(results_all, x="method", y="Precision", fill="method_type", alpha=0.5,
          facet.by = c( "scen"),
          ylab="Precision [TP/TP+FP]",
          palette=custom_palette,legend.title="Method type",
          main="Taxon-taxon correlations") + 
    ggsignif::stat_signif(comparisons = list(c("QMP", "QMP-NR")), test = "t.test", FDR = T) +
    geom_point(aes(fill=method_type), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(filename = "output/taxontaxon/plot_taxontaxon_precision.pdf", device = "pdf", width = 11, height=5)
ggboxplot(results_all, x="method", y="Recall", alpha=0.5,
          fill="method_type", facet.by = c("scen"), 
          ylab="Recall [TP/TP+FN]",
          palette=custom_palette,legend.title="Method type",
          main="Taxon-taxon correlations") + 
    ggsignif::stat_signif(comparisons = list(c("QMP", "QMP-NR")), test = "t.test", FDR = T) +
    geom_point(aes(fill=method_type), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(filename = "output/taxontaxon/plot_taxontaxon_recall.pdf", device = "pdf", width = 11, height=5)
ggboxplot(results_all, x="method", y="false_positive_percent", alpha=0.5,
          fill="method_type", facet.by = c( "scen"), 
          ylab="% False positives",
          palette=custom_palette,legend.title="Method type",
          main="Taxon-taxon correlations") + 
    ggsignif::stat_signif(comparisons = list(c("QMP", "QMP-NR")), test = "t.test", FDR = T) +
    geom_point(aes(fill=method_type), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(filename = "output/taxontaxon/plot_taxontaxon_FP.pdf", device = "pdf", width = 11, height=5)




#### Evaluate taxon-metadata correlations (Figure 2b) ####

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
    
   
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    matrixnum <- c(matrixnum, matrixname)
    scen <- c(scen, scenarioname)
    datatable <- c(datatable, "all")
    true_positive <- c(true_positive, length(which(test_significant & reference_significant & test_sign==reference_sign)))
    false_positive <- c(false_positive, length(which(test_significant & !reference_significant))+length(which(test_significant & reference_significant & test_sign!=reference_sign)))
    true_negative <- c(true_negative, length(which(!test_significant & !reference_significant)))
    false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
    
}

# write results table
results <- tibble(method, spread, scen, matrixnum, datatable, true_positive, false_positive, true_negative, false_negative)
results <- results %>%
    mutate(FDR=100*false_positive/(false_positive+true_positive)) %>% 
    mutate(Recall=100*true_positive/(true_positive+false_negative)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive)) %>% 
    mutate(Accuracy=100*(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(true_positive_percent=100*true_positive/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_positive_percent=100*false_positive/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(true_negative_percent=100*true_negative/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_negative_percent=100*false_negative/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_positive_rate=100*false_positive/(true_negative+false_positive)) %>% 
    dplyr::select(-c(true_positive,true_negative,false_positive,false_negative))

write_tsv(results, "output/taxonmetadata/statistics_taxonmetadata_correlation.tsv", col_names = T)
method_type <- tibble(method=c("RMP", "AST", "CLR", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "VST", "QMP", "QMP-NR"),
                      method_type=c("Traditional transformations", 
                                    rep("Compositional transformations", times=8),
                                    rep("Quantitative transformations", times=2)))

results <- results %>% 
    left_join(method_type, by="method")

results$method_type <- factor(results$method_type, 
                              levels=c("Sequencing", "Traditional transformations", 
                                       "Compositional transformations",
                                       "Quantitative transformations"))
results$method <- factor(results$method, levels=c("RMP", "AST", "CLR", "CSS", "GMPR",
                                                  "UQ", "RLE", "TMM", "VST", "QMP", "QMP-NR"))

results$spread <- factor(results$spread, levels=c("low", "high"))
results$scen <- factor(results$scen, levels=c("Healthy", "Dysbiosis", "Blooming"))
results_all <- results %>% dplyr::filter(datatable=="all")


custom_palette=get_palette("Spectral",11)[c(2,5,10)]
# all data
ggboxplot(results_all, x="method", y="Precision", alpha=0.5,
          fill="method_type", facet.by = c( "scen"), 
          ylab="Precision [TP/TP+FP]",
          palette=custom_palette, legend.title="Method type",
          main="Taxon-metadata correlations") + 
    ggsignif::stat_signif(comparisons = list(c("QMP", "QMP-NR")), test = "t.test", FDR = T) +
    geom_point(aes(fill=method_type), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(filename = "output/taxonmetadata/plot_taxonmetadata_precision.svg", device = "svg", width = 11, height=5)
ggboxplot(results_all, x="method", y="Recall", alpha=0.5,
          fill="method_type", facet.by = c("scen"), 
          ylab="Recall [TP/TP+FN]",
          palette=custom_palette,legend.title="Method type",
          main="Taxon-metadata correlations") + 
    ggsignif::stat_signif(comparisons = list(c("QMP", "QMP-NR")), test = "t.test", FDR = T) +
    geom_point(aes(fill=method_type), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(filename = "output/taxonmetadata/plot_taxonmetadata_recall.svg", device = "svg", width = 11, height=5)
ggboxplot(results_all, x="method", y="false_positive_percent", alpha=0.5,
          fill="method_type", facet.by = c( "scen"), 
          ylab="% False positives",
          palette=custom_palette,legend.title="Method type",
          main="Taxon-metadata correlations") + 
    ggsignif::stat_signif(comparisons = list(c("QMP", "QMP-NR")), test = "t.test", FDR = T) +
    geom_point(aes(fill=method_type), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(filename = "output/taxonmetadata/plot_taxonmetadata_FP.svg", device = "svg", width = 11, height=5)

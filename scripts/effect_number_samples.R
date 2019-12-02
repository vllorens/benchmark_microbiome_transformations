# Test of the effect of the number of samples in the community
# Mon Dec  2 15:58:34 2019 ------------------------------


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
source("R/taxon_correlation.R")
source("R/taxon_metadata_correlation.R")
source("R/GMPR.R")
source("R/taxaToKeep.R")
source("R/metadata.R")
source("R/rarefy_even_sampling_depth.R")

# create matrices (if not there)
system("mkdir -p output/number_samples/")

# remove all previous results (if any)
system("rm -rf output/number_samples/*")



#### Loop to test effect of number of samples ####

# this script simply parses the original file to obtain separated matrices with the original simulated data
set.seed(14102019)

for(numberofsamplestomake in c(20,50,100,200,500,1000)){
    # create matrices (if not there)
    system("mkdir -p data/tax_matrices")
    system("mkdir -p data/seq_matrices")
    system("mkdir -p data/rmp_matrices")
    system("mkdir -p data/qmp_matrices")
    system("mkdir -p ata/qmp2_matrices")
    system("mkdir -p data/metadata_matrices")
    system("mkdir -p data/correlations_taxontaxon")
    system("mkdir -p data/correlations_taxonmetadata")
    
    # remove all previous results (if any)
    system("rm -rf data/tax_matrices/*")
    system("rm -rf data/seq_matrices/*")
    system("rm -rf data/rmp_matrices/*")
    system("rm -rf data/qmp_matrices/*")
    system("rm -rf data/qmp2_matrices/*")
    system("rm -rf data/metadata_matrices/*")
    system("rm -rf data/correlations_taxontaxon/*")
    system("rm -rf data/correlations_taxonmetadata/*")

    # create additional folders for the reference data (on the real matrices)
    system("mkdir -p data/correlations_taxontaxon/reference")
    system("mkdir -p data/correlations_taxonmetadata/reference")
    
    # Load simulated matrices and subsample N random samples from each (N=numberofsamplestomake)
    load("data/raw/20191122_sims_Pedro_v2.4.Rdata")
    
    for(tab in ls(pattern = "Mp")){
        ind_matrix <- get(tab) %>% 
            t() 
        mt <-strsplit(tab, split="_")[[1]][1]
        if(length(grep(mt, pattern="B"))>0){
            scenario <- "Blooming"
        } else if(length(grep(mt, pattern="S"))>0){
            scenario <- "Healthy"
        } else if(length(grep(mt, pattern="D"))>0){
            scenario <- "Dysbiosis"
        }
        ind_matrix <- ind_matrix %>% 
            as.data.frame()
        ind_matrix <- ind_matrix[,sample(colnames(ind_matrix), replace = F, size = numberofsamplestomake)]
        spread <- (apply(ind_matrix,2,sum) %>% max)/(apply(ind_matrix,2,sum) %>% min) 
        if(spread<10){
            spread_cat <- "low"
        } else{
            spread_cat <- "high"
        }
        write.table(ind_matrix, paste0("data/tax_matrices/taxonomy_", mt, "_", spread_cat, "_spread_", scenario, ".tsv"), 
                    col.names = T, row.names=T, quote=F, sep="\t")
    }
    
    # Generate all matrices
    # the sequencing data comes from the real simulated matrix
    for(file in list.files("data/tax_matrices", full.names = T)){
        filename <- basename(file) %>% 
            gsub(x=., pattern="taxonomy", replacement="seqOut_taxonomy")
        tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t")
        seqOut <- data.frame(row.names=rownames(tax_matrix))
        for(col in 1:ncol(tax_matrix)){
            tvect <- sample(rownames(tax_matrix),size=rlnorm(1, meanlog = 10.3, sdlog = 0.3), 
                            prob=tax_matrix[,col], replace=T) %>% table
            seqOut <- cbind(seqOut, tvect[rownames(seqOut)])
        }
        seqOut <- seqOut[,c(F,T)]
        seqOut[is.na(seqOut)] <- 0
        colnames(seqOut) <- colnames(tax_matrix)
        rownames(seqOut) <- rownames(tax_matrix)
        write.table(seqOut, file = paste0("data/seq_matrices/", filename), 
                    quote = FALSE, row.names = TRUE, col.names=T, sep="\t")
    }
    
    # the RMP data comes from the sequencing output
    for(file in list.files("data/seq_matrices", full.names = T)){
        filename <- basename(file) %>% 
            gsub(x=., pattern="seqOut_taxonomy", replacement="RMP_taxonomy")
        seqOut <-  read.table(file, header=T, stringsAsFactors=F, sep="\t")
        RMP = otu_table(seqOut, taxa_are_rows = TRUE)
        RMP = rarefy_even_depth(RMP, trimOTUs = FALSE, 
                                sample.size = min(sample_sums(RMP)), verbose = F) # rarefy
        RMP = as.matrix(RMP) # convert back to matrix
        write.table(RMP, file = paste0("data/rmp_matrices/", filename), 
                    quote = FALSE, row.names = TRUE, col.names=T, sep="\t")
    }
    
    # the QMP data comes from multiplying rmp data by a factor to get numbers proportional to the number of counts
    # i.e. first rarefying to even sampling depth
    for(file in list.files("data/seq_matrices", full.names = T)){
        filename <- basename(file) %>% 
            gsub(x=., pattern="seqOut_taxonomy", replacement="QMP_taxonomy")
        original_file <- paste0("data/tax_matrices/", gsub(filename, pattern="QMP_", replacement=""))
        counts <- read.table(original_file, header=T, stringsAsFactors=F, sep="\t") %>% 
            apply(., 2, sum) %>% 
            as.data.frame()
        seq <- read.table(file, header=T, stringsAsFactors=F, sep="\t")
        QMP <- rarefy_even_sampling_depth(cnv_corrected_abundance_table = seq, cell_counts_table = counts)
        write.table(QMP, file = paste0("data/qmp_matrices/", filename), 
                    quote = FALSE, row.names = TRUE, col.names=T, sep="\t")
    }
    
    # the QMP-NR data comes from multiplying sequencing data by a factor to get numbers proportional to the number of counts
    # i.e. without first rarefying to even sampling depth
    for(file in list.files("data/seq_matrices", full.names = T)){
        filename <- basename(file) %>% 
            gsub(x=., pattern="seqOut_taxonomy", replacement="QMP2_taxonomy")
        original_file <- paste0("data/tax_matrices/", gsub(filename, pattern="QMP2_", replacement=""))
        counts <- read.table(original_file, header=T, stringsAsFactors=F, sep="\t") %>% 
            apply(., 2, sum)
        seq <- read.table(file, header=T, stringsAsFactors=F, sep="\t")
        counts_seq <- apply(seq, 2, sum)
        factors <- counts/counts_seq
        QMP2 <- sweep(seq, MARGIN = 2, factors, '*')
        write.table(QMP2, file = paste0("data/qmp2_matrices/", filename), 
                    quote = FALSE, row.names = TRUE, col.names=T, sep="\t")
    }
    
    # First, define parameters for the metadata
    rho_meta_min <- 0.35       # minimum correlation for the metadata variables
    rho_meta_max <- 0.6        # max correlation for the metadata variables
    num_metadata_noC <- 35     # metadata numeric variables not correlating
    num_metadata_Ct <- 10      # metadata numeric variables correlating w one taxon
    num_metadata_Cc <- 5       # metadata numeric variables correlating with total counts
    cat_metadata_noC <- 35     # metadata categorical variables not correlating
    cat_metadata_Ct <- 10      # metadata categorical variables correlating w one taxon
    cat_metadata_Cc <- 5       # metadata categorical variables correlating with total counts
    
    for(file in list.files("data/tax_matrices", full.names = T)){
        filename <- basename(file)
        taxonomy <- read.table(file, header = T, stringsAsFactors=F, sep="\t")
        metaD <- generate_metadata(taxonomic_matrix = taxonomy, num_noc = num_metadata_noC,
                                   num_ct = num_metadata_Ct, num_cc = num_metadata_Cc,
                                   cat_noc = cat_metadata_noC, cat_ct = cat_metadata_Ct,
                                   cat_cc = cat_metadata_Cc, corr_min = rho_meta_min, 
                                   corr_max = rho_meta_max)
        write.table(metaD, file = paste0("data/metadata_matrices/metadata_", filename),
                    col.names=paste0("variable_", 1:ncol(metaD)), 
                    row.names=paste0("sample_", 1:nrow(metaD)))
    }
    
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
        tax_matrix <- t(zCompositions::cmultRepl(X = t(tax_matrix), output="counts"))
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
    
    # QMP-NR
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
    
    
    # define parameters for significance and initialize vectors for the evaluation data (taxon-taxon)
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
        dplyr::select(-c(true_positive,true_negative,false_positive,false_negative)) %>% 
        drop_na()
    
    
    write_tsv(results, paste0("output/number_samples/statistics_taxontaxon_correlation_", numberofsamplestomake, ".tsv"), col_names = T)
    
    
    # define parameters for significance and initialize vectors for the evaluation data (taxon-metadata)
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
        dplyr::select(-c(true_positive,true_negative,false_positive,false_negative))
    
    write_tsv(results, paste0("output/number_samples/statistics_taxonmetadata_correlation_", numberofsamplestomake, ".tsv"), col_names = T)
}


#### Assess performance of the methods across different sample numbers (taxon-taxon) ####

r20 <- read_tsv("output/number_samples/statistics_taxontaxon_correlation_20.tsv")
r50 <- read_tsv("output/number_samples/statistics_taxontaxon_correlation_50.tsv")
r100 <- read_tsv("output/number_samples/statistics_taxontaxon_correlation_100.tsv")
r200 <- read_tsv("output/number_samples/statistics_taxontaxon_correlation_200.tsv")
r500 <- read_tsv("output/number_samples/statistics_taxontaxon_correlation_500.tsv")
r1000 <- read_tsv("output/number_samples/statistics_taxontaxon_correlation_1000.tsv")
r20 <- r20 %>% mutate(num_samples="20")
r50 <- r50 %>% mutate(num_samples="50")
r100 <- r100 %>% mutate(num_samples="100")
r200 <- r200 %>% mutate(num_samples="200")
r500 <- r500 %>% mutate(num_samples="500")
r1000 <- r1000 %>% mutate(num_samples="1000")
rr <- bind_rows(r20,r50,r100,r200,r500,r1000)
rr$method <- factor(rr$method, levels=c("AST", "CLR", "RMP", "CSS", "GMPR",
                                        "RLE", "TMM", "UQ", "VST",
                                        "QMP", "QMP-NR"))
rr$num_samples <- factor(rr$num_samples, levels=c("20", "50", "100", 
                                                  "200", "500", "1000"))
rr$spread <- factor(rr$spread, levels=c("low", "high"))
rr$scen <- factor(rr$scen, levels=c("Healthy", "Dysbiosis", "Blooming"))

p1 <- ggline(rr %>% dplyr::filter(datatable=="all"), x = "num_samples", y = "Precision", 
             add = c("mean_se"),
             color = "method", palette = "Spectral", facet.by = "scen", xlab="Number of samples") + theme_bw()
ggsave(p1, filename="output/number_samples/plot_samplesize_taxontaxon_precision.pdf", device="pdf", width=11, height=5)
p1 <- ggline(rr %>% dplyr::filter(datatable=="all"), x = "num_samples", y = "Recall", 
             add = c("mean_se"),
             color = "method", palette = "Spectral", facet.by = "scen", xlab="Number of samples") + theme_bw()
ggsave(p1, filename="output/number_samples/plot_samplesize_taxontaxon_recall.pdf", device="pdf", width=11, height=5)
p1 <- ggline(rr %>% dplyr::filter(datatable=="all"), x = "num_samples", y = "false_positive_percent", 
             add = c("mean_se"),
             color = "method", palette = "Spectral", facet.by = "scen", xlab="Number of samples") + theme_bw()
ggsave(p1, filename="output/number_samples/plot_samplesize_taxontaxon_FPS.pdf", device="pdf", width=11, height=5)


#### Assess performance of the methods across different sample numbers (taxon-metadata) ####

r20 <- read_tsv("output/number_samples/statistics_taxonmetadata_correlation_20.tsv")
r50 <- read_tsv("output/number_samples/statistics_taxonmetadata_correlation_50.tsv")
r100 <- read_tsv("output/number_samples/statistics_taxonmetadata_correlation_100.tsv")
r200 <- read_tsv("output/number_samples/statistics_taxonmetadata_correlation_200.tsv")
r500 <- read_tsv("output/number_samples/statistics_taxonmetadata_correlation_500.tsv")
r1000 <- read_tsv("output/number_samples/statistics_taxonmetadata_correlation_1000.tsv")
r20 <- r20 %>% mutate(num_samples="20")
r50 <- r50 %>% mutate(num_samples="50")
r100 <- r100 %>% mutate(num_samples="100")
r200 <- r200 %>% mutate(num_samples="200")
r500 <- r500 %>% mutate(num_samples="500")
r1000 <- r1000 %>% mutate(num_samples="1000")
rr <- bind_rows(r20,r50,r100,r200,r500,r1000)
rr$method <- factor(rr$method, levels=c("AST", "CLR", "RMP", "CSS", "GMPR",
                                        "RLE", "TMM", "UQ", "VST",
                                        "QMP", "QMP-NR"))
rr$num_samples <- factor(rr$num_samples, levels=c("20", "50", "100", 
                                                  "200", "500", "1000"))
rr$spread <- factor(rr$spread, levels=c("low", "high"))
rr$scen <- factor(rr$scen, levels=c("Healthy", "Dysbiosis", "Blooming"))

p1 <- ggline(rr %>% dplyr::filter(datatable=="all"), x = "num_samples", y = "Precision", 
             add = c("mean_se"),
             color = "method", palette = "Spectral", facet.by = "scen", xlab="Number of samples") + theme_bw()
ggsave(p1, filename="output/number_samples/plot_samplesize_taxonmetadata_precision.pdf", device="pdf", width=11, height=5)
p1 <- ggline(rr %>% dplyr::filter(datatable=="all"), x = "num_samples", y = "Recall", 
             add = c("mean_se"),
             color = "method", palette = "Spectral", facet.by = "scen", xlab="Number of samples") + theme_bw()
ggsave(p1, filename="output/number_samples/plot_samplesize_taxonmetadata_recall.pdf", device="pdf", width=11, height=5)
p1 <- ggline(rr %>% dplyr::filter(datatable=="all"), x = "num_samples", y = "false_positive_percent", 
             add = c("mean_se"),
             color = "method", palette = "Spectral", facet.by = "scen", xlab="Number of samples") + theme_bw()
ggsave(p1, filename="output/number_samples/plot_samplesize_taxonmetadata_FPS.pdf", device="pdf", width=11, height=5)



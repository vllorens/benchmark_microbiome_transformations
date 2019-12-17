# Compare QMP and QMP-NR in detecting taxon-metadata associations

# In this script we compare the performance of QMP and QMP-NR in detecting taxon-metadata associations, 
# focusing on invariant taxa (in absolute abundances) and metadata that are associated to cell counts


#### Configure the environment ####

# load packages and functions
suppressPackageStartupMessages(library(CoDaSeq))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggsignif))

# load required functions
source("R/estimateZeros.R")
source("R/taxaToKeep.R")
source("R/taxon_correlation.R")
source("R/taxon_metadata_correlation.R")
source("R/rarefy_even_sampling_depth.R")
source("R/counts_taxa_correlation.R")

# If it has not been run yet, it is necessary to run the script to generate the data (scripts/generate_data.R)
# source("scripts/generate_data.R") # this is needed to generate the original matrices and metadata only

# create folders, if not existing
system("mkdir -p output/qmp_qmpnr")

# clean these folders and remove all previous results (if any) 
system("rm -r output/qmp_qmpnr/*")

set.seed(777)
#### Compare QMP and QMP-NR at different sequencing depths ####

sequencing_avg_depths <- seq(0.000001, 0.00005, by=0.000002)
# Taxa not associated with counts
depths <- c()
fdr_qmp <- c()
fdr_qmpnr <- c()
spread_v <- c()
scenario_v <- c()
fp_qmp <- c()
fp_qmpnr <- c()
matname <- c()
for(file in list.files("data/tax_matrices", full.names = T)){
    # read original file 
    filename <- basename(file)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        round
    original_counts <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(.,2,sum)
    matrix_name <- strsplit(filename, split="_")[[1]][2]
    spread_name <- strsplit(filename, split="_")[[1]][3]
    scenario_name <- strsplit(filename, split="_")[[1]][5] %>% gsub(., pattern="\\.tsv", replacement="")
    counts_original <- apply(tax_matrix, 2, sum) %>% as.data.frame
    # read metadata
    metadatafile <- read.table(paste0("data/metadata_matrices/metadata_", filename), header=T, stringsAsFactors = F)
    # simulate sequencing + rarefying to even sampling depth
    for(depth in sequencing_avg_depths){
        print(depth)
        # sequencing
        seqOut <- data.frame(row.names=rownames(tax_matrix))
        numreads <- round(depth*mean(original_counts))
        for(col in 1:ncol(tax_matrix)){
            tvect <- sample(rownames(tax_matrix),size=rlnorm(1, meanlog = log(numreads), sdlog = 0.3), 
                            prob=tax_matrix[,col], replace=T) %>% table
            seqOut <- cbind(seqOut, tvect[rownames(seqOut)])
        }
        seqOut <- seqOut[,c(F,T)]
        seqOut[is.na(seqOut)] <- 0
        colnames(seqOut) <- colnames(tax_matrix)
        rownames(seqOut) <- rownames(tax_matrix)
        
        # generate QMP and QMP-NR data for this specific sequencing depth
        qmp <- rarefy_even_sampling_depth(cnv_corrected_abundance_table = seqOut, cell_counts_table = counts_original)
        normFactor <- apply(tax_matrix,2,sum)/apply(seqOut,2,sum)
        qmpnr <- sweep(seqOut, MARGIN = 2, normFactor, '*')
        
        # correlation of taxa w number of cell counts in the original data
        taxoncor <- counts_taxa_correlation(as.matrix(tax_matrix), original_counts)
        
        # select taxa that are not correlated to the number of cell counts, and therefore should not be correlated to metadata associated to counts
        not_cor <- which(taxoncor[[2]]>0.05)
        
        if(length(not_cor)>0){    
            #correlations on real data
            file_real <- taxon_metadata_correlation(taxon_matrix = as.matrix(tax_matrix)[not_cor,,drop=F], metadata_matrix = metadatafile[,c(46:50,96:100)])
            #correlations on qmp/qmp-nr data
            file_qmp <- taxon_metadata_correlation(taxon_matrix = as.matrix(qmp)[not_cor,,drop=F], metadata_matrix = metadatafile[,c(46:50,96:100)])
            file_qmpnr <- taxon_metadata_correlation(taxon_matrix = as.matrix(qmpnr)[not_cor,,drop=F], metadata_matrix = metadatafile[,c(46:50,96:100)])
            
            # assess performance
            test_significant_qmp <- c(unlist(file_qmp[[2]] < 0.05))
            test_significant_qmpnr <- c(unlist(file_qmpnr[[2]] < 0.05))
            reference_significant <- c(unlist(file_real[[2]] < 0.05))
            test_sign_qmp <- sign(c(unlist(file_qmp[[1]])))
            test_sign_qmpnr <- sign(c(unlist(file_qmpnr[[1]])))
            reference_sign <- sign(c(unlist(file_real[[1]])))
            test_sign_qmp[is.na(test_sign_qmp)] <- 1
            test_sign_qmpnr[is.na(test_sign_qmpnr)] <- 1
            reference_sign[is.na(reference_sign)] <- 1
            
            tpqmp <- length(which(test_significant_qmp & reference_significant & test_sign_qmp==reference_sign))
            tpqmpnr <- length(which(test_significant_qmpnr & reference_significant & test_sign_qmpnr==reference_sign))
            
            fpqmp <- length(which(test_significant_qmp & !reference_significant))+length(which(test_significant_qmp & reference_significant & test_sign_qmp!=reference_sign))
            fpqmpnr <- length(which(test_significant_qmpnr & !reference_significant))+length(which(test_significant_qmpnr & reference_significant & test_sign_qmpnr!=reference_sign))

            tnqmp <- length(which(!test_significant_qmp & !reference_significant))
            tnqmpnr <- length(which(!test_significant_qmpnr & !reference_significant))
            
            fnqmp <- length(which(!test_significant_qmp & reference_significant))
            fnqmpnr <- length(which(!test_significant_qmpnr & reference_significant))
            
            
            fdr_qmp <- c(fdr_qmp, fpqmp/(fpqmp+tpqmp))
            fdr_qmpnr <- c(fdr_qmpnr, fpqmpnr/(fpqmpnr+tpqmpnr))
            spread_v <- c(spread_v, spread_name)
            scenario_v <- c(scenario_v, scenario_name)
            fp_qmp <- c(fp_qmp, fpqmp/(fpqmp+tnqmp))
            fp_qmpnr <- c(fp_qmpnr, fpqmpnr/(fpqmpnr+tnqmpnr))
            depths <- c(depths, numreads)
            matname <- c(matname, matrix_name)
        }
    }
    result <- tibble(fdr_qmp, fdr_qmpnr, spread_v, scenario_v, fp_qmp, fp_qmpnr, depths, matname)
    write_tsv(result, "output/qmp_qmpnr/qmp_vs_qmpnr_taxonmetadata_depths.tsv", col_names = T)
}


#### Assess performance of both QMP and QMP-NR ####

result <- read_tsv("output/qmp_qmpnr/qmp_vs_qmpnr_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::select(-c(fdr_qmp, fdr_qmpnr)) %>% 
    gather(., key = "method", value="FPR", -c(spread_v, scenario_v, depths, matname))

resultlong$method <- resultlong$method %>% gsub(., pattern="fp_qmpnr", replacement="QMP-NR")
resultlong$method <- resultlong$method %>% gsub(., pattern="fp_qmp", replacement="QMP")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="low", replacement="Low")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="medium", replacement="Medium")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="high", replacement="High")
resultlong$spread_v <- factor(resultlong$spread_v, levels=c("Low", "Medium", "High"))
ggline(resultlong, x="depths", y="FPR", color="method", facet.by="spread_v",
       numeric.x.axis = T, add = "mean_se", palette=get_palette("Spectral", k = 11)[c(2,10)], xlab="Sequencing reads", 
       ylab="FPR [FP/FP+TN]", title = "Associations of invariant taxa with metadata",
       legend.title="Method") + 
    xscale("log10", .format = TRUE) + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave("output/qmp_qmpnr/qmp_vs_qmpNR_invarianttaxa_metadata.pdf", device="pdf", 
       width=11, height=3.5)

resultlong <- resultlong %>% 
    dplyr::filter(scenario_v!="Dysbiosis")
ggline(resultlong, x="depths", y="FPR", color="method", facet.by="spread_v",
       numeric.x.axis = T, add = "mean_se", palette=get_palette("Spectral", k = 11)[c(2,10)], xlab="Sequencing reads", 
       ylab="FPR [FP/FP+TN]", title = "Associations of invariant taxa with metadata",
       legend.title="Method") + 
    xscale("log10", .format = TRUE) + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave("output/qmp_qmpnr/qmp_vs_qmpNR_invarianttaxa_metadata_nodysbiosis.pdf", device="pdf", 
       width=11, height=3.5)



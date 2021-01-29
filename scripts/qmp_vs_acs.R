# Compare QMP and ACS in detecting taxon-metadata associations

# In this script we compare the performance of QMP and ACS in detecting taxon-metadata associations, 
# focusing on invariant taxa (in absolute abundances) and metadata that are associated to cell counts


#### Configure the environment ####

# load packages and functions
suppressPackageStartupMessages(library(CoDaSeq))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(inlmisc))

# load required functions
source("R/estimateZeros.R")
source("R/taxaToKeep.R")
source("R/taxon_correlation.R")
source("R/taxon_metadata_correlation.R")
source("R/rarefy_even_sampling_depth.R")
source("R/counts_taxa_correlation.R")

# If it has not been run yet, it is necessary to run the script to generate the data (scripts/generate_data.R)
source("scripts/generate_data.R") # this is needed to generate the original matrices and metadata only

# create folders, if not existing
system("mkdir -p output/qmp_acs_all")

# clean these folders and remove all previous results (if any) 
system("rm -r output/qmp_acs_all/*")

set.seed(777)

#### Compare QMP and ACS at different sequencing depths ####
sequencing_avg_depths <- seq(9.2, 14.2, length.out=25)

# Taxa not associated with counts
depths <- c()
fdr_qmp <- c()
fdr_acs <- c()
spread_v <- c()
scenario_v <- c()
fp_qmp <- c()
fp_acs <- c()
tp_qmp <- c()
tp_acs <- c()
pr_qmp <- c()
pr_acs <- c()
matname <- c()
valuesused <- c()
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
    filecounts <- gsub(file, pattern="tax_matrices/taxonomy", replacement="counts_estimated/countsEstimated_QMP_taxonomy")
    counts_estimated <- read.table(filecounts,  header=T, stringsAsFactors=F, sep="\t")
    
    # read metadata
    metadatafile <- read.table(paste0("data/metadata_matrices/metadata_", filename), header=T, stringsAsFactors = F)
    #correlations on real data
    file_real <- taxon_metadata_correlation(taxon_matrix = as.matrix(tax_matrix), metadata_matrix = metadatafile) # this is the reference file
    
    
    # simulate sequencing + rarefying to even sampling depth
    for(depth in sequencing_avg_depths){
        print(depth)
        # sequencing
        seqOut <- data.frame(row.names=rownames(tax_matrix))
        for(col in 1:ncol(tax_matrix)){
            tvect <- sample(rownames(tax_matrix),size=rlnorm(1, meanlog = depth, sdlog = 0.3), 
                            prob=tax_matrix[,col], replace=T) %>% table
            seqOut <- cbind(seqOut, tvect[rownames(seqOut)])
        }
        seqOut <- seqOut[,c(F,T)]
        seqOut[is.na(seqOut)] <- 0
        colnames(seqOut) <- colnames(tax_matrix)
        rownames(seqOut) <- rownames(tax_matrix)
        
        # generate QMP and ACS data for this specific sequencing depth
        qmp <- rarefy_even_sampling_depth(cnv_corrected_abundance_table = seqOut, cell_counts_table = counts_estimated)
        normFactor <- counts_estimated[,1]/apply(seqOut,2,sum)
        acs <- sweep(seqOut, MARGIN = 2, normFactor, '*')
        
        #correlations on qmp/ACS data (all taxon-metadata correlations, to determine Precision and Recall)
        data_qmp <- taxon_metadata_correlation(taxon_matrix = as.matrix(qmp), metadata_matrix = metadatafile)
        data_acs <- taxon_metadata_correlation(taxon_matrix = as.matrix(acs), metadata_matrix = metadatafile)

        # correlation of taxa w number of cell counts in the original data
        taxoncor <- counts_taxa_correlation(as.matrix(tax_matrix), original_counts)
        taxoncor[[2]] <- p.adjust(taxoncor[[2]], method="BH")
        
        # select taxa that are not correlated to the number of cell counts, and therefore should not be correlated to metadata associated to counts
        not_cor <- which(taxoncor[[2]]>0.05)
        
        # Assess performance on ALL TAXA
        file_qmp <- data_qmp
        file_qmp[[2]] <- matrix(p.adjust(as.vector(as.matrix(file_qmp[[2]])), method='fdr'),ncol=ncol(file_qmp[[2]]))
        colnames(file_qmp[[2]]) <- colnames(file_qmp[[1]])
        rownames(file_qmp[[2]]) <- rownames(file_qmp[[1]])
        file_acs <- data_acs
        file_acs[[2]] <- matrix(p.adjust(as.vector(as.matrix(file_acs[[2]])), method='fdr'),ncol=ncol(file_acs[[2]]))
        colnames(file_acs[[2]]) <- colnames(file_acs[[1]])
        rownames(file_acs[[2]]) <- rownames(file_acs[[1]])
        file_reference <- file_real
        file_reference[[2]] <- matrix(p.adjust(as.vector(as.matrix(file_reference[[2]])), method='fdr'),ncol=ncol(file_reference[[2]]))
        colnames(file_reference[[2]]) <- colnames(file_reference[[1]])
        rownames(file_reference[[2]]) <- rownames(file_reference[[1]])
        
        # assess performance
        test_significant_qmp_all <- c(unlist(file_qmp[[2]] < 0.05))
        test_significant_acs_all <- c(unlist(file_acs[[2]] < 0.05))
        reference_significant_all <- c(unlist(file_reference[[2]] < 0.05))
        test_sign_qmp_all <- sign(c(unlist(file_qmp[[1]])))
        test_sign_acs_all <- sign(c(unlist(file_acs[[1]])))
        reference_sign_all <- sign(c(unlist(file_reference[[1]])))
        test_sign_qmp_all[is.na(test_sign_qmp_all)] <- 1
        test_sign_acs_all[is.na(test_sign_acs_all)] <- 1
        reference_sign_all[is.na(reference_sign_all)] <- 1
        
        
        tpqmp <- length(which(test_significant_qmp_all & reference_significant_all & test_sign_qmp_all==reference_sign_all))
        tpacs <- length(which(test_significant_acs_all & reference_significant_all & test_sign_acs_all==reference_sign_all))
        
        daqmp <- length(which(test_significant_qmp_all & reference_significant_all & test_sign_qmp_all!=reference_sign_all))
        daacs <- length(which(test_significant_acs_all & reference_significant_all & test_sign_acs_all!=reference_sign_all))
        
        fpqmp <- length(which(test_significant_qmp_all & !reference_significant_all))
        fpacs <- length(which(test_significant_acs_all & !reference_significant_all))
        
        tnqmp <- length(which(!test_significant_qmp_all & !reference_significant_all))
        tnacs <- length(which(!test_significant_acs_all & !reference_significant_all))
        
        fnqmp <- length(which(!test_significant_qmp_all & reference_significant_all))
        fnacs <- length(which(!test_significant_acs_all & reference_significant_all))
        
        fdr_qmp <- c(fdr_qmp, fpqmp/(fpqmp+tpqmp))
        fdr_acs <- c(fdr_acs, fpacs/(fpacs+tpacs))
        spread_v <- c(spread_v, spread_name)
        scenario_v <- c(scenario_v, scenario_name)
        fp_qmp <- c(fp_qmp, (fpqmp+daqmp)/(fpqmp+tnqmp+daqmp))
        fp_acs <- c(fp_acs, (fpacs+daacs)/(fpacs+tnacs+daacs))
        tp_qmp <- c(tp_qmp, tpqmp/(tpqmp+fnqmp+daqmp))
        tp_acs <- c(tp_acs, tpacs/(tpacs+fnacs+daacs))
        pr_qmp <- c(pr_qmp, tpqmp/(tpqmp+fpqmp+daqmp))
        pr_acs <- c(pr_acs, tpacs/(tpacs+fpacs+daacs))
        depths <- c(depths, exp(depth))
        matname <- c(matname, matrix_name)
        valuesused <- c(valuesused, "all")
        
        # Assess performance on INVARIANT TAXA
        if(length(not_cor)>0){  
            file_qmp <- data_qmp
            file_qmp[[1]] <- file_qmp[[1]][not_cor,c(46:50,96:100)]    
            file_qmp[[2]] <- file_qmp[[2]][not_cor,c(46:50,96:100)]
            file_acs <- data_acs
            file_acs[[1]] <- file_acs[[1]][not_cor,c(46:50,96:100)]    
            file_acs[[2]] <- file_acs[[2]][not_cor,c(46:50,96:100)]  
            file_reference <- file_real
            file_reference[[1]] <- file_reference[[1]][not_cor,c(46:50,96:100)]    
            file_reference[[2]] <- file_reference[[2]][not_cor,c(46:50,96:100)]  
            
            # assess performance
            test_significant_qmp <- c(unlist(file_qmp[[2]] < 0.05))
            test_significant_acs <- c(unlist(file_acs[[2]] < 0.05))
            reference_significant <- c(unlist(file_reference[[2]] < 0.05))
            test_sign_qmp <- sign(c(unlist(file_qmp[[1]])))
            test_sign_acs <- sign(c(unlist(file_acs[[1]])))
            reference_sign <- sign(c(unlist(file_reference[[1]])))
            test_sign_qmp[is.na(test_sign_qmp)] <- 1
            test_sign_acs[is.na(test_sign_acs)] <- 1
            reference_sign[is.na(reference_sign)] <- 1
            
            
            tpqmp <- length(which(test_significant_qmp & reference_significant & test_sign_qmp==reference_sign))
            tpacs <- length(which(test_significant_acs & reference_significant & test_sign_acs==reference_sign))
            
            daqmp <- length(which(test_significant_qmp & reference_significant & test_sign_qmp!=reference_sign))
            daacs <- length(which(test_significant_acs & reference_significant & test_sign_acs!=reference_sign))
            
            fpqmp <- length(which(test_significant_qmp & !reference_significant))
            fpacs <- length(which(test_significant_acs & !reference_significant))
            
            tnqmp <- length(which(!test_significant_qmp & !reference_significant))
            tnacs <- length(which(!test_significant_acs & !reference_significant))
            
            fnqmp <- length(which(!test_significant_qmp & reference_significant))
            fnacs <- length(which(!test_significant_acs & reference_significant))
            
            fdr_qmp <- c(fdr_qmp, fpqmp/(fpqmp+tpqmp))
            fdr_acs <- c(fdr_acs, fpacs/(fpacs+tpacs))
            spread_v <- c(spread_v, spread_name)
            scenario_v <- c(scenario_v, scenario_name)
            fp_qmp <- c(fp_qmp, (fpqmp+daqmp)/(fpqmp+tnqmp+daqmp))
            fp_acs <- c(fp_acs, (fpacs+daacs)/(fpacs+tnacs+daacs))
            tp_qmp <- c(tp_qmp, tpqmp/(tpqmp+fnqmp+daqmp))
            tp_acs <- c(tp_acs, tpacs/(tpacs+fnacs+daacs))
            pr_qmp <- c(pr_qmp, tpqmp/(tpqmp+fpqmp+daqmp))
            pr_acs <- c(pr_acs, tpacs/(tpacs+fpacs+daacs))
            depths <- c(depths, exp(depth))
            matname <- c(matname, matrix_name)
            valuesused <- c(valuesused, "invariant")
        }
    }
    
    
    result <- tibble(fdr_qmp, fdr_acs, spread_v, scenario_v, fp_qmp, fp_acs, tp_qmp, tp_acs, pr_qmp, pr_acs, depths, matname, valuesused)
    write_tsv(result, "output/qmp_acs_all/qmp_vs_acs_all_taxonmetadata_depths.tsv", col_names = T)
}

mycolors=c("#83B8D7", "#A50026")

#### Assess performance of both QMP and ACS ####

result <- read_tsv("output/qmp_acs_all/qmp_vs_acs_all_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::filter(valuesused=="all") %>% 
    dplyr::select(-c(fdr_qmp, fdr_acs, tp_qmp, tp_acs, pr_qmp, pr_acs, valuesused)) %>% 
    gather(., key = "method", value="FPR", -c(spread_v, scenario_v, depths, matname))

resultlong$method <- resultlong$method %>% gsub(., pattern="fp_acs", replacement="ACS")
resultlong$method <- resultlong$method %>% gsub(., pattern="fp_qmp", replacement="QMP")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="low", replacement="Low")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="medium", replacement="Medium")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="high", replacement="High")
resultlong$spread_v <- factor(resultlong$spread_v, levels=c("Low", "Medium", "High"))
fprall <- ggline(resultlong, x="depths", y="FPR", color="method", facet.by="spread_v",
       numeric.x.axis = T, add = "mean_se", palette=mycolors, xlab="Sequencing reads", 
       ylab="FPR [FP/FP+TN]", title = "Taxa-metadata associations - FPR", xscale="log10", .format=T,
       legend.title="Method") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(fprall, file="output/qmp_acs_all/qmp_vs_acs_all_alltaxa_fdr_metadata.pdf", device="pdf", 
       width=11, height=3.5)


result <- read_tsv("output/qmp_acs_all/qmp_vs_acs_all_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::filter(valuesused=="all") %>% 
    dplyr::select(-c(fdr_qmp, fdr_acs, tp_qmp, tp_acs, fp_qmp, fp_acs, valuesused)) %>% 
    gather(., key = "method", value="Precision", -c(spread_v, scenario_v, depths, matname))


resultlong$method <- resultlong$method %>% gsub(., pattern="pr_acs", replacement="ACS")
resultlong$method <- resultlong$method %>% gsub(., pattern="pr_qmp", replacement="QMP")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="low", replacement="Low")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="medium", replacement="Medium")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="high", replacement="High")
resultlong$spread_v <- factor(resultlong$spread_v, levels=c("Low", "Medium", "High"))

precisionall <- ggline(resultlong, x="depths", y="Precision", color="method", facet.by="spread_v",
       numeric.x.axis = T, add = "mean_se", palette=mycolors, xlab="Sequencing reads", 
       ylab="Precision [TP/TP+FP]", title = "Taxa-metadata associations - Precision",xscale="log10", .format=T,
       legend.title="Method") + 
  
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(precisionall, file="output/qmp_acs_all/qmp_vs_acs_all_alltaxa_precision_metadata.pdf", device="pdf", 
       width=11, height=3.5)


result <- read_tsv("output/qmp_acs_all/qmp_vs_acs_all_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::filter(valuesused=="all") %>% 
    dplyr::select(-c(fdr_qmp, fdr_acs, pr_qmp, pr_acs, fp_qmp, fp_acs, valuesused)) %>% 
    gather(., key = "method", value="Recall", -c(spread_v, scenario_v, depths, matname))


resultlong$method <- resultlong$method %>% gsub(., pattern="tp_acs", replacement="ACS")
resultlong$method <- resultlong$method %>% gsub(., pattern="tp_qmp", replacement="QMP")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="low", replacement="Low")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="medium", replacement="Medium")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="high", replacement="High")
resultlong$spread_v <- factor(resultlong$spread_v, levels=c("Low", "Medium", "High"))

recallall <- ggline(resultlong, x="depths", y="Recall", color="method", facet.by="spread_v",
       numeric.x.axis = T, add = "mean_se", palette=mycolors, xlab="Sequencing reads", 
       ylab="Recall [TP/TP+FN]", title = "Taxa-metadata associations - Recall", xscale="log10", .format=T,
       legend.title="Method") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(recallall, file="output/qmp_acs_all/qmp_vs_acs_all_alltaxa_recall_metadata.pdf", device="pdf", 
       width=11, height=3.5)


#### Assess performance of both QMP and ACS (invariant taxa only, no dysbiosis) ####

result <- read_tsv("output/qmp_acs_all/qmp_vs_acs_all_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::filter(valuesused=="invariant") %>% 
   # dplyr::filter(scenario_v!="Dysbiosis") %>% 
    dplyr::select(-c(fdr_qmp, fdr_acs, tp_qmp, tp_acs, pr_qmp, pr_acs, valuesused)) %>% 
    gather(., key = "method", value="FPR", -c(spread_v, scenario_v, depths, matname))



resultlong$method <- resultlong$method %>% gsub(., pattern="fp_acs", replacement="ACS")
resultlong$method <- resultlong$method %>% gsub(., pattern="fp_qmp", replacement="QMP")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="low", replacement="Low")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="medium", replacement="Medium")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="high", replacement="High")
resultlong$spread_v <- factor(resultlong$spread_v, levels=c("Low", "Medium", "High"))
fprinv <- ggline(resultlong, x="depths", y="FPR", color="method", facet.by="spread_v",
       numeric.x.axis = T, add = "mean_se", palette=mycolors, xlab="Sequencing reads", 
       ylab="FPR [FP/FP+TN]", title = "Invariant taxa-metadata associations - FPR",
       legend.title="Method", xscale="log10", .format=T) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(fprinv, file="output/qmp_acs_all/qmp_vs_acs_all_invtaxa_fdr_metadata.pdf", device="pdf", 
       width=11, height=3.5)


result <- read_tsv("output/qmp_acs_all/qmp_vs_acs_all_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::filter(valuesused=="invariant") %>% 
    #dplyr::filter(scenario_v!="Dysbiosis") %>% 
    dplyr::select(-c(fdr_qmp, fdr_acs, tp_qmp, tp_acs, fp_qmp, fp_acs, valuesused)) %>% 
    gather(., key = "method", value="Precision", -c(spread_v, scenario_v, depths, matname))


resultlong$method <- resultlong$method %>% gsub(., pattern="pr_acs", replacement="ACS")
resultlong$method <- resultlong$method %>% gsub(., pattern="pr_qmp", replacement="QMP")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="low", replacement="Low")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="medium", replacement="Medium")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="high", replacement="High")
resultlong$spread_v <- factor(resultlong$spread_v, levels=c("Low", "Medium", "High"))

precisioninv <- ggline(resultlong, x="depths", y="Precision", color="method", facet.by="spread_v",
       numeric.x.axis = T, add = "mean_se", palette=mycolors, xlab="Sequencing reads", 
       ylab="Precision [TP/TP+FP]", title = "Invariant taxa-metadata associations - Precision",
       legend.title="Method", xscale="log10", .format=T) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(precisioninv, file="output/qmp_acs_all/qmp_vs_acs_all_invtaxa_precision_metadata.pdf", device="pdf", 
       width=11, height=3.5)


result <- read_tsv("output/qmp_acs_all/qmp_vs_acs_all_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::filter(valuesused=="invariant") %>% 
    #dplyr::filter(scenario_v!="Dysbiosis") %>% 
    dplyr::select(-c(fdr_qmp, fdr_acs, pr_qmp, pr_acs, fp_qmp, fp_acs,valuesused)) %>% 
    gather(., key = "method", value="Recall", -c(spread_v, scenario_v, depths, matname))


resultlong$method <- resultlong$method %>% gsub(., pattern="tp_acs", replacement="ACS")
resultlong$method <- resultlong$method %>% gsub(., pattern="tp_qmp", replacement="QMP")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="low", replacement="Low")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="medium", replacement="Medium")
resultlong$spread_v <- resultlong$spread_v %>% gsub(., pattern="high", replacement="High")
resultlong$spread_v <- factor(resultlong$spread_v, levels=c("Low", "Medium", "High"))

recallinv <- ggline(resultlong, x="depths", y="Recall", color="method", facet.by="spread_v",
       numeric.x.axis = T, add = "mean_se", palette=mycolors, xlab="Sequencing reads", 
       ylab="Recall [TP/TP+FN]", title = "Invariant taxa-metadata associations - Recall",
       legend.title="Method", xscale="log10", .format=T) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(recallinv, file="output/qmp_acs_all/qmp_vs_acs_all_invtaxa_recall_metadata.pdf", device="pdf", 
       width=11, height=3.5)


pub_plot <- ggarrange(recallall, fprall, fprinv, ncol=1, nrow=3, labels="auto",legend="right")
ggsave(pub_plot, file="output/qmp_acs_all/qmp_vs_acs_all_summary.pdf", device="pdf", 
       width=8, height=9,useDingbats=FALSE)



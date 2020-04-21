# Proportionality vs Spearman correlation
# Tue Dec  3 09:53:42 2019 ------------------------------

# This script evaluates the performance of proportionality vs Spearman correlation in capturing taxon-taxon correlations
# It is highly recommended to run this script in a cluster to take maximum advantage of the multicore functionality

#### Configure the environment ####

# load required packages
suppressPackageStartupMessages(library(CoDaSeq))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(furrr))

# this will parallelize the proportionality functions to use all available cores, but it is possible to use only a subset (see function help)
# might not work in some OS
# See also https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html
plan(multicore) 

# load functions
source("R/estimateZeros.R")
source("R/taxon_correlation.R")
source("R/propr_custom.R")

# create matrices (if not there)
system("mkdir -p data/correlations_proportionality")
system("mkdir -p output/proportionality")

# clean these folders and remove all previous results (if any) 
system("rm -r data/correlations_proportionality/*")
system("rm -r output/proportionality/*")

# create additional folders for the reference data (on the real matrices)
system("mkdir -p data/correlations_proportionality/reference")



#### Calculate average abundance of taxa pairs detected in true correlations ####
set.seed(1)

# load data
tax_matrices <- list.files("data/tax_matrices", full.names = T)
cutoff_significance <- 0.05
for(file in tax_matrices){
    # read files
    filename <- basename(file)
    seqname <- gsub(filename, pattern = "taxonomy", replacement="seqOut_taxonomy")
    original_tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t")
    sequencing_matrix <- read.table(paste0("data/seq_matrices/", seqname), header=T, stringsAsFactors=F, sep="\t")
    matrixname <- strsplit(filename, split = "_")[[1]][2]
    spreadname <- strsplit(filename, split = "_")[[1]][3]
    # calculate spearman correlation in the original data
    original_spearman <- taxon_correlation(as.matrix(original_tax_matrix))
    # calculate spearman correlation in original clr-transformed data
    original_nozeros <- estimate0.min(original_tax_matrix)
    original_clr <- t(codaSeq.clr(original_nozeros, samples.by.row = F))
    original_clr_spearman <- taxon_correlation(as.matrix(original_clr))
    # calculate spearman correlation in sequencing clr-transformed data
    sequencing_nozeros <- estimate0.min(sequencing_matrix)
    sequencing_clr <- t(codaSeq.clr(sequencing_nozeros, samples.by.row = F))
    sequencing_clr_spearman <- taxon_correlation(sequencing_clr)
    
    # write output - spearman-real
    outputname_correlation <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_rho_cor_REAL")
    outputname_pval <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_pvalue_cor_REAL")
    write.table(original_spearman[[1]], 
                paste0("data/correlations_proportionality/reference/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(original_spearman[[2]], 
                paste0("data/correlations_proportionality/reference/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output - spearman-real-clr
    outputname_correlation <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_rho_cor_CLR")
    outputname_pval <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_pvalue_cor_CLR")
    write.table(original_clr_spearman[[1]], 
                paste0("data/correlations_proportionality/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(original_clr_spearman[[2]], 
                paste0("data/correlations_proportionality/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output - spearman-seq-clr
    outputname_correlation <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_rho_cor_SeqCLR")
    outputname_pval <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_pvalue_cor_SeqCLR")
    write.table(sequencing_clr_spearman[[1]], 
                paste0("data/correlations_proportionality/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(sequencing_clr_spearman[[2]], 
                paste0("data/correlations_proportionality/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    
    # calculate proportionality in original data
    original_rho <- propr_custom(original_tax_matrix)
    # calculate proportionality in original clr-transformed data
    original_clr_rho <- propr_custom(as.matrix(original_clr))
    # calculate proportionality in sequencing clr-transformed data
    sequencing_clr_rho <- propr_custom(sequencing_clr)
    
    # write output - propr-real
    outputname_correlation <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_rho_propr_REAL")
    outputname_pval <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_pvalue_propr_REAL")
    write.table(original_rho[[1]], 
                paste0("data/correlations_proportionality/reference/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(original_rho[[2]], 
                paste0("data/correlations_proportionality/reference/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output - propr-real-clr
    outputname_correlation <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_rho_propr_CLR")
    outputname_pval <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_pvalue_propr_CLR")
    write.table(original_clr_rho[[1]], 
                paste0("data/correlations_proportionality/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(original_clr_rho[[2]], 
                paste0("data/correlations_proportionality/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
    # write output - propr-seq-clr
    outputname_correlation <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_rho_propr_SeqCLR")
    outputname_pval <- gsub(filename, pattern="taxonomy", replacement="taxontaxon_pvalue_propr_SeqCLR")
    write.table(sequencing_clr_rho[[1]], 
                paste0("data/correlations_proportionality/", outputname_correlation), 
                col.names=T, row.names=T, quote=F, sep="\t")
    write.table(sequencing_clr_rho[[2]], 
                paste0("data/correlations_proportionality/", outputname_pval), 
                col.names=T, row.names=T, quote=F, sep="\t")
}




#### Plot performance of one metric over the other in the original matrices (correlation among them) ####

# Blooming (one example)
file_cor <- "data/correlations_proportionality/reference/taxontaxon_rho_cor_REAL_B1_high_spread_Blooming.tsv"
file_prop <- "data/correlations_proportionality/reference/taxontaxon_rho_propr_REAL_B1_high_spread_Blooming.tsv"
file_cor <- read.table(file_cor, header=T, stringsAsFactors = F,sep="\t")
file_prop <- read.table(file_prop, header=T, stringsAsFactors=F, sep="\t")

values <- tibble(Spearman = c(unlist(file_cor)), Proportionality = c(unlist(file_prop))) %>% 
    drop_na()

# first plot shows direct correlation of values
p1 <- ggscatter(values, x="Spearman", y="Proportionality", fill="#E64B35FF", 
                title="Proportionality vs Spearmann correlation values, blooming",
                add="reg.line", conf.int=T, cor.coef = TRUE, alpha=0.1, shape = 21,
                cor.coeff.args = list(method = "spearman", label.x = -1, label.sep = "\n")) +
    theme_bw() + 
    theme(plot.title = element_text(face = "bold"))
ggsave("output/proportionality/correlation_proportionality_plots_blooming.pdf", device="pdf", width=5, height=5)

# Dysbiosis (one example)
file_cor <- "data/correlations_proportionality/reference/taxontaxon_rho_cor_REAL_D2_high_spread_Dysbiosis.tsv"
file_prop <- "data/correlations_proportionality/reference/taxontaxon_rho_propr_REAL_D2_high_spread_Dysbiosis.tsv"
file_cor <- read.table(file_cor, header=T, stringsAsFactors = F,sep="\t")
file_prop <- read.table(file_prop, header=T, stringsAsFactors=F, sep="\t")

values <- tibble(Spearman = c(unlist(file_cor)), Proportionality = c(unlist(file_prop))) %>% 
    drop_na()

# first plot shows direct correlation of values
p2 <- ggscatter(values, x="Spearman", y="Proportionality", fill="#E64B35FF", 
                title="Proportionality vs Spearmann correlation values, dysbiosis",
                add="reg.line", conf.int=T, cor.coef = TRUE, alpha=0.1, shape = 21,
                cor.coeff.args = list(method = "spearman", label.x = 0, label.sep = "\n")) +
    theme_bw() + 
    theme(plot.title = element_text(face = "bold"))
ggsave("output/proportionality/correlation_proportionality_plots_dysbiosis.pdf", device="pdf", width=5, height=5)

# Healthy (one example)
file_cor <- "data/correlations_proportionality/reference/taxontaxon_rho_cor_REAL_S3_low_spread_Healthy.tsv"
file_prop <- "data/correlations_proportionality/reference/taxontaxon_rho_propr_REAL_S3_low_spread_Healthy.tsv"
file_cor <- read.table(file_cor, header=T, stringsAsFactors = F,sep="\t")
file_prop <- read.table(file_prop, header=T, stringsAsFactors=F, sep="\t")

values <- tibble(Spearman = c(unlist(file_cor)), Proportionality = c(unlist(file_prop))) %>% 
    drop_na()

# first plot shows direct correlation of values
p3 <- ggscatter(values, x="Spearman", y="Proportionality", fill="#E64B35FF", 
                title="Proportionality vs Spearmann correlation values, healthy",
                add="reg.line", conf.int=T, cor.coef = TRUE, alpha=0.1, shape = 21,
                cor.coeff.args = list(method = "spearman", label.x = -1, label.sep = "\n")) +
    theme_bw() + 
    theme(plot.title = element_text(face = "bold"))
ggsave("output/proportionality/correlation_proportionality_plots_healthy.pdf", device="pdf", width=5, height=5)




#### Plot performance of one metric over the other in the original matrices (example showing 40 taxa) ####

# Blooming (example)
postscript("output/proportionality/correlation_proportionality_plots2_blooming.ps", width=6, height=6)
file_cor <- "data/correlations_proportionality/reference/taxontaxon_rho_cor_REAL_B1_high_spread_Blooming.tsv"
file_prop <- "data/correlations_proportionality/reference/taxontaxon_rho_propr_REAL_B1_high_spread_Blooming.tsv"
file_cor <- read.table(file_cor, header=T, stringsAsFactors = F,sep="\t")
file_prop <- read.table(file_prop, header=T, stringsAsFactors=F, sep="\t")

file_prop <- t(file_prop)
merged_file <- file_cor
merged_file[lower.tri(merged_file, diag = F)] <- file_prop[lower.tri(file_prop, diag = F)]
merged_file <- as.matrix(merged_file)

# read pvals
file_pval_cor <- "data/correlations_proportionality/reference/taxontaxon_pvalue_cor_REAL_B1_high_spread_Blooming.tsv"
file_pval_prop <- "data/correlations_proportionality/reference/taxontaxon_pvalue_propr_REAL_B1_high_spread_Blooming.tsv"
file_pval_cor <- read.table(file_pval_cor, header=T, stringsAsFactors = F,sep="\t")
file_pval_prop <- read.table(file_pval_prop, header=T, stringsAsFactors=F, sep="\t")


file_pval_prop <- t(file_pval_prop)
merged_pval_file <- file_pval_cor
merged_pval_file[lower.tri(merged_pval_file, diag = F)] <- file_pval_prop[lower.tri(file_pval_prop, diag = F)]
merged_pval_file <- as.matrix(merged_pval_file)


diag(merged_file) <- 1
diag(merged_pval_file) <- 0.001

corrplot::corrplot(merged_file[1:40, 1:40], type = "full",
                   mar = c(1,2,3,2),
                   title = "Correlation and proportionality, blooming, showing first 40 taxa", 
                   tl.cex = 0.5, tl.col = "black", tl.srt = 40, cl.ratio = 0.1, 
                   p.mat = merged_pval_file[1:40, 1:40], sig.level = 0.05, insig = "blank")
dev.off()

# Dysbiosis (example)
postscript("output/proportionality/correlation_proportionality_plots2_dysbiosis.ps", width=6, height=6)
file_cor <- "data/correlations_proportionality/reference/taxontaxon_rho_cor_REAL_D2_high_spread_Dysbiosis.tsv"
file_prop <- "data/correlations_proportionality/reference/taxontaxon_rho_propr_REAL_D2_high_spread_Dysbiosis.tsv"
file_cor <- read.table(file_cor, header=T, stringsAsFactors = F,sep="\t")
file_prop <- read.table(file_prop, header=T, stringsAsFactors=F, sep="\t")

file_prop <- t(file_prop)
merged_file <- file_cor
merged_file[lower.tri(merged_file, diag = F)] <- file_prop[lower.tri(file_prop, diag = F)]
merged_file <- as.matrix(merged_file)

# read pvals
file_pval_cor <- "data/correlations_proportionality/reference/taxontaxon_pvalue_cor_REAL_D2_high_spread_Dysbiosis.tsv"
file_pval_prop <- "data/correlations_proportionality/reference/taxontaxon_pvalue_propr_REAL_D2_high_spread_Dysbiosis.tsv"
file_pval_cor <- read.table(file_pval_cor, header=T, stringsAsFactors = F,sep="\t")
file_pval_prop <- read.table(file_pval_prop, header=T, stringsAsFactors=F, sep="\t")


file_pval_prop <- t(file_pval_prop)
merged_pval_file <- file_pval_cor
merged_pval_file[lower.tri(merged_pval_file, diag = F)] <- file_pval_prop[lower.tri(file_pval_prop, diag = F)]
merged_pval_file <- as.matrix(merged_pval_file)


diag(merged_file) <- 1
diag(merged_pval_file) <- 0.001

corrplot::corrplot(merged_file[1:40, 1:40], type = "full",
                   mar = c(1,2,3,2),
                   title = "Correlation and proportionality, dysbiosis, showing first 40 taxa", 
                   tl.cex = 0.5, tl.col = "black", tl.srt = 40, cl.ratio = 0.1, 
                   p.mat = merged_pval_file[1:40, 1:40], sig.level = 0.05, insig = "blank")
dev.off()

# Healthy (example)
postscript("output/proportionality/correlation_proportionality_plots2_healthy.ps", width=6, height=6)
file_cor <- "data/correlations_proportionality/reference/taxontaxon_rho_cor_REAL_S3_low_spread_Healthy.tsv"
file_prop <- "data/correlations_proportionality/reference/taxontaxon_rho_propr_REAL_S3_low_spread_Healthy.tsv"
file_cor <- read.table(file_cor, header=T, stringsAsFactors = F,sep="\t")
file_prop <- read.table(file_prop, header=T, stringsAsFactors=F, sep="\t")

file_prop <- t(file_prop)
merged_file <- file_cor
merged_file[lower.tri(merged_file, diag = F)] <- file_prop[lower.tri(file_prop, diag = F)]
merged_file <- as.matrix(merged_file)

# read pvals
file_pval_cor <- "data/correlations_proportionality/reference/taxontaxon_pvalue_cor_REAL_S3_low_spread_Healthy.tsv"
file_pval_prop <- "data/correlations_proportionality/reference/taxontaxon_pvalue_propr_REAL_S3_low_spread_Healthy.tsv"
file_pval_cor <- read.table(file_pval_cor, header=T, stringsAsFactors = F,sep="\t")
file_pval_prop <- read.table(file_pval_prop, header=T, stringsAsFactors=F, sep="\t")


file_pval_prop <- t(file_pval_prop)
merged_pval_file <- file_pval_cor
merged_pval_file[lower.tri(merged_pval_file, diag = F)] <- file_pval_prop[lower.tri(file_pval_prop, diag = F)]
merged_pval_file <- as.matrix(merged_pval_file)


diag(merged_file) <- 1
diag(merged_pval_file) <- 0.001

corrplot::corrplot(merged_file[1:40, 1:40], type = "full",
                   mar = c(1,2,3,2),
                   title = "Correlation and proportionality, healthy, showing first 40 taxa", 
                   tl.cex = 0.5, tl.col = "black", tl.srt = 40, cl.ratio = 0.1, 
                   p.mat = merged_pval_file[1:40, 1:40], sig.level = 0.05, insig = "blank")

dev.off()



#### Calculate precision and recall of the propr vs spearman in the original matrices ####
# define parameters for significance and initialize vectors for the evaluation data
significance_level <- 0.05
true_positive <- c()
false_positive <- c()
true_negative <- c()
false_negative <- c()
method <- c()
spread <- c()
matrixnum <- c()
scenario <- c()
for(file in list.files("data/correlations_proportionality/reference", recursive = F, pattern = "pvalue_propr", full.names = T)){
    # define file name, method used and reference to be compared against
    filename <- basename(file)
    methodname <- strsplit(filename, split = "_")[[1]][3]
    matrixname <- strsplit(filename, split = "_")[[1]][5]
    spreadname <- strsplit(filename, split = "_")[[1]][6]
    scenarioname <- strsplit(filename, split = "_")[[1]][8] %>% gsub(., pattern="\\.tsv", replacement="")
    referencename <- gsub(filename, pattern="propr", replacement="cor")
    rhoname <- gsub(filename, pattern="pvalue", replacement="rho")
    referencerhoname <- gsub(referencename, pattern="pvalue", replacement="rho")
    
    # read files
    test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
    reference_file <- read.table(paste0("data/correlations_proportionality/reference/", referencename),
                                 header=T, stringsAsFactors=F, sep="\t")
    rho_file <- read.table(paste0("data/correlations_proportionality/reference/", rhoname),
                           header=T, stringsAsFactors = F, sep="\t")
    referencerho_file <- read.table(paste0("data/correlations_proportionality/reference/", referencerhoname),
                                    header=T, stringsAsFactors = F, sep = "\t")
    
    # get significance
    test_significant <-  c(unlist(test_file < significance_level))
    reference_significant <- c(unlist(reference_file < significance_level))
    test_sign <- sign(c(unlist(rho_file)))
    reference_sign <- sign(c(unlist(referencerho_file)))
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    true_positive <- c(true_positive, 
                       length(which(test_significant & reference_significant & test_sign==reference_sign)))
    false_positive <- c(false_positive, 
                        length(which(test_significant & !reference_significant))+length(which(test_significant & reference_significant & test_sign!=reference_sign)))
    true_negative <- c(true_negative, 
                       length(which(!test_significant & !reference_significant)))
    false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
    matrixnum <- c(matrixnum, matrixname)
    scenario <- c(scenario, scenarioname)
}

# write results table
results <- tibble(method, spread, matrixnum, scenario,
                  true_positive, false_positive, true_negative, false_negative)
results <- results %>%
    mutate(FDR=100*false_positive/(false_positive+true_positive)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Recall=100*true_positive/(true_positive+false_negative)) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive)) %>% 
    mutate(Accuracy=100*(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(true_positive_percent=100*true_positive/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_positive_percent=100*false_positive/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(true_negative_percent=100*true_negative/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_negative_percent=100*false_negative/(true_positive+true_negative+false_positive+false_negative)) %>% 
    dplyr::select(-c(true_positive,true_negative,false_positive,false_negative))

write_tsv(results, paste0("output/proportionality/statistics_proportionality_spearman_original.tsv"), col_names = T)

# make plot
results <- gather(results, key="metric", value="value", 
                  -c(method, spread, matrixnum, scenario)) %>% 
    drop_na()
results$spread <- factor(results$spread, levels=c("low", "high"))


plot_a <- ggboxplot(results %>% dplyr::filter(metric %in% c("Precision", "Recall")),
                    x="scenario", y="value", fill="scenario", palette=get_palette("Spectral", 6)[c(1,3,6)],
                    add="jitter", 
                    title="Precision and recall of Spearman correlation vs proportionality") + 
    theme_bw()  + theme(plot.title = element_text(face = "bold", size = 10)) +
    facet_grid(metric ~ .) +
    stat_compare_means()
ggsave(plot_a, filename = "output/proportionality/propr_performance_original_matrices.ps", 
       device = "ps", width = 5, height=7)




#### Evaluate performance of spearman and proportionality in the CLR transformed data ####
# define parameters for significance and initialize vectors for the evaluation data
significance_level <- 0.05
true_positive <- c()
false_positive <- c()
true_negative <- c()
false_negative <- c()
method <- c()
spread <- c()
matrixnum <- c()
matrixclass <- c()
scenario <- c()

for(file in list.files("data/correlations_proportionality/", recursive = F, pattern = "pvalue", full.names = T)){
    # define file name, method used and reference to be compared against
    filename <- basename(file)
    methodname <- strsplit(filename, split = "_")[[1]][3]
    matrixtype <- strsplit(filename, split = "_")[[1]][4]
    matrixname <- strsplit(filename, split = "_")[[1]][5]
    spreadname <- strsplit(filename, split = "_")[[1]][6]
    referencename <- gsub(filename, pattern=matrixtype, replacement="REAL")
    rhoname <- gsub(filename, pattern="pvalue", replacement="rho")
    referencerhoname <- gsub(referencename, pattern="pvalue", replacement="rho")
    scenarioname <- strsplit(filename, split = "_")[[1]][8] %>% gsub(., pattern="\\.tsv", replacement="")
    
    # read files
    test_file <- read.table(file, header=T, stringsAsFactors = F,sep="\t")
    reference_file <- read.table(paste0("data/correlations_proportionality/reference/", referencename),
                                 header=T, stringsAsFactors=F, sep="\t")
    rho_file <- read.table(paste0("data/correlations_proportionality/", rhoname),
                           header=T, stringsAsFactors = F, sep="\t")
    referencerho_file <- read.table(paste0("data/correlations_proportionality/reference/", referencerhoname),
                                    header=T, stringsAsFactors = F, sep = "\t")
    
    # get significance
    test_significant <-  c(unlist(test_file < significance_level))
    reference_significant <- c(unlist(reference_file < significance_level))
    test_sign <- sign(c(unlist(rho_file)))
    reference_sign <- sign(c(unlist(referencerho_file)))
    
    # populate evaluation vectors
    method <- c(method, methodname)
    spread <- c(spread, spreadname)
    true_positive <- c(true_positive, 
                       length(which(test_significant & reference_significant & test_sign==reference_sign)))
    false_positive <- c(false_positive, 
                        length(which(test_significant & !reference_significant))+length(which(test_significant & reference_significant & test_sign!=reference_sign)))
    true_negative <- c(true_negative, 
                       length(which(!test_significant & !reference_significant)))
    false_negative <- c(false_negative, length(which(!test_significant & reference_significant)))
    matrixnum <- c(matrixnum, matrixname)
    matrixclass <- c(matrixclass, matrixtype)
    scenario <- c(scenario, scenarioname)
    
}

# write results table
results <- tibble(method, spread, matrixnum, matrixclass, scenario,
                  true_positive, false_positive, true_negative, false_negative)
results <- results %>%
    mutate(FDR=100*false_positive/(false_positive+true_positive)) %>% 
    mutate(Precision=100-FDR) %>% 
    mutate(Recall=100*true_positive/(true_positive+false_negative)) %>% 
    mutate(Specificity=100*true_negative/(true_negative+false_positive)) %>% 
    mutate(Accuracy=100*(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(true_positive_percent=100*true_positive/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_positive_percent=100*false_positive/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(true_negative_percent=100*true_negative/(true_positive+true_negative+false_positive+false_negative)) %>% 
    mutate(false_negative_percent=100*false_negative/(true_positive+true_negative+false_positive+false_negative)) %>% 
    dplyr::select(-c(true_positive,true_negative,false_positive,false_negative))

write_tsv(results, paste0("output/proportionality/statistics_proportionality_spearman_original.tsv"), col_names = T)

# make plot 
results <- gather(results, key="metric", value="value", 
                  -c(method, spread, matrixnum, matrixclass, scenario)) %>% 
    drop_na()
results$spread <- factor(results$spread, levels=c("low",  "high"))

my_comparisons <- list( c("cor", "propr"))
plot_a <- ggboxplot(results %>% dplyr::filter(metric %in% c("Precision", "Recall")),
                    x="scenario", y="value", fill="method", palette=get_palette("Spectral", 6)[c(1,6)], group="matrixnum",
                    add="jitter", facet.by = c("metric","matrixclass"), ylim=c(0,60),
                    title="Precision and recall of Spearman correlation vs proportionality") + 
    theme_bw()  + theme(plot.title = element_text(face = "bold")) +
    stat_compare_means(aes(group = method), label = "p.signif",method = "wilcox.test", paired=T, label.y = 55) 


ggsave(plot_a, filename = "output/proportionality/propr_performance_clr_seqclr.ps", device = "ps", width=10, height=10)



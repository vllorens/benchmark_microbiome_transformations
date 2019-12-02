# Test of the effect of microbial load spread
# Mon Dec  2 15:58:34 2019 ------------------------------


#### Configure the environment ####

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(gdata))


# create folders, if not existing
system("mkdir -p output/microbial_spread/")

# clean these folders and remove all previous results (if any) 
system("rm -r output/microbial_spread/*")

# If it has not been run yet, it is necessary to run the script to generate the data (scripts/generate_data.R)
# source("scripts/generate_data.R")
# Also, the correlations script needs to be run
# source("scripts/correlations.R")


#### Calculate effect of spread of microbial loads in the results ####

# info of microbial load spread
matrixinfo <- read_tsv("data/raw/matrix_stats_3scenarios.tsv")

lowspread <- matrixinfo %>% 
    dplyr::filter(spread < 20) %>% 
    dplyr::select(matrix) %>% 
    mutate(spread="low")

midspread <- matrixinfo %>% 
    dplyr::filter(spread >= 20 & spread < 40) %>% 
    dplyr::select(matrix) %>% 
    mutate(spread="medium")

highspread <- matrixinfo %>% 
    dplyr::filter(spread >= 40) %>% 
    dplyr::select(matrix) %>% 
    mutate(spread="high")

spreads <- bind_rows(lowspread, midspread, highspread)

# read results (taxon taxon)
resultsinfo <- read_tsv("output/taxontaxon/statistics_taxontaxon_correlation.tsv", col_names = T)

resultsinfo <- resultsinfo %>% 
    dplyr::select(-spread) %>% 
    left_join(spreads, by=c("matrixnum"="matrix"))


# plot
resultsinfo$method <- factor(resultsinfo$method, levels=c("AST", "CLR", "RMP", "CSS","GMPR",
                                                          "RLE", "TMM", "UQ", "VST",
                                                          "QMP", "QMP-NR"))

resultsinfo$spread <- factor(resultsinfo$spread, levels=c("low", "medium", "high"))
results_all <- resultsinfo %>% dplyr::filter(datatable=="all")

ggboxplot(results_all, x="spread", y="Precision", alpha=0.5,
          fill="spread", ylim=c(25,110), facet.by = "method", nrow=1,
          ylab="Precision [TP/TP+FP]",
          palette="Spectral",
          main="Taxon-taxon correlations | Effect of microbial load spreads") + 
    geom_signif(comparisons=list(c("low", "medium"),
                                 c("medium", "high"),
                                 c("low", "high")), FDR = T, 
                step_increase = c(0.05,0.05,0.05),tip_length = 0, map_signif_level = T) +
    geom_point(aes(fill=spread), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text() +
    theme(plot.title = element_text(face = "bold")) + theme(axis.title = element_text(size = 14), 
                                                            axis.text.x = element_text(size = 12), 
                                                            axis.text.y = element_text(size = 12), 
                                                            plot.title = element_text(size = 16)) +labs(x = NULL)

ggsave(filename = "output/microbial_spread/plot_taxontaxon_precision_spread.ps", device = "ps", width = 11, height=5)

ggboxplot(results_all, x="spread", y="Recall", alpha=0.5,
          fill="spread", facet.by = "method", nrow=1,
          ylab="Recall [TP/TP+FN]",
          palette="Spectral",
          main="Taxon-taxon correlations | Effect of microbial load spreads") + 
    geom_signif(comparisons=list(c("low", "medium"),
                                 c("medium", "high"),
                                 c("low", "high")), FDR = T, 
                step_increase = c(0.05,0.05,0.05),tip_length = 0, map_signif_level = T) +
    geom_point(aes(fill=spread), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text() +
    theme(plot.title = element_text(face = "bold")) + theme(axis.title = element_text(size = 14), 
                                                            axis.text.x = element_text(size = 12), 
                                                            axis.text.y = element_text(size = 12), 
                                                            plot.title = element_text(size = 16)) +labs(x = NULL)

ggsave(filename = "output/microbial_spread/plot_taxontaxon_recall_spread.ps", device = "ps", width = 11, height=5)

ggboxplot(results_all, x="spread", y="false_positive_percent", alpha=0.5,
          fill="spread", facet.by = "method", nrow=1,
          ylab="% False positives",
          palette="Spectral",
          main="Taxon-taxon correlations | Effect of microbial load spreads") + 
    geom_signif(comparisons=list(c("low", "medium"),
                                 c("medium", "high"),
                                 c("low", "high")), FDR = T, 
                step_increase = c(0.05,0.05,0.05),tip_length = 0, map_signif_level = T) +
    geom_point(aes(fill=spread), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text() +
    theme(plot.title = element_text(face = "bold")) + theme(axis.title = element_text(size = 14), 
                                                            axis.text.x = element_text(size = 12), 
                                                            axis.text.y = element_text(size = 12), 
                                                            plot.title = element_text(size = 16)) +labs(x = NULL)

ggsave(filename = "output/microbial_spread/plot_taxontaxon_FP_spread.ps", device = "ps", width = 11, height=5)





# read results (taxon metadata)
resultsinfo <- read_tsv("output/taxonmetadata/statistics_taxonmetadata_correlation.tsv", col_names = T)

resultsinfo <- resultsinfo %>% 
    dplyr::select(-spread) %>% 
    left_join(spreads, by=c("matrixnum"="matrix"))


# plot
resultsinfo$method <- factor(resultsinfo$method, levels=c("AST", "CLR", "RMP", "CSS","GMPR",
                                                          "RLE", "TMM", "UQ", "VST",
                                                          "QMP", "QMP-NR"))

resultsinfo$spread <- factor(resultsinfo$spread, levels=c("low", "medium", "high"))
results_all <- resultsinfo %>% dplyr::filter(datatable=="all")
results_pos <- resultsinfo %>% dplyr::filter(datatable=="pos")
results_neg <- resultsinfo %>% dplyr::filter(datatable=="neg")
results_not <- resultsinfo %>% dplyr::filter(datatable=="not")

ggboxplot(results_all, x="spread", y="Precision", alpha=0.5,
          fill="spread", facet.by = "method", nrow=1,
          ylab="Precision [TP/TP+FP]",
          palette="Spectral",
          main="Taxon-metadata correlations | Effect of microbial load spreads") + 
    geom_signif(comparisons=list(c("low", "medium"),
                                 c("medium", "high"),
                                 c("low", "high")), FDR = T, 
                step_increase = c(0.05,0.05,0.05),tip_length = 0, map_signif_level = T) +
    geom_point(aes(fill=spread), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text() +
    theme(plot.title = element_text(face = "bold")) + theme(axis.title = element_text(size = 14), 
                                                            axis.text.x = element_text(size = 12), 
                                                            axis.text.y = element_text(size = 12), 
                                                            plot.title = element_text(size = 16)) +labs(x = NULL)

ggsave(filename = "output/microbial_spread/plot_taxonmetadata_precision_spread.ps", device = "ps", width = 11, height=5)

ggboxplot(results_all, x="spread", y="Recall", alpha=0.5,
          fill="spread", facet.by = "method", nrow=1,
          ylab="Recall [TP/TP+FN]",
          palette="Spectral",
          main="Taxon-metadata correlations | Effect of microbial load spreads") + 
    geom_signif(comparisons=list(c("low", "medium"),
                                 c("medium", "high"),
                                 c("low", "high")), FDR = T, 
                step_increase = c(0.05,0.05,0.05),tip_length = 0, map_signif_level = T) +
    geom_point(aes(fill=spread), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text() +
    theme(plot.title = element_text(face = "bold")) + theme(axis.title = element_text(size = 14), 
                                                            axis.text.x = element_text(size = 12), 
                                                            axis.text.y = element_text(size = 12), 
                                                            plot.title = element_text(size = 16)) +labs(x = NULL)

ggsave(filename = "output/microbial_spread/plot_taxonmetadata_recall_spread.ps", device = "ps", width = 11, height=5)

ggboxplot(results_all, x="spread", y="false_positive_percent", alpha=0.5,
          fill="spread", facet.by = "method", nrow=1,
          ylab="% False positives",
          palette="Spectral",
          main="Taxon-metadata correlations | Effect of microbial load spreads") + 
    geom_signif(comparisons=list(c("low", "medium"),
                                 c("medium", "high"),
                                 c("low", "high")), FDR = T, 
                step_increase = c(0.05,0.05,0.05),tip_length = 0, map_signif_level = T) +
    geom_point(aes(fill=spread), size=2, shape=21, colour="grey20",
               position=position_jitter(width=0.2, height=0)) +
    theme_bw()+ rotate_x_text() +
    theme(plot.title = element_text(face = "bold")) + theme(axis.title = element_text(size = 14), 
                                                            axis.text.x = element_text(size = 12), 
                                                            axis.text.y = element_text(size = 12), 
                                                            plot.title = element_text(size = 16)) +labs(x = NULL)

ggsave(filename = "output/microbial_spread/plot_taxonmetadata_FP_spread.ps", device = "ps", width = 11, height=5)


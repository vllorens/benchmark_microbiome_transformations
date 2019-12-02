# Alpha diversity on the simulated matrices
# Mon Dec  2 14:47:38 2019 ------------------------------

# In this script, we calculate different alpha diversity metrics on the original matrices, 
# as well as on the matrices transformed by sequencing, RMP, QMP and QMP-NR

#### Configure the environment ####

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(vegan))

# If it has not been run yet, it is necessary to run the script to generate the data (scripts/generate_data.R)
# source("scripts/generate_data.R")

# create folders, if not existing
system("mkdir -p output/alpha_div")

# clean these folders and remove all previous results (if any) 
system("rm -r output/alpha_div/*")

#### Calculate alpha diversity ####
richness_counts_all <- tibble()
# We first calculate some diversity indices in the simulated matrices
for(file in list.files("data/tax_matrices", full.names = T)){
    # read file 
    filename <- basename(file)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        round
    matrix_name <- strsplit(filename, split="_")[[1]][2]
    spread_name <- strsplit(filename, split="_")[[1]][3]
    scenario_name <- strsplit(filename, split="_")[[1]][5] %>% gsub(., pattern="\\.tsv", replacement="")
    # read total counts file
    original_file <- file
    original_counts <- read.table(original_file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(.,2,sum)
    spread <- range(original_counts)[2]/range(original_counts)[1]
    # make phyloseq object
    abundances <- otu_table(tax_matrix, taxa_are_rows = T)
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>% 
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="Real")
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # read seq file
    seq_matrix <- read.table(paste0("data/seq_matrices/seqOut_", filename), 
                             header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    abundances <- otu_table(seq_matrix, taxa_are_rows = T)
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson",  "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>%
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="Seq")
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # read RMP file
    rmp_matrix <- read.table(paste0("data/rmp_matrices/RMP_", filename), 
                             header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    abundances <- otu_table(rmp_matrix, taxa_are_rows = T)
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson",  "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>%
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="RMP") 
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # read QMP file
    qmp_matrix <- read.table(paste0("data/qmp_matrices/QMP_", filename), 
                             header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% round()
    abundances <- otu_table(qmp_matrix, taxa_are_rows = T)
    richness_counts_table <- estimate_richness(abundances,measures=c("Observed", "Chao1", "Simpson",  "Shannon"))%>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>%
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="QMP") 
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # read QMP-NR file
    qmp_matrix <- read.table(paste0("data/qmp2_matrices/QMP2_", filename), 
                             header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% round()
    abundances <- otu_table(qmp_matrix, taxa_are_rows = T)
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson",  "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>%
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="QMP-NR") 
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
}


#### Assemble results and plot (Figure 1b) ####
richness_counts_all$spread <- factor(richness_counts_all$spread, 
                                     levels=c("low", "high"))
richness_counts_all$Diversity_index <- factor(richness_counts_all$Diversity_index, 
                                              levels=c("Observed", "Chao1", "Simpson", "Shannon"))
richness_counts_all$scenario <- factor(richness_counts_all$scenario,
                                       levels=c("Healthy", "Blooming", "Dysbiosis"))
richness_counts_all$method <- factor(richness_counts_all$method, 
                                     levels=c("Real", "Seq", "RMP", "QMP", "QMP-NR"))
richness_counts_all$matrix <- factor(richness_counts_all$matrix, 
                                     levels=c("S1", "S2", "S3", "S4", "S5", "S6",
                                              "S7", "S8", "S9", "S10",
                                              "D1", "D2", "D3", "D4", "D5", "D6",
                                              "D7", "D8", "D9", "D10", 
                                              "B1", "B2", "B3", "B4", "B5", "B6", 
                                              "B7", "B8", "B9", "B10"))

print(ggscatter(richness_counts_all %>% dplyr::filter(Diversity_index=="Observed"), 
                x="counts", y="Value", group="matrix", facet.by=c("matrix"), nrow=3,
                fill="method", alpha=0.5, shape=21, color="black",
                title="Observed richness vs counts", palette="Spectral",
                ylab="Observed richness", xlab="Actual cell counts",
                add = "reg", add.params = list(color = "method")) + 
          xscale("log10") + 
          theme_bw())
ggsave("output/alpha_div/alphadiv_observed.pdf", device = "pdf")
print(ggscatter(richness_counts_all %>% dplyr::filter(Diversity_index=="Chao1"), 
                x="counts", y="Value", group="matrix", facet.by=c("matrix"), nrow=3,
                fill="method", alpha=0.5, shape=21, color="black",
                title="Chao1 diversity vs counts", palette="Spectral",
                ylab="Chao1 diversity", xlab="Actual cell counts",
                add = "reg", add.params = list(color = "method")) + 
          xscale("log10") + 
          theme_bw())
ggsave("output/alpha_div/alphadiv_chao1.pdf", device = "pdf")
print(ggscatter(richness_counts_all %>% dplyr::filter(Diversity_index=="Simpson"), 
                x="counts", y="Value", group="matrix", facet.by=c("matrix"), nrow=3,
                fill="method", alpha=0.5, shape=21, color="black",
                title="Simpson diversity vs counts", palette="Spectral",
                ylab="Simpson's diversity", xlab="Actual cell counts",
                add = "reg", add.params = list(color = "method")) + 
          xscale("log10") + 
          theme_bw())
ggsave("output/alpha_div/alphadiv_Simpson.pdf", device = "pdf")
print(ggscatter(richness_counts_all %>% dplyr::filter(Diversity_index=="Shannon"), 
                x="counts", y="Value", group="matrix", facet.by=c("matrix"), nrow=3,
                fill="method", alpha=0.5, shape=21, color="black",
                title="Shannon diversity vs counts", palette="Spectral",
                ylab="Shannon's diversity", xlab="Actual cell counts",
                add = "reg", add.params = list(color = "method")) + 
          xscale("log10") + 
          theme_bw())
ggsave("output/alpha_div/alphadiv_Shannon.pdf", device = "pdf")

# example for the figure
richness_ex <- richness_counts_all %>% dplyr::filter(matrix %in% c("S1", "D4", "B1")) %>% 
    dplyr::filter(Diversity_index %in% c("Observed", "Chao1", "Simpson", "Shannon"))
richness_ex$method <- gsub(richness_ex$method, pattern="Seq", replacement="Sequencing")
richness_ex$method <- factor(richness_ex$method, 
                             levels=c("Real", "Sequencing", "RMP", "QMP", "QMP-NR"))
print(ggscatter(richness_ex,
                x="counts", y="Value", group="matrix", ylab="Alpha diversity index value",
                fill="method", alpha=0.9, shape=21, color="black",
                scales="free", palette="Spectral", xlab="Absolute cell counts",
                add = "reg", add.params = list(color = "method"),
                facet.by=c("Diversity_index", "scenario"), 
                title="Alpha diversity indices in simulated scenarios") + 
          xscale("log10") + 
          theme_bw() + theme(plot.title = element_text(face = "bold")))
ggsave("output/alpha_div/alphadiv_examples_fig1b.pdf",device = "pdf" , width=7, height=7)

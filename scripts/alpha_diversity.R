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
suppressPackageStartupMessages(library(CoDaSeq))
suppressPackageStartupMessages(library(metagenomeSeq))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))

# load functions
source("R/GMPR.R")

# If it has not been run yet, it is necessary to run the script to generate the data (scripts/generate_data.R)
# source("scripts/generate_data.R")

# create folders, if not existing
system("mkdir -p output/alpha_div")

# clean these folders and remove all previous results (if any) 
system("rm -r output/alpha_div/*")

set.seed(777)

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
    # # CLR transformation (omitted as negative values don't make sense)
    # tax_matrix <- t(zCompositions::cmultRepl(X = t(seq_matrix), output="counts")) # estimates zeros
    # tax_matrix <- t(codaSeq.clr(tax_matrix, samples.by.row=F)) %>% round()
    # abundances <- otu_table(tax_matrix, taxa_are_rows = T)
    # richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1")) %>% 
    #     as_tibble() %>% 
    #     mutate(samples=colnames(abundances)) %>% 
    #     mutate(counts=original_counts) %>% 
    #     gather(key = "Diversity_index", value = "Value", -samples, -counts) %>% 
    #     mutate(matrix=matrix_name) %>% 
    #     mutate(spread=spread_name) %>% 
    #     mutate(scenario=scenario_name) %>% 
    #     mutate(method="CLR")
    # richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # GMPR transformation
    gmpr.size.factor <- GMPR(comm=seq_matrix)
    tax_matrix <- sweep(seq_matrix, MARGIN = 2, gmpr.size.factor, '/')
    abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>% 
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="GMPR")
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # CSS - metagenomeSeq transformation
    tax_matrix <- newMRexperiment(seq_matrix)
    p <- cumNormStatFast(tax_matrix)
    tax_matrix <- cumNorm(tax_matrix, p = p)
    tax_matrix <- MRcounts(tax_matrix, norm=T)
    abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>% 
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="CSS")
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # TMM - edgeR transformation
    tax_matrix <- DGEList(counts=seq_matrix)
    tax_matrix <- calcNormFactors(tax_matrix, method = "TMM")
    tax_matrix <- cpm(tax_matrix)
    abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>% 
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="TMM")
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # UQ - edgeR transformation
    tax_matrix <- DGEList(counts=seq_matrix)
    tax_matrix <- calcNormFactors(tax_matrix, method = "upperquartile")
    tax_matrix <- cpm(tax_matrix)
    abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>% 
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="UQ")
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # RLE - edgeR transformation
    tax_matrix <- DGEList(counts=seq_matrix)
    tax_matrix <- calcNormFactors(tax_matrix, method = "RLE")
    tax_matrix <- cpm(tax_matrix)
    abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts) %>% 
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="RLE")
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # # AST transformation (omitted as only values between 0 and 1 returned)
    # tax_matrix = apply(seq_matrix, 2, function(x){x/sum(x)}) # TSS
    # tax_matrix <- asin(sqrt(tax_matrix)) # AST
    # abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
    # richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
    #     as_tibble() %>% 
    #     mutate(samples=colnames(abundances)) %>% 
    #     mutate(counts=original_counts) %>% 
    #     gather(key = "Diversity_index", value = "Value", -samples, -counts) %>% 
    #     mutate(matrix=matrix_name) %>% 
    #     mutate(spread=spread_name) %>% 
    #     mutate(scenario=scenario_name) %>% 
    #     mutate(method="AST")
    # richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # # VST - DESeq2 (omitted as negative values don't make sense)
    # tax_matrix <- varianceStabilizingTransformation(seq_matrix)
    # abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
    # richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1")) %>% 
    #     as_tibble() %>% 
    #     mutate(samples=colnames(abundances)) %>% 
    #     mutate(counts=original_counts) %>% 
    #     gather(key = "Diversity_index", value = "Value", -samples, -counts) %>% 
    #     mutate(matrix=matrix_name) %>% 
    #     mutate(spread=spread_name) %>% 
    #     mutate(scenario=scenario_name) %>% 
    #     mutate(method="VST")
    # richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
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
# first separate the "real" alpha diversity to use as reference
richness_counts_all_2 <- spread(richness_counts_all, key="method", value="Value") %>% 
    dplyr::filter(Diversity_index!="se.chao1")
richness_counts_all_2 <- gather(richness_counts_all_2, key="method", value="Value", 
                              -c("samples", "counts", "Diversity_index", "matrix",
                                 "spread", "scenario", "Real"))
richness_counts_all_2 <- richness_counts_all_2 %>% drop_na()

# fit to the real data
fit_summary <- tibble()
for(mat in unique(richness_counts_all_2$matrix)){
    temp_matrix <- richness_counts_all_2 %>% 
        dplyr::filter(matrix==mat)
    for(met in unique(temp_matrix$method)){
        temp_matrix2 <- temp_matrix %>% 
            dplyr::filter(method==met)
        for(div in unique(temp_matrix2$Diversity_index)){
            temp_matrix3 <- temp_matrix2 %>% 
                dplyr::filter(Diversity_index==div)
            tib_temp <- tibble(Diversity_index=div, matrix=mat, method=met,
                               scenario=unique(temp_matrix3$scenario),
                               fit=coef(lm(Real ~ Value, data=temp_matrix3))[2],
                               rho=cor(temp_matrix3$Real, temp_matrix3$Value, method="spearman"))
            fit_summary <- bind_rows(fit_summary,tib_temp)
        }
    }
}

method_type <- tibble(method=c("Seq", "RMP", "CSS", "GMPR",
                               "UQ", "RLE", "TMM","QMP", "QMP-NR"),
                      method_type=c("Sequencing", "Traditional transformations", 
                                    rep("Compositional transformations", times=5),
                                    rep("Quantitative transformations", times=2)))
fit_summary <- fit_summary %>% 
    left_join(method_type, by="method")

fit_summary$method <- factor(fit_summary$method, 
                                     levels=c("Seq", "RMP", "CSS", "GMPR",
                                              "UQ", "RLE", "TMM", "QMP", "QMP-NR"))
fit_summary$method_type <- factor(fit_summary$method_type, 
                             levels=c("Sequencing", "Traditional transformations", 
                                      "Compositional transformations",
                                      "Quantitative transformations"))
fit_summary$scenario <- factor(fit_summary$scenario,
                                       levels=c("Healthy", "Blooming", "Dysbiosis"))
fit_summary$Diversity_index <- factor(fit_summary$Diversity_index, 
                                              levels=c("Observed", "Chao1", "Simpson", "Shannon"))


fig1b <- ggerrorplot(fit_summary, x="method", y="rho", color="method_type", 
          facet.by = c("Diversity_index", "scenario"), desc_stat="mean_ci", size=0.3,
          fill="method_type", 
          error.plot = "pointrange", palette=get_palette("Spectral",11)[c(2,5,8,10)], legend.title="Method type",
          title="Correlation of alpha diversity metrics with the real values") + 
    theme_bw() + 
    geom_hline(yintercept = 1) + rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
    panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
    plot.title = element_text(size = 14, 
        face = "bold"), legend.title = element_text(face = "bold")) +
    labs(x = "Method")

fig1b
ggsave(fig1b, filename = "output/alpha_div/correlation_alpha_div_real.ps", device="ps", height=8, width=8)



# For supplementary
method_type <- tibble(method=c("Real", "Seq", "RMP", "CSS", "GMPR",
                               "UQ", "RLE", "TMM","QMP", "QMP-NR"),
                      method_type=c("No transformation", "Sequencing", 
                                    "Traditional transformations", 
                                    rep("Compositional transformations", times=5),
                                    rep("Quantitative transformations", times=2)))
richness_counts_all <- richness_counts_all %>% 
    left_join(method_type, by="method")
richness_counts_all$spread <- factor(richness_counts_all$spread, 
                                     levels=c("low", "high"))
richness_counts_all$Diversity_index <- factor(richness_counts_all$Diversity_index, 
                                              levels=c("Observed", "Chao1", "Simpson", "Shannon"))
richness_counts_all$scenario <- factor(richness_counts_all$scenario,
                                       levels=c("Healthy", "Blooming", "Dysbiosis"))
richness_counts_all$method <- factor(richness_counts_all$method, 
                                     levels=c("Real", "Seq", "RMP", "CSS", "GMPR",
                                              "UQ", "RLE", "TMM", "QMP", "QMP-NR"))
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
                ylab="Observed richness", xlab="Microbial load",
                add = "reg", add.params = list(color = "method")) + 
          xscale("log10") + 
          theme_bw())+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))
ggsave("output/alpha_div/alphadiv_observed.ps", device = "ps", width=11, height=5)
print(ggscatter(richness_counts_all %>% dplyr::filter(Diversity_index=="Chao1"), 
                x="counts", y="Value", group="matrix", facet.by=c("matrix"), nrow=3,
                fill="method", alpha=0.5, shape=21, color="black",
                title="Chao1 diversity vs counts", palette="Spectral",
                ylab="Chao1 diversity", xlab="Microbial load",
                add = "reg", add.params = list(color = "method")) + 
          xscale("log10") + 
          theme_bw())+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))
ggsave("output/alpha_div/alphadiv_chao1.ps", device = "ps", width=11, height=5)
print(ggscatter(richness_counts_all %>% dplyr::filter(Diversity_index=="Simpson"), 
                x="counts", y="Value", group="matrix", facet.by=c("matrix"), nrow=3,
                fill="method", alpha=0.5, shape=21, color="black",
                title="Simpson diversity vs counts", palette="Spectral",
                ylab="Simpson's diversity", xlab="Microbial load",
                add = "reg", add.params = list(color = "method")) + 
          xscale("log10") + 
          theme_bw())+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))
ggsave("output/alpha_div/alphadiv_Simpson.ps", device = "ps", width=11, height=5)
print(ggscatter(richness_counts_all %>% dplyr::filter(Diversity_index=="Shannon"), 
                x="counts", y="Value", group="matrix", facet.by=c("matrix"), nrow=3,
                fill="method", alpha=0.5, shape=21, color="black",
                title="Shannon diversity vs counts", palette="Spectral",
                ylab="Shannon's diversity",xlab="Microbial load",
                add = "reg", add.params = list(color = "method")) + 
          xscale("log10") + 
          theme_bw())+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold"))
ggsave("output/alpha_div/alphadiv_Shannon.ps", device = "ps", width=11, height=5)

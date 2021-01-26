# Alpha diversity on the simulated matrices
# Mon Dec  2 14:47:38 2019 ------------------------------

# In this script, we calculate different alpha diversity metrics on the original matrices, 
# as well as on the matrices transformed by sequencing, RMP, QMP and ACS

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
suppressPackageStartupMessages(library(ggtext))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(inlmisc))

# load functions
source("R/GMPR.R")
source("R/plot_comparisons.R")

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
    seq_matrix <- read.table(paste0("data/seq_matrices/seqOut_", filename), 
                             header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    seq_counts <- apply(seq_matrix, 2, sum)
    sampling_depth <- 100*seq_counts/original_counts
    # make phyloseq object
    abundances <- otu_table(tax_matrix, taxa_are_rows = T)
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        mutate(sam_depth=sampling_depth) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
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
        mutate(sam_depth=sampling_depth) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
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
        mutate(sam_depth=sampling_depth) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
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
        mutate(sam_depth=sampling_depth) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
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
        mutate(sam_depth=sampling_depth) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="TMM")
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # UQ - edgeR transformation
    tax_matrix <- DGEList(counts=seq_matrix)
    tax_matrix <- calcNormFactors(tax_matrix, method = "upperquartile", p=0.75)
    if(length(which(is.infinite(tax_matrix$samples$norm.factors)))==0){
        #whenever the distrubution is too uneven after sequencing, normalization factors cannot be calculated, so we skip the analysis in those cases
        tax_matrix <- cpm(tax_matrix)
        abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
        richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
            as_tibble() %>% 
            mutate(samples=colnames(abundances)) %>% 
            mutate(counts=original_counts) %>% 
            mutate(sam_depth=sampling_depth) %>% 
            gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
            mutate(matrix=matrix_name) %>% 
            mutate(spread=spread_name) %>% 
            mutate(scenario=scenario_name) %>% 
            mutate(method="UQ")
        richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    }
    # RLE - edgeR transformation
    tax_matrix <- DGEList(counts=seq_matrix)
    tax_matrix <- calcNormFactors(tax_matrix, method = "RLE")
    tax_matrix <- cpm(tax_matrix)
    abundances <- otu_table(tax_matrix, taxa_are_rows = T) %>% round()
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson", "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        mutate(sam_depth=sampling_depth) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
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
        mutate(sam_depth=sampling_depth) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
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
        mutate(sam_depth=sampling_depth) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="QMP") 
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
    # read ACS file
    acs_matrix <- read.table(paste0("data/acs_matrices/ACS_", filename), 
                             header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% round()
    abundances <- otu_table(acs_matrix, taxa_are_rows = T)
    richness_counts_table <- estimate_richness(abundances, measures=c("Observed", "Chao1", "Simpson",  "Shannon")) %>% 
        as_tibble() %>% 
        mutate(samples=colnames(abundances)) %>% 
        mutate(counts=original_counts) %>% 
        mutate(sam_depth=sampling_depth) %>% 
        gather(key = "Diversity_index", value = "Value", -samples, -counts, -sam_depth) %>% 
        mutate(matrix=matrix_name) %>% 
        mutate(spread=spread_name) %>% 
        mutate(scenario=scenario_name) %>% 
        mutate(method="ACS") 
    richness_counts_all <- bind_rows(richness_counts_all, richness_counts_table)
}


#### Assemble results ####
# first separate the "real" alpha diversity to use as reference
richness_counts_all_2 <- spread(richness_counts_all, key="method", value="Value") %>% 
    dplyr::filter(Diversity_index!="se.chao1")
richness_counts_all_2 <- gather(richness_counts_all_2, key="method", value="Value", 
                              -c("samples", "counts", "sam_depth", "Diversity_index", "matrix",
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
                               fit=coef(lm(Value ~ Real, data=temp_matrix3))[2],
                               rsq=summary(lm(Value ~ Real, data=temp_matrix3))[[9]],
                               rho=cor(temp_matrix3$Real, temp_matrix3$Value, method="pearson"),
                               rho_p=cor.test(temp_matrix3$Real, temp_matrix3$Value, method="pearson")[[3]])
            
            fit_summary <- bind_rows(fit_summary,tib_temp)
        }
    }
}

method_type <- tibble(method=c("Seq", "RMP", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "ACS", "QMP"),
                      method_type=c("Sequencing", "Relative transformations", 
                                    rep("Compositional transformations", times=5),
                                    rep("Quantitative transformations", times=2)))
fit_summary <- fit_summary %>% 
    left_join(method_type, by="method")

fit_summary$method <- factor(fit_summary$method, 
                                     levels=c("Seq", "RMP", "CSS", "GMPR",
                                              "UQ", "RLE", "TMM", "ACS", "QMP"))
fit_summary$method_type <- factor(fit_summary$method_type, 
                             levels=c("Sequencing", "Relative transformations", 
                                      "Compositional transformations",
                                      "Quantitative transformations"))
fit_summary$scenario <- factor(fit_summary$scenario,
                                       levels=c("Healthy", "Blooming", "Dysbiosis"))
fit_summary$Diversity_index <- factor(fit_summary$Diversity_index, 
                                              levels=c("Observed", "Chao1", "Simpson", "Shannon"))

fit_summary <- fit_summary %>% 
    mutate(correlation_category="High correlation (R>0.8)")
fit_summary[fit_summary$rho<0.8 & fit_summary$rho_p<0.05, "correlation_category"] <- "Moderate correlation (R>0.5)"
fit_summary[fit_summary$rho<0.5 & fit_summary$rho_p<0.05, "correlation_category"] <- "Mild correlation (R<0.5)"
fit_summary[fit_summary$rho<0 & fit_summary$rho_p<0.05, "correlation_category"] <- "Negative correlation (R<0)"
fit_summary[fit_summary$rho_p>=0.05, "correlation_category"] <- "Non-significant correlation"


## statistics on method vs method comparison
kruskal_rho <- fit_summary %>% 
    group_by(scenario, Diversity_index) %>% 
    rstatix::kruskal_test(rho ~ method)

dunn_rho <- fit_summary %>% 
    group_by(scenario, Diversity_index) %>% 
    rstatix::dunn_test(rho ~ method)


## summary figure observed richness (Figure 2)
toplot <- fit_summary %>% 
    dplyr::select(Diversity_index, method, scenario, correlation_category) %>% 
    table %>%
    as_tibble() %>% 
    group_by(method, scenario, Diversity_index) %>% mutate(percentage=100*n/sum(n)) %>% 
    ungroup 


toplot$scenario <- factor(toplot$scenario,
                          levels=c("Healthy", "Blooming", "Dysbiosis"))
toplot$Diversity_index <- factor(toplot$Diversity_index, 
                                 levels=c("Observed", "Chao1", "Simpson", "Shannon"))
toplot$correlation_category <- factor(toplot$correlation_category,
                                      levels=c("High correlation (R>0.8)",
                                               "Moderate correlation (R>0.5)",
                                               "Mild correlation (R<0.5)",
                                               "Non-significant correlation",
                                               "Negative correlation (R<0)"))

fig2a <- ggbarplot(toplot %>% dplyr::filter(Diversity_index=="Observed"), 
                   x="method", y="percentage", color="white", 
                   fill="correlation_category", 
                   palette=rev(inlmisc::GetTolColors(5, scheme = "sunset")), 
                   title="Observed richness: correlation between transformed values and synthetic communities",
                   legend.title="Correlation category", xlab="Method",
                   ylab="Percentage",
                   facet.by="scenario")+
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(face = "bold"))

legend_2a <- get_legend(fig2a)
fig2a <- fig2a + rremove("legend")


dunn1 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Healthy", Diversity_index=="Observed"), title="Healthy succession")
dunn2 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Blooming", Diversity_index=="Observed"), title="Blooming")
dunn3 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Dysbiosis", Diversity_index=="Observed"), title="Dysbiosis")
fig2b <- ggarrange(dunn1, dunn2, dunn3, nrow=1, ncol=3, widths=c(1,1,1), legend="top")


legend_2b <- ggplot() + theme_void()
plot1 <- ggarrange(fig2a, fig2b, nrow=2, align="hv", heights=c(1,1.3), labels = c("a", "b"))
plot2 <- ggarrange(legend_2a, legend_2b, nrow=2, heights=c(1,1.3),  align="hv")
plot_fin <- ggarrange(plot1, plot2, nrow=1, ncol=2, widths=c(3.5,1))


ggsave(plot_fin, filename = "output/alpha_div/correlation_observed_richness_real_transformed.pdf", device="pdf", height=8, width=9, useDingbats=FALSE)



#### Additional figures (some of these go to supp material, "as is" or with reorganized panels) ####
richness_counts_all_2 <- richness_counts_all_2 %>% 
    left_join(method_type, by="method")

richness_counts_all_2$method <- factor(richness_counts_all_2$method, 
                                       levels=c("Seq", "RMP", "CSS", "GMPR",
                                                "UQ", "RLE", "TMM", "ACS", "QMP"))
richness_counts_all_2$method_type <- factor(richness_counts_all_2$method_type, 
                                            levels=c("Sequencing", "Relative transformations", 
                                                     "Compositional transformations",
                                                     "Quantitative transformations"))

richness_dysbiosis_real <- tibble()
for(file in list.files("data/tax_matrices/", pattern = "D", full.names = T)){
    mat <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        round
    matname <- strsplit(file, split="_")[[1]][3]
    total_counts <- apply(mat,2,sum)
    healthy <- which(total_counts >= median(total_counts))
    dysbiosis <- which(total_counts < median(total_counts))
    mat_healthy <- mat[,healthy]
    mat_healthy <- otu_table(mat_healthy, taxa_are_rows = T)
    rich_healthy <- estimate_richness(mat_healthy,measures=c("Observed")) %>% 
        as.tibble() %>% 
        mutate(cat="Healthy") %>% 
        mutate(mat=matname)
    mat_dysbiosis <- mat[,dysbiosis]
    mat_dysbiosis <- otu_table(mat_dysbiosis, taxa_are_rows = T)
    rich_dysbiosis <- estimate_richness(mat_dysbiosis,measures=c("Observed")) %>% 
        as.tibble() %>% 
        mutate(cat="Dysbiosis") %>% 
        mutate(mat=matname)
    richness_dysbiosis_real <- bind_rows(richness_dysbiosis_real, rich_healthy, rich_dysbiosis)
}

richness_dysbiosis_real$mat <- factor(richness_dysbiosis_real$mat,
                                      levels=c("D1", "D2", "D3", "D4", "D5",
                                               "D6", "D7", "D8", "D9", "D10"))

pdf("output/alpha_div/richness_healthy_vs_dysbiotic.pdf", height=6, width=9)
ggboxplot(richness_dysbiosis_real, x = "cat", y="Observed", facet.by="mat", nrow=2, fill="cat", 
          palette=inlmisc::GetTolColors(2, scheme = "sunset"),
          title="Richness in healthy vs dysbiotic samples", legend="right", legend.title="Category",
          xlab="Category", ylab="Observed richness", ylim=c(150,320)) + 
    theme_bw() + 
    stat_compare_means(comparisons = list(c("Healthy", "Dysbiosis")), label = "p.signif") + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(face = "bold"))
dev.off()

pdf("output/alpha_div/plots_diversity_correlation_slope_rsquared.pdf", height=8, width=11)
ggerrorplot(fit_summary, x="method", y="rho", color="method_type",  scales="free_y",
                     facet.by = c("Diversity_index", "scenario"), desc_stat="mean_ci", size=0.3,
                     fill="method_type", palette=inlmisc::GetTolColors(4, scheme = "sunset"),
                     error.plot = "pointrange",legend.title="Method type",
                     title="Pearson correlation coefficient: Method alpha diversity vs Real alpha diversity") + 
    theme_bw() + stat_compare_means(label.x = 3, label.y.npc = "bottom") + 
    geom_hline(yintercept = 1) + rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(face = "bold")) +
    labs(x = "Method")


ggerrorplot(fit_summary, x="method", y="rsq", color="method_type",  scales="free_y",
            facet.by = c("Diversity_index", "scenario"), desc_stat="mean_ci", size=0.3,
            fill="method_type",  palette=inlmisc::GetTolColors(4, scheme = "sunset"),
            error.plot = "pointrange",legend.title="Method type",
            title="R-squared value of the fit: Method alpha diversity vs Real alpha diversity") + 
    theme_bw() + stat_compare_means(label.x = 3, label.y.npc = "bottom") + 
    geom_hline(yintercept = 1) + rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(face = "bold")) +
    labs(x = "Method")

ggerrorplot(fit_summary, x="method", y="fit", color="method_type",  scales="free_y",
            facet.by = c("Diversity_index", "scenario"), desc_stat="mean_ci", size=0.3,
            fill="method_type",  palette=inlmisc::GetTolColors(4, scheme = "sunset"),
            error.plot = "pointrange", legend.title="Method type",
            title="Slope of the fit: Method alpha diversity vs Real alpha diversity") + 
    theme_bw() + stat_compare_means(label.x = 3, label.y.npc = "bottom") +
    geom_hline(yintercept = 1) + rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(face = "bold")) +
    labs(x = "Method")
dev.off()



## summary chao1
p1 <- ggbarplot(toplot %>% dplyr::filter(Diversity_index=="Chao1"), 
                x="method", y="percentage", fill="correlation_category", 
                color="white",
                palette=rev(inlmisc::GetTolColors(5, scheme = "sunset")),
                title="Chao1 richness: correlation between transformed and real values",
                legend.title="Correlation category", xlab="Method",
                ylab="Percentage",
                facet.by="scenario")+
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(face = "bold"))


legend_p1 <- get_legend(p1)
p1 <- p1 + rremove("legend")

p2 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Healthy", Diversity_index=="Chao1"), title="Healthy succession")
p3 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Blooming", Diversity_index=="Chao1"), title="Blooming")
p4 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Dysbiosis", Diversity_index=="Chao1"), title="Dysbiosis")
plotb <- ggarrange(p2, p3, p4, nrow=1, ncol=3, common.legend = F, legend="top")


plot1 <- ggarrange(p1, plotb, nrow=2, align="hv", heights=c(1,1.3), labels = c("a", "b"))
plot2 <- ggarrange(legend_p1, legend_2b, nrow=2, heights=c(1,1.3),  align="hv")
plot_fin <- ggarrange(plot1, plot2, nrow=1, ncol=2, widths=c(3.5,1))
ggsave(plot_fin, filename = "output/alpha_div/correlation_chao1_richness_real_transformed.pdf", device="pdf", height=8, width=9, useDingbats=FALSE)


## summary simpson
p1 <- ggbarplot(toplot %>% dplyr::filter(Diversity_index=="Simpson"), 
                x="method", y="percentage", fill="correlation_category", 
                color="white",
                palette=rev(inlmisc::GetTolColors(5, scheme = "sunset")),
                title="Simpson diversity: correlation between transformed and real values",
                legend.title="Correlation category", xlab="Method",
                ylab="Percentage",
                facet.by="scenario")+
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(face = "bold"))


legend_p1 <- get_legend(p1)
p1 <- p1 + rremove("legend")

p2 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Healthy", Diversity_index=="Simpson"), title="Healthy succession")
p3 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Blooming", Diversity_index=="Simpson"), title="Blooming")
p4 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Dysbiosis", Diversity_index=="Simpson"), title="Dysbiosis")
plotb <- ggarrange(p2, p3, p4, nrow=1, ncol=3, common.legend = F, legend="top")

plot1 <- ggarrange(p1, plotb, nrow=2, align="hv", heights=c(1,1.3), labels = c("a", "b"))
plot2 <- ggarrange(legend_p1, legend_2b, nrow=2, heights=c(1,1.3),  align="hv")
plot_fin <- ggarrange(plot1, plot2, nrow=1, ncol=2, widths=c(3.5,1))
ggsave(plot_fin, filename = "output/alpha_div/correlation_simpson_diversity_real_transformed.pdf", device="pdf", height=8, width=9, useDingbats=FALSE)


## summary shannon
p1 <- ggbarplot(toplot %>% dplyr::filter(Diversity_index=="Shannon"), 
                x="method", y="percentage", fill="correlation_category", 
                color="white",
                palette=rev(inlmisc::GetTolColors(5, scheme = "sunset")),
                title="Shannon diversity: correlation between transformed and real values",
                legend.title="Correlation category", xlab="Method",
                ylab="Percentage",
                facet.by="scenario")+
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(face = "bold"))


legend_p1 <- get_legend(p1)
p1 <- p1 + rremove("legend")

p2 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Healthy", Diversity_index=="Shannon"), title="Healthy succession")
p3 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Blooming", Diversity_index=="Shannon"), title="Blooming")
p4 <- plot_comparisons(dunn=dunn_rho %>% filter(scenario=="Dysbiosis", Diversity_index=="Shannon"), title="Dysbiosis")
plotb <- ggarrange(p2, p3, p4, nrow=1, ncol=3, common.legend = F, legend="top")

plot1 <- ggarrange(p1, plotb, nrow=2, align="hv", heights=c(1,1.3), labels = c("a", "b"))
plot2 <- ggarrange(legend_p1, legend_2b, nrow=2, heights=c(1,1.3),  align="hv")
plot_fin <- ggarrange(plot1, plot2, nrow=1, ncol=2, widths=c(3.5,1))
ggsave(plot_fin, filename = "output/alpha_div/correlation_shannon_diversity_real_transformed.pdf", device="pdf", height=8, width=9, useDingbats=FALSE)



## observed richness vs microbial load - example
method_type <- tibble(method=c("Real", "Seq", "RMP", "CSS", "GMPR",
                               "UQ", "RLE", "TMM", "ACS", "QMP"),
                      method_type=c("No transformation", "Sequencing", 
                                    "Traditional transformations", 
                                    rep("Compositional transformations", times=5),
                                    rep("Transformations incorporating microbial loads", times=2)))
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
                                              "UQ", "RLE", "TMM", "ACS", "QMP"))
richness_counts_all$matrix <- factor(richness_counts_all$matrix, 
                                     levels=c("S1", "S2", "S3", "S4", "S5", "S6",
                                              "S7", "S8", "S9", "S10",
                                              "D1", "D2", "D3", "D4", "D5", "D6",
                                              "D7", "D8", "D9", "D10", 
                                              "B1", "B2", "B3", "B4", "B5", "B6", 
                                              "B7", "B8", "B9", "B10"))

cor_richness_counts <- richness_counts_all %>% 
    dplyr::filter(Diversity_index!="se.chao1") %>% 
    group_by(method, matrix, Diversity_index) %>% 
    mutate(cor_counts_diversity=cor(counts, Value)) %>% 
    mutate(corp_counts_diversity=cor.test(counts, Value)[[3]]) %>% 
    ungroup %>% 
    dplyr::select(Diversity_index, matrix, scenario, method, cor_counts_diversity, corp_counts_diversity) %>% 
    distinct() %>% 
    mutate(correlation_category="Non-significant correlation")


cor_richness_counts[cor_richness_counts$cor_counts_diversity>0.8 & 
                        cor_richness_counts$corp_counts_diversity<0.05,"correlation_category"] <- "High correlation (R>0.8)"
cor_richness_counts[cor_richness_counts$cor_counts_diversity>0.5 & 
                        cor_richness_counts$corp_counts_diversity<0.05,"correlation_category"] <- "Moderate correlation (R>0.5)"
cor_richness_counts[cor_richness_counts$cor_counts_diversity<0.5 & 
                        cor_richness_counts$corp_counts_diversity<0.05,"correlation_category"] <- "Mild correlation (R<0.5)"
cor_richness_counts[cor_richness_counts$cor_counts_diversity<0 & 
                        cor_richness_counts$corp_counts_diversity<0.05,"correlation_category"] <- "Negative correlation (R<0)"

toplot <- cor_richness_counts %>% 
    dplyr::select(Diversity_index, method, scenario, correlation_category) %>% 
    table %>%
    as_tibble() %>% 
    group_by(method, scenario, Diversity_index) %>% 
    mutate(percentage=100*n/sum(n)) %>% 
    ungroup 


toplot$scenario <- factor(toplot$scenario,
                          levels=c("Healthy", "Blooming", "Dysbiosis"))
toplot$Diversity_index <- factor(toplot$Diversity_index, 
                                 levels=c("Observed", "Chao1", "Simpson", "Shannon"))
toplot$correlation_category <- factor(toplot$correlation_category,
                                      levels=c("High correlation (R>0.8)",
                                               "Moderate correlation (R>0.5)",
                                               "Mild correlation (R<0.5)",
                                               "Non-significant correlation",
                                               "Negative correlation (R<0)"))

toplot$method <- factor(toplot$method, levels=c("Real", "Seq", "RMP", "CSS", "GMPR",
                                               "UQ", "RLE", "TMM", "ACS", "QMP"))


ps1 <- ggbarplot(toplot, 
          x="method", y="percentage", fill="correlation_category", 
          color="white",
          palette=rev(inlmisc::GetTolColors(5, scheme = "sunset"))[2:5],
          title="Correlation between alpha diversity and microbial loads in real and transformed data",
          legend.title="Correlation category", xlab="Real data | Transformation method",
          ylab="Percentage",
          facet.by=c("Diversity_index", "scenario"))+
    theme_bw()+ rotate_x_text(45) +
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(face = "bold"))

ggsave(ps1, filename = "output/alpha_div/correlation_alphdiv_microbialload.pdf", device="pdf", height=10, width=11)


richness_example <- richness_counts_all_2 %>% dplyr::filter(Diversity_index=="Observed", matrix %in% c("S1", "D1", "B1"))
richness_example$matrix <- gsub(richness_example$matrix, pattern="S1", replacement="Matrix S1 (Healthy succession)")
richness_example$matrix <- gsub(richness_example$matrix, pattern="B1", replacement="Matrix B1 (Blooming)")
richness_example$matrix <- gsub(richness_example$matrix, pattern="D1", replacement="Matrix D1 (Dysbiosis)")

richness_example$matrix <- factor(richness_example$matrix, levels=c("Matrix S1 (Healthy succession)", "Matrix B1 (Blooming)", "Matrix D1 (Dysbiosis)"))
ggscatter(richness_example, 
                x="counts", y="Value", group="matrix", facet.by=c("matrix", "method"),
                fill="method_type", alpha=0.5, shape=21, color="black",
                title="Observed richness vs sampling depth", 
                palette=(inlmisc::GetTolColors(4, scheme = "sunset")),
                ylab="Observed richness", xlab="Cell densities", legend.title="Method type",
                add = "reg", add.params = list(color = "method_type")) + 
          xscale("log10") + 
          theme_bw() + rotate_x_text(45) +
    geom_point(aes(x=counts, y=Real), fill="gray", alpha=0.5) +
    geom_smooth(aes(x=counts, y=Real),method = "lm", se = FALSE, color="gray60")+
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), 
          legend.title = element_text(face = "bold")) + xscale(.scale="none", .format=T)
ggsave("output/alpha_div/alphadiv_observed_example.pdf", device = "pdf", width=16, height=6, useDingbats=FALSE)




dunn_rho <- dunn_rho %>% dplyr::filter(!(Diversity_index=="Observed" & scenario=="Dysbiosis"))
write_tsv(dunn_rho, path="output/alpha_div/dunn_tests.txt")

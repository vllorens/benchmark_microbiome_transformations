# Create simulations
# Mon Dec  2 14:08:31 2019 ------------------------------

# In this script, first we obtain individual simulated taxonomy matrices. 
# From the data object with all matrices, we write txt files with the original matries and we also produce the sequencing output, the RMP matrices, and the QMP matrices.
# Finally, we also generate metadata for all the original taxonomy matrices

#### Configure the environment ####

# load required packages
suppressPackageStartupMessages(library(CoDaSeq))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(inlmisc))

# load required functions
source("R/metadata.R")
source("R/rarefy_even_sampling_depth.R")
source("R/counts_taxa_correlation.R")

set.seed(14102019)

# create folders, if not existing
system("mkdir -p data/tax_matrices")
system("mkdir -p data/seq_matrices")
system("mkdir -p data/rmp_matrices")
system("mkdir -p data/qmp_matrices")
system("mkdir -p data/acs_matrices")
system("mkdir -p data/metadata_matrices")
system("mkdir -p data/counts_estimated")
system("mkdir -p output/rawdataplots")

# clean these folders and remove all previous results (if any) 
system("rm -r data/tax_matrices/*")
system("rm -r data/seq_matrices/*")
system("rm -r data/rmp_matrices/*")
system("rm -r data/qmp_matrices/*")
system("rm -r data/acs_matrices/*")
system("rm -r data/counts_estimated/*")
system("rm -r data/metadata_matrices/*")
system("rm -r output/rawdataplots/*")

# remove all objects in the environment
rm(list = setdiff(ls(), lsf.str()))


#### Generate original matrices ####

# Load simulated matrices
load("data/raw/20200707_sims_Pedro_v5.2.Rdata")
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
    # include a random subset of 200 samples per matrix (of a total of 1000 simulated)
    ind_matrix <- ind_matrix[,sample(colnames(ind_matrix), replace = F, size = 200)]
    spread <- (apply(ind_matrix,2,sum) %>% max)/(apply(ind_matrix,2,sum) %>% min) 
    print(scenario)
    print(spread)
    print((apply(ind_matrix,2,sum) %>% max))
    print((apply(ind_matrix,2,sum) %>% min))
    if(spread<20){
        spread_cat <- "low"
    } else if(spread >= 20 & spread < 40){
        spread_cat <- "medium"
    }else{
        spread_cat <- "high"
    }
    write.table(ind_matrix, paste0("data/tax_matrices/taxonomy_", mt, "_", spread_cat, "_spread_", scenario, ".tsv"), 
                col.names = T, row.names=T, quote=F, sep="\t")
}

#### Generate sequencing data, RMP, QMP and ACS matrices ####

# the sequencing data comes from the real simulated matrix
# sequencing matrices include an average of 30000 sequencing reads/counts per sample
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
# we rarefy to the # of reads of the sample with the lowest sequencing depth
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

# the QMP data comes from rarefaction to even sampling depth the sequencing data and then scaling to total counts
# i.e. first rarefying to even sampling depth

# first generate random counts (as estimated from flow cytometry) with a correlation ~90% with the actual counts from the simulation
cor_counts_min <- 0.85
cor_counts_max <- 0.95

# function to generate correlated data
# function from https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
complement <- function(y, rho, x) { 
    if (missing(x)) x <- rnorm(length(y)) 
    y.perp <- residuals(lm(x ~ y))
    rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}


for(file in list.files("data/seq_matrices", full.names = T)){
    filename <- basename(file) %>% 
        gsub(x=., pattern="seqOut_taxonomy", replacement="QMP_taxonomy")
    original_file <- paste0("data/tax_matrices/", gsub(filename, pattern="QMP_", replacement=""))
    # original counts
    counts <- read.table(original_file, header=T, stringsAsFactors=F, sep="\t") %>% 
        apply(., 2, sum) %>% 
        as.data.frame()
    # generate correlated counts
    isvalid <- F
    while(isvalid==F){
        corr_counts <- runif(1, cor_counts_min, cor_counts_max)
        logcounts <- log(counts[,1])
        counts_fc <- complement(y=logcounts, rho=corr_counts)
        # scale to ensure that range is the same
        counts_fc_scaled <- (counts_fc - min(counts_fc))/(max(counts_fc)-min(counts_fc)) * (max(logcounts) - min(logcounts)) + min(logcounts) 
        # ensure all counts are positive (sometimes negative numbers may be generated)
        counts_fc_scaled <- abs(counts_fc_scaled)
        # validate that correlation still holds
        correlation = cor.test(exp(counts_fc_scaled),counts[,1], exact = F)$estimate
        correlation2 <- cor.test(counts_fc_scaled,logcounts, exact = F)$estimate
        if(correlation>cor_counts_min & correlation2>cor_counts_min){
            isvalid <- T
        }
    }
    counts_fc_scaled <- exp(counts_fc_scaled)
    seqcounts <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        apply(., 2, sum) %>% 
        as.data.frame()
    counts_fc_scaled <- as.data.frame(counts_fc_scaled)
    rownames(counts_fc_scaled) <- rownames(counts)
    write.table(counts_fc_scaled, file = paste0("data/counts_estimated/countsEstimated_", filename), 
                quote = FALSE, row.names = TRUE, col.names=T, sep="\t")
}


# then scale the counts
for(file in list.files("data/seq_matrices", full.names = T)){
    filename <- basename(file) %>% 
        gsub(x=., pattern="seqOut_taxonomy", replacement="QMP_taxonomy")
    # read estimated counts
    counts <- read.table(paste0("data/counts_estimated/countsEstimated_", filename), header=T, stringsAsFactors=F, sep="\t") 
    
    seq <- read.table(file, header=T, stringsAsFactors=F, sep="\t")
    QMP <- rarefy_even_sampling_depth(cnv_corrected_abundance_table = seq, cell_counts_table = counts)
    write.table(QMP, file = paste0("data/qmp_matrices/", filename), 
                quote = FALSE, row.names = TRUE, col.names=T, sep="\t")
}

# the ACS data comes from multiplying sequencing data by a factor to get numbers proportional to the number of counts
# i.e. without first rarefying to even sampling depth
for(file in list.files("data/seq_matrices", full.names = T)){
    filename <- basename(file) %>% 
        gsub(x=., pattern="seqOut_taxonomy", replacement="QMP_taxonomy")
    # read estimated counts
    counts <- read.table(paste0("data/counts_estimated/countsEstimated_", filename), header=T, stringsAsFactors=F, sep="\t") 
    filename <- basename(file) %>% 
        gsub(x=., pattern="seqOut_taxonomy", replacement="ACS_taxonomy")
    seq <- read.table(file, header=T, stringsAsFactors=F, sep="\t")
    counts_seq <- apply(seq, 2, sum)
    factors <- counts/counts_seq
    ACS <- sweep(seq, MARGIN = 2, factors[,1], '*')
    write.table(ACS, file = paste0("data/acs_matrices/", filename), 
                quote = FALSE, row.names = TRUE, col.names=T, sep="\t")
}


#### Generate metadata ####
# First, define parameters for the metadata
rho_meta_min <- 0.3       # minimum correlation for the metadata variables
rho_meta_max <- 0.6        # max correlation for the metadata variables
num_metadata_noC <- 35     # metadata numeric variables not correlating
num_metadata_Ct <- 10      # metadata numeric variables correlating w one taxon
num_metadata_Cc <- 5       # metadata numeric variables correlating with total counts
cat_metadata_noC <- 35     # metadata categorical variables not correlating
cat_metadata_Ct <- 10      # metadata categorical variables correlating w one taxon
cat_metadata_Cc <- 5       # metadata categorical variables correlating with total counts

# Make metadata tables for all of the original matrices
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


#### Calculate matrix basic stats ####

# Matrix stats --> a reduced version of this with some columns removed is table S1
matrix_stats <- tibble()
stats <- as_tibble(stats_table)
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
    # read estimated counts file
    est_counts <- read.table(paste0("data/counts_estimated/countsEstimated_QMP_", filename), header=T, stringsAsFactors=F, sep="\t")
    correlation_counts <- cor(original_counts, est_counts[,1])
    # correlation of taxa w number of cell counts
    taxoncor <- counts_taxa_correlation(tax_matrix, original_counts)
    taxoncor[[2]] <- p.adjust(taxoncor[[2]], method="BH")
    
    not_cor <- length(which(taxoncor[[2]]>0.05))
    neg_cor <- length(which(taxoncor[[2]]<0.05 & taxoncor[[1]]<0))
    pos_cor <- length(which(taxoncor[[2]]<0.05 & taxoncor[[1]]>0))
    # specialtaxa
    sptax <- stats %>% 
        dplyr::filter(grepl(paste0(matrix_name, "_"), matrix)) %>% 
        pull(bloomer_or_dysbiotic_taxon)
    cor_special_effsize <- taxoncor[[1]][sptax]
    cor_special_pvalue <- taxoncor[[2]][sptax]
    # specialtaxa_2_onlydysbiosis
    if(scenario_name=="Dysbiosis"){
        sptax2 <- stats %>% 
            dplyr::filter(grepl(paste0(matrix_name, "_"), matrix)) %>% 
            pull(flat_taxa) %>% strsplit(., split=",") %>% unlist()
        falseflat <- names(which(taxoncor[[2]][sptax2]<0.05))
        sptax2 <- setdiff(sptax2, falseflat)
        cor_special_effsize2 <- taxoncor[[1]][sptax2] %>% mean()
        cor_special_pvalue2 <- taxoncor[[2]][sptax2] %>% mean()
        sptax2 <- paste(sptax2, collapse=",")
    } else{
        sptax2 <- NA
        cor_special_effsize2 <- NA
        cor_special_pvalue2 <- NA
    }
    
    # richness (number of species) and evenness (pielou's)
    richness_obs <- vegan::specnumber(t(tax_matrix))
    richness_shannnon <- vegan::diversity(t(tax_matrix), index = "shannon")
    evenness_pielou <- richness_shannnon/log(richness_obs)
    
    # read sequencing data file
    seq_matrix <- read.table(paste0("data/seq_matrices/seqOut_",filename), header=T, stringsAsFactors=F, sep="\t")
    seq_counts <- seq_matrix %>% 
        as.matrix() %>% 
        apply(.,2,sum)
    
    # sampling depth
    sampling_depth <- 100*seq_counts/original_counts
    
    # put everythiing together
    matrix_stats_current <- tibble(matrix=matrix_name, spread_type=spread_name,
                                   scenario=scenario_name, spread=spread, 
                                   correlation_counts=correlation_counts,
                                   
                                   non_correlated_taxa=not_cor,
                                   positively_correlated_taxa=pos_cor,
                                   negatively_correlated_taxa=neg_cor,
                                   
                                   special_taxon=sptax,
                                   flat_taxon_dysbiosis=sptax2,
                                   
                                   correlation_special_counts=cor_special_effsize,
                                   correlation_special_counts_p=cor_special_pvalue,
                                   correlation_flat_counts_dysbiosis=cor_special_effsize2,
                                   correlation_flat_counts_dysbiosis_p=cor_special_pvalue2,
                                   
                                   cell_density_median=formatC(median(original_counts), format="e", digits = 2),
                                   cell_density_range=paste(formatC(range(original_counts), format="e", digits = 2), collapse=" - "),
                                   
                                   richness_median=median(richness_obs),
                                   richness_range=paste(range(richness_obs), collapse=" - "),
                                   
                                   evenness_median=round(median(evenness_pielou),2),
                                   evenness_range=paste(round(range(evenness_pielou),2), collapse=" - "),
                                   
                                   seq_depth_median=median(seq_counts),
                                   seq_depth_range=paste(range(seq_counts), collapse=" - "),
                                   
                                   sampling_depth_median=median(sampling_depth),
                                   sampling_depth_range=paste(formatC(range(sampling_depth), format="e", digits = 2), collapse=" - ")
                                   )
    matrix_stats <- bind_rows(matrix_stats, matrix_stats_current)
}

write_tsv(matrix_stats, "data/raw/matrix_stats_3scenarios.tsv", col_names = T)

#### Calculate sampling depth for all samples ####
sequencing_stats <- tibble()
# Sampling depth is the number of reads divided by the total cell counts
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
    
    # read sequencing data file
    seq_matrix <- read.table(paste0("data/seq_matrices/seqOut_",filename), header=T, stringsAsFactors=F, sep="\t")
    seq_counts <- seq_matrix %>% 
        as.matrix() %>% 
        apply(.,2,sum)
    # sampling depth in %
    sampling_depth <- seq_counts*100/original_counts
    tib <- tibble(seq_counts, sampling_depth, matrix_name, spread_name, scenario_name)
    sequencing_stats <- bind_rows(sequencing_stats, tib)
}

gghistogram(sequencing_stats, x = "seq_counts", fill = "lightgray", add="median", 
            facet.by="scenario_name", title="Sequencing depth", xlab="Sequencing counts",
            ylab="Frequency") + theme_bw()

gghistogram(sequencing_stats, x = "sampling_depth", fill = "lightgray", add="median", 
            facet.by="scenario_name", title="Sampling depth", xlab="Sampling depth",
            ylab="Frequency") + theme_bw()

ggboxplot(sequencing_stats, x = "scenario_name", y="sampling_depth", fill = "lightgray", 
          title="Sampling depth", xlab="Scenario",
            ylab="Sampling depth") + theme_bw()

write_tsv(sequencing_stats, "data/raw/sequencing_stats_3scenarios.tsv", col_names = T)


# rarest (non-zero) taxa
rarest_abundance <- tibble()
for(file in list.files("data/tax_matrices", full.names = T)){
    # read file 
    filename <- basename(file)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        round
    matrix_name <- strsplit(filename, split="_")[[1]][2]
    spread_name <- strsplit(filename, split="_")[[1]][3]
    scenario_name <- strsplit(filename, split="_")[[1]][5] %>% gsub(., pattern="\\.tsv", replacement="")
    # rarest non-zero taxa
    tax_matrix[tax_matrix==0] <- NA
    rarest <- apply(tax_matrix, 2, which.min)
    rarest_original <- c()
    for(i in 1:ncol(tax_matrix)){
        rarest_original <- c(rarest_original, tax_matrix[rarest[i],i])
    }
    original_counts <- read.table(original_file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix() %>% 
        apply(.,2,sum)
    rarest_original_percentage <- rarest_original*100/original_counts
    # read sequencing data file
    seq_matrix <- read.table(paste0("data/seq_matrices/seqOut_",filename), header=T, stringsAsFactors=F, sep="\t")
    rarest_seq <- c()
    for(i in 1:ncol(tax_matrix)){
        rarest_seq <- c(rarest_seq, seq_matrix[rarest[i],i])
    }
    tib <- tibble(rarest_original, rarest_original_percentage, rarest_seq, matrix_name, scenario_name, spread_name)
    rarest_abundance <- bind_rows(rarest_abundance, tib)
}
rarest_abundance <- rarest_abundance %>% dplyr::filter(rarest_original>0)
gghistogram(rarest_abundance, x = "rarest_original", fill = "lightgray", add="median", 
            facet.by="scenario_name", title="Rarest taxon (non-zero) abundance", xlab="Rarest taxon abundance",
            ylab="Frequency") + theme_bw()
gghistogram(rarest_abundance, x = "rarest_original_percentage", fill = "lightgray", add="median", 
            facet.by="scenario_name", title="Rarest taxon (non-zero) abundance", xlab="Rarest taxon abundance (% over total community)",
            ylab="Frequency") + theme_bw()


#### random samples from succession scenario ####
tax_matrix <- read.table("data/tax_matrices/taxonomy_S1_low_spread_Healthy.tsv", header=T, stringsAsFactors=F, sep="\t") %>% 
    as.matrix() %>% 
    round

set.seed(42)
tax_matrix <- tax_matrix[,sample(1:ncol(tax_matrix),size=1, replace = F)] %>%
    as.tibble %>%
    mutate(taxon=rownames(tax_matrix)) %>% 
    gather(key = "sample", value="counts", -taxon)

tax_matrix$taxon <- factor(tax_matrix$taxon, levels=unique(unlist(tax_matrix[order(tax_matrix$counts),"taxon"])))


ggbarplot(tax_matrix, x="taxon", y="counts",  title="Taxon abundances", xlab="Taxon", ylab="Abundance") + theme_bw()


#### Plot original datasets (Figure 1a) ####
# plot top 10 taxa per sample (relative + absolute) AND special taxon
for(file in list.files("data/tax_matrices", full.names = T)){
    # read file 
    filename <- basename(file)
    tax_matrix <- read.table(file, header=T, stringsAsFactors=F, sep="\t") %>% 
        as.matrix()
    matrix_name <- strsplit(filename, split="_")[[1]][2]
    spread_name <- strsplit(filename, split="_")[[1]][3]
    scenario_name <- strsplit(filename, split="_")[[1]][5] %>% gsub(., pattern="\\.tsv", replacement="")
    
    tax_tibble <- as_tibble(tax_matrix) %>% 
        mutate(taxa=rownames(tax_matrix))
    
    ## choose top 10 and discard the rest
    tax_long <- tax_tibble %>% 
        gather(key = "sample", value = "counts", -taxa)
    top10 <- tax_long %>%
        group_by(taxa) %>%
        summarize(counts = sum(counts, na.rm=T)) %>%
        arrange(desc(counts)) %>%
        pull(taxa) %>%
        head(.,10)
    
    ## add special taxa in dysbiosis and in blooming
    if(scenario_name=="Dysbiosis"){
        opportunist <- c()
        opportunist <- matrix_stats %>% 
            filter(matrix==matrix_name) %>% 
            pull(special_taxon)
        top10sp <- unique(c(top10, opportunist))
        if(length(top10sp)!=10){ # keep the taxa plotted to 10 max
            top10sp <- unique(c(top10[1:(10-(length(top10sp)-length(top10)))], 
                                opportunist))
        }
    }
    if(scenario_name=="Blooming"){
        bloomer <- c()
        bloomer <- matrix_stats %>% 
            filter(matrix==matrix_name) %>% 
            pull(special_taxon)
        top10sp <- unique(c(top10, bloomer))
        if(length(top10sp)!=10){ # keep the taxa plotted to 10 max
            top10sp <- unique(c(top10[1:(10-(length(top10sp)-length(top10)))], 
                                bloomer))
        }
    }
    if(scenario_name=="Healthy"){ # No special taxa
        top10sp <- top10
    }
    
    ## regroup and plot
    simplifiedData <- tax_long
    simplifiedData[!simplifiedData$taxa %in% top10sp,"taxa"] <- NA
    
    simplifiedData <- simplifiedData %>%
        group_by(taxa, sample) %>%
        summarize_all(sum)
    
    simplifiedData$taxa <- as.factor(simplifiedData$taxa)
    simplifiedData$sample <- as.factor(simplifiedData$sample)
    
    simplifiedData$sample <- reorder(simplifiedData$sample, simplifiedData$counts)
    
    if(scenario_name=="Blooming"){
        simplifiedData$taxa <- gsub(simplifiedData$taxa, pattern=bloomer, 
                                    replacement=paste0(bloomer, "_Bloomer"))
        top10sp <- gsub(top10sp, pattern=bloomer, 
                                    replacement=paste0(bloomer, "_Bloomer"))
        top10sp <- sort(top10sp)
        legendtitle <- "Top 10 and special taxa"
        simplifiedData$taxa <- factor(simplifiedData$taxa)
        
    }
    if(scenario_name=="Dysbiosis"){
        simplifiedData$taxa <- gsub(simplifiedData$taxa, pattern=opportunist, 
                                    replacement=paste0(opportunist, "_Opportunist"))
        top10sp <- gsub(top10sp, pattern=opportunist, 
                      replacement=paste0(opportunist, "_Opportunist"))
        top10sp <- sort(top10sp)
        legendtitle <- "Top 10 and special taxa"
        simplifiedData$taxa <- factor(simplifiedData$taxa)
    }
    if(scenario_name=="Healthy"){ # No special taxa
        top10sp <- sort(top10sp)
        legendtitle <- "Top 10 taxa"
        simplifiedData$taxa <- factor(simplifiedData$taxa)
        
    }
    
    
    p1 <- ggbarplot(simplifiedData, x = "sample", y="counts", color="taxa", fill="taxa",
                    legend="right",font.legend = c(8, "plain", "black"),
                    legend.title=legendtitle, main=paste0("Absolute counts - simulation - Matrix: ", matrix_name, " - Spread: ", spread_name),
                    font.main = c(14,"bold", "black"), font.x = c(12, "bold"), 
                    font.y=c(12,"bold"), xlab="Sample", ylab="Absolute cell counts") + 
        scale_colour_manual(values=rev(inlmisc::GetColors(10, scheme = "sunset")), na.value="grey", labels = c(levels(simplifiedData$taxa), "Other")) +
        scale_fill_manual(values=rev(inlmisc::GetColors(10, scheme = "sunset")),  na.value="grey", labels = c(levels(simplifiedData$taxa), "Other")) +
        theme_bw() + 
        theme(axis.text.x=element_blank(),
              axis.ticks.x = element_blank()) + 
        theme(panel.grid.major = element_line(colour = "gray97"), 
              panel.grid.minor = element_line(colour = "gray97")) + 
        theme(axis.title = element_text(face = "bold"), 
              plot.title = element_text(size = 14, 
                                        face = "bold"), 
              legend.title = element_text(face = "bold"))    


    ggsave(plot = p1, filename=paste0("output/rawdataplots/toptaxa_", matrix_name, 
                             "_", spread_name, "_spread_",  scenario_name, ".pdf"), 
           device = "pdf", width = 7, height = 7)
}


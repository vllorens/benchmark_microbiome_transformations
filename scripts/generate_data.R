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

# load required functions
source("R/metadata.R")
source("R/rarefy_even_sampling_depth.R")

set.seed(14102019)

# create folders, if not existing
system("mkdir -p data/tax_matrices")
system("mkdir -p data/seq_matrices")
system("mkdir -p data/rmp_matrices")
system("mkdir -p data/qmp_matrices")
system("mkdir -p data/qmp2_matrices")
system("mkdir -p data/metadata_matrices")
system("mkdir -p output/rawdataplots")

# clean these folders and remove all previous results (if any) 
system("rm -r data/tax_matrices/*")
system("rm -r data/seq_matrices/*")
system("rm -r data/rmp_matrices/*")
system("rm -r data/qmp_matrices/*")
system("rm -r data/qmp2_matrices/*")
system("rm -r data/metadata_matrices/*")
system("rm -r output/rawdataplots")

# remove all objects in the environment
rm(list = setdiff(ls(), lsf.str()))


#### Generate original matrices ####

# Load simulated matrices
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
    # include a random subset of 200 samples per matrix (of a total of 1000 simulated)
    ind_matrix <- ind_matrix[,sample(colnames(ind_matrix), replace = F, size = 200)]
    spread <- (apply(ind_matrix,2,sum) %>% max)/(apply(ind_matrix,2,sum) %>% min) 
    print(scenario)
    print(spread)
    print((apply(ind_matrix,2,sum) %>% max))
    print((apply(ind_matrix,2,sum) %>% min))
    if(spread<10){
        spread_cat <- "low"
    } else{
        spread_cat <- "high"
    }
    write.table(ind_matrix, paste0("data/tax_matrices/taxonomy_", mt, "_", spread_cat, "_spread_", scenario, ".tsv"), 
                col.names = T, row.names=T, quote=F, sep="\t")
}

#### Generate sequencing data, RMP, QMP and QMP-NR matrices ####

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

#### Generate metadata ####
# First, define parameters for the metadata
rho_meta_min <- 0.35       # minimum correlation for the metadata variables
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


#### Plot original datasets (Figure 1a) ####
# plot top 10 taxa per sample (relative + absolute)
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
    
    simplifiedData <- tax_long
    simplifiedData[!simplifiedData$taxa %in% top10,"taxa"] <- NA
    
    simplifiedData <- simplifiedData %>%
        group_by(taxa, sample) %>%
        summarize_all(sum)
    
    simplifiedData$taxa <- as.factor(simplifiedData$taxa)
    simplifiedData$sample <- as.factor(simplifiedData$sample)
    
    simplifiedData$sample <- reorder(simplifiedData$sample, simplifiedData$counts)
    
    p1 <- ggbarplot(simplifiedData, x = "sample", y="counts", color="taxa", fill="taxa",
                    legend="right",font.legend = c(8, "plain", "black"),
                    legend.title="Top 10 taxa", main=paste0("Absolute counts - simulation - Matrix: ", matrix_name, " - Spread: ", spread_name),
                    font.main = c(14,"bold", "black"), font.x = c(12, "bold"), 
                    font.y=c(12,"bold"), xlab="Sample", ylab="Absolute cell counts") + 
        scale_colour_brewer(palette="Spectral", na.value="grey") +
        scale_fill_brewer(palette="Spectral", na.value="grey") +
        theme_bw() + 
        theme(axis.text.x=element_blank(),
              axis.ticks.x = element_blank()) + 
        theme(plot.title = element_text(face = "bold"))
    
    ggsave(plot = p1, paste0("output/rawdataplots/toptaxa_", matrix_name, 
                             "_", spread_name, "_spread_",  scenario_name, ".pdf"), 
           device = "pdf", width = 7, height = 7)
}


#### Calculate matrix basic stats ####

# Matrix stats
matrix_stats <- tibble()
stats <- as.tibble(stats_table)
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
    # correlation of taxa w number of cell counts
    taxoncor <- counts_taxa_correlation(tax_matrix, original_counts)
    taxoncor[[2]] <- p.adjust(taxoncor[[2]], method="BH")
    
    not_cor <- length(which(taxoncor[[2]]>0.05))
    neg_cor <- length(which(taxoncor[[2]]<0.05 & taxoncor[[1]]<0))
    pos_cor <- length(which(taxoncor[[2]]<0.05 & taxoncor[[1]]>0))
    # specialtaxa
    sptax <- stats %>% 
        dplyr::filter(grepl(paste0(matrix_name, "_"), matrix)) %>% 
        pull(dysbiotic_taxon)
    cor_special_effsize <- taxoncor[[1]][sptax]
    cor_special_pvalue <- taxoncor[[2]][sptax]
    # specialtaxa_2_onlydysbiosis
    if(scenario_name=="Dysbiosis"){
        sptax2 <- stats %>% 
            dplyr::filter(grepl(paste0(matrix_name, "_"), matrix)) %>% 
            pull(flat_taxon)
        cor_special_effsize2 <- taxoncor[[1]][sptax2]
        cor_special_pvalue2 <- taxoncor[[2]][sptax2]
    } else{
        sptax2 <- NA
        cor_special_effsize2 <- NA
        cor_special_pvalue2 <- NA
    }
    
    # put everythiing together
    matrix_stats_current <- tibble(matrix=matrix_name, spread_type=spread_name,
                                   scenario=scenario_name, spread=spread, 
                                   non_correlated_taxa=not_cor,
                                   positively_correlated_taxa=pos_cor,
                                   negatively_correlated_taxa=neg_cor,
                                   special_taxon=sptax,
                                   flat_taxon_dysbiosis=sptax2,
                                   correlation_special_counts=cor_special_effsize,
                                   correlation_special_counts_p=cor_special_pvalue,
                                   correlation_flat_counts_dysbiosis=cor_special_effsize2,
                                   correlation_flat_counts_dysbiosis_p=cor_special_pvalue2)
    matrix_stats <- bind_rows(matrix_stats, matrix_stats_current)
}

write_tsv(matrix_stats, "data/raw/matrix_stats_3scenarios.tsv", col_names = T)


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
# source("scripts/generate_data.R") # this is needed to generate the original matrices and metadata only

# create folders, if not existing
system("mkdir -p output/qmp_acs")

# clean these folders and remove all previous results (if any) 
system("rm -r output/qmp_acs/*")

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

stats_matrix <- read_tsv("data/raw/matrix_stats_3scenarios.tsv")

for(file in list.files("data/tax_matrices", full.names = T, pattern="Dysbiosis")){
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
    
    # instead of reading the metadata, make just one metadata variable (healthy vs diseased) that correlates with original cell counts, but not perfectly
    val <- 0
    chisq <- 1
    wilcox <- 1
    sptax <- stats_matrix %>% 
        dplyr::filter(matrix==matrix_name) %>% 
        pull(special_taxon)
    # Accessory function to calculate correlated metadata
    # function from https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
    complement <- function(y, rho, x) { 
        if (missing(x)) x <- rnorm(length(y)) 
        y.perp <- residuals(lm(x ~ y))
        rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
    }
    while(!(val < 0.9 & val > 0.7 & chisq < 0.01 & wilcox < 0.01)){ 
        rho_val <- runif(1, 0.7, 0.9) # get target rho correlation value
            
        z <- complement(y=original_counts, rho=rho_val)
            
        correlation = cor.test(original_counts,z, method = "spearman", exact = F) # spearman
        val <- abs(correlation$estimate)
        
        if (val < 0.9 & val > 0.7){ # make sure rho is within the specified range.
            numCats <- 2
            valCats <- c("groupA", "groupB")
            proportions <- runif(n = numCats,0.4,0.6) # to make sure that more or less half of the data falls within each category
            proportions <- cumsum(proportions)/sum(proportions) # to define the quantiles
                
            speciesZ <- z
            species_quantiles <- rev(quantile(speciesZ, probs=proportions))
            metaLine <- rep(NA, times=length(speciesZ))
            for(q in 1:length(species_quantiles)){
                metaLine[speciesZ<=species_quantiles[q]] <- valCats[q]
            }
            dftest <- data.frame(values=original_counts, categories=metaLine)
            chisq <- kruskal.test(values ~ categories, data=dftest)[[3]]
            dfsptax <- data.frame(values=c(unlist(tax_matrix[sptax,])), categories=metaLine)
            wilcox <- wilcox.test(values ~ categories, data=dfsptax)[[3]]
             
        } 
        
    }
    
    groupClass <- metaLine
    
    # check association between counts and the disease status
    counts_status <- tibble(original_counts, groupClass) %>% 
        rstatix::wilcox_test(original_counts ~ groupClass)

    group_class <- tibble(sample=colnames(tax_matrix), groupClass)
    
    # associations between taxa and disease status on real data
    realmatrix <- tax_matrix %>% 
        as_tibble() %>% 
        mutate(taxa=rownames(tax_matrix)) %>% 
        gather(key="sample", value="counts", -taxa) %>% 
        left_join(group_class, by="sample") %>% 
        group_by(taxa) 
    
    real_assoc_p <- realmatrix %>% 
        rstatix::wilcox_test(counts ~ groupClass) %>% 
        mutate(p_adj = p.adjust(p, method="BH"))
    
    real_assoc_sign <- realmatrix %>% 
        select(-sample) %>% 
        group_by(taxa, groupClass) %>% 
        summarize(counts=mean(counts)) %>% 
        spread(key="groupClass", value="counts") %>% 
        ungroup %>% 
        mutate(diff=groupA-groupB) %>% 
        mutate(sign=sign(diff))
    
    # plot example
    orig_data <- tibble(counts_original=counts_original[,1],  group=groupClass )
    orig_data$group <- factor(orig_data$group, levels = c("groupA", "groupB"))
    p1 <- ggboxplot(orig_data, x="group", y="counts_original", 
                    title="Disease associated to lower microbial loads (dysbiosis)",
                    ylab="Microbial loads (cells/gram of stool)", xlab="",
                    palette=inlmisc::GetTolColors(5, scheme = "sunset")[c(2,5)],
                    fill="group", add="jitter") + 
        stat_compare_means(comparisons = list(c("groupA", "groupB"))) +
        theme_bw() +
        theme(panel.grid.major = element_line(colour = "gray97"), 
              panel.grid.minor = element_line(colour = "gray97")) + 
        theme(axis.title = element_text(face = "bold"), 
              plot.title = element_text(size = 14, 
                                        face = "bold"), legend.title = element_text(face = "bold"))
    ggsave(p1, file=paste0("output/qmp_acs/example_", matrix_name, "_counts.pdf"), device="pdf", 
           width=5, height=3.5, useDingbats=FALSE)
    opportunist_taxa <- stats_matrix %>% dplyr::filter(matrix==matrix_name) %>% pull(special_taxon)
    flat_taxa <- stats_matrix %>% dplyr::filter(matrix==matrix_name) %>% pull(flat_taxon_dysbiosis) 
    flat_taxa <- strsplit(flat_taxa, split=",")[[1]][1]
    normal_taxa <- sample(setdiff(rownames(tax_matrix), c(opportunist_taxa, flat_taxa)),1)
    print(paste0("Matrix ", matrix_name))
    print(paste0("Normal taxa: ", normal_taxa))
    print(paste0("Opportunist taxa: ", opportunist_taxa))
    print(paste0("Flat taxa: ", flat_taxa))
    opp_matrix <- tax_matrix[opportunist_taxa,, drop=F] %>% 
        as_tibble %>% 
        mutate(taxa=opportunist_taxa) %>% 
        gather(key="sample", value="counts", -taxa) %>% 
        mutate(type="opportunist")
    flat_matrix <- tax_matrix[flat_taxa,, drop=F] %>% 
        as_tibble %>% 
        mutate(taxa=flat_taxa) %>% 
        gather(key="sample", value="counts", -taxa) %>% 
        mutate(type="flat")
    normal_matrix <- tax_matrix[normal_taxa,, drop=F] %>% 
        as_tibble %>% 
        mutate(taxa=normal_taxa) %>% 
        gather(key="sample", value="counts", -taxa) %>% 
        mutate(type="normal")
    toplot <- bind_rows(opp_matrix, flat_matrix, normal_matrix) %>% 
        left_join(group_class, by="sample")
    toplot$groupClass <- factor(toplot$groupClass, levels = c("groupA", "groupB"))
    toplot$type <- factor(toplot$type, levels = c("normal", "opportunist", "flat"))
    p2 <- ggboxplot(toplot, x="groupClass", y="counts", 
              fill="groupClass", add="jitter",
              palette=inlmisc::GetTolColors(5, scheme = "sunset")[c(2,5)],
              facet.by = "type", scales="free_y", 
              xlab="", ylab="Taxon counts") + 
        stat_compare_means(comparisons = list(c("groupA", "groupB"))) +
        theme_bw() +
        theme(panel.grid.major = element_line(colour = "gray97"), 
              panel.grid.minor = element_line(colour = "gray97")) + 
        theme(axis.title = element_text(face = "bold"), 
              plot.title = element_text(size = 14, 
                                        face = "bold"), legend.title = element_text(face = "bold"))
    ggsave(p2, file=paste0("output/qmp_acs/example_", matrix_name, "_taxa.pdf"), device="pdf", 
           width=11, height=3.5, useDingbats=FALSE)
    
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
        
        # remove taxa that has zero counts throughout all qmp samples - to avoid an error in the comparisons
        taxa_to_remove <- which(apply(qmp,1,sum)==0)
        if(length(taxa_to_remove)>0){
            qmp <- qmp[-taxa_to_remove,]
            acs <- acs[-taxa_to_remove,]
        }
        
        #correlations on qmp/ACS data (all taxon-metadata correlations, to determine Precision and Recall)
        data_qmp_p <- qmp %>% 
            as_tibble() %>% 
            mutate(taxa=rownames(qmp)) %>% 
            gather(key="sample", value="counts", -taxa) %>% 
            left_join(group_class, by="sample") %>% 
            group_by(taxa) %>% 
            rstatix::wilcox_test(counts ~ groupClass) %>% 
            mutate(p_adj = p.adjust(p, method="BH"))
        
        data_qmp_sign <- qmp %>% 
            as_tibble() %>% 
            mutate(taxa=rownames(qmp)) %>% 
            gather(key="sample", value="counts", -taxa) %>% 
            left_join(group_class, by="sample") %>% 
            select(-sample) %>% 
            group_by(taxa, groupClass) %>% 
            summarize(counts=mean(counts)) %>% 
            spread(key="groupClass", value="counts") %>% 
            ungroup %>% 
            mutate(diff=groupA-groupB) %>% 
            mutate(sign=sign(diff))
        
        data_acs_p <- acs %>% 
            as_tibble() %>% 
            mutate(taxa=rownames(acs)) %>% 
            gather(key="sample", value="counts", -taxa) %>% 
            left_join(group_class, by="sample") %>% 
            group_by(taxa) %>% 
            rstatix::wilcox_test(counts ~ groupClass) %>% 
            mutate(p_adj = p.adjust(p, method="BH"))
        
        data_acs_sign <- acs %>% 
            as_tibble() %>% 
            mutate(taxa=rownames(acs)) %>% 
            gather(key="sample", value="counts", -taxa) %>% 
            left_join(group_class, by="sample") %>% 
            select(-sample) %>% 
            group_by(taxa, groupClass) %>% 
            summarize(counts=mean(counts)) %>% 
            spread(key="groupClass", value="counts") %>% 
            ungroup %>% 
            mutate(diff=groupA-groupB) %>% 
            mutate(sign=sign(diff))
        
        
        # first remove taxa from the original matrix as well 
        real_assoc_p_filtered <- real_assoc_p %>% dplyr::filter(taxa %in% data_qmp_p$taxa)
        real_assoc_sign_filtered <- real_assoc_sign %>% dplyr::filter(taxa %in% data_qmp_sign$taxa)
        
        # assess performance
        test_significant_qmp_all <- data_qmp_p$p_adj < 0.05
        test_significant_acs_all <- data_acs_p$p_adj < 0.05
        reference_significant_all <- real_assoc_p_filtered$p_adj < 0.05
        test_sign_qmp_all <- data_qmp_sign$sign
        test_sign_acs_all <- data_acs_sign$sign
        reference_sign_all <- real_assoc_sign_filtered$sign
        
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
        
        # Assess performance specifically on taxa associated to disease
        # taxa associated to disease
        sign_disease <- table(reference_sign_all) %>% which.min %>% names() %>% as.numeric()
        taxa_disease <- real_assoc_sign_filtered %>% 
            dplyr::filter(sign==sign_disease) %>% 
            pull(taxa)
        taxa_disease <- real_assoc_p_filtered %>% 
            dplyr::filter(taxa %in% taxa_disease, p_adj<0.05) %>% 
            pull(taxa)
        test_significant_qmp_all <- data_qmp_p %>% 
            dplyr::filter(taxa %in% taxa_disease) %>% 
            pull(p_adj) < 0.05
        test_significant_acs_all <- data_acs_p %>% 
            dplyr::filter(taxa %in% taxa_disease) %>% 
            pull(p_adj) < 0.05
        reference_significant_all <- real_assoc_p_filtered %>% 
            dplyr::filter(taxa %in% taxa_disease) %>% 
            pull(p_adj) < 0.05
        test_sign_qmp_all <- data_qmp_sign %>% 
            dplyr::filter(taxa %in% taxa_disease) %>% 
            pull(sign)
        test_sign_acs_all <- data_acs_sign %>% 
            dplyr::filter(taxa %in% taxa_disease) %>% 
            pull(sign)
        reference_sign_all <- real_assoc_sign_filtered %>% 
            dplyr::filter(taxa %in% taxa_disease) %>% 
            pull(sign)
        
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
        valuesused <- c(valuesused, "disease")
        
        # Assess performance specifically on invariant taxa
        # taxa associated to disease
        taxa_flat <-  real_assoc_p_filtered %>% 
            dplyr::filter(p_adj>0.05) %>% 
            pull(taxa)
        test_significant_qmp_all <- data_qmp_p %>% 
            dplyr::filter(taxa %in% taxa_flat) %>% 
            pull(p_adj) < 0.05
        test_significant_acs_all <- data_acs_p %>% 
            dplyr::filter(taxa %in% taxa_flat) %>% 
            pull(p_adj) < 0.05
        reference_significant_all <- real_assoc_p_filtered %>% 
            dplyr::filter(taxa %in% taxa_flat) %>% 
            pull(p_adj) < 0.05
        test_sign_qmp_all <- data_qmp_sign %>% 
            dplyr::filter(taxa %in% taxa_flat) %>% 
            pull(sign)
        test_sign_acs_all <- data_acs_sign %>% 
            dplyr::filter(taxa %in% taxa_flat) %>% 
            pull(sign)
        reference_sign_all <- real_assoc_sign_filtered %>% 
            dplyr::filter(taxa %in% taxa_flat) %>% 
            pull(sign)
        # test_sign_qmp_all[is.na(test_sign_qmp_all)] <- 1
        # test_sign_acs_all[is.na(test_sign_acs_all)] <- 1
        # reference_sign_all[is.na(reference_sign_all)] <- 1
        
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
        valuesused <- c(valuesused, "flat")
    }
    
    result <- tibble(fdr_qmp, fdr_acs, spread_v, scenario_v, fp_qmp, fp_acs, tp_qmp, tp_acs, pr_qmp, pr_acs, depths, matname, valuesused)
    write_tsv(result, "output/qmp_acs/qmp_vs_acs_taxonmetadata_depths.tsv", col_names = T)
}


#### Assess performance of both QMP and ACS ####

result <- read_tsv("output/qmp_acs/qmp_vs_acs_taxonmetadata_depths.tsv", col_names=T)
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
                 numeric.x.axis = T, add = "mean_se", palette=inlmisc::GetTolColors(5, scheme = "sunset")[c(2,5)], xlab="Sequencing reads", 
                 ylab="FPR [FP/FP+TN]", title = "Taxa-metadata associations - FPR", xscale="log10", .format=T,
                 legend.title="Method") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(fprall, file="output/qmp_acs/qmp_vs_acs_alltaxa_fdr_metadata.pdf", device="pdf", 
       width=8, height=3,useDingbats=FALSE)


result <- read_tsv("output/qmp_acs/qmp_vs_acs_taxonmetadata_depths.tsv", col_names=T)
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
                       numeric.x.axis = T, add = "mean_se", palette=inlmisc::GetTolColors(5, scheme = "sunset")[c(2,5)], xlab="Sequencing reads", 
                       ylab="Precision [TP/TP+FP]", title = "Taxa-metadata associations - Precision",xscale="log10", .format=T,
                       legend.title="Method") + 
    
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(precisionall, file="output/qmp_acs/qmp_vs_acs_alltaxa_precision_metadata.pdf", device="pdf", 
       width=8, height=3,useDingbats=FALSE)


result <- read_tsv("output/qmp_acs/qmp_vs_acs_taxonmetadata_depths.tsv", col_names=T)
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
                    numeric.x.axis = T, add = "mean_se", palette=inlmisc::GetTolColors(5, scheme = "sunset")[c(2,5)], xlab="Sequencing reads", 
                    ylab="Recall [TP/TP+FN]", title = "Taxa-metadata associations - Recall", xscale="log10", .format=T,
                    legend.title="Method") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(recallall, file="output/qmp_acs/qmp_vs_acs_alltaxa_recall_metadata.pdf", device="pdf", 
       width=8, height=3,useDingbats=FALSE)


#### Assess performance of both QMP and ACS (invariant taxa only, no dysbiosis) ####

result <- read_tsv("output/qmp_acs/qmp_vs_acs_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::filter(valuesused=="flat") %>% 
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
                 numeric.x.axis = T, add = "mean_se", palette=inlmisc::GetTolColors(5, scheme = "sunset")[c(2,5)],
                 xlab="Sequencing reads", 
                 ylab="FPR [FP/FP+TN]", title = "Invariant taxa-metadata associations - FPR",
                 legend.title="Method", xscale="log10", .format=T) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(fprinv, file="output/qmp_acs/qmp_vs_acs_invtaxa_fdr_metadata.pdf", device="pdf", 
       width=8, height=3, useDingbats=F)


result <- read_tsv("output/qmp_acs/qmp_vs_acs_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::filter(valuesused=="flat") %>% 
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
                       numeric.x.axis = T, add = "mean_se", palette=inlmisc::GetTolColors(5, scheme = "sunset")[c(2,5)], xlab="Sequencing reads", 
                       ylab="Precision [TP/TP+FP]", title = "Invariant taxa-metadata associations - Precision",
                       legend.title="Method", xscale="log10", .format=T) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(precisioninv, file="output/qmp_acs/qmp_vs_acs_invtaxa_precision_metadata.pdf", device="pdf", 
       width=8, height=3, useDingbats=F)


result <- read_tsv("output/qmp_acs/qmp_vs_acs_taxonmetadata_depths.tsv", col_names=T)
resultlong <- result %>% 
    dplyr::filter(valuesused=="disease") %>% 
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
                    numeric.x.axis = T, add = "mean_se", palette=inlmisc::GetTolColors(5, scheme = "sunset")[c(2,5)],
                    xlab="Sequencing reads", 
                    ylab="Recall [TP/TP+FN]", title = "Opportunist taxa-metadata associations - Recall",
                    legend.title="Method", xscale="log10", .format=T) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "gray97"), 
          panel.grid.minor = element_line(colour = "gray97")) + 
    theme(axis.title = element_text(face = "bold"), 
          plot.title = element_text(size = 14, 
                                    face = "bold"), legend.title = element_text(face = "bold"))
ggsave(recallinv, file="output/qmp_acs/qmp_vs_acs_opportunist_recall_metadata.pdf", device="pdf", 
       width=8, height=3,useDingbats=FALSE)


pub_plot <- ggarrange(recallall, fprall, fprinv, ncol=1, nrow=3, labels="auto",legend="right")
ggsave(pub_plot, file="output/qmp_acs/qmp_vs_acs_summary.pdf", device="pdf", 
       width=8, height=9,useDingbats=FALSE)



#' Generate metadata
#' 
#' This function generates the metadata for a single taxonomy matrix
#' 
#' @param taxonomic_matrix matrix with rows as taxa and columns as samples
#' @param num_noc number of metadata features, numeric and non-correlated with taxonomic composition
#' @param num_ct number of metadata features, numeric and correlated with some taxon abundance
#' @param num_cc number of metadata features, numeric and correlated with total counts
#' @param cat_noc number of metadata features, categoric and non-correlated with taxonomic composition
#' @param cat_ct number of metadata features, categoric and correlated with some taxon abundance
#' @param cat_cc number of metadata features, categoric and correlated with total counts
#' @param corr_min minimum correlation value (absolute value) for the correlated features
#' @param corr_max maximum correlation value (absolute value) for the correlated features
#' 
#' @return metadata data frame

generate_metadata <- function(taxonomic_matrix=NULL, num_noc=0, num_ct=0, num_cc=0, 
                              cat_noc=0, cat_ct=0, cat_cc=0, corr_min, corr_max){

    # Evaluate arguments
    if(num_noc<0) stop("Number of metadata features must be non-negative")
    if(num_ct<0) stop("Number of metadata features must be non-negative")
    if(num_cc<0) stop("Number of metadata features must be non-negative")
    if(cat_noc<0) stop("Number of metadata features must be non-negative")
    if(cat_ct<0) stop("Number of metadata features must be non-negative")
    if(cat_cc<0) stop("Number of metadata features must be non-negative")
    
    N <- num_noc + num_ct + num_cc + cat_noc + cat_ct + cat_cc
    if(N==0) stop("Please include at least one metadata column of any type")
    
    if(corr_min>corr_max) stop("corr_min value needs to be smaller than corr_max value")
    
    if(!is.matrix(taxonomic_matrix) & !is.data.frame(taxonomic_matrix))
        stop("Please provide a taxonomic matrix")
    
    if(is.data.frame(taxonomic_matrix))
        taxonomic_matrix <- as.matrix(taxonomic_matrix)
    
    # Accessory function to calculate correlated metadata
    # function from https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
    complement <- function(y, rho, x) { 
        if (missing(x)) x <- rnorm(length(y)) 
        y.perp <- residuals(lm(x ~ y))
        rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
    }
    
    # Print message showing number of samples and taxons so that user can transpose taxonomy matrix if necessary
    print(paste0("The taxonomy matrix contains ", nrow(taxonomic_matrix), " taxa and ", 
                 ncol(taxonomic_matrix), " samples. If this is not correct, TRANSPOSE!"))
    
    # Define metaData dimensions
    metaData = matrix(nrow = ncol(taxonomic_matrix), ncol = N)
    a=1
    
    # Populate the matrix
    # NUMERIC, NON-CORRELATED
    if(num_noc>0){
        for (i in 1:num_noc){
            typeDist <- sample(c(1:3),1) # different types of distributions can be simulated for the uncorrelated: unif, normal, nbiomial
            if(typeDist==1){ # uniform
                parameters <- sample(1:100,2)
                metaData[,a] = runif(ncol(taxonomic_matrix), min=min(parameters), max=max(parameters))
            } else if(typeDist==2){ # gaussian
                parameters <- sample(1:100,2)
                metaData[,a] = rnorm(ncol(taxonomic_matrix), mean = max(parameters), sd=min(parameters))
            } else if(typeDist==3){ #neg binomial
                metaData[,a] <- rnbinom(ncol(taxonomic_matrix), size = runif(1, min=1, max=100), prob=runif(1, min=0.1, max=0.9))
            }
            a <- a+1
        }
    }
    # NUMERIC, CORRELATED TO SOME TAXON
    if(num_cc>0){
        for (i in 1:num_cc){
            val <- 0
            isvalid <- F
            while(!(val < corr_max & val > corr_min & isvalid == T)){
                index = sample(1:nrow(taxonomic_matrix),1) # get the index of the taxon correlating
                isvalid <- length(which(taxonomic_matrix[index,]==0))<0.5*ncol(taxonomic_matrix) # ensure there are not too many zeros
                rho_val <- runif(1, corr_min, corr_max)
                rho_sign <- sample(c(-1,1),1)
                rho_val <- rho_val*rho_sign
                
                z <- complement(y=taxonomic_matrix[index,], rho=rho_val)
                
                correlation = cor.test(taxonomic_matrix[index,],z, method = "spearman", exact = F) # spearman
                val <- abs(correlation$estimate)
                if(is.na(val)){
                    val=0
                }
                if (val < corr_max & val > corr_min){ # make sure rho is within the specified range.
                    metaData[,a] = z
                } 
            }    
            a <- a+1
        }
    }
    # NUMERIC, CORRELATED TO TOTAL COUNTS
    if(num_ct>0){
        for (i in 1:num_ct){
            total_counts <- apply(taxonomic_matrix, 2,sum)
            val <- 0
            while(!(val < corr_max & val > corr_min)){
                rho_val <- runif(1, corr_min, corr_max)
                rho_sign <- sample(c(-1,1),1)
                rho_val <- rho_val*rho_sign
                
                z <- complement(y=total_counts, rho=rho_val)
                
                correlation = cor.test(total_counts,z, method = "spearman", exact = F) # spearman
                val <- abs(correlation$estimate)
                
                if (val < corr_max & val > corr_min){ # make sure rho is within the specified range.
                    metaData[,a] = z
                } 
            }    
            a <- a+1
        }
    }
    metaData <- as.data.frame(metaData)
    # CATEGORIC, NON-CORRELATED
    if(cat_noc>0){
        for (i in 1:cat_noc){
            numCats <- sample(c(2,4,6,8),1) # select how many categories in this variable
            valCats <- letters[1:numCats]
            proportions <- runif(n = numCats,0.2,1) # select proportions of each category
            metaData[,a] <- sample(valCats, size = ncol(taxonomic_matrix), prob = proportions, replace=T)
            a <- a+1
        }
    }
    # CATEGORIC, CORRELATED TO SOME TAXON
    # in order to avoid perfect correlations, we first generate a correlated numeric 
    # variable, and then the categorical variable maps to the correlated numeric
    if(cat_ct>0){
        for (i in 1:cat_ct){
            val <- 0
            chisq <- 1
            isvalid <- F
            while(!(val < corr_max & val > corr_min & chisq < 0.01 & isvalid == T)){
                index = sample(1:nrow(taxonomic_matrix),1) # get the index of the taxon correlating
                isvalid <- length(which(taxonomic_matrix[index,]==0))<0.5*ncol(taxonomic_matrix) # ensure there are not too many zeros
                rho_val <- runif(1, corr_min, corr_max) # get target rho correlatio value
                
                z <- complement(y=taxonomic_matrix[index,], rho=rho_val)
                
                correlation <- cor.test(taxonomic_matrix[index,],z, method = "spearman", exact = F) # spearman
                val <- abs(correlation$estimate)
                if(is.na(val)){
                    val <- 0
                }
                if (val < corr_max & val > corr_min){ # make sure rho is within the specified range.
                    numCats <- sample(c(2,4,6,8),1)
                    valCats <- letters[1:numCats]
                    proportions <- runif(n = numCats,0,1)
                    proportions <- cumsum(proportions)/sum(proportions) # to define the quantiles
                    
                    speciesZ <- z
                    species_quantiles <- rev(quantile(speciesZ, probs=proportions))
                    metaLine <- rep(NA, times=length(speciesZ))
                    for(q in 1:length(species_quantiles)){
                        metaLine[speciesZ<=species_quantiles[q]] <- valCats[q]
                    }
                    dftest <- data.frame(values=taxonomic_matrix[index,], categories=metaLine)
                    chisq <- kruskal.test(values ~ categories, data=dftest)[[3]]
                    metaData[,a] <- metaLine 
                } 
            }   
            a=a+1
        }
    }
    # CATEGORIC, CORRELATED TO TOTAL COUNTS
    # in order to avoid perfect correlations, we first generate a correlated numeric
    # variable, and then the categorical variable maps to the correlated numeric
    if(cat_cc>0){
        for (i in 1:cat_cc){
            total_counts <- apply(taxonomic_matrix, 2,sum)
            val <- 0
            chisq <- 1
            while(!(val < corr_max & val > corr_min & chisq < 0.01)){ 
                rho_val <- runif(1, corr_min, corr_max) # get target rho correlatio value
                
                z <- complement(y=total_counts, rho=rho_val)
                
                correlation = cor.test(total_counts,z, method = "spearman", exact = F) # spearman
                val <- abs(correlation$estimate)
                
                if (val < corr_max & val > corr_min){ # make sure rho is within the specified range.
                    numCats <- sample(c(2,4,6,8),1)
                    valCats <- letters[1:numCats]
                    proportions <- runif(n = numCats,0,1)
                    proportions <- cumsum(proportions)/sum(proportions) # to define the quantiles
                    
                    speciesZ <- z
                    species_quantiles <- rev(quantile(speciesZ, probs=proportions))
                    metaLine <- rep(NA, times=length(speciesZ))
                    for(q in 1:length(species_quantiles)){
                        metaLine[speciesZ<=species_quantiles[q]] <- valCats[q]
                    }
                    dftest <- data.frame(values=total_counts, categories=metaLine)
                    chisq <- kruskal.test(values ~ categories, data=dftest)[[3]]
                    metaData[,a] <- metaLine 
                } 
            }    
            a=a+1
        }
    }
    return(metaData)
}


# Microbiota simulations

This repository contains the scripts to generate and analyze the simulations from the paper "xxx xxx xxx". 


## How to use this repository

Clone or download this repository to your local computer. 

The different folders and files contained here are:

* **scripts**: contains the scripts to generate the simulated data and run the analyses included in the paper
* **data/raw**: contains the data object including the raw simulations of the three types of communities, from which, the original matrices are generated
* **R**: contains dditional helper functions, called from the scripts
* **README.md**: this file
* **microbiota_simultions.Rproj**: the R project from which to run all the analyses (for easier portability)
* **.gitignore**: all the files not to include here

## Steps to reproduce the analyses of the paper

### 1. Generate data

Run the script **generate_data.R**. In this script, first we obtain individual simulated taxonomy matrices. Three types of populations are simulated:

* Healthy, where there is an ecological succession in which higher microbial load correlates to higher richness.
* Dysbiosis, where in a fraction of the samples there is a drop in microbial load. Most taxa decrease or dissapear, some remain stable and others increase with lower cell counts.
* Blooming, where a single taxon blooms and overtakes the community.

From the data object with all matrices, we write txt files with the original matries and we also produce the sequencing output, the RMP matrices, and the QMP matrices. Finally, we also generate metadata for all the original taxonomy matrices, and calculate some basic statistics about the generated matrices.

### 2. Calculate alpha diversity on the simulated matrices

Run the script **alpha_diversity.R**. In this script, we calculate different alpha diversity metrics on the original matrices, as well as on the matrices transformed by sequencing, RMP, QMP and QMP-NR.

### 3. Calculate taxon-taxon and taxon-metadata correlations

Run the script **correlations.R**. In this script, we use different data transformation methods to determine taxon-taxon and taxon-metadata associations in all our populations. For comparisons, we will only use Spearman correlations for taxon-taxon and either Spearman correlation or Kruskal-wallis test for taxon-metadata associations (Spearman for numeric and Kruskal-Wallis for factors in the metadata)

The following methods will be tested:

* REAL: Original data (as reference). Not transformed.
* CLR: Sequencing data. CLR transformed.
* CSS: Sequencing data. Normalized using cumulative sum scaling (MetagenomeSeq implemetnation).
* AST: Sequencing data. Normalized using total sum scaling + arcsine squared-root transformation.
* GMPR: Sequencing data. Size factors for normalization are calculated using the Geometric Mean of Pairwise Ratios (implementation using GMPR package).
* TMM: Sequencing data. Normalization using trimmed mean of M-values (EdgeR implementation).
* RLE: Sequencing data. Normalization using relative log-expression (EdgeR implementation).
* UQ: Sequencing data. Normalization using upper quartile (EdgeR implementation).
* VST: Sequencing data. Using variance stabilizing transformation from DESeq2.
* RMP: RMP data. Rarefied to equal sequencing depth.
* QMP: QMP data. Normalized by rarefying to even sampling depth and scaling to total counts.
* QMP-NR: QMP, non-rarefied data. Normalized by scaling to total counts.

In all of these methods, we will only evaluate those taxa with prevalence > 50%. We will use the sequencing matrix to determine this threshold on the zeros.


### 4. Comparison QMP vs QMP-NR

Run the script **qmp_vs_qmpnr.R**. In this script we compare the performance of QMP and QMP-NR in detecting taxon-metadata associations, focusing on invariant taxa (in absolute abundances) and metadata that are associated to cell counts, to find which method performs better in terms of lower false positive rates.


### 5. Test the effect of increasing microbial load spreads

Run the script **effect_microbial_spread.R**. This script tests the performance of the different methods tested in **section 3**, but instead of focusing on the different scenarios modeled, we focus on the different spreads in microbial load that are observed in the generated matrices.

### 6. Test the effect of the number of samples, sequencing depth and taxon prevalence threshold

Run the following scripts: 

* **effect_number_samples.R**
* **effect_sequencing_depth.R**
* **effect_taxon_prevalence.R**

These three scripts repeat the correlation analysis from **section 3**, but screening throughout the different parameters to evaluate their influcence on the different transformations. These scripts require generating the data and calculating all of the correlations several times, so they may take a long time. You may consider running them in the command line.

### 7. Evaluate proportionality vs Spearman correlation

Run the script **proportionality_vs_spearman.R**. This script evaluates the performance of proportionality vs Spearman correlation in capturing taxon-taxon correlations. It is highly recommended to run this script in a cluster to take maximum advantage of the multicore functionality. Otherwise this will take a long time to run.







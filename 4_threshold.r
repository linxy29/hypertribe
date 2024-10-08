# Checkpoint!
# What we have now: 
# - An excel sheet of counts of ref and alt alleles for each RNA edit site, across all 9 samples

# What do we want to know: 
# - Which edit sites are TRUE edit sites (i.e. our fusion protein generated a mismatch there)?
# - Which edit sites are NOISY/FALSE edit sites (i.e. also present in control; may have baseline alt allele count)?

# How do we do this? 
# - We need to make some sort of comparison of the edit frequencies between the fusion protein samples, and the control samples. 
# - If they are different enough, we can say that those edit sites are indeed TRUE edit sites, caused by our fusion protein. 
# - If they are not different enough between fusion protein and control samples, we likely say that those edit sites are probably FALSE edit sites (i.e. not generated by our fusion protein)

# ---> This script can help tell us if the edit frequencies are "different enough" between our FP samples and controls.

# We first begin by trying to model our edit frequencies after a certain statistical distribution:

#################################
# 1. Beta-binomial distribution
#################################
# The beta-binomial distribution is a family of discrete probability distributions on a finite support of non-negative integers arising when the probability of success in each of a fixed or known number of Bernoulli (Bernoulli distribution is the discrete probability distribution of a random variable which takes the value 1 with probability p and the value 0 with probability) trials is either unknown or random.

# The beta-binomial distribution is the binomial distribution in which the probability of success at each of n trials is not fixed but randomly drawn from a beta distribution.

# Has 3 parameters: n (no. of trials), alpha and beta (a parameter is any measured quantity of a statistical population that summarises or describes an aspect of the population, such as a mean or a standard deviation), where alpha > 0 and beta > 0

### Why beta-binomial distribution?
# Single parameter distributions, such as Poisson and Binomial, imply that the variance is determined by the mean.

# However in many cases and especially in analysis of biological data, this mean-variance relationship fails mainly due to presence of overdispersion, where the data have a higher variance than anticipated under the simple model. 

# In biology, overdispersion occurs since there is more than one source of variation: technical variation due to error measurements coming from the experiment design and biological variation between the subjects of interest, e.g. different cells on different conditions.

# By using beta-binomial, we allow each edit site to be independent of each other, thus having a different probability of success (i.e. A-G edit). We can then account for overdispersion when the probabilities of success are thought to vary from edit to edit.

# A good read that explains the use of beta-binomial distribution to model overdispersion 
# https://rpubs.com/cakapourani/beta-binomial

# We then want to estimate the parameters of our beta-binomial distribution:

##################################
# 2. Maximum likelihood estimation
##################################
# MLE is a method of estimating the parameters of a probability distribution by maximizing a likelihood function, so that under the assumed statistical model the observed data is most probable. (i.e. given a set of observations/data, what parameters best explains the data --> returns the highest probability for the observed data?). Intuitively, it selects the parameter values that make the observed data most probable.

# The point in the parameter space that maximizes the likelihood function is called the maximum likelihood estimate.

# Often convenient to work with the natural logarithm of the likelihood function, called the log-likelihood (easier to work with logs, especially for exponential distributions)



# Libraries required
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("bbmle")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("statmod")
library(bbmle)
library(statmod)
library(openxlsx)
library(parallel)
library(stringr)
library(VGAM)
library(DEXSeq)
library(tidyverse)

# Read in excel file of snp counts
setwd("/home/linxy29/data/TRIBE/HGC20240511005-0002/extract_RNAedit")
df_snp_counts <- read.xlsx( "snp_counts_dedupped.xlsx", sheet = 1 )

# Using pre-defined threshold to define the RNA editing site detection function with customizable thresholds and column names
detect_rna_editing_site <- function(df, 
                                    A_count_gDNA_col, 
                                    G_count_gDNA_col, 
                                    A_count_RNA_col, 
                                    G_count_RNA_col,
                                    gDNA_coverage_threshold = 4, 
                                    gDNA_editing_threshold = 0.005, 
                                    A_count_threshold = 0.8, 
                                    RNA_editing_threshold = 0) {
  
  # Create a new column in the dataframe to store the results
  is_editing_site <- with(df, {
    # Access the specified columns dynamically
    A_count_gDNA <- df[[A_count_gDNA_col]]
    G_count_gDNA <- df[[G_count_gDNA_col]]
    A_count_RNA <- df[[A_count_RNA_col]]
    G_count_RNA <- df[[G_count_RNA_col]]
    
    # Calculate the total gDNA coverage
    gDNA_coverage <- A_count_gDNA + G_count_gDNA
    
    # Apply the criteria step-by-step with the customizable thresholds
    is_coverage_sufficient <- gDNA_coverage > gDNA_coverage_threshold
    is_gDNA_editing_low <- (G_count_gDNA / gDNA_coverage) < gDNA_editing_threshold
    is_A_count_dominant <- (A_count_gDNA / gDNA_coverage) >= A_count_threshold
    is_RNA_editing_present <- G_count_RNA > RNA_editing_threshold
    
    # Return TRUE if all conditions are satisfied, FALSE otherwise
    is_editing_site <- is_coverage_sufficient & is_gDNA_editing_low & is_A_count_dominant & is_RNA_editing_present
    return(is_editing_site)
  })
  
  # Return the modified dataframe with the RNA editing site classification
  return(is_editing_site)
}

df_snp_counts$A1O5_editing <- detect_rna_editing_site(df_snp_counts, A_count_gDNA_col = "A1_dedupped.bam.ref.count", G_count_gDNA_col = "A1_dedupped.bam.alt.count", 
                                             A_count_RNA_col = "AO5_dedupped.bam.ref.count", G_count_RNA_col = "AO5_dedupped.bam.alt.count", 
                                             gDNA_coverage_threshold = 4, gDNA_editing_threshold = 0.005, A_count_threshold = 0.8, RNA_editing_threshold = 0)
df_snp_counts$A3O7_editing <- detect_rna_editing_site(df_snp_counts, A_count_gDNA_col = "A3_dedupped.bam.ref.count", G_count_gDNA_col = "A3_dedupped.bam.alt.count", 
                                                      A_count_RNA_col = "AO7_dedupped.bam.ref.count", G_count_RNA_col = "AO7_dedupped.bam.alt.count", 
                                                      gDNA_coverage_threshold = 4, gDNA_editing_threshold = 0.005, A_count_threshold = 0.8, RNA_editing_threshold = 0)

                              
write.xlsx(df_snp_counts, file = "RNAedit.xlsx")

# Generate BedGraph file

convert2bedgraph <- function(df, sample, filter = NA, file) {
  if (!is.na(filter)){
    df = df[df[[filter]] == TRUE,]
  }
  # Adding 1 to the start position to make it zero-based as required in some genome browsers
  df$end <- df$start + 1 
  
  # Dynamically detect the columns for the ref and alt counts based on the sample name
  ref_col <- grep(paste0(sample, ".*ref\\.count"), colnames(df), value = TRUE)
  alt_col <- grep(paste0(sample, ".*alt\\.count"), colnames(df), value = TRUE)
  
  # If the columns are found, proceed with the calculations
  if (length(ref_col) == 1 && length(alt_col) == 1) {
    
    # Calculating a score, for example, the percentage of alternative reads
    df$score <- with(df, ifelse(df[[ref_col]] + df[[alt_col]] > 0, 
                                df[[alt_col]] / (df[[ref_col]] + df[[alt_col]]) * 100, 
                                0))
    
    # Prepare the data frame for output
    output_df <- df[c("seqnames", "start", "end", "score", "gene.symbol", "strand", "annotation")]
    
    # Write to bedgraph file
    write.table(output_df, 
                file = file, 
                quote = FALSE, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = FALSE)
    
  } else {
    stop("Could not find the ref or alt count columns for the specified sample.")
  }
}

## Get the bed file
convert2bedgraph(df = df_snp_counts, sample = "A1", file = "A1.bedgraph")
convert2bedgraph(df = df_snp_counts, sample = "A3", file = "A3.bedgraph")
convert2bedgraph(df = df_snp_counts, sample = "AO5", file = "AO5.bedgraph")
convert2bedgraph(df = df_snp_counts, sample = "AO7", file = "AO7.bedgraph")
convert2bedgraph(df = df_snp_counts, sample = "AO5", filter = "A1O5_editing", file = "A1O5_editing.bedgraph")
convert2bedgraph(df = df_snp_counts, sample = "AO7", filter = "A3O7_editing", file = "A3O7_editing.bedgraph")

## Get the summary file
summarize_bedgraph <- function(df) {
  
  # Create an empty list to store results for each gene
  results <- list()
  
  # Loop over each unique gene.symbol in the data
  for (gene in unique(df$gene.symbol)) {
    # Subset the data for this specific gene
    gene_data <- df[df$gene.symbol == gene, ]
    
    # Count the number of editing sites for this gene
    num_edit_sites <- nrow(gene_data)
    
    # Calculate the average editing percentage (score)
    avg_editing <- mean(gene_data$score)
    
    # Concatenate the edit_percentage and read info
    edit_percent_read_str <- paste0(round(gene_data$score, 1), "%_", gene_data$end - gene_data$start, "r", collapse = ",")
    
    # Concatenate the gene feature information
    gene_features <- paste(gene_data$annotation, collapse = ",")
    
    # Create a unique identifier string for the gene (seqnames and start-end positions)
    site_identifiers <- paste(gene_data$seqnames, gene_data$start, sep = "_", collapse = ";")
    
    # Store the result in a list
    results[[gene]] <- data.frame(
      Gene_name = gene,
      Num_edit_sites = num_edit_sites,
      Avg_editing = round(avg_editing, 1),
      Features = gene_features,
      Edit_percent_read_str = edit_percent_read_str,
      Identifier_str = site_identifiers
    )
  }
  
  # Combine all results into a single data frame
  summary_df <- do.call(rbind, results)
  
  summary_df = summary_df[order(summary_df$Avg_editing, summary_df$Num_edit_sites, decreasing = TRUE),]
  
  return(summary_df)
}

df <- read.table("A3O7_editing.bedgraph", header = FALSE, sep = "\t", col.names = c("seqnames", "start", "end", "score", "gene.symbol", "strand", "annotation"))
summary_result <- summarize_bedgraph(df)
write.xlsx(df, file = "A3O7_editing_summary.xlsx")


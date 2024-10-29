## This file contains script using thresholds to find out RNA editing site.
## The criteria used to detect RNA editing site (https://github.com/rosbashlab/HyperTRIBE/blob/master/CODE/find_rnaeditsites.pl#L132):
## 1. Total coverage in gDNA (Totalcount) must be greater than 9 reads to ensure sufficient data quality.
## 2. Proportion of G counts in gDNA must be less than 0.005 (i.e., less than 0.5% of total reads), indicating that the G base is rare in the genome.
## 3. Proportion of A counts in gDNA must be greater than or equal to 0.8 (i.e., A dominates the position in gDNA).
## 4. G counts in RNA must be present (greater than 0) to indicate potential RNA editing.
## The criteria used to further filter RNA editing site: (https://github.com/rosbashlab/HyperTRIBE/blob/master/CODE/filter_by_threshold_without_header.pl)
## the editing sites have at least 5% editing and at least a coverage of 20 reads 

# Libraries required
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("bbmle")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("statmod")
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

# Create a subfolder called 'step4' if it doesn't exist
output_dir <- "step4"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Using pre-defined threshold to define the RNA editing site detection function with customizable thresholds and column names
detect_rna_editing_site <- function(df, 
                                    A_count_gDNA_col, 
                                    G_count_gDNA_col, 
                                    A_count_RNA_col, 
                                    G_count_RNA_col,
                                    gDNA_coverage_threshold = 9, 
                                    gDNA_editing_threshold = 0.005, 
                                    A_count_threshold = 0.8, 
                                    RNA_editing_threshold = 0.05,
                                    RNA_coverage_threshold = 20) {
  
  # Create a new column in the dataframe to store the results
  is_editing_site <- with(df, {
    # Access the specified columns dynamically
    A_count_gDNA <- df[[A_count_gDNA_col]]
    G_count_gDNA <- df[[G_count_gDNA_col]]
    A_count_RNA <- df[[A_count_RNA_col]]
    G_count_RNA <- df[[G_count_RNA_col]]
    
    # Calculate the total gDNA coverage
    gDNA_coverage <- A_count_gDNA + G_count_gDNA
    RNA_coverage <- A_count_RNA + G_count_RNA
    
    # Apply the criteria step-by-step with the customizable thresholds
    is_gDNA_coverage_sufficient <- gDNA_coverage > gDNA_coverage_threshold
    is_gDNA_editing_low <- (G_count_gDNA / gDNA_coverage) < gDNA_editing_threshold
    is_A_count_dominant <- (A_count_gDNA / gDNA_coverage) >= A_count_threshold
    is_RNA_editing_present <- G_count_RNA > 0
    
    # Further selection on RNA editing site
    is_RNA_coverage_sufficient <- RNA_coverage > RNA_coverage_threshold
    is_RNA_editing_sufficient <- (G_count_RNA / RNA_coverage) > RNA_editing_threshold
    
    # Return TRUE if all conditions are satisfied, FALSE otherwise
    is_editing_site <- is_gDNA_coverage_sufficient & is_gDNA_editing_low & is_A_count_dominant & is_RNA_editing_present & is_RNA_coverage_sufficient & is_RNA_editing_sufficient
    return(is_editing_site)
  })
  
  # Return the modified dataframe with the RNA editing site classification
  return(is_editing_site)
}

df_snp_counts$A1O5_editing <- detect_rna_editing_site(df_snp_counts, A_count_gDNA_col = "A1_dedupped.bam.ref.count", G_count_gDNA_col = "A1_dedupped.bam.alt.count", 
                                             A_count_RNA_col = "AO5_dedupped.bam.ref.count", G_count_RNA_col = "AO5_dedupped.bam.alt.count", 
                                             gDNA_coverage_threshold = 9, gDNA_editing_threshold = 0.005, A_count_threshold = 0.8, RNA_editing_threshold = 0.05, RNA_coverage_threshold = 20)
df_snp_counts$A3O7_editing <- detect_rna_editing_site(df_snp_counts, A_count_gDNA_col = "A3_dedupped.bam.ref.count", G_count_gDNA_col = "A3_dedupped.bam.alt.count", 
                                                      A_count_RNA_col = "AO7_dedupped.bam.ref.count", G_count_RNA_col = "AO7_dedupped.bam.alt.count", 
                                                      gDNA_coverage_threshold = 9, gDNA_editing_threshold = 0.005, A_count_threshold = 0.8, RNA_editing_threshold = 0.05, RNA_coverage_threshold = 20)

                              
write.xlsx(df_snp_counts, file = file.path(output_dir, "RNAedit.xlsx"))

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
    
    # Calculate the total coverage (ref count + alt count)
    df$coverage <- df[[ref_col]] + df[[alt_col]]
    
    # Calculating the score as the percentage of alternative reads
    df$score <- with(df, ifelse(df$coverage > 0, 
                                df[[alt_col]] / df$coverage * 100, 
                                0))
    
    # Prepare the data frame for output
    output_df <- df[c("seqnames", "start", "end", "score", "gene.symbol", "strand", "annotation", "coverage")]
    
    # Write to bedgraph file
    write.table(output_df, 
                file = file.path(output_dir, file), 
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
    
    # Calculate the average coverage
    avg_coverage <- mean(gene_data$coverage)
    
    # Concatenate the edit_percentage and read info
    edit_percent_read_str <- paste0(round(gene_data$score, 1), "%_", gene_data$coverage, "r", collapse = ",")
    
    # Concatenate the gene feature information
    gene_features <- paste(gene_data$annotation, collapse = ",")
    
    # Create a unique identifier string for the gene (seqnames and start-end positions)
    site_identifiers <- paste(gene_data$seqnames, gene_data$start, sep = "_", collapse = ";")
    
    # Store the result in a list
    results[[gene]] <- data.frame(
      Gene_name = gene,
      Num_edit_sites = num_edit_sites,
      Avg_editing = round(avg_editing, 1),
      Avg_coverage = round(avg_coverage, 1),
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

df <- read.table(file.path(output_dir, "A1O5_editing.bedgraph"), header = FALSE, sep = "\t", col.names = c("seqnames", "start", "end", "score", "gene.symbol", "strand", "annotation", "coverage"))
summary_result <- summarize_bedgraph(df)
write.xlsx(summary_result, file = file.path(output_dir, "A1O5_editing_summary_gene.xlsx"))

df <- read.table(file.path(output_dir, "A3O7_editing.bedgraph"), header = FALSE, sep = "\t", col.names = c("seqnames", "start", "end", "score", "gene.symbol", "strand", "annotation", "coverage"))
summary_result <- summarize_bedgraph(df)
write.xlsx(summary_result, file = file.path(output_dir, "A3O7_editing_summary_gene.xlsx"))

# Write the content to the file

parameters <- "gDNA_coverage_threshold = 9, gDNA_editing_threshold = 0.005, A_count_threshold = 0.8, RNA_editing_threshold = 0.05, RNA_coverage_threshold = 20"
writeLines(parameters, con = file.path(output_dir, "parameters.txt"))
## ----setup------------------------------------------------------------------------------------------------------------
.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/", "file_paths_and_colours.R"))
library(dplyr)
library(data.table)
library(MatrixEQTL)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges
)

## ----expression_data--------------------------------------------------------------------------------------------------
expression_data <- fread("../RNA_counts/n579_TPM_combat_adjusted_counts_log2.txt")
# Set the row names using the first column and remove it from the data frame
rownames(expression_data) <- expression_data$V1
expression_data$V1 <- NULL
column_order <- colnames(expression_data)


## ----covariates_data--------------------------------------------------------------------------------------------------
covar <- read.delim("../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes_imputed_final.txt")
meta <- read.delim("../metadata/n579_metadata_matched_for_eQTL.txt")

# Assume 'batch' is a factor variable representing batch IDs
batch <- factor(meta$flowcell_and_lane)  

# Create the model matrix for the batch variable
#batch_matrix <- model.matrix(~ batch - 1)
# Convert to a data frame or matrix as needed
#batch_matrix <- as.data.frame(batch_matrix)
# Include these columns in your covariate matrix for MatrixEQTL
#meta$batch <- batch_matrix


table(meta$TYPE)
meta <- meta[, c(1, 9, 44)] # RIn 44, disease status 61, batch 52-59
#meta_df <- cbind(batch_matrix, meta)
#df  <- merge(covar, meta_df, by = "IID") # combine covariates and metadata so that RIN is now included
df  <- merge(covar, meta, by = "IID") # combine covariates and metadata so that RIN is now included
# Remove FID and pheno column 
df$FID <- NULL
df$Pheno <- NULL

# Further define the sample types by A-T-S status 
df$TYPE.ATS <- ifelse(df$TYPE == "LBD" & df$Braak <= 3 & df$Thal < 2, "LBD_S",
                             ifelse(df$TYPE == "LBD" & df$Braak > 3 & df$Thal >= 2, "LBD_ATS",
                                    ifelse(df$TYPE == "LBD" & df$Braak <= 3 & df$Thal >= 2, "LBD_AS", # high amyloid 
                                      ifelse(df$TYPE == "AD", "AD",
                                           ifelse(df$TYPE == "PA", "PA",
                                                  ifelse(df$TYPE == "CONTROL", "CONTROL", "LBD_TS") # high Braak
                                           )))))

# check
table(df$TYPE)
table(df$TYPE.ATS)
df$TYPE <- NULL
df$TYPE.ATS <- factor(df$TYPE.ATS, levels=c("CONTROL", "PA", "AD", "LBD_S", "LBD_TS", "LBD_AS", "LBD_ATS"))
table(df$TYPE.ATS) # inspect


## ----genotype_data----------------------------------------------------------------------------------------------------
CWOW_snps <- fread("../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC_genotype_matrix.txt")
colnames(CWOW_snps)

# Rename last column to snpid
setnames(CWOW_snps, old = colnames(CWOW_snps)[ncol(CWOW_snps)], new = "snpid")

# Remove allele suffix if present
CWOW_snps[, snpid := sub("_.*", "", snpid)]

# Move snpid column to the front
CWOW_snps <- CWOW_snps[, c("snpid", setdiff(names(CWOW_snps), "snpid")), with = FALSE]

# Optional sanity checks
stopifnot(!anyDuplicated(CWOW_snps$snpid))

cat("Rows (SNPs):", nrow(CWOW_snps), "\n")
cat("Columns (including snpid):", ncol(CWOW_snps), "\n")

# Write SNP matrix for MatrixEQTL
fwrite(
  CWOW_snps,
  "../snp_array/Filtered_n579_CWOW_genotype_data_no_allele_info.txt",
  sep = "\t"
)

# Save SNP IDs separately
fwrite(
  CWOW_snps[, .(snpid)],
  "../snp_array/CWOW_final_imputed_snp_ids_GRCh38.txt",
  sep = "\t"
)

## ----subset-----------------------------------------------------------------------------------------------------------
# Next we need the sample IDs (aka the IIDs) to be the columns and the rows to be the covariates
# This section reformats the data frame & subsets by disease type 
snp_ids <- read.delim("../snp_array/CWOW_final_imputed_snp_ids_GRCh38.txt")

# Define  subsets
subsets <- list(
  "LBD_control" = c("CONTROL", "LBD_S", "LBD_AS", "LBD_ATS", "LBD_TS"), # all LBD types
  "AD_control" = c("CONTROL", "AD"),
  "PA_control" = c("CONTROL", "PA"),
  "LBD_S_control" = c("CONTROL", "LBD_S"),
  "LBD_AS_control" = c("CONTROL", "LBD_AS"),
  "LBD_ATS_control" = c("CONTROL", "LBD_ATS")
)

# Iterate over each subset
for (subset_name in names(subsets)) {
  
  # Subset the data based on disease type
  df_subset <- df %>% 
    filter(TYPE.ATS %in% subsets[[subset_name]])
  
  # Filter the expression data to include only IIDs present in the subset
  expression_data_subset <- expression_data[, .SD, .SDcols = df_subset$IID]
  # Write the output to a file
  exp_output_file <- paste0("../RNA_counts/cwow_", subset_name, "_imputed.txt")
  write.table(expression_data_subset, exp_output_file, quote = FALSE, sep = "\t", row.names = TRUE)
  column_order <- colnames(expression_data_subset)
  
  df_subset$Pheno <- ifelse(df_subset$TYPE.ATS == "CONTROL", 1, 2) # 1=unaffected (control), 2=affected (case)
  df_subset$TYPE.ATS <- NULL
  # Reformat the data 
  df_converted <- df_subset %>%
    mutate(across(-IID, as.character))
  
  df_long <- df_converted %>%
    pivot_longer(cols = -IID, names_to = "Variable", values_to = "Value")
  
  df_wide <- df_long %>%
    pivot_wider(names_from = IID, values_from = Value)
  
  # Reorder based on expression data
  df1_reordered <- df_wide[, column_order]
  
  # Convert to numeric
  data_numeric <- data.frame(lapply(df1_reordered, function(x) {
    as.numeric(as.character(x))
  }))
  
  # Combine covariates and phenotype data
  covars <- cbind(df_wide$Variable, data_numeric)
  
  # Calculate variance for each row
 row_var <- apply(covars[, -1], 1, var) # ignore the first column of the variables

  # Filter out rows with zero variance
  cvrt_filtered <- covars[row_var != 0, ]

  # Check the filtered covariate matrix
  apply(cvrt_filtered[, -1], 1, var)
  names(cvrt_filtered)[1] <- "variable"
  
  # Write the output to a file
  cov_output_file <- paste0("../snp_array/covariates_and_phenotype_files/cwow_covariates_", subset_name, "_imputed.txt")
  write.table(cvrt_filtered, cov_output_file, quote = FALSE, sep = "\t", row.names = FALSE)
  
  CWOW_snps_reordered <- CWOW_snps[, ..column_order] # reorder to match the order of the expression data  
  rownames(CWOW_snps_reordered) <- snp_ids$x
  snp_output_file <- paste0("../snp_array/cwow_", subset_name, "_imputed.txt")
  fwrite(CWOW_snps_reordered, snp_output_file, sep = " ")
}


## ----SNP_and_gene_data------------------------------------------------------------------------------------------------
# Read in gene information. These are the genes used in the bulk RNAseq differential expression analysis.
genes <- fread("../RNA_counts/genes_filtered_protein_coding.txt")
gene_info <- genes[,c(1:4)] # Keep only columns 1 - 4
colnames(gene_info) <- c("geneid","chr", "left", "right") # Rename the columns 

# Read in SNP annotation
# The SNP annotation was download from the dbSNP NCBI website 
# The vcf file was processed using a python script to obtain chromosome, position, SNP id, and gene name
snp_anno <- fread("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/snp_array/reference/SNP_positions_with_gene.tsv", 
                  header = TRUE,  fill = TRUE)
# Make an alleles column 
snp_anno$Alleles <- paste(snp_anno$REF, snp_anno$ALT, sep = "/")

# Convert snp annotation to GRanges object
gr <- GRanges(
  seqnames = paste0("chr", snp_anno$`#CHROM`),
  ranges = IRanges(start = snp_anno$POS, end = snp_anno$POS),
  names = snp_anno$ID,
  alleles = snp_anno$Alleles, 
  Gene = snp_anno$GENE
  )

# Lift SNP annotation to GRCh38
chain <- import.chain("../snp_array/reference/hg19ToHg38.over.chain")
lifted <- liftOver(gr, chain)
lifted_gr <- unlist(lifted)

GRCh38_snp_annotation <- as.data.frame(lifted_gr)
write.table(GRCh38_snp_annotation,"../snp_array/reference/NCBI_SNP_positions_with_gene_hg19ToHg38_lift_over.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#------ End
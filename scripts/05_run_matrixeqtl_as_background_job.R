# Setup
setwd("/research/labs/neurology/fryer/m239830/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")
library(dplyr)
library(data.table)
library(MatrixEQTL)
library(tidyverse)

saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}


# Read inputs for Matrix eQTL 
#Inputs include:\
#genotype - columns are samples, rows are snps. 0 is homozygous reference, 1 is heterozygous, 2 is homozygous alternate 
#expression - columns are samples, rows are genes. Values are TPM expression 
#covariates - columns are samples, rows are phenotype information such as sex, age, and RIN
#gene info - gene_id, chromosome, start, and end positions
#snp annotation - snpids and the corresponding gene_id 

# Read in gene information. These are the genes used in the bulk RNAseq differential expression analysis.
genes <- fread("../RNA_counts/genes_filtered.txt")
gene_info <- genes[,c(1:4)] # Keep only columns 1 - 4
colnames(gene_info) <- c("geneid","chr", "left", "right") # Rename the columns 
# Read in SNP annotation
# The SNP annotation was download from the illumina website 
snp_anno <- fread("../snp_array/illumina_OminExpress24/InfiniumOmniExpress-24v1-3_A1.annotated.txt", 
                  header = TRUE,  fill = TRUE)
snp_anno$gene <- sub(",.*", "", snp_anno$`Gene(s)`) # create a gene column 
snp_anno <- snp_anno[,c(1:3)] # Keep only columns 1 - 3
colnames(snp_anno) <- c("snpid","chr", "pos") # rename the columns 

# Thresholds, Only associations significant at this level will be saved
pvOutputThreshold = 1e-8 
# The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. 
# For larger datasets the threshold should be lower. 
# Setting the threshold to a high value for a large dataset may cause excessively large output files.
useModel = modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS

snps = SlicedData$new()
snps$fileDelimiter = " "      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snps$LoadFile("../snp_array/Filtered_n579_CWOW_genotype_data_no_allele_info.txt")

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile("../RNA_counts/n579_TPM_combat_adjusted_counts.txt")

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
cvrt$LoadFile("../snp_array/covariates_and_phenotype_files/covariates.txt")


# Order check 
# Make sure that the gene expression (gene) and genetic datasets (snps) have the samples in the same order!

stopifnot(all(colnames(snps) == colnames(gene)))


# Multiple testing 
# To adjust for the number of tests, we could divide 0.05 by the number of tests to be performed. 

snps_rows <- as.numeric(nrow(snps))
gene_rows <- as.numeric(nrow(gene))

# How many tests
snps_rows*gene_rows 

lower.bound <- nrow(snps)
upper.bound <- nrow(snps) * gene_rows
threshold <- 0.05/lower.bound


# Matrix eQTL 
# Each significant gene-SNP association is recorded in a separate line in the output file and in the returned object me.
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name =  file.path("../MatrixEQTL/n579_all_eQTL"),
  pvOutputThreshold = pvOutputThreshold, 
  useModel = modelLINEAR, ## test using linear models 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
## Plot the histogram of all p-values
plot(me)
saveToPDF("../MatrixEQTL/qqplot_pvalues_eqtl.pdf", width = 7, height = 6)
# Save as RDS 
saveRDS(me, "../MatrixEQTL/n579_eQTL.Rds")

# Run matrix eqtl and loop through the disease groups
# I.e only run matrix eQTL for the AD and control samples

# Setup
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/", "file_paths_and_colours.R"))

#-------------------------------
# Read inputs for Matrix eQTL 
#-------------------------------
# Inputs include:\
# genotype - columns are samples, rows are snps. 0 is homozygous reference, 1 is heterozygous, 2 is homozygous alternate 
# expression - columns are samples, rows are genes. Values are TPM expression 
# covariates - columns are samples, rows are phenotype information such as sex, age, and RIN
# gene info - gene_id, chromosome, start, and end positions
# snp annotation - snpids and the corresponding gene_id 

#------------
# Gene information 
#------------
# Read in gene information. These are the genes used in the bulk RNAseq differential expression analysis.
genes <- fread("../RNA_counts/genes_filtered.txt")
gene_info <- genes[,c(1:4)] # Keep only columns 1 - 4
colnames(gene_info) <- c("geneid","chr", "left", "right") # Rename the columns 

#------------
# SNP annotation
#------------
# The SNP annotation was download from the illumina website and lifted over to GRCh38
snp_anno <- fread("../snp_array/reference/NCBI_SNP_positions_with_gene_hg19ToHg38_lift_over.tsv", 
                  header = TRUE,  fill = TRUE)

snp_anno <- snp_anno[,c(6,1,2)] # Keep 
colnames(snp_anno) <- c("snpid","chr", "pos") # rename the columns 
#------------
# Matrix eQTL thresholds
#------------
# Thresholds, Only associations significant at this level will be saved
pvOutputThreshold = 1e-5
# The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. 
# For larger datasets the threshold should be lower. 
# Setting the threshold to a high value for a large dataset may cause excessively large output files.
useModel = modelLINEAR_CROSS # modelANOVA or modelLINEAR or modelLINEAR_CROSS

#------------
# Read in the matrix eQTL inputs 
#------------

# Begin the loop
# Define the subsets
subsets <- c(
  "LBD_control", # all LBD samples
  "AD_control", 
  "PA_control",
  "LBD_S_control",
  "LBD_AS_control",
  "LBD_ATS_control"
)

# Loop through each subset
for (subset in subsets) {
  
  # Define input file paths based on the subset
  snp_file <- paste0("../snp_array/Filtered_n580_CWOW_genotype_data_", subset, "_imputed.txt")
  gene_file <- paste0("../RNA_counts/n580_TPM_combat_adjusted_counts_log2_", subset, "_imputed.txt")
  cvrt_file <- paste0("../snp_array/covariates_and_phenotype_files/covariates_", subset, "_imputed.txt")
  
  # Load SNP data
  snps <- SlicedData$new()
  snps$fileDelimiter <- " "      
  snps$fileOmitCharacters <- "NA"
  snps$fileSkipRows <- 1
  snps$fileSkipColumns <- 1
  snps$fileSliceSize <- 2000
  snps$LoadFile(snp_file)
  
  # Load gene expression data
  gene <- SlicedData$new()
  gene$fileDelimiter <- "\t"
  gene$fileOmitCharacters <- "NA"
  gene$fileSkipRows <- 1
  gene$fileSkipColumns <- 1
  gene$fileSliceSize <- 2000
  gene$LoadFile(gene_file)
  
  # Load covariate data
  cvrt <- SlicedData$new()
  cvrt$fileDelimiter <- "\t"
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows <- 1
  cvrt$fileSkipColumns <- 1
  cvrt$fileSliceSize <- 2000
  cvrt$LoadFile(cvrt_file)
  
  # Ensure the samples are in the same order between SNPs and gene expression data and covariates
  stopifnot(all(colnames(snps) == colnames(gene)))
  stopifnot(all(colnames(cvrt) == colnames(gene)))
  
  # Calculate number of tests and threshold
  # To adjust for the number of tests, we could divide 0.05 by the number of tests to be performed. 
  snps_rows <- as.numeric(nrow(snps))
  gene_rows <- as.numeric(nrow(gene))
  
  # How many tests
  tests <- snps_rows*gene_rows 
  
  lower.bound <- nrow(snps)
  upper.bound <- nrow(snps) * gene_rows
  threshold <- 0.05/lower.bound
  
  # Define output file paths
  output_file_trans <- file.path("../MatrixEQTL/", paste0("trans_eQTL_", subset, "_imputed_modelLINEAR_CROSS"))
  output_file_cis <- file.path("../MatrixEQTL/", paste0("cis_eQTL_", subset, "_imputed_modelLINEAR_CROSS"))
  
  # Run the eQTL analysis
  dis_me <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_trans,
    pvOutputThreshold = 1e-5,
    useModel = useModel, 
    pvOutputThreshold.cis = 1, 
    output_file_name.cis = output_file_cis,
    snpspos = snp_anno, 
    genepos = gene_info[, c("geneid", "chr", "left", "right")],
    cisDist = 1e6, #1e6 1 million
    verbose = TRUE,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )
  
  # Create results table 
  results <- data.frame(
    Description = c(
      "Subset",
      "Analysis time (seconds)",
      "Number of tests (snps_rows * gene_rows)",
      "Multiple testing threshold",
      "Detected local eQTLs",
      "Detected distant eQTLs"
    ),
    Value = c(
      subset,
      dis_me$time.in.sec,
      tests,
      threshold,
      nrow(dis_me$cis$eqtls),
      nrow(dis_me$trans$eqtls)
    )
  )
  
  # Save the results table as a CSV file named according to the subset
  output_table <- paste0("../MatrixEQTL/summary_", subset, "_imputed_modelLINEAR_CROSS.txt")
  write.table(results, file = output_table, row.names = FALSE, quote = FALSE)
  
}
# End 
#-------------------------------------------------------------------------------
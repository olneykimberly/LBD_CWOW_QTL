#!/usr/bin/env Rscript

# ===============================================================
# MatrixEQTL interaction eQTL pipeline
# Loops through subset-specific comparisons
# Uses subset-specific SNP / expression / covariate files
# Includes full input validation before each run
#
# modelLINEAR_CROSS tests genotype x Pheno interaction
# so Pheno MUST be the LAST covariate row.
# ===============================================================

.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))

suppressPackageStartupMessages({
  library(data.table)
  library(MatrixEQTL)
  library(ggplot2)
})

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")


dir.create("../MatrixEQTL/results", recursive = TRUE, showWarnings = FALSE)
dir.create("../MatrixEQTL/results/logs", recursive = TRUE, showWarnings = FALSE)

# ===============================================================
# Comparisons to run
# LBD_TS_control removed as requested
# ===============================================================

subsets <- c(
#  "LBD_control",
  "AD_control"
#  "PA_control",
#  "LBD_S_control",
#  "LBD_AS_control",
#  "LBD_ATS_control"
)

# ===============================================================
# Global settings
# ===============================================================

pvOutputThreshold <- 1e-8
useModel <- modelLINEAR_CROSS

# ===============================================================
# Loop through subsets
# ===============================================================

summary_list <- list()

for (subset_name in subsets) {
  
  cat("\n====================================================\n")
  cat("Running MatrixEQTL for subset:", subset_name, "\n")
  cat("====================================================\n\n")
  
  # -------------------------------------------------------------
  # File paths
  # -------------------------------------------------------------
  snp_file <- paste0("../MatrixEQTL/inputs/cwow_", subset_name, "_snps.txt")
  expr_file <- paste0("../MatrixEQTL/inputs/cwow_", subset_name, "_expression_invnorm.txt")
  covar_file <- paste0("../MatrixEQTL/inputs/cwow_", subset_name, "_covariates_invnorm.txt")
  
  out_prefix <- paste0("../MatrixEQTL/results/cwow_", subset_name, "_eQTL")
  log_file <- paste0("../MatrixEQTL/results/logs/cwow_", subset_name, "_input_summary.txt")

  }
  # -------------------------------------------------------------
  # Load MatrixEQTL sliced data
  # -------------------------------------------------------------
  cat("Loading SlicedData objects...\n")
  
  snps <- SlicedData$new()
  snps$fileDelimiter <- " "
  snps$fileOmitCharacters <- "NA"
  snps$fileSkipRows <- 1
  snps$fileSkipColumns <- 1
  snps$fileSliceSize <- 2000
  snps$LoadFile(snp_file)
  
  gene <- SlicedData$new()
  gene$fileDelimiter <- "\t"
  gene$fileOmitCharacters <- "NA"
  gene$fileSkipRows <- 1
  gene$fileSkipColumns <- 1
  gene$fileSliceSize <- 2000
  gene$LoadFile(expr_file)
  
  cvrt <- SlicedData$new()
  cvrt$fileDelimiter <- "\t"
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows <- 1
  cvrt$fileSkipColumns <- 1
  cvrt$fileSliceSize <- 2000
  cvrt$LoadFile(covar_file)
  
  # Confirm sample order again inside MatrixEQTL objects
  stopifnot(all(colnames(snps) == colnames(gene)))
  stopifnot(all(colnames(gene) == colnames(cvrt)))
  
  # -------------------------------------------------------------
  # Multiple testing context
  # -------------------------------------------------------------
  snps_rows <- as.numeric(nrow(snps))
  gene_rows <- as.numeric(nrow(gene))
  threshold <- 0.05 / snps_rows
  
  cat("Bonferroni lower-bound threshold:", threshold, "\n")
  cat("Running MatrixEQTL...\n\n")
  
  # -------------------------------------------------------------
  # Run MatrixEQTL
  # -------------------------------------------------------------
  me <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = paste0(out_prefix, ".txt"),
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )
  
  # -------------------------------------------------------------
  # Save results
  # -------------------------------------------------------------
  cat("\nAnalysis complete for", subset_name, "\n")
  cat("Runtime:", me$time.in.sec, "seconds\n\n")
  
  pdf(paste0(out_prefix, "_qqplot.pdf"), width = 7, height = 6)
  plot(me)
  dev.off()
  
  saveRDS(
    me,
    paste0(out_prefix, ".Rds")
  )
  
  summary_list[[subset_name]] <- data.table(
    subset = subset_name,
    samples = n_samples,
    genes = n_genes,
    snps = n_snps,
    covariates = n_covars,
    tests = n_tests,
    runtime_sec = me$time.in.sec
  )
  
  # Clean up between runs
  rm(snps, gene, cvrt, me)
  gc()


# ===============================================================
# Save run summary
# ===============================================================

if (length(summary_list) > 0) {
  run_summary <- rbindlist(summary_list, fill = TRUE)
  fwrite(
    run_summary,
    "../MatrixEQTL/results/MatrixEQTL_run_summary.tsv",
    sep = "\t"
  )
  
  cat("\n====================================================\n")
  cat("All completed runs summary:\n")
  print(run_summary)
  cat("====================================================\n")
} else {
  cat("\nNo subsets were run.\n")
}
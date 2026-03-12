#!/usr/bin/env Rscript

# ======================================================================
# Validate and write MatrixEQTL-ready subset-specific SNP, expression,
# and covariate files with identical sample columns in identical order.
#
# Inputs:
#   1) Full SNP matrix (rows = SNPs, cols = samples)
#   2) Subset-specific expression files already created
#   3) Subset-specific covariate files already created
#
# Outputs for each subset:
#   - SNP matrix reordered to match expression/covariates
#   - expression matrix reordered/validated
#   - covariate matrix reordered/validated
#   - sample manifest
#
# This script assumes:
#   - first column of SNP/expression files is row ID
#   - first column of covariate file is "variable"
# ======================================================================

.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))

suppressPackageStartupMessages({
  library(data.table)
})

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")

dir.create("../MatrixEQTL/inputs", recursive = TRUE, showWarnings = FALSE)
dir.create("../MatrixEQTL/inputs/manifests", recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# Master full SNP matrix
# ----------------------------------------------------------------------
# Update this path if needed. This should be the full 579-sample SNP matrix
# already prepared for eQTL and using rsIDs / final SNP IDs as row names.
full_snp_file <- "../snp_array/Filtered_n579_CWOW_genotype_data_no_allele_info.txt"

# ----------------------------------------------------------------------
# Subsets to validate and write
# ----------------------------------------------------------------------
subsets <- c(
  "AD_control",
  "PA_control",
  "LBD_S_control",
  "LBD_AS_control",
  "LBD_ATS_control",
  "LBD_TS_control",
  "LBD_control"
)

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
read_matrix_with_ids <- function(file, id_col_name = NULL) {
  dt <- fread(file)
  first_col <- colnames(dt)[1]
  if (!is.null(id_col_name)) setnames(dt, first_col, id_col_name)
  dt
}

get_sample_cols <- function(dt) {
  colnames(dt)[-1]
}

assert_no_duplicates <- function(x, label) {
  dups <- unique(x[duplicated(x)])
  if (length(dups) > 0) {
    stop(sprintf(
      "%s contains duplicated sample IDs. First duplicates: %s",
      label,
      paste(head(dups, 10), collapse = ", ")
    ))
  }
}

# ----------------------------------------------------------------------
# Load full SNP matrix once
# ----------------------------------------------------------------------
cat("Reading full SNP matrix...\n")
snp_full <- read_matrix_with_ids(full_snp_file, id_col_name = "snpid")

full_snp_samples <- get_sample_cols(snp_full)
assert_no_duplicates(full_snp_samples, "Full SNP matrix")

cat("Full SNP matrix dimensions:\n")
cat("  SNPs:", nrow(snp_full), "\n")
cat("  Samples:", length(full_snp_samples), "\n\n")

# ----------------------------------------------------------------------
# Process each subset
# ----------------------------------------------------------------------
for (subset_name in subsets) {
  cat("============================================================\n")
  cat("Processing subset:", subset_name, "\n")
  
  expr_file <- paste0("../RNA_counts/eQTL_ready/subsets/cwow_", subset_name, "_expression_invnorm.txt")
  covar_file <- paste0("../snp_array/covariates_and_phenotype_files/eQTL_ready/cwow_covariates_", subset_name, "_invnorm.txt")
  
  out_expr_file <- paste0("../MatrixEQTL/inputs/cwow_", subset_name, "_expression_invnorm.txt")
  out_covar_file <- paste0("../MatrixEQTL/inputs/cwow_", subset_name, "_covariates_invnorm.txt")
  out_snp_file <- paste0("../MatrixEQTL/inputs/cwow_", subset_name, "_snps.txt")
  out_manifest_file <- paste0("../MatrixEQTL/inputs/manifests/cwow_", subset_name, "_sample_order.txt")
  
  # Check file existence
  if (!file.exists(expr_file)) stop("Missing expression file: ", expr_file)
  if (!file.exists(covar_file)) stop("Missing covariate file: ", covar_file)
  
  # Read expression and covariates
  expr_dt <- read_matrix_with_ids(expr_file, id_col_name = "geneid")
  covar_dt <- read_matrix_with_ids(covar_file, id_col_name = "variable")
  
  expr_samples <- get_sample_cols(expr_dt)
  covar_samples <- get_sample_cols(covar_dt)
  
  assert_no_duplicates(expr_samples, paste0("Expression file for ", subset_name))
  assert_no_duplicates(covar_samples, paste0("Covariate file for ", subset_name))
  
  # Check expression vs covariates sample identity
  if (!setequal(expr_samples, covar_samples)) {
    only_expr <- setdiff(expr_samples, covar_samples)
    only_covar <- setdiff(covar_samples, expr_samples)
    
    stop(
      paste0(
        "Sample mismatch between expression and covariates for subset ", subset_name, ".\n",
        "Only in expression: ", paste(head(only_expr, 20), collapse = ", "), "\n",
        "Only in covariates: ", paste(head(only_covar, 20), collapse = ", ")
      )
    )
  }
  
  # Define master sample order from expression
  sample_order <- expr_samples
  
  # Reorder covariates to expression order
  covar_dt <- covar_dt[, c("variable", sample_order), with = FALSE]
  
  # Match SNPs to subset samples
  if (!all(sample_order %in% full_snp_samples)) {
    missing_in_snp <- setdiff(sample_order, full_snp_samples)
    stop(
      paste0(
        "Subset ", subset_name, " has samples missing from full SNP matrix:\n",
        paste(head(missing_in_snp, 20), collapse = ", ")
      )
    )
  }
  
  snp_subset_dt <- snp_full[, c("snpid", sample_order), with = FALSE]
  
  # Final order checks
  stopifnot(identical(get_sample_cols(expr_dt), sample_order))
  stopifnot(identical(get_sample_cols(covar_dt), sample_order))
  stopifnot(identical(get_sample_cols(snp_subset_dt), sample_order))
  
  # Write manifest
  manifest_dt <- data.table(IID = sample_order, order = seq_along(sample_order))
  fwrite(manifest_dt, out_manifest_file, sep = "\t")
  
  # Write validated outputs
  write.table(
    expr_dt,
    file = out_expr_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  write.table(
    covar_dt,
    file = out_covar_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  write.table(
    snp_subset_dt,
    file = out_snp_file,
    sep = " ",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  cat("Validated and wrote:\n")
  cat("  Expression:", out_expr_file, "\n")
  cat("  Covariates:", out_covar_file, "\n")
  cat("  SNPs:", out_snp_file, "\n")
  cat("  Samples:", length(sample_order), "\n")
  cat("  Genes:", nrow(expr_dt), "\n")
  cat("  SNPs:", nrow(snp_subset_dt), "\n\n")
}

cat("All subset-specific MatrixEQTL inputs written successfully.\n")
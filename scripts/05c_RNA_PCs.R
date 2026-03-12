## ----setup------------------------------------------------------------------------------------------------------------
.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/", "file_paths_and_colours.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

dir.create("../RNA_counts/eQTL_ready/subsets", recursive = TRUE, showWarnings = FALSE)
dir.create("../snp_array/covariates_and_phenotype_files/eQTL_ready", recursive = TRUE, showWarnings = FALSE)

## ----expression_data--------------------------------------------------------------------------------------------------
# Read inverse-normal transformed, filtered expression matrix
# first column = gene ID
# remaining columns = sample IIDs
expression_data <- fread("../RNA_counts/eQTL_ready/n579_TPM_combat_adjusted_counts_log2_filtered_invnorm.txt")

rownames(expression_data) <- expression_data[[1]]
expression_data[[1]] <- NULL

expression_sample_order <- colnames(expression_data)

## ----covariates_data--------------------------------------------------------------------------------------------------
# Genotype-side covariates
covar <- fread("../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes_imputed_final.txt")

# Matched metadata
meta <- fread("../metadata/n579_metadata_matched_for_eQTL.txt")

# Hidden covariates from RNA PCA
hidden_covars <- fread("../RNA_counts/eQTL_ready/n579_expression_hidden_covariates.txt")

# Clean RIN column
meta_sub <- copy(meta)

# Keep only metadata fields needed here
meta_sub <- meta_sub[, .(FID, IID, TYPE, RIN)]

# Merge all covariate sources
df <- merge(covar, meta_sub, by = c("FID", "IID"), all.x = TRUE)
df <- merge(df, hidden_covars, by = c("FID", "IID"), all.x = TRUE)

# Define finer disease groups by A/T/S status
df[, TYPE.ATS := fifelse(TYPE == "LBD" & Braak <= 3 & Thal < 2, "LBD_S",
                  fifelse(TYPE == "LBD" & Braak > 3 & Thal >= 2, "LBD_ATS",
                  fifelse(TYPE == "LBD" & Braak <= 3 & Thal >= 2, "LBD_AS",
                  fifelse(TYPE == "LBD" & Braak > 3 & Thal < 2, "LBD_TS",
                  fifelse(TYPE == "AD", "AD",
                  fifelse(TYPE == "PA", "PA",
                  fifelse(TYPE == "CONTROL", "CONTROL", NA_character_)))))))]

df[, TYPE.ATS := factor(TYPE.ATS,
                        levels = c("CONTROL", "PA", "AD", "LBD_S", "LBD_TS", "LBD_AS", "LBD_ATS"))]

cat("\nTYPE counts:\n")
print(table(df$TYPE, useNA = "ifany"))

cat("\nTYPE.ATS counts:\n")
print(table(df$TYPE.ATS, useNA = "ifany"))

# Build final per-sample covariate table
hidden_cols <- grep("^exprPC[0-9]+$", colnames(df), value = TRUE)

sample_covars <- df[, c(
  "FID", "IID", "TYPE", "TYPE.ATS",
  "Age", "Sex", "Braak", "Thal", "CingLB",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "RIN",
  hidden_cols
), with = FALSE]


# Keep only samples present in the expression matrix and reorder to match expression columns
sample_covars <- sample_covars[IID %in% colnames(expression_data)]
sample_covars <- sample_covars[match(colnames(expression_data), sample_covars$IID)]

stopifnot(all(sample_covars$IID == colnames(expression_data)))

## ----subset-----------------------------------------------------------------------------------------------------------
subsets <- list(
  "LBD_control"     = c("CONTROL", "LBD_S", "LBD_AS", "LBD_ATS", "LBD_TS"),
  "AD_control"      = c("CONTROL", "AD"),
  "PA_control"      = c("CONTROL", "PA"),
  "LBD_S_control"   = c("CONTROL", "LBD_S"),
  "LBD_AS_control"  = c("CONTROL", "LBD_AS"),
  "LBD_ATS_control" = c("CONTROL", "LBD_ATS"),
  "LBD_TS_control"  = c("CONTROL", "LBD_TS")
)

for (subset_name in names(subsets)) {
  
  cat("\nProcessing subset:", subset_name, "\n")
  
  df_subset <- sample_covars %>%
    filter(TYPE.ATS %in% subsets[[subset_name]])
  
  # Expression subset
  expression_data_subset <- expression_data[, df_subset$IID, with = FALSE]
  stopifnot(all(colnames(expression_data_subset) == df_subset$IID))
  
  # Write expression file
  exp_output_file <- paste0("../RNA_counts/eQTL_ready/subsets/cwow_", subset_name, "_expression_invnorm.txt")
  write.table(
    expression_data_subset,
    exp_output_file,
    quote = FALSE,
    sep = "\t",
    row.names = TRUE,
    col.names = NA
  )
  
  # Binary phenotype for subset
  df_subset <- as.data.table(df_subset)
  df_subset[, Pheno := fifelse(TYPE.ATS == "CONTROL", 1, 2)]
  
  # Drop grouping columns not intended as covariates
  df_subset[, c("TYPE", "TYPE.ATS") := NULL]
  
  # MatrixEQTL covariate file format: rows = covariates, columns = samples
  covar_long <- melt(
    df_subset[, !c("FID"), with = FALSE],
    id.vars = "IID",
    variable.name = "variable",
    value.name = "value"
  )
  
  covar_wide <- dcast(covar_long, variable ~ IID, value.var = "value")
  
  # Remove zero-variance covariates
  covar_mat <- as.data.frame(covar_wide)
  row_var <- apply(covar_mat[, -1, drop = FALSE], 1, function(x) var(as.numeric(x), na.rm = TRUE))
  covar_mat <- covar_mat[row_var != 0, , drop = FALSE]
  names(covar_mat)[1] <- "variable"
  
  # Write covariate file
  cov_output_file <- paste0("../snp_array/covariates_and_phenotype_files/eQTL_ready/cwow_covariates_", subset_name, "_invnorm.txt")
  write.table(
    covar_mat,
    cov_output_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  # Write sample manifest
  manifest_output_file <- paste0("../snp_array/covariates_and_phenotype_files/eQTL_ready/cwow_samples_", subset_name, "_invnorm.txt")
  fwrite(
    df_subset[, .(FID, IID, Pheno)],
    manifest_output_file,
    sep = "\t"
  )
  
  cat("  Samples:", ncol(expression_data_subset), "\n")
  cat("  Genes:", nrow(expression_data_subset), "\n")
  cat("  Covariates:", nrow(covar_mat), "\n")
}

## ----optional_full_set------------------------------------------------------------------------------------------------
# Full matched expression matrix
write.table(
  expression_data,
  "../RNA_counts/eQTL_ready/n579_expression_invnorm_full.txt",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

full_covars <- copy(sample_covars)
full_covars[, Pheno := fifelse(TYPE.ATS == "CONTROL", 1, 2)]
full_covars[, c("TYPE", "TYPE.ATS") := NULL]

full_covar_long <- melt(
  full_covars[, !c("FID"), with = FALSE],
  id.vars = "IID",
  variable.name = "variable",
  value.name = "value"
)

full_covar_wide <- dcast(full_covar_long, variable ~ IID, value.var = "value")
full_covar_mat <- as.data.frame(full_covar_wide)

row_var_full <- apply(full_covar_mat[, -1, drop = FALSE], 1, function(x) var(as.numeric(x), na.rm = TRUE))
full_covar_mat <- full_covar_mat[row_var_full != 0, , drop = FALSE]
names(full_covar_mat)[1] <- "variable"

write.table(
  full_covar_mat,
  "../snp_array/covariates_and_phenotype_files/eQTL_ready/n579_covariates_invnorm_full.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

cat("\nDone. Wrote subset-specific expression and covariate files for eQTL.\n")
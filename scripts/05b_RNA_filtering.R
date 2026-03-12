#!/usr/bin/env Rscript

.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(ggplot2)
  library(pheatmap)
})

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")

dir.create("../results/eQTL_prep", recursive = TRUE, showWarnings = FALSE)
dir.create("../RNA_counts/eQTL_ready", recursive = TRUE, showWarnings = FALSE)

#-----------------------------
# input files
#-----------------------------
expr_file <- "../RNA_counts/n579_TPM_combat_adjusted_counts_log2.txt"
meta_file <- "../metadata/n579_metadata.txt"

#-----------------------------
# read metadata
#-----------------------------
meta <- fread(meta_file)

#-----------------------------
# read expression matrix
#-----------------------------
expr <- fread(expr_file)

# assume first column is gene ID
gene_ids <- expr[[1]]
expr_df <- as.data.frame(expr)
rownames(expr_df) <- gene_ids
expr_df[[1]] <- NULL

# keep only metadata samples and reorder to match metadata IID
expr_df <- expr_df[, meta$IID, drop = FALSE]

stopifnot(all(colnames(expr_df) == meta$IID))

expr_mat <- as.matrix(expr_df)
mode(expr_mat) <- "numeric"

cat("Initial expression matrix dimensions:\n")
print(dim(expr_mat))

#-----------------------------
# 1. Filter low-expression genes
#-----------------------------
# Since data are log2(TPM), convert threshold into log2 space.
# Keep genes with TPM > 0.1 in at least 20% of samples.
# log2(0.1 + small offset) is roughly -3.32; since your transform may vary,
# use a simple threshold on the log2 TPM values.
#
# If TPM was log2(TPM), then TPM > 0.1 corresponds to > log2(0.1) = -3.32.
# If TPM was log2(TPM + 1), nearly all values are >= 0, so use > 0 for TPM > 0.
#
# A practical conservative filter for log2 TPM data:
expr_threshold <- 0
min_prop <- 0.20
min_n <- ceiling(ncol(expr_mat) * min_prop)

keep_genes <- rowSums(expr_mat > expr_threshold, na.rm = TRUE) >= min_n
expr_filt <- expr_mat[keep_genes, , drop = FALSE]

cat("\nAfter expression filtering:\n")
print(dim(expr_filt))
cat("Genes retained:", nrow(expr_filt), "\n")

write.table(
  rownames(expr_filt),
  file = "../RNA_counts/eQTL_ready/n579_genes_retained_for_eQTL.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

#-----------------------------
# 2. Sample PCA on filtered expression
#-----------------------------
# center genes across samples
expr_scaled <- t(scale(t(expr_filt), center = TRUE, scale = TRUE))
expr_scaled[is.na(expr_scaled)] <- 0

pca <- prcomp(t(expr_scaled), center = FALSE, scale. = FALSE)

pca_df <- data.frame(
  IID = colnames(expr_filt),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  PC4 = pca$x[, 4],
  PC5 = pca$x[, 5]
)

pca_df <- merge(pca_df, as.data.frame(meta), by = "IID", all.x = TRUE)

# Scree / variance explained
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
pca_var <- data.frame(
  PC = seq_along(var_explained),
  VarianceExplained = var_explained,
  CumulativeVariance = cumsum(var_explained)
)

p1 <- ggplot(pca_var[1:20, ], aes(x = PC, y = VarianceExplained)) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw(base_size = 12) +
  labs(
    title = "Expression PCA scree plot",
    x = "PC",
    y = "Proportion variance explained"
  )
p1
ggsave("../results/eQTL_prep/expression_PCA_scree_plot.pdf", p1, width = 5.5, height = 4)

# PC1 vs PC2
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(TYPE))) +
  geom_point(size = 2, alpha = 0.85) +
  theme_bw(base_size = 12) +
  labs(
    title = "Expression PCA: PC1 vs PC2",
    color = "TYPE"
  )
p2
ggsave("../results/eQTL_prep/expression_PCA_PC1_PC2_by_TYPE.pdf", p2, width = 6, height = 5)

# By batch if available
if ("batch" %in% colnames(pca_df)) {
  p3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(batch))) +
    geom_point(size = 2, alpha = 0.85) +
    theme_bw(base_size = 12) +
    labs(
      title = "Expression PCA: PC1 vs PC2 by batch",
      color = "batch"
    )
  
  ggsave("../results/eQTL_prep/expression_PCA_PC1_PC2_by_batch.pdf", p3, width = 6, height = 5)
}

#-----------------------------
# 3. Correlate RNA PCs with covariates
#-----------------------------
# identify available PC columns in pca_df
pc_names <- grep("^PC[0-9]+$", colnames(pca_df), value = TRUE)

assoc_list <- lapply(pc_names, function(pc) {
  out <- list(PC = pc)
  
  if ("Age" %in% colnames(pca_df)) {
    fit_age <- lm(as.formula(paste(pc, "~ Age")), data = pca_df)
    out$Age_p <- summary(fit_age)$coefficients[2, 4]
  }
  
  if ("Sex" %in% colnames(pca_df)) {
    fit_sex <- lm(as.formula(paste(pc, "~ Sex")), data = pca_df)
    out$Sex_p <- summary(fit_sex)$coefficients[2, 4]
  }
  
  if ("Pheno" %in% colnames(pca_df)) {
    fit_pheno <- lm(as.formula(paste(pc, "~ Pheno")), data = pca_df)
    out$Pheno_p <- summary(fit_pheno)$coefficients[2, 4]
  }
  
  if ("batch" %in% colnames(pca_df)) {
    fit_batch <- lm(as.formula(paste(pc, "~ as.factor(batch)")), data = pca_df)
    out$batch_p <- anova(fit_batch)[1, "Pr(>F)"]
  }
  
  as.data.table(out)
})

pc_assoc <- rbindlist(assoc_list, fill = TRUE)
fwrite(pc_assoc, "../results/eQTL_prep/expression_PC_covariate_associations.tsv", sep = "\t")


#-----------------------------
# 4. Build hidden covariates from expression PCs
#-----------------------------
# Use the top expression PCs as hidden covariates.
# A common starting point for this sample size is 10-20 hidden factors.
n_hidden <- 10

hidden_covars <- data.table(
  FID = meta$FID,
  IID = meta$IID
)

for (i in seq_len(n_hidden)) {
  hidden_covars[[paste0("exprPC", i)]] <- pca$x[, i]
}

fwrite(
  hidden_covars,
  "../RNA_counts/eQTL_ready/n579_expression_hidden_covariates.txt",
  sep = "\t"
)

#-----------------------------
# 5. Inverse-normal transform each gene
#-----------------------------
invnorm <- function(x) {
  ok <- is.finite(x)
  r <- rank(x[ok], ties.method = "average")
  z <- qnorm((r - 0.5) / length(r))
  out <- rep(NA_real_, length(x))
  out[ok] <- z
  out
}

expr_invnorm <- t(apply(expr_filt, 1, invnorm))
colnames(expr_invnorm) <- colnames(expr_filt)
rownames(expr_invnorm) <- rownames(expr_filt)

cat("\nInverse-normal transformed matrix dimensions:\n")
print(dim(expr_invnorm))

#-----------------------------
# 6. Save eQTL-ready outputs
#-----------------------------
write.table(
  expr_filt,
  file = "../RNA_counts/eQTL_ready/n579_TPM_combat_adjusted_counts_log2_filtered.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  expr_invnorm,
  file = "../RNA_counts/eQTL_ready/n579_TPM_combat_adjusted_counts_log2_filtered_invnorm.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

fwrite(
  pca_var,
  "../results/eQTL_prep/expression_PCA_variance_explained.tsv",
  sep = "\t"
)

fwrite(
  pca_df,
  "../results/eQTL_prep/expression_PCA_sample_scores.tsv",
  sep = "\t"
)

cat("\nSaved files:\n")
cat(" - filtered expression matrix\n")
cat(" - inverse-normal transformed expression matrix\n")
cat(" - expression hidden covariates\n")
cat(" - expression PCA scores\n")
cat(" - expression PCA variance explained\n")
cat(" - PC/covariate association table\n")
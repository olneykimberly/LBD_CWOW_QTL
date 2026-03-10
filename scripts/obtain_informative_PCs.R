library(data.table)
library(ggplot2)
.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")

eigenval_file <- "../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC.eigenval"
eigenvec_file <- "../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC.eigenvec"
pheno_file    <- "../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes.txt"

eigval <- fread(eigenval_file, header = FALSE)
setnames(eigval, "V1", "eigenvalue")
eigval[, PC := paste0("PC", seq_len(nrow(eigval)))]
eigval[, PC_num := seq_len(nrow(eigval))]
eigval[, prop_var := eigenvalue / sum(eigenvalue)]
eigval[, cum_var := cumsum(prop_var)]

pcs <- fread(eigenvec_file, header = FALSE)
setnames(pcs, c("FID", "IID", paste0("PC", seq_len(ncol(pcs) - 2))))

pheno <- fread(pheno_file)

dat <- merge(
  pcs,
  pheno[, .(FID, IID, Age, Sex, Pheno, Braak.NFT, Thal.amyloid, Cing.LB)],
  by = c("FID", "IID"),
  all.x = TRUE
)

p1 <- ggplot(eigval, aes(x = PC_num, y = prop_var)) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw(base_size = 12) +
  labs(
    title = "Scree plot",
    x = "Principal component",
    y = "Proportion of variance explained"
  )

print(p1)

p2 <- ggplot(eigval, aes(x = PC_num, y = cum_var)) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw(base_size = 12) +
  labs(
    title = "Cumulative variance explained",
    x = "Principal component",
    y = "Cumulative proportion of variance explained"
  )

print(p2)

pc_names <- paste0("PC", seq_len(min(20, ncol(pcs) - 2)))

test_results <- rbindlist(lapply(pc_names, function(pc) {
  out <- list(PC = pc)

  fit_age   <- lm(as.formula(paste("Age ~", pc)), data = dat)
  fit_sex   <- lm(as.formula(paste("Sex ~", pc)), data = dat)
  fit_pheno <- lm(as.formula(paste("Pheno ~", pc)), data = dat)
  fit_braak <- lm(as.formula(paste("`Braak.NFT` ~", pc)), data = dat)
  fit_thal  <- lm(as.formula(paste("`Thal.amyloid` ~", pc)), data = dat)
  fit_cing  <- lm(as.formula(paste("`Cing.LB` ~", pc)), data = dat)

  out$Age_p   <- summary(fit_age)$coefficients[2, 4]
  out$Sex_p   <- summary(fit_sex)$coefficients[2, 4]
  out$Pheno_p <- summary(fit_pheno)$coefficients[2, 4]
  out$Braak_p <- summary(fit_braak)$coefficients[2, 4]
  out$Thal_p  <- summary(fit_thal)$coefficients[2, 4]
  out$CingLB_p <- summary(fit_cing)$coefficients[2, 4]

  as.data.table(out)
}))

print(eigval)
print(test_results)
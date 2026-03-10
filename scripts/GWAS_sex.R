.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggrastr)

# 1. Load and DEDUPLICATE
females <- fread("../snp_array/associations_with_imputed_snps/GWAS_Females_Only.assoc.logistic") %>%
  group_by(SNP) %>% filter(P == min(P, na.rm = TRUE)) %>% slice(1) %>% ungroup()
males <- fread("../snp_array/associations_with_imputed_snps/GWAS_Males_Only.assoc.logistic") %>%
  group_by(SNP) %>% filter(P == min(P, na.rm = TRUE)) %>% slice(1) %>% ungroup()
snp_anno <- fread("../snp_array/reference/NCBI_SNP_positions_with_gene_hg19ToHg38_lift_over.tsv")
snp_anno[, `:=`(gene = Gene, SNP = names)] 

# 2. Merge, Calculate, and THIN the data immediately
combined <- inner_join(females, males, by="SNP", suffix=c("_F", "_M")) %>%
  filter(CHR_F != 23, !is.na(P_F), !is.na(P_M)) %>%
  left_join(snp_anno[, .(SNP, gene)], by="SNP") %>%
  mutate(
    BETA_F = log(OR_F), BETA_M = log(OR_M),
    SE_F = abs(BETA_F) / abs(STAT_F), SE_M = abs(BETA_M) / abs(STAT_M),
    z_diff = (BETA_M - BETA_F) / sqrt(SE_F^2 + SE_M^2),
    abs_z = abs(z_diff)
  ) %>%
  # PERFORMANCE FIX: Only keep significant hits, top Z-scores, and 1% of the background noise
  filter(P_F < 0.001 | P_M < 0.001 | runif(n()) < 0.01)

# 3. Setup Labels
top_10 <- combined %>% filter(!is.na(gene) & gene != "") %>%
  group_by(gene) %>% summarise(max_z = max(abs_z)) %>%
  arrange(desc(max_z)) %>% slice(1:10)

combined_labeled <- combined %>%
  mutate(label = ifelse(gene %in% top_10$gene & abs_z == ave(abs_z, gene, FUN = max), gene, NA))

# 4. Plotting (Miami)
pdf("../results/gwas/Miami_Plot_Fast.pdf", width = 12, height = 6)
ggplot(combined_labeled, aes(x = BP_F, y = -log10(P_F))) +
  rasterise(geom_point(aes(color = as.factor(CHR_F)), size = 0.4), dpi = 300) +
  rasterise(geom_point(aes(y = log10(P_M), color = as.factor(CHR_F)), size = 0.4), dpi = 300) +
  geom_text_repel(aes(label = label), size = 3, force = 50, clip = "off") +
  scale_y_continuous(labels = abs) +
  facet_wrap(~CHR_F, scales = "free_x", nrow = 1) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
dev.off()

# 5. Plotting (Z-score)
pdf("../results/gwas/Zscore_Plot_Fast.pdf", width = 12, height = 6)
ggplot(combined_labeled, aes(x = BP_F, y = z_diff)) +
  rasterise(geom_point(aes(color = as.factor(CHR_F)), size = 0.4), dpi = 300) +
  geom_text_repel(aes(label = label), size = 3, force = 50, clip = "off") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~CHR_F, scales = "free_x", nrow = 1) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
dev.off()
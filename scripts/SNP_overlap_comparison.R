library(data.table)  # for fast data loading
library(ggvenn)      # for creating a Venn diagram
#install.packages("ggvenn")
getwd()
# Read the .bim files
bim1 <- fread("snp_array/Ross_TOPMED_imputed/CWOW_TOPMED_Rsq08_QC_maf01.bim", header = FALSE)
bim2 <- fread("snp_array/TOPMED_imput/final_filtered_data.bim", header = FALSE)

# Assign column names
colnames(bim1) <- c("CHR", "SNP", "cM", "BP", "A1", "A2")
colnames(bim2) <- c("CHR", "SNP", "cM", "BP", "A1", "A2")


# Merge by SNP ID, Position (BP), Allele 1 (A1), and Allele 2 (A2)
merged_bim <- merge(bim1, bim2, by = c("CHR", "BP", "A1", "A2"), all = TRUE, suffixes = c("_bim1", "_bim2"))

# Identify common and unique SNPs
common_snps <- merged_bim[complete.cases(merged_bim), ]
unique_to_bim1 <- merged_bim[is.na(merged_bim$SNP_bim2), ]
unique_to_bim2 <- merged_bim[is.na(merged_bim$SNP_bim1), ]


# Prepare the data for the Venn diagram
snp_sets <- list(
  BIM1 = merged_bim$SNP_bim1,
  BIM2 = merged_bim$SNP_bim2
)

# Create the Venn diagram
ggvenn(snp_sets, fill_color = c("#F8766D", "#00BFC4"), show_percentage = FALSE)

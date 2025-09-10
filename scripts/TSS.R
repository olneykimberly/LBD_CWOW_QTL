# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")

# --- 1. Load the data ---

# Load the cis eQTL results
# The file is a tab-separated value (.tsv), so we use read_tsv()
eqtl_results <- read_tsv("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/MatrixEQTL/cis_eQTL_LBD_ATS_control_imputed_modelLINEAR_CROSS")
#fread
# Load SNP coordinates
snp_coords <- read_tsv("../snp_array/reference/NCBI_SNP_positions_with_gene_hg19ToHg38_lift_over.tsv")

# Load gene coordinates
# The gene file is a tab-separated value (.txt) but read_tsv() is robust enough
gene_coords <- fread("../RNA_counts/genes_filtered.txt")

spec(snp_coords)
# --- 2. Merge all data frames and prepare for calculation ---

# Rename columns for merging
# Find the column named "names" and rename it to "SNP"
names(snp_coords)[names(snp_coords) == "names"] <- "SNP"

# Check the new column names to confirm the change
colnames(snp_coords)
names(gene_coords)[names(gene_coords) == "gene_id"] <- "gene"

# Join the tables
merged_data <- eqtl_results %>%
  left_join(snp_coords, by = "SNP") %>%
  left_join(gene_coords, by = "gene", suffix = c("_snp", "_gene"))

glimpse(gene_coords)
glimpse(eqtl_results)
# --- 3. Determine the TSS and calculate the signed distance ---

# Calculate the Transcription Start Site (TSS) position
# based on the gene's strand (+ or -)
merged_data <- merged_data %>%
  mutate(tss = ifelse(strand_gene == "+", start_gene, end_gene))

# Calculate the signed distance in kilobases (Kb)
merged_data <- merged_data %>%
  mutate(distance_kb = (start_snp - tss) / 1000)

# --- 4. Plot the histogram using ggplot2 ---
# `geom_histogram` creates the bars, and `geom_vline` adds the vertical line

merged_data_q01 <- merged_data %>%
  filter(FDR <= 0.01)

ggplot(merged_data_q01 , aes(x = distance_kb)) +
  geom_histogram(bins = 100, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = 0), color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Distribution of cis-eQTLs Relative to TSS",
    x = "Distance from TSS (Kb)",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 12)
  )


# --- 3. Determine the Transcription Start Site (TSS) for each gene ---

# The TSS is the 'start' position for genes on the positive (+) strand
# and the 'end' position for genes on the negative (-) strand.
gene_coords <- gene_coords %>%
  mutate(tss = ifelse(strand == "+", start, end))

# --- 4. Merge all data frames ---

# First, merge eQTL results with SNP coordinates based on the SNP column
merged_data <- left_join(eqtl_results, snp_coords, by = "SNP")

# Then, merge the combined data with gene coordinates based on the gene column
merged_data <- left_join(merged_data, gene_coords, by = "gene", suffix = c("_snp", "_gene"))

# --- 5. Calculate the distance and filter ---

# Calculate the distance between the SNP and the gene's TSS
# A '5kb' distance is 5000 base pairs.
filtered_results <- merged_data %>%
  mutate(distance_to_tss = abs(start_snp - tss)) %>%
  filter(distance_to_tss <= 5000)  %>%
  filter(FDR <= 0.05)

# --- 6. Output the final filtered results ---

# You can write the filtered data to a new file
write_tsv(filtered_results, "cis_eQTL_within_500b_of_TSS.tsv")

# To view the first few rows of the filtered results
head(filtered_results)
glimpse(snp_coords)


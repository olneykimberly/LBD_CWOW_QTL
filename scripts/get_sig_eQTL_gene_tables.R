# Setup
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/", "file_paths_and_colours.R"))

genes <- fread("../RNA_counts/genes_filtered.txt")
gene_info <- genes[,c(1:4, 12, 13)] # Keep only columns 1 - 4
colnames(gene_info) <- c("gene","chr", "start", "end", "gene_type", "gene_name") # Rename the columns 


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
  cis <- fread(paste0("../MatrixEQTL/cis_eQTL_", subset,"_imputed_modelLINEAR_CROSS"), header = TRUE,  fill = TRUE)
  sig_cis <- subset(cis, FDR < 0.001)
  sig_cis_gene_info <- merge(sig_cis, gene_info, by = "gene")
  # Save the results 
  output_table <- paste0("../MatrixEQTL/cis_eQTL_", subset,"_imputed_modelLINEAR_CROSS_FDR0.001.txt")
  write.table(sig_cis_gene_info, file = output_table, row.names = FALSE, quote = FALSE, sep = "\t")
  
  sig_cis_up <- subset(sig_cis_gene_info, `t-stat` > 5)
  sig_cis_up_unique_gene <- sig_cis_up[!duplicated(sig_cis_up$gene), ]
  output_table <- paste0("../MatrixEQTL/cis_eQTL_", subset,"_imputed_modelLINEAR_CROSS_FDR0.001_abststat5_positive_tstat.txt")
  write.table(sig_cis_up_unique_gene, file = output_table, row.names = FALSE, quote = FALSE, sep = "\t")
  
  sig_cis_down <- subset(sig_cis_gene_info, `t-stat` < -5)
  sig_cis_down_unique_gene <- sig_cis_down[!duplicated(sig_cis_down$gene), ]
  output_table <- paste0("../MatrixEQTL/cis_eQTL_", subset,"_imputed_modelLINEAR_CROSS_FDR0.001_abststat5_negative_tstat.txt")
  write.table(sig_cis_down_unique_gene, file = output_table, row.names = FALSE, quote = FALSE, sep = "\t")
}

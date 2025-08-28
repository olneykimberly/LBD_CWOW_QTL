## ----setup------------------------------------------------------------------------------------------------------------
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/", "file_paths_and_colours.R"))
library(dplyr)
library(data.table)
library(MatrixEQTL)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

## ----SNP_and_gene_data------------------------------------------------------------------------------------------------
# Read in gene information. These are the genes used in the bulk RNAseq differential expression analysis.
genes <- fread("../RNA_counts/genes_filtered.txt")
gene_info <- genes[,c(1:4)] # Keep only columns 1 - 4
colnames(gene_info) <- c("geneid","chr", "left", "right") # Rename the columns 

# Read in SNP annotation
# The SNP annotation was download from the dbSNP NCBI website 
# The vcf file was processed using a python script to obtain chromosome, position, SNP id, and gene name
snp_anno <- fread("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/snp_array/reference/filtered_snp_positions.tsv", 
                  header = TRUE,  fill = TRUE)
# Make an alleles column 
snp_anno$Alleles <- paste(snp_anno$REF, snp_anno$ALT, sep = "/")

# Convert snp annotation to GRanges object
gr <- GRanges(
  seqnames = paste0("chr", snp_anno$`#CHROM`),
  ranges = IRanges(start = snp_anno$POS, end = snp_anno$POS),
  names = snp_anno$ID,
  alleles = snp_anno$Alleles, 
  Gene = snp_anno$GENE
)

# Lift SNP annotation to GRCh38
chain <- import.chain("../snp_array/reference/hg19ToHg38.over.chain")
lifted <- liftOver(gr, chain)
lifted_gr <- unlist(lifted)

GRCh38_snp_annotation <- as.data.frame(lifted_gr)
write.table(GRCh38_snp_annotation,"../snp_array/reference/NCBI_SNP_positions_with_gene_hg19ToHg38_lift_over.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#------ End
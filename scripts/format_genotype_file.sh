#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=150G
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/format_genotype.job.%j
#SBATCH --job-name=format_genotype # Job name
#SBATCH --mail-type=END     # Mail events (NONE, BEGIN, END, FAIL, ALL)

# Go to working directory
cd /tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/
conda activate QTL

# Make logs directory if it does not already exist
mkdir -p logs

# Input/output files
INPUT_FILE="snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC.clean_genotype.raw"
STEP1_OUT="snp_array/final_gwas_dataset/output_file_additional_filtering.txt"
STEP2_OUT="snp_array/final_gwas_dataset/transposed_data_additional_filters.txt"
FINAL_OUT="snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC_n579.genotype_formatted.txt"

echo "Starting job at: $(date)"
echo "Working directory: $(pwd)"

# Replace spaces with tabs
awk '{gsub(/ /,"\t"); print}' "${INPUT_FILE}" > "${STEP1_OUT}"

# Transpose the file
datamash transpose < "${STEP1_OUT}" > "${STEP2_OUT}"

# Remove line 1 and lines 3-6 from the transposed file
sed -e '1d; 3,6d' "${STEP2_OUT}" > "${FINAL_OUT}"

echo "Finished job at: $(date)"
echo "Final output written to: ${FINAL_OUT}"

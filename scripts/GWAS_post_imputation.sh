#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=68G
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=//tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/GWAS.job.%j
#SBATCH --job-name=GWAS # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)

conda activate QTL

## GWAS
# Genome-wide association studies (GWAS) test genetic variants across the genome
# to identify those associated with quantitative or binary traits.
# Here we test linear associations with Braak NFT stage, Thal amyloid phase,
# cingulate Lewy body counts, and age.
# Models adjust for PC1-PC5, age, and sex as appropriate.
cd /tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts

mkdir -p ../snp_array/associations_with_imputed_snps
mkdir -p ../snp_array/plink_logs

FINAL_BFILE="../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC"
PHENO_COVAR="../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes_imputed_final.txt"
OUTDIR="../snp_array/associations_with_imputed_snps"

# Quantitative GWAS: Braak NFT stage
plink \
  --bfile ${FINAL_BFILE} \
  --linear \
  --pheno ${PHENO_COVAR} \
  --pheno-name Braak \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Sex \
  --hide-covar \
  --out ${OUTDIR}/Braak_association_imputed_covar_adj

# Quantitative GWAS: Thal amyloid phase
plink \
  --bfile ${FINAL_BFILE} \
  --linear \
  --pheno ${PHENO_COVAR} \
  --pheno-name Thal \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Sex \
  --hide-covar \
  --out ${OUTDIR}/Thal_association_imputed_covar_adj

# Quantitative GWAS: Cingulate Lewy body counts
plink \
  --bfile ${FINAL_BFILE} \
  --linear \
  --pheno ${PHENO_COVAR} \
  --pheno-name CingLB \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Sex \
  --hide-covar \
  --out ${OUTDIR}/CingLB_association_imputed_covar_adj

# Quantitative GWAS: Age
# Adjust for ancestry PCs, sex, and disease status
plink \
  --bfile ${FINAL_BFILE} \
  --linear \
  --pheno ${PHENO_COVAR} \
  --pheno-name Age \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Sex,Pheno \
  --hide-covar \
  --out ${OUTDIR}/Age_association_imputed_covar_adj


## Binary GWAS
# Disease status and sex are coded as:
# Sex: 1 = male, 2 = female
# Pheno: 1 = controls/PA, 2 = AD/LBD

# GWAS of disease status: controls/PA (1) vs AD/LBD (2)
plink \
  --bfile ${FINAL_BFILE} \
  --logistic \
  --pheno ${PHENO_COVAR} \
  --pheno-name Pheno \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Sex \
  --hide-covar \
  --out ${OUTDIR}/GWAS_DiseaseStatus_Control_vs_Case

# GWAS of sex
# Adjust for ancestry PCs, age, and disease status
plink \
  --bfile ${FINAL_BFILE} \
  --logistic \
  --pheno ${PHENO_COVAR} \
  --pheno-name Sex \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Pheno \
  --hide-covar \
  --out ${OUTDIR}/GWAS_Sex

# Female-only disease GWAS
plink \
  --bfile ${FINAL_BFILE} \
  --logistic \
  --pheno ${PHENO_COVAR} \
  --pheno-name Pheno \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age \
  --filter-females \
  --hide-covar \
  --out ${OUTDIR}/GWAS_Females_Only

# Male-only disease GWAS
plink \
  --bfile ${FINAL_BFILE} \
  --logistic \
  --pheno ${PHENO_COVAR} \
  --pheno-name Pheno \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age \
  --filter-males \
  --hide-covar \
  --out ${OUTDIR}/GWAS_Males_Only

  # GWAS of age-related SNPs

# ANALYSIS 1: Longevity / healthy aging in controls
plink \
  --bfile ${FINAL_BFILE} \
  --linear \
  --keep ../snp_array/covariates_and_phenotype_files/controls_and_PA.txt \
  --pheno ${PHENO_COVAR} \
  --pheno-name Age \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Sex \
  --hide-covar \
  --out ${OUTDIR}/Longevity_GWAS_Controls_Only

# ANALYSIS 2: Age-associated SNPs in AD/LBD cases
plink \
  --bfile ${FINAL_BFILE} \
  --linear \
  --keep ../snp_array/covariates_and_phenotype_files/cases_AD_LBD.txt \
  --pheno ${PHENO_COVAR} \
  --pheno-name Age \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Sex \
  --hide-covar \
  --out ${OUTDIR}/Age_of_Onset_GWAS_Cases_Only

# clean up logs
mv ${OUTDIR}/*.log ../snp_array/plink_logs 2>/dev/null
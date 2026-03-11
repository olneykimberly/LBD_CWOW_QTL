#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=68G
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/GWAS_control_vs_disease.job.%j
#SBATCH --job-name=GWAS # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)

## GWAS
# Genome-wide association studies (GWAS) test genetic variants across the genome
# to identify those associated with quantitative or binary traits.
# Here we test linear associations with Braak NFT stage, Thal amyloid phase,
# cingulate Lewy body counts, and age.
# Models adjust for PC1-PC5, age, and sex as appropriate.
cd /tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts
conda activate QTL

mkdir -p ../snp_array/associations_with_imputed_snps
mkdir -p ../snp_array/plink_logs

# controls+PA vs AD
# tests SNP associations with AD risk relative to neurologically normal/presymptomatic individuals

# controls+PA vs LBD
# tests SNP associations with LBD risk relative to neurologically normal/presymptomatic individuals

# AD vs LBD
# exploratory analysis testing whether SNPs differ between the two disease groups directly

# controls only vs AD/LBD = cleaner neuropathology-free reference
# controls + PA vs AD/LBD = broader nonclinical/preclinical reference

FINAL_BFILE="../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC"
PHENO_COVAR="../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes_imputed_final.txt"
OUTDIR="../snp_array/associations_with_imputed_snps"

# controls + PA vs AD
plink \
  --bfile ${FINAL_BFILE} \
  --logistic \
  --ci 0.95 \
  --pheno ../snp_array/covariates_and_phenotype_files/Pheno_controlsPA_vs_AD.txt \
  --pheno-name Pheno_ctrlPA_vs_AD \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Sex \
  --hide-covar \
  --out ${OUTDIR}/GWAS_ControlsPA_vs_AD

# controls + PA vs LBD
plink \
  --bfile ${FINAL_BFILE} \
  --logistic \
  --ci 0.95 \
  --pheno ../snp_array/covariates_and_phenotype_files/Pheno_controlsPA_vs_LBD.txt \
  --pheno-name Pheno_ctrlPA_vs_LBD \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Sex \
  --hide-covar \
  --out ${OUTDIR}/GWAS_ControlsPA_vs_LBD

# AD vs LBD (exploratory)
plink \
  --bfile ${FINAL_BFILE} \
  --logistic \
  --ci 0.95 \
  --pheno ../snp_array/covariates_and_phenotype_files/Pheno_AD_vs_LBD.txt \
  --pheno-name Pheno_AD_vs_LBD \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Sex \
  --hide-covar \
  --out ${OUTDIR}/GWAS_AD_vs_LBD


# controls only vs AD
plink \
  --bfile ${FINAL_BFILE} \
  --logistic \
  --ci 0.95 \
  --pheno ../snp_array/covariates_and_phenotype_files/Pheno_controlsOnly_vs_AD.txt \
  --pheno-name Pheno_ctrl_vs_AD \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Sex \
  --hide-covar \
  --out ${OUTDIR}/GWAS_ControlsOnly_vs_AD

# controls only vs LBD
plink \
  --bfile ${FINAL_BFILE} \
  --logistic \
  --ci 0.95 \
  --pheno ../snp_array/covariates_and_phenotype_files/Pheno_controlsOnly_vs_LBD.txt \
  --pheno-name Pheno_ctrl_vs_LBD \
  --covar ${PHENO_COVAR} \
  --covar-name PC1,PC2,PC3,PC4,PC5,Age,Sex \
  --hide-covar \
  --out ${OUTDIR}/GWAS_ControlsOnly_vs_LBD
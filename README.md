# LBD_CWOW_QTL
Expression quantitative trait loci (eQTL)

SNP array data available on Synapse: syn51238188

The n=632 samples part of the Lewy body dementia multi-institutional Center Without Walls (CWOW) initiative were genotyped on the Infinium OmniExpress 24 Illumina SNP Array with a genotyping rate of 99.82%.
Quality control was performed to remove duplicates, indels and multiallelic variants. Imputation was performed using TopMed Imputation server (TOPMed Imputation Server (nih.gov)) at hg38 with the imputation Rsq threshold set to 0.8.
CWOW_TOPMED_Rsq08_QC files were generated when HWE p < 1 x 10-5 was applied and contain 22,701,528 variants.
CWOW_TOPMED_Rsq08_QC_maf01 were generated when a MAF 1% filter was applied and contain 8,212,690 variants.

The original files (632_CWOW) used for imputation are also provided and they contain 714238 variants.

Identity by descent analysis (IBD) was performed on the original (un-imputed) files and we identified 3 pairs of individuals that are potentially related PI_HAT > 0.25. For any future analysis we recommend excluding an individual from each of the following pairs.

MC10078-MC05057
MC13264-MC17200
MC07261-MC09314

## Set up conda environment
This workflow uses conda. For information on how to install conda [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

To create the environment:
```
conda env create -n QTL --file QTL.yml

# To activate this environment, use
#
#     $ conda activate QTL
#
# To deactivate an active environment, use
#
#     $ conda deactivate QTL

```
### Step 1: Add metadata information to .fam
Add sex information to the .fam file, check for sex discrepancies with PLINK, and remove individuals from related pairs. 
```
R 01_check_meta.Rmd
plink --bfile 632_CWOW --check-sex --out sex_check
```

### Step 2: Remove excluded IIDs
Remove individuals that don't have corresponding RNA data, remove a individual from each sample pair. A list of samples to exclude was created when running the R script 01_check_meta.Rmd
```
plink --bfile yourdata_with_sex --remove exclude.txt --make-bed --out yourdata_filtered
```


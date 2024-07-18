# CWOW LBD QTL study
Expression quantitative trait loci (eQTL)

SNP array data available on Synapse: syn51238188\
Bulk RNAseq expressiond data available on Synapse: syn52394100

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
cd snp_array # move to the snp folder
plink --bfile 632_CWOW --remove ../metadata/exclude_fam_IID.txt --make-bed --out Filtered_n598_CWOW

plink --bfile CWOW_TOPMED_Rsq08_QC --remove ../metadata/exclude_fam_IID_for_imputed_files.txt --make-bed --out Filtered_n598_CWOW_TOPMED_Rsq08_QC

plink --bfile CWOW_TOPMED_Rsq08_QC_maf01 --remove ../metadata/exclude_fam_IID_for_imputed_files.txt --make-bed --out Filtered_n598_CWOW_TOPMED_Rsq08_QC_maf01

# inspect filtering
wc -l Filtered_n598_CWOW_TOPMED_Rsq08_QC.fam # there should be 598
```

### Step 3: Create genotype file
Convert PLINK files 
```
plink --bfile Filtered_n598_CWOW --recodeA --out Filtered_n598_CWOW_genotype    
```

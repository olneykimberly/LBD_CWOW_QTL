# CWOW LBD QTL study
Expression quantitative trait loci (eQTL)

SNP array data available on Synapse: syn51238188
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
plink --bfile yourdata_with_sex --remove exclude.txt --make-bed --out yourdata_filtered
```


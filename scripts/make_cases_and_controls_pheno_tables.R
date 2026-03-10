# make_cases_and_controls_pheno_tables.R
.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))
library(data.table)

df <- fread("../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes_imputed_final.txt")

# Controls and PA only (Pheno == 1)
controls <- df[Pheno == 1, .(FID, IID)]
fwrite(controls,
       "../snp_array/covariates_and_phenotype_files/controls_and_PA.txt",
       sep = "\t", col.names = FALSE)

# Cases: AD and LBD (Pheno == 2)
cases <- df[Pheno == 2, .(FID, IID)]
fwrite(cases,
       "../snp_array/covariates_and_phenotype_files/cases_AD_LBD.txt",
       sep = "\t", col.names = FALSE)
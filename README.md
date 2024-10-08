# CWOW LBD QTL study
Expression quantitative trait loci (eQTL)

SNP array data available on Synapse: syn51238188\
Bulk RNAseq expressiond data available on Synapse: syn52394100

## Set up conda environment
This workflow uses conda. For information on how to install conda may be found  [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

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

# Additionally, the liftover tool will need to be obtained
# See https://genome.ucsc.edu/cgi-bin/hgLiftOver for instructions. 
```
## Quality control genotype data
### Step 1: Match RNAseq samples to samples with SNP array data 
The output will be a filtered metadata that contains the information on the individuals that have both bulk RNAseq expression data and SNP array data. 
```
cd scripts # All commands below assume you are in the scripts folder. 
R 01_check_meta.Rmd
```

### Step 2: Remove excluded IIDs
Remove individuals that don't have corresponding RNA data and remove an individual from each related sample pair. A list of samples to exclude was created when running the R script 01_check_meta.Rmd
```
plink --bfile ../snp_array/632_CWOW \
      --remove ../metadata/exclude_fam_IID.txt \
      --make-bed --out ../snp_array/Filtered_n598_CWOW

# inspect filtering
wc -l  ../snp_array/Filtered_n598_CWOW.fam # there should be 598
```

### Step 3: Quality control 
First, run Principal component analysis (PCA) with plink to determine identify outliers. 
PCA with ancestry information will also be looked at, but for now, we will look at the PCA of only the CWOW samples. 
```
plink --bfile ../snp_array/Filtered_n598_CWOW \
      --pca 20 \
      --out ../snp_array/Filtered_n598_CWOW_pca_results

# Plot the PCA
R 02_PLINK_population_PCA.Rmd
```
Clean up the pca results after creating the PCA plot. These files are no longer needed. 
```
rm ../snp_array/Filtered_n598_CWOW_pca_results*
# move plink logs
mv ../snp_array/*.log plink_logs
```
The above PCA only contains individuals in the CWOW dataset. The identification of individuals of divergent ancestry can be achieved by combining the genotypes of the the CWOW population with genotypes of a reference dataset consisting of individuals from known ethnicities (for instance individuals from the Hapmap or 1000 genomes study). Below describes how to download the HapMap III data and merge with the CWOW data. 

#### Reference HAPMAP data 
Following the tutorial outlined in plinkQC HapMap [here](https://meyer-lab-cshl.github.io/plinkQC/articles/HapMap.html)
```
# First create a reference folder in which the hapmap data will be stored
cd ../snp_array/
mkdir reference 
cd reference 

# Download the hapmap genotype files from NCBI
wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2

# unzip 
bunzip2 hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
bunzip2 hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2

# Run plink to create a bed file 
plink --file hapmap3_r2_b36_fwd.consensus.qc.poly \
      --make-bed --out HapMapIII_NCBI36
```

Hapmap chromosome data is encoded numerically, with chrX represented by chr23, and chrY as chr24. In order to match to data encoded by chrX and chrY, we will have to rename these hapmap chromosomes. 
Be sure to be in the snp_array/reference folder 
```
# Assumes your in the snp_array/reference/ folder 
awk '{print "chr" $1, $4 -1, $4, $2 }' HapMapIII_NCBI36.bim | \
     sed 's/chr23/chrX/' | \
     sed 's/chr24/chrY/' > HapMapIII_NCBI36.tolift
     
awk '{print "chr" $1, $4 -1, $4, $2 }' ../Filtered_n598_CWOW.bim | \
     sed 's/chr23/chrX/' | \
     sed 's/chr24/chrY/' > ../Filtered_n598_CWOW.tolift
     
```

The genome build of HapMap III data is NCBI36. In order to update the HapMap III data to GRCh38, we use the UCSC liftOver tool. The liftOver tool takes information in a format similar to the PLINK .bim format, the UCSC bed format and a liftover chain, containing the mapping information between the old genome (target) and new genome (query). It returns the updated annotation and a file with unmappable variants. 

The liftOver tool will need to be first downloaded, which is freely available for academic use. It can not be obtained via conda. [Instructions](https://genome.ucsc.edu/cgi-bin/hgLiftOver) on how to download liftOver. 

Additionally the appropriate chain file will need to be download. The file names reflect the assembly conversion data contained within
in the format <db1>To<Db2>.over.chain.gz. For example, a file named hg19ToHg38.over.chain.gz file contains the liftOver data needed to
convert hg19 coordinates to hg38. 

```
# Download chain, be sure to be in the reference folder 
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# unzip 
gunzip hg18ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
 
# liftOver tool had to be downloaded first, see: https://genome-store.ucsc.edu/products/
# Lift over HapMap data
liftOver HapMapIII_NCBI36.tolift \
         hg18ToHg38.over.chain \
         HapMapIII_CGRCh38 HapMapIII_NCBI36.unMapped
# Lift over CWOW data 
liftOver ../Filtered_n598_CWOW.tolift \
        hg19ToHg38.over.chain \
        ../CWOW_n598_GRCh38 ../CWOW_n598_GRCh38.unMapped

# Ectract mapped variants
awk '{print $4}' HapMapIII_CGRCh38 > HapMapIII_CGRCh38.snps
# Ectract updated positions
awk '{print $4, $3}' HapMapIII_CGRCh38 > HapMapIII_CGRCh38.pos

# Update the hapmap reference 
# Extract the mappable variants from the old build and update their position
plink --bfile HapMapIII_NCBI36 \
      --extract HapMapIII_CGRCh38.snps \
      --update-map HapMapIII_CGRCh38.pos \
      --make-bed --out HapMapIII_CGRCh38_updated

## Repeat for CWOW data 
# Ectract mapped variants
awk '{print $4}' ../CWOW_n598_GRCh38 > ../CWOW_n598_GRCh38.snps
# Ectract updated positions
awk '{print $4, $3}' ../CWOW_n598_GRCh38 > ../CWOW_n598_GRCh38.pos

# Update the CWOW data 
# Extract the mappable variants from the old build and update their position
plink --bfile ../Filtered_n598_CWOW \
      --extract ../CWOW_n598_GRCh38.snps \
      --update-map ../CWOW_n598_GRCh38.pos \
      --make-bed --out ../CWOW_n598_GRCh38_updated
```

clean up by moving log files to the plink_logs folder and removing intermediate files 
```
# remove intermediate HapMap files
rm HapMapIII_NCBI36.hh HapMapIII_NCBI36.bed HapMapIII_NCBI36.fam HapMapIII_NCBI36.bim HapMapIII_NCBI36.tolift HapMapIII_CGRCh38 HapMapIII_NCBI36.unMapped HapMapIII_CGRCh38.snps HapMapIII_CGRCh38.pos
mv *.log ../plink_logs

# remove intermediate CWOW files
rm ../Filtered_n598_CWOW* ../CWOW_n598_CGRCh38 ../CWOW_n598_CGRCh38.unMapped ../CWOW_n598_CGRCh38.snps ../CWOW_n598_CGRCh38.pos
mv ../*.log ../plink_logs
```

### strand and allele flip correction
Strand flips and allele flips are common issues encountered in SNP array data when aligning genotypes between different datasets or comparing to a reference panel. These issues arise because SNPs can be represented in different ways depending on the strand or the order of the alleles.

A strand flip occurs when the SNP is read from the opposite DNA strand. For example, an A/T SNP on one strand will appear as a T/A SNP on the complementary strand. If one dataset refers to the forward strand and another to the reverse strand, the SNP may appear as if it has different alleles, leading to discrepancies unless corrected.
An allele flip refers to cases where the alleles of a SNP are swapped (e.g., A/G versus G/A). Although technically the same SNP, the order of alleles might differ between datasets or when compared to a reference. Allele flips can create inconsistency in analysis, especially when comparing allele frequencies or performing meta-analyses.

We will use PLINK’s --flip option to correct strand flips. The --flip command can be used after generating a strand report with --flip-scan, which detects potential strand flips based on allele frequencies.
We can use PLINK’s --a1-allele option to ensure allele coding is consistent with a reference.
Finally, we will compare allele frequencies between our CWOW dataset and a reference panel using PLINK’s --freq command. Large discrepancies might indicate allele mismatches or strand issues.
We may also want to remove SNPs that don’t match the reference dataset using --extract or --exclude options in PLINK.
```
cd ../ # move out of the reference folder and into the snp_array folder 

# Strand flip scan 
plink --bfile CWOW_n598_CGRCh38_updated --flip-scan --out flip_report

# If a variant ID appears multiple times with different information (e.g., different allele flip types) the variant ID can only appear once 
awk '!seen[$1]++' flip_report.flipscan > flip_report_unique.flipscan

# Flip correction 
plink --bfile CWOW_n598_CGRCh38_updated \
      --flip flip_report_unique.flipscan \
      --make-bed --out CWOW_n598_CGRCh38_strand_corrected

# Reference alleles
awk '{print $2, $5}' reference/HapMapIII_CGRCh38_updated.bim > reference/HapMapIII_reference_alleles.txt

# Keep only variants that are in the CWOW data
awk 'NR==FNR{a[$2]; next} $1 in a' CWOW_n598_CGRCh38_strand_corrected.bim \
    reference/HapMapIII_reference_alleles.txt > reference/HapMapIII_filtered_reference_alleles.txt

# Correct allele flips by comparing CWOW data to HapMap data 
plink --bfile CWOW_n598_CGRCh38_strand_corrected \
      --bmerge reference/HapMapIII_CGRCh38_updated.bed \
      reference/HapMapIII_CGRCh38_updated.bim \
      reference/HapMapIII_CGRCh38_updated.fam \
      --merge-mode 6 --out flip_check
# The command above will generate a file flip_check.missnp listing SNPs that have potential strand flips or mismatches.

# Fix strand flips using the list of mismatched SNPs flip_check.missnp
plink --bfile CWOW_n598_CGRCh38_strand_corrected \
      --flip flip_check.missnp \
      --make-bed --out CWOW_flipped

# Merge CWOW data with HapMap to ensure allele flipping worked 
plink --bfile CWOW_flipped \
      --bmerge reference/HapMapIII_CGRCh38_updated.bed \
      reference/HapMapIII_CGRCh38_updated.bim \
      reference/HapMapIII_CGRCh38_updated.fam \
      --make-bed --out merged_data
```
Check output: Next we will ensure that the strand alignment of our CWOW target dataset matches that of the reference dataset. We can compare allele frequencies and manually check or use software like Genotype Harmonizer.
We will also ensure that SNP identifiers match between our CWOW dataset and the reference. Tools like liftOver or dbSNP can help verify the correct SNP coordinates and annotations.

Clean up 
```
mv *.log plink_logs
mv flip_report* plink_logs
mv flip_check* plink_logs
rm CWOW_n598_CGRCh38_updated*
rm CWOW_n598_CGRCh38_strand_corrected.*
rm merged_data*
```
After the above steps have been completed the HapMap III & CWOW data sets can be used for inferring study ancestry as described below. 

#### Match CWOW genotypes with hapmap reference data
In order to compute joint principal components of the reference hapmap and the CWOW study population, we’ll need to combine the two datasets. The plink –merge function enables this merge, but requires the variants in the datasets to be matching by chromosome, position and alleles. The following sections show how to extract the relevant data from the reference and study data set and how to filter matching variants, as described [here](https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html) in the plinkQC tutorial. 

Filter reference and study data for non A-T or G-C SNPs as these SNPs are more difficult to align and only a subset of SNPs is required for the analysis, we will remove them from both the reference and study data set.
```
# Filter hapmap data and create bed file 
awk 'BEGIN {OFS="\t"}  \
    ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") \
    {print $2}' reference/HapMapIII_CGRCh38_updated.bim > reference/HapMapIII_CGRCh38.ac_gt_snps

plink --bfile reference/HapMapIII_CGRCh38_updated \
      --exclude reference/HapMapIII_CGRCh38.ac_gt_snps \
      --make-bed --out reference/HapMapIII_CGRCh38.ac_gt_snps.no_ac_gt_snps

# Filter CWOW data and create bed file 
awk 'BEGIN {OFS="\t"} \
    ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") \
    {print $2}' CWOW_flipped.bim  > CWOW_flipped.ac_gt_snps

plink --bfile CWOW_flipped \
      --exclude CWOW_flipped.ac_gt_snps \
      --make-bed --out CWOW_flipped.no_ac_gt_snps
```

#### Prune
Conduct principle component analysis on genetic variants that are pruned for variants in linkage disequilibrium (LD) with an 𝑟2>0.2 in a 50kb window. The LD-pruned data set is generated below, using plink–indep-pairwise to compute the LD-variants. Additionally exclude range is used to remove genomic ranges of known high-LD structure. This file is available in file.path(find.package('plinkQC'),'extdata','high-LD-regions.txt').
```
# Create prune in file 
plink --bfile CWOW_flipped.no_ac_gt_snps \
      --exclude range reference/high-LD-regions-hg38-GRCh38.txt \
      --indep-pairwise 50 5 0.2 \
      --out CWOW_flipped.no_ac_gt_snps.no_ac_gt_snps

# Create bed file of the newly pruned data
plink --bfile CWOW_flipped.no_ac_gt_snps \
      --extract CWOW_flipped.no_ac_gt_snps.no_ac_gt_snps.prune.in \
      --make-bed \
      --out CWOW_flipped.no_ac_gt_snps.pruned

# Filter hapmap reference data for the same SNP set as in the CWOW study
plink --bfile reference/HapMapIII_CGRCh38_updated \
      --extract CWOW_flipped.no_ac_gt_snps.no_ac_gt_snps.prune.in \
      --make-bed \
      --out reference/HapMapIII_CGRCh38_updated.pruned
```

Check and correct chromosome mismatch. The following section uses an awk to check that the variant IDs of the reference data have the same chromosome ID as the study data.  Merging the files via PLINK will only work for variants with perfectly matching attributes. For simplicity, and not crucial to the final task of inferring ancestory, we will ignore XY-encoded sex chromosomes (via sed -n '/^[XY]/!p').
```
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    CWOW_flipped.no_ac_gt_snps.pruned.bim reference/HapMapIII_CGRCh38_updated.pruned.bim | \
    sed -n '/^[XY]/!p' > reference/HapMapIII_CGRCh38_updated.toUpdateChr

plink --bfile reference/HapMapIII_CGRCh38_updated.pruned \
      --update-chr reference/HapMapIII_CGRCh38_updated.toUpdateChr 1 2 \
      --make-bed \
      --out reference/HapMapIII_CGRCh38_updated.updateChr
```

Position mismatch - Similar to the chromosome matching, we use awk to find variants with mis-matching chromosomal positions
```
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    CWOW_flipped.no_ac_gt_snps.pruned.bim reference/HapMapIII_CGRCh38_updated.pruned.bim > \
    reference/HapMapIII_CGRCh38_updated.toUpdatePos
```

Unlike chromosomal and base-pair annotation, mismatching allele-annotations will not only prevent the plink –merge, but also mean that it is likely that actually a different genotype was measured. Initially, we can use the following awk to check if non-matching allele codes are a simple case of allele flips.
```
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    CWOW_flipped.no_ac_gt_snps.pruned.bim reference/HapMapIII_CGRCh38_updated.pruned.bim > \
    reference/HapMapIII_CGRCh38_updated.toFlip
```

We use plink to update the mismatching positions and possible allele-flips identified above.
Any alleles that do not match after allele flipping, are identified and removed from the reference dataset.
```
plink --bfile reference/HapMapIII_CGRCh38_updated.updateChr \
      --update-map reference/HapMapIII_CGRCh38_updated.toUpdatePos 1 2 \
      --flip reference/HapMapIII_CGRCh38_updated.toFlip \
      --make-bed \
      --out reference/HapMapIII_CGRCh38_updated.flipped

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
     CWOW_flipped.no_ac_gt_snps.pruned.bim reference/HapMapIII_CGRCh38_updated.flipped.bim > \
     reference/HapMapIII_CGRCh38_updated.mismatch

plink --bfile reference/HapMapIII_CGRCh38_updated.flipped \
      --exclude reference/HapMapIII_CGRCh38_updated.mismatch \
      --make-bed \
      --out reference/HapMapIII_CGRCh38_updated.clean
```
### Merge
Merge study genotypes and reference data
The matching study and reference dataset can now be merged into a combined dataset with plink –bmerge. If all steps outlined above were conducted successfully, no mismatch errors should occur.
```
# Merge 
plink --bfile CWOW_flipped.no_ac_gt_snps.pruned \
      --bmerge reference/HapMapIII_CGRCh38_updated.clean.bed \
               reference/HapMapIII_CGRCh38_updated.clean.bim \
               reference/HapMapIII_CGRCh38_updated.clean.fam \
      --make-bed \
      --out merge_HAP_CWOW

# PCA with merged hapmap and CWOW data 
plink --bfile merge_HAP_CWOW \
      --pca \
      --out merge_HAP_CWOW_pca
```
We can now use the merge_HAP_CWOW_pca.eigenvec file to estimate the ancestry of the CWOW study samples. Identifying individuals of divergent ancestry is implemented in check_ancestry in the R script described below. 
Before proceeding, clean up the snp_array folder by moving all of the log files to the plink_logs folder
```
mv *log plink_logs
mv reference/*log plink_logs
rm reference/*ac_gt_snps* reference/HapMapIII_CGRCh38_updated.pruned* reference/*updateChr* reference/*flipped* reference/*.mismatch reference/*toUpdateChr reference/*toUpdatePos reference/*toFlip
rm *.no_ac_gt_snps*
```
#### Per individual plink QC checks 
The R script 02_PLINK_QC.Rmd will run the PlinkQC protocol which will implement  three main functions:
1) The per-individual quality control (perIndividualQC)
2) The per-marker quality control (perMarkerQC)
3) The generation of the new, quality control dataset (cleanData)

The script will output a clean dataset after removing outlier samples and markers.
```
# move to the scripts folder
cd ../scripts
R 02_PLINK_QC.Rmd
```
An overview of the results may be viewed [here](https://rpubs.com/olneykimberly/PlinkQC_LBD_CWOW_SNP_array)
The output is cleaned data that removed individuals and markers that didn't pass the quality control checks as described in the 02_PLINK_QC.Rmd 

There are now 580 individuals that passed QC and 615,293 SNPs. 
Total number of samples remaining post QC: 580
Total SNPs remaining post QC: 615,293 before QC was 714,238
```
# clean up 
cd ../snp_array/
mv *no_failIDs* post_plinkQC/
mv *fail* post_plinkQC/
mv *.log plink_logs/
mv *.remove* post_plinkQC/
mv *.prune* post_plinkQC/
rm CWOW_flipped.sexcheck CWOW_flipped.lmiss CWOW_flipped.imiss CWOW_flipped.het
```

### Step 4: Create PCA with the filtered CWOW data
Create PCA and update metadata file to reflect that there are now 580 individuals 
```
cd ../scripts
# PCA 
plink --bfile ../snp_array/CWOW_flipped.clean \
      --pca 20 --out ../snp_array/CWOW_flipped.clean_pca_results

# update metadata
R 03_update_meta.Rmd
```
### Step 5: Create genotype file
```
cd ../snp_array
plink --bfile CWOW_flipped.clean \
      --recodeA --out CWOW_flipped.clean_genotype    
```
### Step 6: Update metadata file, counts data, and genotype file 
```
R 04_process_genotype.Rmd
```
### Step 7: Run Matrix eQTL 
Firstly, the files will need some additional formatting. Then we can run the eQTl analysis, with or without including an interaction term with disease. 
```
# Format the inputs 
R 05a_format_inputs_for_MatrixEQTL.Rmd

# Run matrix eQTL, note that it takes several hours to run. Submit as a background job. 
R 05b_run_MatrixEQTL.Rmd

# There is also a loop script that will loop through each disease type. 
05_run_matrixeqlt_loop.R
```

### Step 8: Plot the eQTL results 
```
# eQTL plots for top SNP to gene associations
R 06_plot_eQTLs.Rmd
```

## GWAS 
Genome-wide association studies (GWAS) test thousands of genetic variants across many genomes to find those statistically associated with a specific trait or disease. Here we will determine linear assocaitions between variants and disease traits of Braak NFT stage, Thal amyloid phase, and the counts of Lewy bodies in the cingulate cortex. 
```
cd ../ # In the main project folder
mkdir GWAS
cd GWAS

plink --bfile ../snp_array/CWOW_flipped.clean \
      --linear \
      --out cingLBD_association \
      --pheno ../snp_array/covariates_and_phenotype_files/CingLB_phenotypes.txt
      
plink --bfile ../snp_array/CWOW_flipped.clean \
      --linear \
      --out Braak_association \
      --pheno ../snp_array/covariates_and_phenotype_files/Braak_phenotypes.txt
      
plink --bfile ../snp_array/CWOW_flipped.clean \
      --linear \
      --out Thal_association \
      --pheno ../snp_array/covariates_and_phenotype_files/Thal_phenotypes.txt
      
# clean up 
mv *.log ../snp_array/plink_logs
```

Make Manhattan plots from GWAS association tests
```
# Manhattan plots from GWAS 
R 07_GWAS.Rmd
```

## Format files for imputing 
The above was completed for the clean unimputed genotypes. Now that we have clean genotypes from unrelated individuals, we will impute genotypes following the TOPMed impute protocol. https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/

If you use the Imputation TOPMed Server, cite:
Das S, Forer L, Schönherr S, Sidore C, Locke AE, Kwong A, Vrieze S, Chew EY, Levy S, McGue M, Schlessinger D, Stambolian D, Loh PR, Iacono WG, Swaroop A, Scott LJ, Cucca F, Kronenberg F, Boehnke M, Abecasis GR, Fuchsberger C. Next-generation genotype imputation service and methods. Nature Genetics 48, 1284–1287 (2016).

### Register and install packages
In the snp_array/reference folder. See https://www.chg.ox.ac.uk/~wrayner/tools/ for details. 
```
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
unzip HRC-1000G-check-bim-v4.3.0.zip
```

### Obtain the 1000 Genomes Phase 3 reference panel and lift over to GRCh38
```
wget https://www.chg.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz
gzip 1000GP_Phase3_combined.legend.gz
```
The 1000G reference panel is for hg19 and will need to be lifted over to GRCh38
```
# Format legened into a .bed
awk 'NR > 1 {print "chr"$2"\t"$3-1"\t"$3"\t"$1"\t0\t+"}' \
            1000GP_Phase3_combined.legend > 1000GP_Phase3_combined.bed
            
liftOver 1000GP_Phase3_combined.bed \ # input
         hg19ToHg38.over.chain \ # chain 
         1000GP_Phase3_combined_lifted_to_GRCh38.bed \ # output lifted to GRCh38
         1000GP_Phase3_combined_unlifted.bed # ouptut of unlifted variants

# Opitional, for checking the file 
# Format lifted .bed into legend format
awk '{print $4"\t"$1"\t"$3"\t"substr($4, index($4, ":")+1)}' \
        1000GP_Phase3_combined_lifted_to_GRCh38.bed > 1000GP_Phase3_combined_lifted_to_GRCh38.legend
```

Then in order to keep the 1000GP legend information, such as population ancestry frequency information, the lifted .bed will need to be merged with the unlifted 1000GP_Phase3_combined.legend.
```
# Step 1 - merge
python 01_merge_legend_with_lifted_bed.py
# inputs 1000GP_Phase3_combined.legend and 1000GP_Phase3_combined_lifted_to_GRCh38.bed
# Output 1000GP_Phase3_combined_lifted_to_GRCh38_merged_output.txt

# Step 2 - format
python 02_format_merged_legend.py 
# input 1000GP_Phase3_combined_lifted_to_GRCh38_merged_output.txt 
# output GRCh38_1000GP_Phase3_combined.legend
```

### Check CWOW SNPs against the 1000G reference panel
Execute HRC-1000G-check-bim.pl script which will generate the Run-plink.sh file 
The script will check plink .bim files against 1000G for strand, id names, positions, alleles, ref/alt assignment. 
The output will be Run-plink.sh. 
```
# Create a frequency file 
plink --freq \
      --bfile ../../CWOW_n598_GRCh38_updated.clean.bim \
      --out ../../CWOW.clean.frequency.frq

perl HRC-1000G-check-bim.pl \
     -b ../../CWOW_n598_GRCh38_updated.clean.bim \
     -f ../../CWOW.clean.frequency.frq \
     -r GRCh38_1000GP_Phase3_combined.legend -g -p EUR
# 1000G population will default to ALL if not specified. Most samples in CWOW are self reported white. 
# Writing plink commands to: Run-plink.sh
```
### Run-plink.sh to correct strand and allele flips against the 1000GP
Strand flips and allele flips are common issues encountered in SNP array data when aligning genotypes between different datasets or comparing to a reference panel. These issues arise because SNPs can be represented in different ways depending on the strand or the order of the alleles.

A strand flip occurs when the SNP is read from the opposite DNA strand. For example, an A/T SNP on one strand will appear as a T/A SNP on the complementary strand. If one dataset refers to the forward strand and another to the reverse strand, the SNP may appear as if it has different alleles, leading to discrepancies unless corrected. An allele flip refers to cases where the alleles of a SNP are swapped (e.g., A/G versus G/A). Although technically the same SNP, the order of alleles might differ between datasets or when compared to a reference. Allele flips can create inconsistency in analysis, especially when comparing allele frequencies or performing meta-analyses.

We will use PLINK’s --flip option to correct strand flips. The --flip command can be used after generating a strand report with --flip-scan, which detects potential strand flips based on allele frequencies. We can use PLINK’s --a1-allele option to ensure allele coding is consistent with a reference. Finally, we will compare allele frequencies between our CWOW dataset and a reference panel using PLINK’s --freq command. Large discrepancies might indicate allele mismatches or strand issues. We may also want to remove SNPs that don’t match the reference dataset using --extract or --exclude options in PLINK.
```
# Exclude individuals
plink --bfile CWOW_n598_GRCh38_updated.clean \
      --exclude Exclude-CWOW_n598_GRCh38_updated.clean-1000G.txt \
      --make-bed --out TEMP1 --allow-extra-chr

# Update chromosomes
plink --bfile CWOW_n598_GRCh38_updated.clean \
      --update-map Chromosome-CWOW_n598_GRCh38_updated.clean-1000G.txt \
      --update-chr --make-bed --out TEMP2 --allow-extra-chr

# Update positions       
plink --bfile TEMP2 \
      --update-map Position-CWOW_n598_GRCh38_updated.clean-1000G.txt 
      --make-bed --out TEMP3 --allow-extra-chr
      
# Strand flip
plink --bfile TEMP2 \
      --flip Strand-Flip-CWOW_n598_GRCh38_updated.clean-1000G.txt \
      --make-bed --out TEMP --allow-extra-chr

# Allele flip      
plink --bfile TEMP \
      --a2-allele Force-Allele1-CWOW_n598_GRCh38_updated.clean-1000G.txt \
      --make-bed --out CWOW_n598_GRCh38_updated.clean-updated --allow-extra-chr
```

### Prepare CWOW plink files for impuation 
TOPMed Imputation Server accepts VCF files compressed with bgzip. CWOW variants are in plink format and will need to be converted to VCF.
```
# Create VCF
plink --bfile CWOW_n598_GRCh38_updated.clean-updated  \
      --keep-allele-order \
      --recode vcf \
      --out CWOW
# clean up
mv *.log plink_logs/

# Sort VCF by genomic position
bcftools sort CWOW_no_allele_flip.vcf\
        -o CWOW_sorted.vcf

# zip 
bgzip CWOW_sorted.vcf

# index 
bcftools index CWOW_sorted.vcf.gz

# obtain chromosome names 
bcftools index -s CWOW_sorted.vcf.gz\
         | cut -f1 | sort | uniq > chrom_names.txt

# rename chromosome IDs to include "chr"
bcftools annotate \
         --rename-chrs rename_map.txt\
         -o CWOW_chr_renamed_file.vcf.gz -O z CWOW_sorted.vcf.gz
# index
bcftools index CWOW_chr_renamed_file.vcf.gz

# filter to exclude duplications, indels, and retain only biallelic SNPs
bcftools norm -d all -o CWOW_no_duplicates.vcf CWOW_chr_renamed_file.vcf.gz
bcftools view -v snps -o CWOW_no_duplicates_no_indels.vcf CWOW_no_duplicates.vcf
bcftools view -m2 -M2 -o CWOW_final.vcf CWOW_no_duplicates_no_indels.vcf

# compress & index
bgzip CWOW_final.vcf
bcftools index CWOW_final.vcf.gz

# Create separate VCFs for each chromosome
for chr in {1..22} X Y; do
    bcftools view -r chr$chr CWOW_chr_renamed_file.vcf.gz -o chr$chr.vcf
done

# Compress and index the VCF files
for chr in {1..22} X Y; do
    bgzip chr$chr.vcf
    tabix -p vcf chr$chr.vcf.gz
done
```

## Imputation via TOPMed
Go to https://imputation.biodatacatalyst.nhlbi.nih.gov/#! and create an account. Login and then select Run Genotype Imputation (Minimac4) 2.0.0-beta3. \
Input SNPs: 614955\
Number of samples: 580\
Chromosomes: 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 X\
Name: CWOW \
reference panel: TopMed r3\
Array build: GRCh38\
rsq Filter:0.2\
Phasing:eagle\
Population:all\
Mode:impuation\
Submit job


### Post imputation 
Index 
```
for chr in {1..22} X; do
    tabix -p vcf chr$chr.dose.vcf.gz
done
```

merge
```
bcftools concat -a -O z -o CWOW_TOPMED_imputed.vcf.gz \
chr1.dose.vcf.gz \
chr2.dose.vcf.gz \
chr3.dose.vcf.gz \
chr4.dose.vcf.gz \
chr5.dose.vcf.gz \
chr6.dose.vcf.gz \
chr7.dose.vcf.gz \
chr8.dose.vcf.gz \
chr9.dose.vcf.gz \
chr10.dose.vcf.gz \
chr11.dose.vcf.gz \
chr12.dose.vcf.gz \
chr13.dose.vcf.gz \
chr14.dose.vcf.gz \
chr15.dose.vcf.gz \
chr16.dose.vcf.gz \
chr17.dose.vcf.gz \
chr18.dose.vcf.gz \
chr19.dose.vcf.gz \
chr20.dose.vcf.gz \
chr21.dose.vcf.gz \
chr22.dose.vcf.gz \
chrX.dose.vcf.gz
```

remove of duplicates, indels and multiallelic variants 
```
plink --bfile CWOW_TOPMED_imputed --list-duplicate-vars --out duplicates
plink --bfile CWOW_TOPMED_imputed --exclude duplicates.dupvar --make-bed --out data_no_dups
plink --bfile data_no_dups --snps-only just-acgt --make-bed --out data_no_indels
plink --bfile data_no_indels --biallelic-only strict --make-bed --out clean_data
plink --bfile data_no_indels --biallelic-only strict --make-bed --out clean_data
```

HWE  HWE p < 1 x 10-5
```
plink --bfile clean_data --hardy --out hwe_results
# To filter for controls (recommended for case/control studies), add the midp option:
# plink --bfile clean_data --hwe 1e-5 midp --make-bed --out hwe_filtered_data

plink --bfile clean_data --hwe 1e-5 --make-bed --out hwe_filtered_data
```

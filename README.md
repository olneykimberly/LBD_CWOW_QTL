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

### Step 1: Match RNAseq samples to samples with SNP array data 
The output will be a filtered metadata that contains the information on the individuals that have both bulk RNAseq expression data and SNP array data. 
```
cd scripts # All commands below assume you are in the scripts folder. 
R 01_check_meta.Rmd
```

### Step 2: Remove excluded IIDs
Remove individuals that don't have corresponding RNA data and remove an individual from each related sample pair. A list of samples to exclude was created when running the R script 01_check_meta.Rmd
```
plink --bfile ../snp_array/632_CWOW --remove ../metadata/exclude_fam_IID.txt --make-bed --out ../snp_array/Filtered_n598_CWOW

# inspect filtering
wc -l  ../snp_array/Filtered_n598_CWOW.fam # there should be 598
```

### Step 3: Quality control 
First, run Principal component analysis (PCA) with plink to determine population ancestry and identify outliers. 
```
plink --bfile ../snp_array/Filtered_n598_CWOW --pca 20 --out ../snp_array/Filtered_n598_CWOW_pca_results

# Plot the PCA. 
R 02_PLINK_population_PCA.Rmd
```
The above PCA only contains individuals in the CWOW dataset. The identification of individuals of divergent ancestry can be achieved by combining the genotypes of the the CWOW population with genotypes of a reference dataset consisting of individuals from known ethnicities (for instance individuals from the Hapmap or 1000 genomes study). Below describes how to download the HapMap III data and merge with the CWOW data. 

#### Reference HAPMAP data 
Following the tutorial outlined here: https://meyer-lab-cshl.github.io/plinkQC/articles/HapMap.html
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
plink --file hapmap3_r2_b36_fwd.consensus.qc.poly --make-bed --out HapMapIII_NCBI36
```

Hapmap chromosome data is encoded numerically, with chrX represented by chr23, and chrY as chr24. In order to match to data encoded by chrX and chrY, we will have to rename these hapmap chromosomes. 
Be sure to be in the snp_array/reference folder 
```
# Assumes your in the snp_array/reference folder 
awk '{print "chr" $1, $4 -1, $4, $2 }' HapMapIII_NCBI36.bim | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > HapMapIII_NCBI36.tolift
```

The genome build of HapMap III data is NCBI36. In order to update the HapMap III data to GRCh38, we use the UCSC liftOver tool. The liftOver tool takes information in a format similar to the PLINK .bim format, the UCSC bed format and a liftover chain, containing the mapping information between the old genome (target) and new genome (query). It returns the updated annotation and a file with unmappable variants. 

The liftOver tool will need to be first downloaded, which is freely available for academic use. It can not be obtained via conda. See https://genome.ucsc.edu/cgi-bin/hgLiftOver for instructions on how to download liftOver. 

Additionally the appropriate chain file will need to be download. The file names reflect the assembly conversion data contained within
in the format <db1>To<Db2>.over.chain.gz. For example, a file named hg18ToHg38.over.chain.gz file contains the liftOver data needed to
convert hg18 coordinates to hg38. 

```
# Download chain    
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz
# unzip 
gunzip hg18ToHg38.over.chain.gz
 
# liftOver tool had to be downloaded first, see: https://genome-store.ucsc.edu/products/
liftOver HapMapIII_NCBI36.tolift hg18ToHg38.over.chain HapMapIII_CGRCh38 HapMapIII_NCBI36.unMapped

# ectract mapped variants
awk '{print $4}' HapMapIII_CGRCh38 > HapMapIII_CGRCh38.snps
# ectract updated positions
awk '{print $4, $3}' HapMapIII_CGRCh38 > HapMapIII_CGRCh38.pos

# update the hapmap reference by extracting the mappable variants from the old build and update their position. 
plink --bfile HapMapIII_NCBI36 --extract HapMapIII_CGRCh38.snps --update-map HapMapIII_CGRCh38.pos --make-bed --out HapMapIII_CGRCh38logs
```
After the above steps, the HapMap III dataset can be used for inferring study ancestry as described below. 

#### Match CWOW genotypes with hapmap reference data
In order to compute joint principal components of the reference hapmap and the CWOW study population, we’ll need to combine the two datasets. The plink –merge function enables this merge, but requires the variants in the datasets to be matching by chromosome, position and alleles. The following sections show how to extract the relevant data from the reference and study dataset and how to filter matching variants, as described here: https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html

Filter reference and study data for non A-T or G-C SNPs. These SNPs are more difficult to align and only a subset of SNPs is required for the analysis, we will remove them from both the reference and study data set.
```
cd ../ # Move out of the reference folder and into the snp_array folder 
# Filter CWOW data and create bed file 
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA")  {print $2}' Filtered_n598_CWOW.bim  > Filtered_n598_CWOW.ac_gt_snps
plink --bfile  Filtered_n598_CWOW --exclude Filtered_n598_CWOW.ac_gt_snps --make-bed --out Filtered_n598_CWOW.no_ac_gt_snps
      
# Filter hapmap data and create bed file 
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA")  {print $2}' reference/HapMapIII_CGRCh38.bim > reference/HapMapIII_CGRCh38.ac_gt_snps
plink --bfile  reference/HapMapIII_CGRCh38 --exclude reference/HapMapIII_CGRCh38.ac_gt_snps --make-bed --out reference/HapMapIII_CGRCh38.ac_gt_snps.no_ac_gt_snps
```

#### Prune
Conduct principle component analysis on genetic variants that are pruned for variants in linkage disequilibrium (LD) with an 𝑟2>0.2 in a 50kb window. The LD-pruned dataset is generated below, using plink–indep-pairwise to compute the LD-variants. additionally exclude range is used to remove genomic ranges of known high-LD structure. This file was originally provided by [6] and is available in file.path(find.package('plinkQC'),'extdata','high-LD-regions.txt').
```
# Create prune in file 
plink --bfile  Filtered_n598_CWOW.no_ac_gt_snps \
      --exclude range  reference/high-LD-regions-hg38-GRCh38.txt \
      --indep-pairwise 50 5 0.2 \
      --out Filtered_n598_CWOW.no_ac_gt_snps.no_ac_gt_snps

# Create bed file of the newly pruned data
plink --bfile  Filtered_n598_CWOW.no_ac_gt_snps \
      --extract Filtered_n598_CWOW.no_ac_gt_snps.no_ac_gt_snps.prune.in \
      --make-bed \
      --out Filtered_n598_CWOW.no_ac_gt_snps.pruned

# Filter hapmap reference data for the same SNP set as in the CWOW study
plink --bfile  reference/HapMapIII_CGRCh38 \
      --extract Filtered_n598_CWOW.no_ac_gt_snps.no_ac_gt_snps.prune.in \
      --make-bed \
      --out reference/HapMapIII_CGRCh38.pruned
```

Check and correct chromosome mismatch. The following section uses an awk to check that the variant IDs of the reference data have the same chromosome ID as the study data.  Merging the files via PLINK will only work for variants with perfectly matching attributes. For simplicity, and not crucial to the final task of inferring ancestory, we will ignore XY-encoded sex chromosomes (via sed -n '/^[XY]/!p').
```
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    Filtered_n598_CWOW.no_ac_gt_snps.pruned.bim reference/HapMapIII_CGRCh38.pruned.bim | \
    sed -n '/^[XY]/!p' > reference/HapMapIII_CGRCh38.toUpdateChr

plink --bfile reference/HapMapIII_CGRCh38.pruned \
      --update-chr reference/HapMapIII_CGRCh38.toUpdateChr 1 2 \
      --make-bed \
      --out reference/HapMapIII_CGRCh38.updateChr
```

Position mismatch - Similar to the chromosome matching, we use an awk-script to find variants with mis-matching chromosomal positions
```
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    Filtered_n598_CWOW.no_ac_gt_snps.pruned.bim reference/HapMapIII_CGRCh38.pruned.bim > \
    reference/HapMapIII_CGRCh38.toUpdatePos
```

Unlike chromosomal and base-pair annotation, mismatching allele-annotations will not only prevent the plink –merge, but also mean that it is likely that actually a different genotype was measured. Initially, we can use the following awk-script to check if non-matching allele codes are a simple case of allele flips.
```
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    Filtered_n598_CWOW.no_ac_gt_snps.pruned.bim reference/HapMapIII_CGRCh38.pruned.bim > \
    reference/HapMapIII_CGRCh38.toFlip
```
We use plink to update the mismatching positions and possible allele-flips identified above.
Any alleles that do not match after allele flipping, are identified and removed from the reference dataset.
```
plink --bfile reference/HapMapIII_CGRCh38.updateChr \
      --update-map reference/HapMapIII_CGRCh38.toUpdatePos 1 2 \
      --flip reference/HapMapIII_CGRCh38.toFlip \
      --make-bed \
      --out reference/HapMapIII_CGRCh38.flipped

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
     Filtered_n598_CWOW.no_ac_gt_snps.pruned.bim reference/HapMapIII_CGRCh38.flipped.bim > \
     reference/HapMapIII_CGRCh38.mismatch

plink --bfile reference/HapMapIII_CGRCh38.flipped \
      --exclude reference/HapMapIII_CGRCh38.mismatch \
      --make-bed \
      --out reference/HapMapIII_CGRCh38.clean
```

Merge study genotypes and reference data
The matching study and reference dataset can now be merged into a combined dataset with plink –bmerge. If all steps outlined above were conducted successfully, no mismatch errors should occur.
```
# Merge 
plink --bfile Filtered_n598_CWOW.no_ac_gt_snps.pruned  \
      --bmerge reference/HapMapIII_CGRCh38.clean.bed reference/HapMapIII_CGRCh38.clean.bim reference/HapMapIII_CGRCh38.clean.fam  \
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
mkdir plink_logs
mv *log plink_logs
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
An overview of the results may be viewed here: https://rpubs.com/olneykimberly/PlinkQC_LBD_CWOW_SNP_array 
The output is cleaned data that removed individuals and markers that didn't pass the quality control checks as described in the 02_PLINK_QC.Rmd 

There are now 579 individuals that passed QC and 616,183 SNPs. Lets rename the files to reflect that there are now only 579 individuals. 
```
cd ../snp_array
mv Filtered_n598_CWOW.clean.hh  Filtered_n579_CWOW.clean.hh
mv Filtered_n598_CWOW.clean.bed Filtered_n579_CWOW.clean.bed
mv Filtered_n598_CWOW.clean.fam Filtered_n579_CWOW.clean.fam
mv Filtered_n598_CWOW.clean.bim Filtered_n579_CWOW.clean.bim
```

### Step 4: Create PCA with the filtered CWOW data
Create PCA and update metadata file to reflect that there are now 579 individuals 
```
# PCA 
plink --bfile ../snp_array/Filtered_n579_CWOW.clean --pca 20 --out ../snp_array/Filtered_n579_CWOW_pca_results

# update metadata
R 03_update_meta.Rmd
```
### Step 5: Create genotype file
```
plink --bfile Filtered_n579_CWOW.clean --recodeA --out Filtered_n579_CWOW.clean_genotype    
```
### Step 6: Update metadata file, counts data, and genotype file 
```
R 04_process_genotype.Rmd
```



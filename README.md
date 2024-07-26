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



# Reference HAPMAP data 
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
plink --bfile HapMapIII_NCBI36 --extract HapMapIII_CGRCh38.snps --update-map HapMapIII_CGRCh38.pos --make-bed --out HapMapIII_CGRCh38
```

After the above steps, the HapMap III dataset can be used for inferring study ancestry as described below. 

Filter reference and study data for non A-T or G-C SNPs
```
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    Filtered_n598_CWOW.bim  > \
    Filtered_n598_CWOW.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    reference/HapMapIII_CGRCh37.bim > \
    HapMapIII_CGRCh37.ac_gt_snps
   
plink --bfile  reference/HapMapIII_CGRCh37 \
      --exclude HapMapIII_CGRCh37.bim.ac_gt_snps \
      --make-bed \
      --out HapMapIII_CGRCh37.ac_gt_snps.no_ac_gt_snps

plink --bfile  Filtered_n598_CWOW \
      --exclude Filtered_n598_CWOW.ac_gt_snps \
      --make-bed \
      --out Filtered_n598_CWOW.no_ac_gt_snps
```

# Prune
```
plink --bfile  Filtered_n598_CWOW.no_ac_gt_snps \
      --exclude range  /research/labs/neurology/fryer/m239830/home/R/x86_64-pc-linux-gnu-library/4.3/plinkQC/extdata/high-LD-regions-hg38-GRCh38.txt \
      --indep-pairwise 50 5 0.2 \
      --out Filtered_n598_CWOW.no_ac_gt_snps.no_ac_gt_snps

plink --bfile  Filtered_n598_CWOW.no_ac_gt_snps \
      --extract Filtered_n598_CWOW.no_ac_gt_snps.no_ac_gt_snps.prune.in \
      --make-bed \
      --out Filtered_n598_CWOW.no_ac_gt_snps.pruned
mv  $qcdir/$name.pruned.log $qcdir/plink_log/$name.pruned.log
```

```
plink --bfile  reference/HapMapIII_CGRCh37 \
      --extract Filtered_n598_CWOW.no_ac_gt_snps.no_ac_gt_snps.prune.in \
      --make-bed \
      --out HapMapIII_CGRCh37.pruned
```

```
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    Filtered_n598_CWOW.no_ac_gt_snps.pruned.bim HapMapIII_CGRCh37.pruned.bim | \
    sed -n '/^[XY]/!p' > HapMapIII_CGRCh37.toUpdateChr

plink --bfile HapMapIII_CGRCh37.pruned \
      --update-chr HapMapIII_CGRCh37.toUpdateChr 1 2 \
      --make-bed \
      --out HapMapIII_CGRCh37.updateChr
mv $qcdir/$refname.updateChr.log $qcdir/plink_log/$refname.updateChr.log
```

```
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    Filtered_n598_CWOW.no_ac_gt_snps.pruned.bim HapMapIII_CGRCh37.pruned.bim > \
    HapMapIII_CGRCh37.toUpdatePos
    
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    Filtered_n598_CWOW.no_ac_gt_snps.pruned.bim HapMapIII_CGRCh37.pruned.bim > \
    HapMapIII_CGRCh37.toFlip
    
plink --bfile HapMapIII_CGRCh37.updateChr \
      --update-map HapMapIII_CGRCh37.toUpdatePos 1 2 \
      --flip HapMapIII_CGRCh37.toFlip \
      --make-bed \
      --out HapMapIII_CGRCh37.flipped
mv $qcdir/$refname.flipped.log $qcdir/plink_log/$refname.flipped.log
```

```
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
     Filtered_n598_CWOW.no_ac_gt_snps.pruned.bim HapMapIII_CGRCh37.flipped.bim > \
    HapMapIII_CGRCh37.mismatch

plink --bfile HapMapIII_CGRCh37.flipped \
      --exclude HapMapIII_CGRCh37.mismatch \
      --make-bed \
      --out HapMapIII_CGRCh37.clean
mv $qcdir/$refname.clean.log $qcdir/plink_log/$refname.clean.log
```
merge
```
plink --bfile Filtered_n598_CWOW.no_ac_gt_snps.pruned  \
      --bmerge HapMapIII_CGRCh37.clean.bed HapMapIII_CGRCh37.clean.bim \
         HapMapIII_CGRCh37.clean.fam  \
      --make-bed \
      --out merge_HAP_CWOW
mv $qcdir/$name.merge.$refname.log $qcdir/plink_log

plink --bfile merge_HAP_CWOW \
      --pca \
      --out merge_HAP_CWOW_pca
mv $qcdir/$name.$reference.log $qcdir/plink_log
```


The R script 02_PLINK_QC.Rmd will run the PlinkQC protocol which will implement  three main functions:
1) The per-individual quality control (perIndividualQC)
2) The per-marker quality control (perMarkerQC)
3) The generation of the new, quality control dataset (cleanData)

The script will output a clean dataset after removing outlier samples and markers.
```
R 02_PLINK_population_PCA.Rmd
```
An overview of the results may be viewed here: https://rpubs.com/olneykimberly/PlinkQC_LBD_CWOW_SNP_array 

### Step 5: Create genotype file
Convert PLINK files 
```
plink --bfile Filtered_n598_CWOW --recodeA --out Filtered_n598_CWOW_genotype    
```


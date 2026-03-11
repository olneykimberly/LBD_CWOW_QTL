# CWOW LBD QTL study
SNP array data available on Synapse: syn51238188\
Bulk RNAseq expressiond data available on Synapse: syn52394100\

Abstract:\
Lewy body disease (LBD) and Alzheimer's disease (AD) are two of the most common neurodegenerative disorders. LBD is primarily marked by the accumulation of α-synuclein aggregates (Lewy bodies), often accompanied by amyloid-β (Aβ) plaques and neurofibrillary tangles (NFTs), which are the hallmark of AD. The APOE ε4 allele is associated with an increased risk of developing AD and LBD. Additionally, variants in SNCA and LRRK2 are linked to a higher risk of developing Lewy body pathology in AD patients. Building on these established associations, we tested the hypothesis that genomic variants significantly contribute to gene dysregulation across neuropathologies. We employed an integrative genomic approach, combining bulk RNA sequencing and SNP array data from n = 580 neuropathologically defined cases, to better understand the molecular mechanisms underlying these neuropathologies. Through eQTL (expression quantitative trait loci) analysis, we identified LBD-associated cis-eQTLs that were involved in pathways related to thiamine metabolism, which is critical for maintaining neuronal energy balance. Additionally, top LBD-associated cis-eQTL were enriched in pathways related to programmed cell death, suggesting a role for apoptosis and cellular regulation in the pathogenesis of LBD. Another significant finding was the identification of APOE as the top hit in a genome-wide association study (GWAS) analyzing the Thal amyloid phase in the cingulate cortex, further supporting its known role in amyloid deposition.
Conclusion: In conclusion, our eQTL analysis underscores the importance of regulatory variants in shaping the transcriptional landscape of these diseases, thereby providing a deeper understanding of the genetic mechanisms underlying neurodegeneration.

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
# output ../metadata/filtered_metadata.txt
```

### Step 2: Remove excluded IIDs
Remove individuals that don't have corresponding RNA data and remove an individual from each related sample pair. A list of samples to exclude was created when running the R script 01_check_meta.Rmd
```
plink --bfile ../snp_array/632_CWOW \
      --remove ../metadata/exclude_fam_IID.txt \
      --make-bed --out ../snp_array/Filtered_n598_CWOW

# inspect filtering
wc -l  ../snp_array/Filtered_n598_CWOW.fam # there should be 598 

# Filtered n598 files were moved to preprocessing_snp_array_before_imputation
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
# output is results/pca/CWOW_pca_dim1&2.pdf and dim2&3.pds
```

Clean up the pca files after creating the PCA plot. These files are no longer needed. 
```
rm ../snp_array/Filtered_n598_CWOW_pca_results*
# move plink logs
mv ../snp_array/*.log plink_logs
```
The above PCA only contains individuals in the CWOW dataset. The identification of individuals of divergent ancestry can be achieved by combining the genotypes of the the CWOW population with genotypes of a reference dataset consisting of individuals from known ethnicity (for instance individuals from the Hapmap or 1000 genomes study). Below describes how to download the HapMap III data and merge with the CWOW data. 

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
      
# Note - CWOW_n598_GRCh38_updated files were moved to preprocessing_snp_array_before_imputation/ folder 
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
plink --bfile ../snp_array/post_plinkQC/CWOW_flipped.clean \
      --pca 20 --out ../snp_array/CWOW_flipped.clean_pca_results

# update metadata
R 03_update_meta.Rmd
```
### Step 5: Create genotype file
```
cd ../snp_array
plink --bfile preprocessing_snp_array_before_imputation/CWOW_flipped.clean_pca_results \
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

# Format files for imputing 
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

The output is imputed variants per chromosome as vcf files.

### Processing TOPMed Imputed Genotypes and Filtering High-Confidence Variants
Index vcf files
```
for chr in {1..22} X; do
    tabix -p vcf chr$chr.dose.vcf.gz
done
```

Imputed genotype data from the TOPMed imputation server were provided as per-chromosome dosage VCF files (`chr*.dose.vcf.gz`) together with imputation quality summary files (`chr*.info.gz`). These files are then converted into PLINK2 format, merged into a genome-wide dataset, and filtered to retain only high-confidence imputed variants using the Minimac imputation accuracy metric (`R2`).
```
# Create sex file for chromosome X import
# PLINK2 requires sample sex information when importing chromosome X 
# This is used to correctly handle genotype ploidy.
# where `1 = male` and `2 = female`.

awk 'NR==1{next} {print $1, $2, $8}' \
../../covariates_and_phenotype_files/covariates_and_phenotypes.txt \
> ../../covariates_and_phenotype_files/sex_file.txt
```

### Convert imputed VCF files to PLINK2 format
Each imputed chromosome VCF was converted into PLINK2 format while preserving imputed genotype dosage values.
```
for chr in {1..22}; do
  plink2 \
    --vcf chr${chr}.dose.vcf.gz dosage=DS \
    --make-pgen \
    --out chr${chr}_TOPMED
done

plink2 \
  --vcf chrX.dose.vcf.gz dosage=DS \
  --update-sex ../../covariates_and_phenotype_files/sex_file.txt \
  --make-pgen \
  --out chrX_TOPMED
```

### Merge chromosomes into a genome-wide dataset
All chromosome datasets were merged into a single genome-wide PLINK2 dataset. Variant IDs were reassigned using the format `CHR:POS:REF:ALT` to ensure uniqueness across split multiallelic variants present in the imputed data.
```
for chr in {1..22} X; do
  echo chr${chr}_TOPMED
done > pmerge_list.txt

plink2 \
  --pmerge-list pmerge_list.txt \
  --set-all-var-ids @:#:\$r:\$a \
  --new-id-max-allele-len 96 \
  --make-pgen \
  --out CWOW_TOPMED_allchr
```

### Identify high-confidence imputed variants
Imputation introduces uncertainty in genotype estimation. The Minimac imputation accuracy metric (`R2`) estimates the correlation between imputed and true genotypes. To retain high-confidence variants for downstream GWAS and eQTL analyses, variants were filtered to `R2 ≥ 0.8`.
Because variant IDs were reassigned during merging (`CHR:POS:REF:ALT`), the INFO-derived SNP list must use the same format.
```
# make directory to store the r2 variant list per chromosome
mkdir -p info_r2_0.8_lists
# for each chr get the variants with r2 >= 0.8
for chr in {1..22} X; do
  zcat chr${chr}.info.gz | \
  awk -F'\t' '
  BEGIN{OFS=""}
  !/^#/ {
    match($8, /R2=([0-9.]+)/, a)
    if (a[1] >= 0.8) print $1, ":", $2, ":", $4, ":", $5
  }' > info_r2_0.8_lists/chr${chr}.R2_0.8.CHRPOSREFALT.snplist
done

# The number of high-confidence variants retained per chromosome can be inspected with:
for chr in {1..22} X; do
  wc -l info_r2_0.8_lists/chr${chr}.R2_0.8.CHRPOSREFALT.snplist
done
```

### Combine high-confidence variants across chromosomes
This produces a genome-wide list of high-confidence imputed variants.
```
cat info_r2_0.8_lists/chr*.R2_0.8.CHRPOSREFALT.snplist \
> info_r2_0.8_lists/allchr.R2_0.8.CHRPOSREFALT.snplist

sort -u info_r2_0.8_lists/allchr.R2_0.8.CHRPOSREFALT.snplist \
> info_r2_0.8_lists/allchr.R2_0.8.unique.CHRPOSREFALT.snplist

wc -l info_r2_0.8_lists/allchr.R2_0.8.unique.CHRPOSREFALT.snplist
# 20,484,716

# remove chr as PLINK won't have that in the name 
sed 's/^chr//' info_r2_0.8_lists/allchr.R2_0.8.unique.CHRPOSREFALT.snplist \
> info_r2_0.8_lists/allchr.R2_0.8.unique.CHRPOSREFALT.nochr.snplist
```

### Apply the INFO filter to the merged dataset
The resulting filtered dataset contains genome-wide imputed genotypes restricted to variants with high imputation accuracy (`R2 ≥ 0.8`), suitable for downstream quality control and association analyses.
```
plink2 \
  --pfile CWOW_TOPMED_allchr \
  --extract info_r2_0.8_lists/allchr.R2_0.8.unique.CHRPOSREFALT.snplist \
  --make-pgen \
  --out CWOW_TOPMED_allchr_R2_0.8
```

### Before moving on to filtering, make a mapping key of the rsIDs and positions
Because the variants are now named as CHR:POS:REF:ALT we will have a key that matches this to the rsIDs
```
for chr in {1..22} X; do
  awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1 ":" $2 ":" $4 ":" $5, $3}' chr${chr}_TOPMED.pvar
done > CHRPOSREFALT_to_rsid.txt
```

### Filter to retain only biallelic variants
-max-alleles 2 keeps only biallelic variants
--snps-only just-acgt keeps only standard A/C/G/T SNPs, removing indels and other non-SNP alleles
writes a new filtered PLINK2 dataset
```
 plink2 \
  --pfile CWOW_TOPMED_allchr_R2_0.8 \
  --max-alleles 2 \
  --snps-only just-acgt \
  --make-pgen \
  --out CWOW_TOPMED_allchr_R2_0.8_biallelic_snps
```
19,182,512 variants remaining after filtering 

### Additional filtering and Create genotype file
| Filter                     | Threshold                             |
| -------------------------- | ------------------------------------- |
| Missingness                | `geno ≤ 0.05`                         |
| Hardy–Weinberg equilibrium | `HWE ≥ 1e-6` (controls)               |
| Minor allele count         | `MAC ≥ 20` (better for small samples) |
```
# create a final clean genotype file for all downstream data analysis
plink2 \
  --pfile CWOW_TOPMED_allchr_R2_0.8_biallelic_snps \
  --geno 0.05 \
  --mac 20 \
  --make-pgen \
  --out CWOW_TOPMED_allchr_R2_0.8_biallelic_snps_geno0.05_mac20
```
6,793,783 variants remaining after mac and geno filters.

### HWE filtering (controls)
Hardy–Weinberg Equilibrium (HWE) is a population genetics principle describing the expected relationship between allele frequencies and genotype frequencies in a randomly mating population.\
In case–control GWAS, true disease-associated variants can deviate from HWE in cases. For example, if allele A increases disease risk, then cases will contain excess AA or Aa genotypes relative to the general population. That deviation from equilibrium is real biology, not a genotyping error. Thus only apply HWE to the controls. 
```
# Obtain list of control samples which is the controls and PA samples
awk 'NR>1 && $9==1 {print $1 "_" $2}' \
../../covariates_and_phenotype_files/covariates_and_phenotypes.txt \
> controls.keep

plink2 \
  --pfile CWOW_TOPMED_allchr_R2_0.8_biallelic_snps_geno0.05_mac20 \
  --keep controls.keep \
  --hwe 1e-6 midp \
  --write-snplist \
  --out controls_hwe_pass
# 6,793,778 variants remaining after main filters. (only 5 variants)

plink2 \
  --pfile CWOW_TOPMED_allchr_R2_0.8_biallelic_snps_geno0.05_mac20 \
  --extract controls_hwe_pass.snplist \
  --make-pgen \
  --out CWOW_TOPMED_allchr_R2_0.8_biallelic_snps_geno0.05_mac20_hwe1e6
```

### Re-run PLINKQC on the final imputed dataset
We will need to merge with HapMap to obtain ancestry information. Some variants are not consistnt between HapMap and our CWOW dataset, so they need to be either flipped or removed. 
```
cd TOPMED_imput/cwow_merge_with_hapmap/
# Flipping the problematic SNPs
plink --bfile TOPMED_imput/gwas_filtered_data \
      --flip TOPMED_imputmerge_HAP_CWOW_filtered-merge.missnp \
      --make-bed \
      --out TOPMED_imput/gwas_filtered_data_flipped

# Try the merge again
plink --bfile TOPMED_imput/gwas_filtered_data_flipped \
      --bmerge reference/HapMapIII_CGRCh38_updated.clean.bed \
               reference/HapMapIII_CGRCh38_updated.clean.bim \
               reference/HapMapIII_CGRCh38_updated.clean.fam \
      --make-bed \
      --out TOPMED_imputmerge_HAP_CWOW_filtered_try2

# Remove problematic SNPs from CWOW dataset
plink --bfile TOPMED_imput/gwas_filtered_data_flipped \
      --exclude TOPMED_imputmerge_HAP_CWOW_filtered_try2-merge.missnp \
      --make-bed \
      --out TOPMED_imput/gwas_filtered_data_flipped_clean

# Remove same SNPs from HapMap reference
plink --bfile reference/HapMapIII_CGRCh38_updated.clean \
      --exclude TOPMED_imputmerge_HAP_CWOW_filtered_try2-merge.missnp \
      --make-bed \
      --out reference/HapMapIII_CGRCh38_updated.clean_no_missnp

# Merge cleaned datasets
plink --bfile TOPMED_imput/gwas_filtered_data_flipped_clean \
      --bmerge reference/HapMapIII_CGRCh38_updated.clean_no_missnp.bed \
               reference/HapMapIII_CGRCh38_updated.clean_no_missnp.bim \
               reference/HapMapIII_CGRCh38_updated.clean_no_missnp.fam \
      --make-bed \
      --out TOPMED_imputmerge_HAP_CWOW_filtered_final
```

### Add phenotype information to the imputed PLINK files
Fix IDs so PLINK uses real FID/IID instead of MC00092_MC00092
```
# First obtain the list of FID and IIDs
# This creates an ID update file from our CWOW covariate table:
awk 'BEGIN{OFS="\t"} NR>1 {print 0, $1 "_" $2, $1, $2}' \
../../covariates_and_phenotype_files/covariates_and_phenotypes.txt \
> update_ids.txt

# Now apply to our PLINK file
plink2 \
  --pfile CWOW_TOPMED_allchr_R2_0.8_biallelic_snps_geno0.05_mac20_hwe1e6 \
  --update-ids update_ids.txt \
  --make-pgen \
  --out CWOW_qc_ids

# check 
head CWOW_qc_ids.psam
```

### Create sex and pheno type files for downstream use
```
# sex only file
awk 'BEGIN{OFS="\t"; print "#FID","IID","SEX"} NR>1 {print $1,$2,$8}' \
../../covariates_and_phenotype_files/covariates_and_phenotypes.txt \
> ../../covariates_and_phenotype_files/sex_file_plink.txt

# pheno only file
awk 'BEGIN{OFS="\t"; print "#FID","IID","Pheno"} NR>1 {print $1,$2,$9}' \
../../covariates_and_phenotype_files/covariates_and_phenotypes.txt \
> ../../covariates_and_phenotype_files/pheno_file_plink.txt

# covariates file for later GWAS/PCA plotting
awk 'BEGIN{OFS="\t"; print "#FID","IID","Age","Sex","Braak.NFT","Thal.amyloid","Cing.LB"} NR>1 {print $1,$2,$13,$8,$10,$11,$12}' \
../../covariates_and_phenotype_files/covariates_and_phenotypes.txt \
> ../../covariates_and_phenotype_files/covar_file_plink.txt
```

### Add the sex and pheno files, and split PAR on chrX
Our CWOW data is for hg38/TOPMed. Thus we will use --split-par b38 before sex-checking chrX. PLINK 2 explicitly recommends splitting PAR to avoid chrX handling problems
```
plink2 \
  --pfile CWOW_qc_ids \
  --update-sex sex_file_plink.txt \
  --pheno pheno_file_plink.txt \
  --split-par b38 \
  --make-pgen \
  --out CWOW_qc_ids_sex_pheno_splitPAR

# check
head CWOW_qc_ids_sex_pheno_splitPAR.psam
```
At this point, sex should no longer be NA, and the phenotype should be loaded for downstream PLINK runs. PLINK 2 loads phenotypes from --pheno, and --make-pgen writes a new fileset after applying these operations.

### Make PLINK1 files for PLINKQC 
plinkQC is built around PLINK QC outputs and is most natural with a PLINK1 fileset.
```
# convert PLINK2 to PLINK1 outputs
plink2 \
  --pfile CWOW_qc_ids_sex_pheno_splitPAR \
  --make-bed \
  --out CWOW_qc_ids_sex_pheno_splitPAR

# check
head CWOW_qc_ids_sex_pheno_splitPAR.fam
```
Should now have proper FID IID populated SEX and phenotype in column 6

### PLINKQC 
indir = /snp_array/TOPMED_imput/imputed_vcf\
name  = CWOW_qc_ids_sex_pheno_splitPAR\
This input will be used for check_sex(), check_het_and_miss(), and check_relatedness()\ 
A sample failed heterozygosiry check and we will remove that sample before merging with HapMap for ancestry checking. 
```
# create file of the failed sample
echo -e "MC12181\tMC12181" > remove_het_fail.txt

# create a new PLINK file that has that sample removed
plink \
  --bfile CWOW_qc_ids_sex_pheno_splitPAR \
  --remove remove_het_fail.txt \
  --make-bed \
  --out CWOW_qc_ids_sex_pheno_splitPAR_noHetFail
```
There are now 579 individuals retained.

### Merge with HapMap
To check for population ancestry, we will now use this clean dataset and restrict to the autosomes only. 
```
# get autosomes only
plink \
  --bfile CWOW_qc_ids_sex_pheno_splitPAR_noHetFail \
  --autosome \
  --make-bed \
  --out CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome
```
rewrite the CWOW BIM IDs from CHR:POS:REF:ALT back to rsIDs where available, then merge with HapMap.
```
# Start from the autosome-only CWOW file after removing the heterozygosity outlier
# Current file:
#   CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome

# Make a unique CHR:POS:REF:ALT -> rsID mapping.
# Keep only rows with a real rsID, and only keep rsIDs that appear once.
awk 'BEGIN{FS=OFS="\t"}
     $2 != "." {count[$2]++; key[NR]=$1; val[NR]=$2}
     END{
       for(i=1;i<=NR;i++){
         if(val[i] != "." && count[val[i]]==1) print key[i], val[i]
       }
     }' CHRPOSREFALT_to_rsid.txt > CHRPOSREFALT_to_rsid.unique.txt

# Rewrite the BIM file so column 2 uses rsIDs where a unique mapping exists
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {map[$1]=$2; next}
     {if($2 in map) $2=map[$2]; print}' \
CHRPOSREFALT_to_rsid.unique.txt \
CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome.bim \
> CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome_rsid.bim

# Copy over the BED and FAM so this becomes a valid PLINK dataset
cp CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome.bed CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome_rsid.bed
cp CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome.fam CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome_rsid.fam

# Check that rsIDs are now present
head CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome_rsid.bim

# Find overlap SNPs between CWOW and HapMap by rsID
awk '{print $2}' CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome_rsid.bim | sort -u > cwow.rsids
awk '{print $2}' ../../reference/HapMapIII_CGRCh38_updated.clean.bim | sort -u > hapmap.rsids
comm -12 cwow.rsids hapmap.rsids > overlap.rsids

# Subset both datasets to the overlapping autosomal SNPs
plink \
  --bfile CWOW_qc_ids_sex_pheno_splitPAR_noHetFail_autosome_rsid \
  --extract overlap.rsids \
  --make-bed \
  --out CWOW_overlap_rsid

plink \
  --bfile ../../reference/HapMapIII_CGRCh38_updated.clean \
  --autosome \
  --extract overlap.rsids \
  --make-bed \
  --out HapMap_overlap_rsid

# Try the merge
plink \
  --bfile CWOW_overlap_rsid \
  --bmerge HapMap_overlap_rsid.bed HapMap_overlap_rsid.bim HapMap_overlap_rsid.fam \
  --make-bed \
  --out merge_HAP_CWOW_pca

# If PLINK reports strand mismatches, flip the problematic CWOW SNPs
plink \
  --bfile CWOW_overlap_rsid \
  --flip merge_HAP_CWOW_pca-merge.missnp \
  --make-bed \
  --out CWOW_overlap_rsid_flipped

# Try the merge again
plink \
  --bfile CWOW_overlap_rsid_flipped \
  --bmerge HapMap_overlap_rsid.bed HapMap_overlap_rsid.bim HapMap_overlap_rsid.fam \
  --make-bed \
  --out merge_HAP_CWOW_pca_try2

# If there are still problematic SNPs, remove them from both datasets
plink \
  --bfile CWOW_overlap_rsid_flipped \
  --exclude merge_HAP_CWOW_pca_try2-merge.missnp \
  --make-bed \
  --out CWOW_overlap_rsid_flipped_clean

plink \
  --bfile HapMap_overlap_rsid \
  --exclude merge_HAP_CWOW_pca_try2-merge.missnp \
  --make-bed \
  --out HapMap_overlap_rsid_clean

# Final cleaned merge for ancestry PCA
plink \
  --bfile CWOW_overlap_rsid_flipped_clean \
  --bmerge HapMap_overlap_rsid_clean.bed HapMap_overlap_rsid_clean.bim HapMap_overlap_rsid_clean.fam \
  --make-bed \
  --out merge_HAP_CWOW_pca_final
```

### PLINKQC on post imputed SNPs
plinkQC was used as a diagnostic/QC visualization step, while the actual filtering had already been applied upstream in your TOPMed/PLINK pipeline.
```
R 02_PLINK_QC_post_imputation.Rmd
```

# Final PLINK file for downstream
imputed_vcf/CWOW_qc_ids_sex_pheno_splitPAR_noHetFail\
Rename and move the final PLINK files
```
mkdir -p ../../final_gwas_dataset

cp CWOW_qc_ids_sex_pheno_splitPAR_noHetFail.bed ../../final_gwas_dataset/CWOW_TOPMED_final_postQC.bed
cp CWOW_qc_ids_sex_pheno_splitPAR_noHetFail.bim ../../final_gwas_dataset/CWOW_TOPMED_final_postQC.bim
cp CWOW_qc_ids_sex_pheno_splitPAR_noHetFail.fam ../../final_gwas_dataset/CWOW_TOPMED_final_postQC.fam

# Compute fresh PCs for GWAS
# Use autosomes only and LD pruning.
plink \
  --bfile ../../final_gwas_dataset/CWOW_TOPMED_final_postQC \
  --autosome \
  --indep-pairwise 200 50 0.2 \
  --out ../../final_gwas_dataset/CWOW_TOPMED_final_postQC_pruned

plink \
  --bfile ../../final_gwas_dataset/CWOW_TOPMED_final_postQC \
  --autosome \
  --extract ../../final_gwas_dataset/CWOW_TOPMED_final_postQC_pruned.prune.in \
  --pca 20 \
  --out ../../final_gwas_dataset/CWOW_TOPMED_final_postQC
```

### Update pheno table with new PC infomration
```
r post_imputation_pca_pheno_update.R
r obtain_informative_PCs.R
r make_cases_and_controls_pheno_tables.R
```
final PLINK prefix:\
../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC

final phenotype/covariate table:\
../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes_imputed_final.txt

## GWAS 
Genome-wide association studies (GWAS) test thousands of genetic variants across many genomes to find those statistically associated with a specific trait or disease. Here we will determine linear assocaitions between variants and disease traits of Braak NFT stage, Thal amyloid phase, and the counts of Lewy bodies in the cingulate cortex. Model adjustes for covariates PC1-5 to account for population ancestry differences and age and sex. 
```
sh GWAS_post_imputation.sh
sh GWAS_post_imputation_control_vs_disease.sh
```

Make Manhattan plots from GWAS association tests
```
# Manhattan plots from GWAS 
R 07_GWAS_imputed_linear.Rmd
R 07_GWAS_imputed_age.Rmd
R 07_GWAS_imputed_logistic.Rmd
R 08_GWAS_imputed_sex.R 
```

### Update metadata file, counts data, and genotype file 
regenerate the additive genotype export from
```
# Regenerate the additive genotype file from the final 579-sample PLINK set
cd /LBD_CWOW_QTL/scripts/

plink \
  --bfile ../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC \
  --recodeA \
  --out ../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC.clean_genotype

cd /tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/snp_array/final_gwas_dataset

r create_n579_metadata.R
```

# Process files for eQTL
```
R chr_to_rsid.R
plink \
  --bfile ../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC \
  --update-name ../snp_array/reference/update_rsid_map.txt \
  --make-bed \
  --out ../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC_rsID
# Now our .bim will contain rsIDs.


sh format_genotype_file.sh
R process_files_for_eQTL.Rmd
```




# --- HERE
```
awk '{gsub(/ /,"\t"); print}' gwas_filtered_data.clean_genotype.raw > output_file_additional_filtering.txt
salloc --mem=200G --cpus-per-task=4 --time=2:00:00

datamash transpose < output_file_additional_filtering.txt > transposed_data_additional_filters.txt
sed -e '1d; 3,6d' transposed_data.txt > ../gwas_filtered_data.genotype_formatted.txt

R 04_process_impute_genotype.Rmd
```

# Get SNP annotation
```
sh get_GRCh38_SNP_annotations.sh # downloads the dbSNP NCBI All_20180418.vcf.gz
tabix -p vcf All_20180418.vcf.gz # index 
bcftools view -m2 -M2 All_20180418.vcf.gz > All_20180418_biallelic_output.vcf.gz # get only biallelic sites 
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


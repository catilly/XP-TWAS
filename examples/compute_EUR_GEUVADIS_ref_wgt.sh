#!/bin/bash

#####
# This script is designed to compute weights using the Geuvadis gene expression and genotype data within European individuals.
# It uses lasso, top1, enet, and blup to do so.
# The resulting weights can be used to impute gene expression scores into another population given new genotype data
# OR used as references to obtain xpen and xpen2 weights in another population.
#
# Example: sh compute_EUR_GEUVADIS_ref_wgt.sh /bcb/agusevlab/cali/GEUVADIS/GeneQuant_filtered/GeneQuantRPKM.filt.chr21.txt /bcb/agusevlab/cali/GENOS/1000G.21 /bcb/agusevlab/cali/FUSION/Tests
#####

QUANT_DATA=$1 # Full path to quantification data table, assumed to be from a particular chromosome
# We used Geuvadis quantification data (GD462.GeneQuantRPKM.50FN.samplename.resk10.txt),
# which can be downloaded from https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/;
# and broke it up by chromosome into 22 files
ORIG_PLINK=$2 # Full path to the genotypes (in plink format) that correspond with the Geuvadis gene expression data;
# Should be from the same chromosome as QUANT_DATA
OUT_PATH=$3 # Full path to directory where output will be stored

GCTA=<PATH TO GCTA>
PLINK=<PATH TO PLINK>
GEMMA=<PATH TO GEMMA>
R_LOC=<PATH TO DIRECTORY WHERE R AND RSCRIPT LIVE>
LDREF=<PATH TO LD REFERENCE>
# LDREF can be obtained from the FUSION website or https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2;
# it's used to restrict input SNPs to reference IDs only
GENERATE_PHENOS=<path/to/generate_phenos.R>
COMPUTE_WEIGHTS=<path.to/FUSION.compute_weights.CL.R>
# generate_phenos.R and FUSION.compute_weights.CL.R are included in this repo

#####

mkdir -p tmp

# Get chromosome number from input file
CHR=$(echo $QUANT_DATA| cut -d'.' -f 3)
CHR_NUM=`echo $CHR | tr -d -c 0-9`

# Create SNP lists and phenotype files for EUR and YRI
$R_LOC/R --slave --args $QUANT_DATA $LDREF/1000G.EUR.${CHR_NUM}.bim $OUT_PATH < ${GENERATE_PHENOS}

# Get list of TargetIDs from input file
ID_LIST=$(awk 'NR>1{print $1}' ${QUANT_DATA})

# Loop over genes
for i in ${ID_LIST}; do

  # Create relevant plinks
  for PHENO in ${OUT_PATH}/pop1.pheno.${CHR}.TargetID_${i} ${OUT_PATH}/pop2.pheno.${CHR}.TargetID_${i}; do
    $PLINK --bfile $ORIG_PLINK --extract ${OUT_PATH}/locus.SNP.list.${CHR}.TargetID_${i}.txt --keep ${PHENO}.txt --pheno ${PHENO}.txt --make-bed --allow-no-sex --out ${PHENO}
  done

  # Compute weights
  # Note that we have --hsq_set 0.05 to prevent gene skipping
  $R_LOC/Rscript $COMPUTE_WEIGHTS --bfile ${OUT_PATH}/pop1.pheno.${CHR}.TargetID_${i} --tmp tmp/EUR.${i}.tmp --out ${OUT_PATH}/EUR.${i} --save_hsq --PATH_plink $PLINK --PATH_gemma $GEMMA --PATH_gcta $GCTA --models lasso,top1,enet,blup --hsq_set 0.05

done

# Cleanup
rm ${OUT_PATH}/*pheno.${CHR}.TargetID_*
rm ${OUT_PATH}/locus.SNP.list.${CHR}.TargetID_*

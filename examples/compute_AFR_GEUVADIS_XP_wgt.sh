#!/bin/bash

#####
# This script is designed to compute weights using the Geuvadis gene expression and genotype data within African individuals.
# It is designed to be used AFTER running compute_EUR_GEUVADIS_ref_wgt.sh
# because, in addition to using lasso, top1, enet, and blup to compute weights,
# it also incorporates European reference weights in implementing xpen and xpen2.
# The resulting weights can be used to impute gene expression scores into another population given new genotype data.
#
# Example: sh compute_AFR_GEUVADIS_XP_wgt.sh /bcb/agusevlab/cali/GEUVADIS/GeneQuant_filtered/GeneQuantRPKM.filt.chr21.txt /bcb/agusevlab/cali/GENOS/1000G.21 /bcb/agusevlab/cali/FUSION/Tests /bcb/agusevlab/cali/FUSION/WEIGHTS/NTR_YFS_lookup
#####

QUANT_DATA=$1 # Full path to quantification data table, assumed to be from a particular chromosome
# We used Geuvadis quantification data (GD462.GeneQuantRPKM.50FN.samplename.resk10.txt),
# which can be downloaded from https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/;
# and broke it up by chromosome into 22 files
ORIG_PLINK=$2 # Full path to the genotypes (in plink format) that correspond with the Geuvadis gene expression data;
# Should be from the same chromosome as QUANT_DATA
OUT_PATH=$3 # Full path to directory where output will be stored
LOOKUP_DIR=${4:-none} # Optional; directory with files, for each chromosome, that contain paths to additional reference weights for each TargetID
# The first column of the lookup file should contain gene IDs (e.g., ENSG00000120437)
# and the second column should contain paths to references separated by commas
# (e.g., </path/to/reference_1_weights_for_this_GeneID.wgt.RDat,/path/to/reference_2_weights_for_this_GeneID.wgt.RDat>).
# If not specified, only the EUR Geuvadis reference data in DEFAULT_REFS_DIR will be used

GCTA=<PATH TO GCTA>
PLINK=<PATH TO PLINK>
GEMMA=<PATH TO GEMMA>
R_LOC=<PATH TO DIRECTORY WHERE R AND RSCRIPT LIVE>
LDREF=<PATH TO LD REFERENCE>
# LDREF CAN BE OBTAINED FROM FUSION WEBSITE or https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
# THIS IS USED TO RESTRICT INPUT SNPS TO REFERENCE IDS ONLY
GENERATE_PHENOS=</path/to/generate_phenos.R>
COMPUTE_WEIGHTS=</path/to/FUSION.compute_weights.CL.R>
# generate_phenos.R and FUSION.compute_weights.CL.R are included in this repo
DEFAULT_REFS_DIR=<PATH TO DIRECTORY WHERE EUR wgt.RDat FILES ARE STORED>

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

  if [ ${LOOKUP_DIR} == none ]
  then
    FINAL_REFS=${DEFAULT_REFS_DIR}/EUR.${i}.wgt.RDat
  else
    # Chop off final digit in transcript ID, then look up reference weight file paths and store in REFS
    LOOKUP_KEY=$(echo "$i" | awk -F'.' '{print $1}')
    LOOKUP_FILE=${LOOKUP_DIR}/lookup.chr${CHR_NUM}.txt
    REFS=$(awk -vLOOKUPVAL=${LOOKUP_KEY} '$1 == LOOKUPVAL { print $2 }' < ${LOOKUP_FILE})
    # A WORKAROUND -- fill in the EUR Geuvadis reference weights, even if this transcript is not in our lookup table
    if [ -z "${REFS}" ]; then
      FINAL_REFS=${DEFAULT_REFS_DIR}/EUR.${i}.wgt.RDat
    else
      FINAL_REFS=${DEFAULT_REFS_DIR}/EUR.${i}.wgt.RDat,${REFS}
    fi
  fi

  # Compute weights
  # Note that we have --hsq_set 0.05 to prevent gene skipping
  $R_LOC/Rscript $COMPUTE_WEIGHTS --bfile ${OUT_PATH}/pop2.pheno.${CHR}.TargetID_${i} --tmp tmp/AFR.${i}.tmp --out ${OUT_PATH}/AFR.${i} --save_hsq --PATH_plink $PLINK --PATH_gemma $GEMMA --PATH_gcta $GCTA --models lasso,top1,enet,blup,xpen,xpen2 --ref_wgt ${FINAL_REFS} --ref_cutoff 0.05 --hsq_set 0.05

done

# Cleanup
rm ${OUT_PATH}/*pheno.${CHR}.TargetID_*
rm ${OUT_PATH}/locus.SNP.list.${CHR}.TargetID_*

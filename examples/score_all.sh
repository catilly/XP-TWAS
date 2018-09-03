#!/bin/bash

#####
# This script uses FUSION weights to predict gene expression scores into a new population;
# these scores can then be used to find gene-phenotype associations.
#
# Example: sh score_all.sh /bcb/agusevlab/cali/FUSION/WEIGHTS/AFR.GEUVADIS /bcb/agusevlab/DATA/GWAS/UKBB/ukbPCs.AFR /bcb/agusevlab/cali/FUSION/SCORES/ukbPCs.AFR_using_GEUV_NTR_YFS xpen2
#####

WGT_DIR=$1 # Full path of directory to look for FUSION wgt.RDat files
GENOS=$2 # Genotypes (in PLINK format) in which we'll make predictions
OUT_PATH=$3 # Full path to directory where output will be stored
FORCE_MODEL=${4:-best} # Optional; forces the scoring script to use a specific model (e.g., blup, lasso, enet, xpen, xpen2, top1)
# If not specified, the model with the best CV performance will be used

PLINK=<PATH TO PLINK>
R_LOC=<PATH TO DIRECTORY WHERE R AND RSCRIPT LIVE>
SCORE_SCRIPT=</path/to/FUSION.score.CL.R> # Included in this repo

#####

mkdir -p ${OUT_PATH}/err
mkdir -p ${OUT_PATH}/out

for file in ${WGT_DIR}/*.wgt.RDat; do
  BASENAME=$(basename $file .wgt.RDat)
  if [ ${FORCE_MODEL} == best ]
  then
    qsub -b y -cwd -o ${OUT_PATH}/out -e ${OUT_PATH}/err ${R_LOC}/Rscript ${SCORE_SCRIPT} --wgt $file --bfile $GENOS --out ${OUT_PATH}/${BASENAME} --tmp ${OUT_PATH}/${BASENAME}.tmp --PATH_plink ${PLINK}
  else
    qsub -b y -cwd -o ${OUT_PATH}/out -e ${OUT_PATH}/err ${R_LOC}/Rscript ${SCORE_SCRIPT} --wgt $file --bfile $GENOS --out ${OUT_PATH}/${BASENAME} --tmp ${OUT_PATH}/${BASENAME}.${FORCE_MODEL}.tmp --PATH_plink ${PLINK} --force_model ${FORCE_MODEL}
  fi
done
# XP-TWAS

## Introduction

This repository is a extension of the [FUSION/TWAS software](https://github.com/gusevlab/fusion_twas) (Gusev et al. 2016). It enables the use of two new methods, cross-population elastic net (XPEN) and two-step cross-population elastic net (XPEN2), in (1) computing weights that relate genotype to gene expression, and (2) imputing gene expression scores for another population given new genotype data.

The resulting gene expression scores can be used to test for gene-phenotype associations; the XP methods, in particular, are designed to increase discovery power in understudied populations.

## Description of Contents

This directory contains the following sub-directories and files:
- **FUSION.compute_weights.CL.R**: Extension of the original FUSION.compute_weights.R script by Gusev et al. 2016 that enables the use of XPEN and XPEN2 to compute weights modeling the relationship between genotype and gene expression.
- **FUSION.score.CL.R**: Extension of the original make_score.R script by Gusev et al. 2016 that can handle XPEN and XPEN2 input weights to produce individual-level gene expression scores given new genotype data.
- **examples/**: Some examples to show how these R utilities are used in practice.
	- **compute_EUR_GEUVADIS_ref_wgt.sh**: An example bash script that computes weights in European Geuvadis data.
	- **compute_AFR_GEUVADIS_XP_wgt.sh**: An example bash script that computes weights in African Geuvadis data. It is designed for use AFTER weights in the Geuvadis Europeans have already been computed; it uses the European weights as references in implementing XPEN and XPEN2.
	- **generate_phenos.R**: A helper R script that is called in the two compute_wgt.sh scripts above in order to generate phenotype files for each population.
	- **score_all.sh**: An example bash script that uses weight files to compute gene expression scores in new genotype data.
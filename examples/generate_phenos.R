### GENERATE_PHENOS ###
# This R script is used in the compute_*_wgt.sh scripts in the examples directory.
#
# It loops through genes in a specified chromosome and generates
# 1. a SNP list that represents the intersection of
# common SNPs (>1% MAF in either population) lying the region +-500kb around the gene TSS
# AND an LD reference file,
# 2. a phenotype file for the EUR individuals in the 1000 Genomes data,
# 3. and a phenotype file for the YRI individuals in the 1000 Genomes data.

### CHANGE THESE PATHS AS NECESSARY ###
wd <- "/bcb/agusevlab/cali" # Working directory
1000G_loc <- "./GENOS/1000G." # Path to 1000 Genomes genotype data in plink format.
# Note: 1000G data should be split by chromosome, and should include 1000G.POP, which indicates which individuals are members of each population

setwd(wd)
library(snpStats)
library(dplyr)
source("simulate_expression_prediction_v2.R")

arg = commandArgs(trailingOnly=T)
QUANT.DATA = toString(arg[1])
LD.REF = toString(arg[2])
PAR.OUTDIR = toString(arg[3])

set.seed( 1 )

### Load in quantification data (assumed to be all from the same chromosome!)
quant <- read.table(QUANT.DATA, header = TRUE)
### Load in LDREF data
ldref <- read.table(LD.REF)
ldref.snps <- as.character(ldref[, 2])
### Load in data (only from the relevant chromosome to speed things up) ###
plink <- read.plink(paste0(1000G_loc, quant[1,]$Chr, ".bed"),
                    paste0(1000G_loc, quant[1,]$Chr, ".bim"),
                    paste0(1000G_loc, quant[1,]$Chr, ".fam"))
pop = read.table(paste0(1000G_loc, "POP")) # Read in population data

# Initialize parameters and vectors to store corresponding scores
window.size = 1000000 # Specify window size
rarity = "common"

for (i in 1:nrow(quant)){
  
  # Grab the relevant expression values for
  E <- t(quant[i, -c(1:4)]); colnames(E) <- "E"
  # Extract SNPs that are near the gene of interest
  X = extract.X.from.plink(plink, quant[i,]$Coord - window.size/2, quant[i,]$Coord + window.size/2)
  
  # If there are zero or one SNPs near this locus, skip
  if ((ncol(X) < 2) || (is.null(ncol(X))) || (ncol(X) > 10000)){
    cat("Too few SNPs selected -- skipping gene\n")
    next
  }
  
  # Check if we want only the variants with MAF >1% in either population
  if (rarity == "common"){
    MAF1 = apply( X[which(pop[,2] != "YRI"),], 2, function(x) min( mean(x-1)/2, 1-mean(x-1)/2 ) )
    MAF2 = apply( X[which(pop[,2] == "YRI"),], 2, function(x) min( mean(x-1)/2, 1-mean(x-1)/2 ) )
    SNP.subset <- which(MAF1 > 0.01 | MAF2 > 0.01) # Note that we define common as >1% MAF here
    
    X <- X[, SNP.subset, drop = FALSE]
    
    # Check again if there are enough SNPs
    if ((ncol(X) < 2) || (is.null(ncol(X))) || (ncol(X) > 10000)){
      cat("Too few SNPs selected -- skipping gene\n")
      next
    }
  }

  # Write locus SNP list
  locus.SNP.list = intersect(colnames(X), ldref.snps)
  write.table( locus.SNP.list ,
               file = file.path(PAR.OUTDIR,
                                paste("locus.SNP.list.chr", quant[1,]$Chr,
                                      ".TargetID_", quant[i,]$TargetID,
                                      ".txt", sep="")),
               row.names = FALSE,
               col.names = FALSE,
               quote = F )
  
  # Merge the expression and genotype data; format
  df = transform(merge(X, E, by = "row.names", all.x = TRUE), row.names=Row.names, Row.names=NULL)
  df <- df %>% select(E, everything())
  
  # Subset by population
  df1 = df[which(pop[,2] != "YRI"),]
  df2 = df[which(pop[,2] == "YRI"),]
  
  # Pop 1 phenotypes (used for getting reference weights / variable selection)
  pop1.pheno = data.frame(row.names(df1), row.names(df1), df1[,1])
  write.table( pop1.pheno ,
               file = file.path(PAR.OUTDIR, paste("pop1.pheno.chr", quant[1,]$Chr,
                                                  ".TargetID_", quant[i,]$TargetID,
                                                  ".txt", sep="")),
               row.names = FALSE,
               col.names = FALSE,
               quote = F )
  
  # Pop 2 phenotypes (used for ethnic-specific weights)
  pop2.pheno = data.frame(row.names(df2), row.names(df2), df2[,1])
  write.table( pop2.pheno ,
               file = file.path(PAR.OUTDIR, paste("pop2.pheno.chr", quant[1,]$Chr,
                                                  ".TargetID_", quant[i,]$TargetID,
                                                  ".txt", sep="")),
               row.names = FALSE,
               col.names = FALSE,
               quote = F )
}


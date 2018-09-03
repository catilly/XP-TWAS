# --- FUSION.score.CL.R
# Builds on the original FUSION individual-levle scoring script, make_score.R (Gusev et al. 2016).
# In addition to the capabilities of the original script,
# this script can perform individual-level prediction in a PLINK file using XPEN and XPEN2.

# Note that multiple score files are generated when XPEN2 is used for scoring: these represent
# the predictions from each set of reference weights individually,
# only the combined contribution of the reference weights,
# only the contribution of the target population weights,
# and the final combination of reference and target weights.
# ---

suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('methods'))

allele.qc = function(a1,a2,ref1,ref2) {
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)
  
  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip
  
  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;
  
  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  
  return(snp)
}


option_list = list(
  make_option("--wgt", action="store", default=NA, type='character',
              help="Path to FUSION wgt.RDat file [required]"),
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--force_model", action="store", default=NA, type='character',
              help="Name of model to use for scoring. Must be present in the wgt.RDat file. If not specified, model with the best cv performance is used [optional]\n
              Available models:\n
              top1:\tTop eQTL (standard marginal eQTL Z-scores always computed and stored)\n
              blup:\t Best Unbiased Linear Predictor (dual of ridge regression)\n
              bslmm:\t Bayesian Sparse Linear Model (spike/slab MCMC)\n
              lasso:\t LASSO regression (with heritability used as lambda)\n
              enet:\t Elastic-net regression (with mixing parameter of 0.5)\n
              xpen:\t Cross-population elastic-net\n
              xpen2:\t Two-step cross-population elastic-net\n")
  )

opt = parse_args(OptionParser(option_list=option_list))

### Load in weights

load( opt$wgt )

### Load in genotypes

# Write snp list using the wgt file
snp.file = paste(opt$tmp,".snps",sep='')
snp.list = snps[,2]
write.table(snp.list,quote=F,row.names=F,col.names=F,file=snp.file)

# Restrict PLINK to relevant SNPs only, then read PLINK
geno.file = paste(opt$tmp,".geno",sep='')
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --extract ",snp.file," --make-bed --out ",geno.file,sep='')
system(arg)
genos = read_plink(geno.file,impute="avg")
# important : genotypes are standardized and scaled here:
genos$bed = scale(genos$bed)

# --- CLEANUP AND I/O CHECKS
cleanup = function() {
  arg = paste("rm -f " , opt$tmp , "*", sep='')
  system(arg)
}

if ( system( paste(opt$PATH_plink,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
  cat( "ERROR: plink could not be executed, set with --PATH_plink\n" , sep='', file=stderr() )
  cleanup()
  q()
}

if ( !(opt$force_model %in% colnames(wgt.matrix)) ){
  cat( "ERROR: model specified with --force_model not found in wgt.RDat\n" , sep='', file=stderr() )
  cleanup()
  q()
}

#---

### Weight matrix preprocessing

# Get weights for the model to be used for scoring
if (is.na(opt$force_model)){
  best = which.min(cv.performance[2,])
  
  # Check if the "best" model actually has valid overall weights;
  # sometimes one of the CV folds will appear good, but the overall fit doesn't converge
  # (observed occasionally with lasso)
  while (all(is.na(wgt.matrix[,best]))){
    cv.performance = cv.performance[, -best]
    wgt.matrix = wgt.matrix[, -best]
    best = which.min(cv.performance[2,])
  }
} else {
  best = which(colnames(wgt.matrix) == opt$force_model)
  names(best) = opt$force_model
}

# Replace NAs with zeros now
wgt.matrix[is.na(wgt.matrix)] <- 0

# Match up SNP ordering
m = match( genos$bim[,2] , snps[,2])
m.keep = !is.na(m)
snps = snps[m[m.keep],] # Now the reference snps are matched with genos
wgt.matrix = wgt.matrix[m[m.keep], , drop = FALSE] # and so the are reference weights
cur.bim = genos$bim[m.keep, ]
# Flip weights for mismatching alleles
qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

### Use weights to predict into genotypes

if ( names(best) == "xpen2" ) {

  W = length(ref.coefs)
  
  ref.preds = matrix(nrow = nrow(genos$bed), ncol = W)
  for (w in 1:W){
    cur.ref.colname = grep(paste0("ref.", w), colnames(wgt.matrix), value=TRUE)
    cur.ref.pred = genos$bed %*% wgt.matrix[,cur.ref.colname]
    ref.preds[, w] = cur.ref.pred
    
    # Write intermediate prediction
    # NOTE: these are NOT scaled
    write.table( cbind( genos$fam[,1] , genos$fam[,2] , cur.ref.pred) ,
                 quote=F ,
                 row.names=F , col.names=F ,
                 sep='\t' ,
                 file=paste0(opt$out, ".ref.", w, ".profile"))
    
  }
  target.pred = genos$bed %*% wgt.matrix[,"xpen2"]
  scaled.ref.preds = scale(ref.preds); scaled.ref.preds[is.na(scaled.ref.preds)] <- 0 # (prevents NA propagation if ref preds are all the same)
  combined.ref.pred = scaled.ref.preds %*% ref.coefs
  final.pred = target.pred + combined.ref.pred
  
  # Write prediction using only population-specific xpen2 weights
  write.table( cbind( genos$fam[,1] , genos$fam[,2] , target.pred) ,
               quote=F ,
               row.names=F , col.names=F ,
               sep='\t' ,
               file=paste0(opt$out, ".xpen2.target.profile"))
  
  # Write prediction using only the contribution from the reference populations
  # NOTE: these are scaled and weighted using the reference coefficients
  write.table( cbind( genos$fam[,1] , genos$fam[,2] , combined.ref.pred) ,
               quote=F ,
               row.names=F , col.names=F ,
               sep='\t' ,
               file=paste0(opt$out, ".xpen2.refs.profile"))
  
  # Write prediction using only population-specific xpen2 weights
  write.table( cbind( genos$fam[,1] , genos$fam[,2] , final.pred) ,
               quote=F ,
               row.names=F , col.names=F ,
               sep='\t' ,
               file=paste0(opt$out, ".xpen2.combined.profile"))
  
} else {

  best.wgt = wgt.matrix[,best]
  if ( names(best) == "top1" ){
    best.wgt[ - which.max( best.wgt^2 ) ] = 0
  }
  final.pred = genos$bed %*% best.wgt
  write.table( cbind( genos$fam[,1] , genos$fam[,2] , final.pred) ,
               quote=F ,
               row.names=F , col.names=F ,
               sep='\t' ,
               file=paste0(opt$out, ".", names(best), ".profile"))
}

# Cleanup
cleanup()
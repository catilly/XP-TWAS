# --- FUSION.compute_weights.CL.R
# Builds on the original FUSION weight-computing script, FUSION.compute_weights.R (Gusev et al. 2016).
# In addition to the capabilities of the original,
# this script can compute weights using XPEN and XPEN2 given one or more reference weight files.
# ---

# --- Change log
# 6/30/18: Implemented XPEN and XPEN2.
# 7/4/2018: Implemented XPEN2 with multiple reference populations.

# --- TODO
# * Make sure BLUP/BSLMM weights are being scaled properly based on MAF
# * Implement the ref_model option for multiple reference populations in XPEN2; currently always uses best model according to CV results
# * Allow user to specify different reference weights for XPEN and XPEN2
# * Store reference weights in a matrix rather than a list -- would speed things up


suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
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
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--ref_wgt", action="store", default=NA, type='character',
              help="Comma-separated list of paths to files containing weights previously computed from reference populations. 
              Only the first set of weights will be used for xpen; all sets will be used for xpen2. [required for xpen, xpen2]"),
  make_option("--ref_cutoff", action="store", default=0.05, type='double',
              help="P-value cutoff or number of top SNPs.
               Used in xpen to select variables from a reference population [default: %default]"),
  # make_option("--ref_model", action="store", default=NA, type='character',
  #             help="Comma-separated list of models that will be used for the first step of xpen2 [optional; model with the best cv metrics used otherwise]\n
  #             Available models:\n
  #             top1:\tTop eQTL (standard marginal eQTL Z-scores always computed and stored)\n
  #             blup:\t Best Unbiased Linear Predictor (dual of ridge regression)\n
  #             bslmm:\t Bayesian Sparse Linear Model (spike/slab MCMC)\n
  #             lasso:\t LASSO regression (with heritability used as lambda)\n
  #             enet:\t Elastic-net regression (with mixing parameter of 0.5)\n"),
  make_option("--pheno", action="store", default=NA, type='character',
              help="Path to molecular phenotype file (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gcta", action="store", default="gcta_nr_robust", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gemma", action="store", default="gemma", type='character',
              help="Path to plink executable [%default]"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
  make_option("--hsq_p", action="store", default=0.01, type='double',
              help="Minimum heritability p-value for which to compute weights [default: %default]"),
  make_option("--hsq_set", action="store", default=NA, type='double',
              help="Skip heritability estimation and set hsq estimate to this value [optional]"),
  make_option("--crossval", action="store", default=5, type='double',
              help="How many folds of cross-validation, 0 to skip [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--rn", action="store_true", default=FALSE,
              help="Rank-normalize the phenotype after all QC: [default: %default]"),
  make_option("--save_hsq", action="store_true", default=FALSE,
              help="Save heritability results even if weights are not computed [default: %default]"),			  
  make_option("--models", action="store", default="blup,lasso,top1,enet", type='character',
              help="Comma-separated list of prediction models [default: %default]\n
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
models = unique( c(unlist(strsplit(opt$models,',')),"top1") )
M = length(models)

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

# --- PREDICTION MODELS

# BSLMM
weights.bslmm = function( input , bv_type , snp , out=NA ) {
  if ( is.na(out) ) out = paste(input,".BSLMM",sep='')
  
  arg = paste( opt$PATH_gemma , " -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -bfile " , input , " -bslmm " , bv_type , " -o " , out , sep='' )
  system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
  eff = read.table( paste(out,".param.txt",sep=''),head=T,as.is=T)
  eff.wgt = rep(NA,length(snp))
  m = match( snp , eff$rs )
  m.keep = !is.na(m)
  m = m[m.keep]
  eff.wgt[m.keep] = (eff$alpha + eff$beta * eff$gamma)[m]
  return( eff.wgt )
}

# PLINK: LASSO
weights.lasso = function( input , hsq , snp , out=NA ) {
  if ( is.na(out) ) out = paste(input,".LASSO",sep='')
  
  arg = paste( opt$PATH_plink , " --allow-no-sex --bfile " , input , " --lasso " , hsq , " --out " , out , sep='' )
  system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT )
  if ( !file.exists(paste(out,".lasso",sep='')) ) {
    cat( paste(out,".lasso",sep='') , " LASSO output did not exist\n" )
    eff.wgt = rep(NA,length(snp))
  } else {
    eff = read.table( paste(out,".lasso",sep=''),head=T,as.is=T)
    eff.wgt = rep(0,length(snp))
    m = match( snp , eff$SNP )
    m.keep = !is.na(m)	
    m = m[m.keep]
    eff.wgt[m.keep] = eff$EFFECT[m]
  }
  return( eff.wgt )
}

# Marginal Z-scores (used for top1)
weights.marginal = function( genos , pheno , beta=F ) {
  if ( beta ) eff.wgt = t( genos ) %*% (pheno) / ( nrow(pheno) - 1)
  else eff.wgt = t( genos ) %*% (pheno) / sqrt( nrow(pheno) - 1 )
  return( eff.wgt )
}

# Elastic Net
weights.enet = function( genos , pheno , alpha=0.5 ) {
  eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
  # remove monomorphics
  sds = apply( genos  , 2 , sd )
  keep = sds != 0 & !is.na(sds)
  enet = cv.glmnet( genos[,keep] , pheno , alpha=alpha , nfold=5 , intercept=T , standardize=F )
  eff.wgt[ keep ] = coef( enet , s = "lambda.min")[2:(sum(keep)+1)]
  return( eff.wgt )
}

# CL: Cross-Population Elastic Net
weights.xpen = function( genos , pheno , ref , cutoff , alpha=0.5 ) {
  eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
  # remove monomorphics
  sds = apply( genos  , 2 , sd )
  keep = sds != 0 & !is.na(sds)
  
  # XP-EN
  # Get C1 SNPs based on reference weights and cutoff
  if (cutoff < 1){
    # If using a p-value cutoff, convert to z-score and get SNPs that fall beyond that score
    z.cutoff = qnorm(1 - cutoff/2)
    C1 = names(ref[abs(ref) > z.cutoff])
  } else {
    # Else, take the n SNPs with the largest effect sizes (in magnitude)
    C1 = names(head(sort(abs(ref), decreasing = TRUE), n = round(cutoff)))
  }
  
  cur.best.cvm = 1000000
  best.model = NA
  # First check if C1 is length 0 to save some computational time
  if (length(C1) == 0){
    best.model = cv.glmnet(genos[,keep], pheno, alpha=alpha, nfold=5, intercept=T, standardize=F)
  } else {
    for (scaling.factor in c(0.01, 0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000)){
      # Check which SNPs are in C1
      # We shrink the SNPs that aren't in C1 more by a factor of scaling.factor
      penalty = colnames(genos[,keep]) %in% C1
      penalty[penalty == 0] <- scaling.factor
      # Fit EN with alpha = 0.5
      model = cv.glmnet(genos[,keep], pheno, alpha=alpha, nfold=5, intercept=T, standardize=F,
                        penalty.factor = penalty)
      # Get model with the best scaling.factor
      if (min(model$cvm) < cur.best.cvm){
        best.model = model
        cur.best.cvm = min(model$cvm)
      }
    }
  }
  
  eff.wgt[ keep ] = coef( best.model , s = "lambda.min")[2:(sum(keep)+1)]
  return( eff.wgt )
}

# CL: Two-Step Cross-Population Elastic Net
weights.xpen2 = function( genos , pheno , refs , alpha=0.5 ) {
  eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
  # remove monomorphics
  sds = apply( genos  , 2 , sd )
  keep = sds != 0 & !is.na(sds)
  
  # Two-step XP-EN
  # Use reference weights (assumed to be in a list) to predict into sample
  W = length(refs)
  ref.preds = matrix(nrow = nrow(genos), ncol = W)
  for (w in 1:W){
    ref.preds[, w] = scale(genos %*% refs[[w]])
  }
  ref.preds[is.na(ref.preds)] <- 0 # Prevent NAs from propagating if the reference weights are bad
  new.train = cbind(genos[,keep], ref.preds)
  
  # Incorporate Step 1 output columns as unpenalized predictors in EN trained in pop2
  penalty = rep(1, ncol(new.train)); penalty[(ncol(new.train)-W+1):ncol(new.train)] = 0
  xpen2 = cv.glmnet(new.train, pheno, alpha=alpha, nfold=5, intercept=T, standardize=F,
                    penalty.factor = penalty)
  eff.wgt[ keep ] = coef( xpen2 , s = "lambda.min")[2:(sum(keep)+1)]
  ref.coefs = coef( xpen2 , s = "lambda.min")[(sum(keep)+2):(sum(keep)+W+1)]
  return( list(eff.wgt, ref.coefs) )
}

# --- CLEANUP
cleanup = function() {
  if ( ! opt$noclean ) {
    arg = paste("rm -f " , opt$tmp , "*", sep='')
    system(arg)
  }
}

# Perform i/o checks here:
files = paste(opt$bfile,c(".bed",".bim",".fam"),sep='')
if ( !is.na(opt$pheno) ) files = c(files,opt$pheno)
if ( !is.na(opt$covar) ) files = c(files,opt$covar)

for ( f in files ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
    cleanup()
    q()
  }
}

if ( system( paste(opt$PATH_plink,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
  cat( "ERROR: plink could not be executed, set with --PATH_plink\n" , sep='', file=stderr() )
  cleanup()
  q()
}

if ( !is.na(opt$hsq_set) && system( opt$PATH_gcta , ignore.stdout=T,ignore.stderr=T ) != 0 ){
  cat( "ERROR: gcta could not be executed, set with --PATH_gcta\n" , sep='', file=stderr() )
  cleanup()
  q()
}

if ( sum(models=="bslmm" | models=="blup") != 0 && system( paste(opt$PATH_gemma,"-h") , ignore.stdout=T,ignore.stderr=T ) != 0 ){
  cat( "ERROR: gemma could not be executed, set with --PATH_gemma or remove 'bslmm' and 'blup' from models\n" , sep='', file=stderr() )
  cleanup()
  q()
}

# CL: Check for reference weights if xpen or xpen2 is among the models
if ( sum(models=="xpen" | models=="xpen2") != 0 && is.na(opt$ref_wgt) ){
  cat( "ERROR: no reference weights found, set with --ref_wgt or remove 'xpen' and 'xpen2' from models\n" , sep='', file=stderr() )
  cleanup()
  q()
}
# ---

fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)

# Make/fetch the phenotype file
if ( !is.na(opt$pheno) ) {
  pheno.file = opt$pheno
  pheno = read.table(pheno.file,as.is=T)
  # Match up data
  m = match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
  m.keep = !is.na(m)
  fam = fam[m.keep,]
  m = m[m.keep]
  pheno = pheno[m,]
} else {
  pheno.file = paste(opt$tmp,".pheno",sep='')
  pheno = fam[,c(1,2,6)]
  write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

if ( opt$rn ) {
  library('GenABEL')
  library(preprocessCore)
  pheno[,3] = rntransform( pheno[,3] )
  write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

# Load in the covariates if needed
if ( !is.na(opt$covar) ) {
  covar = ( read.table(opt$covar,as.is=T,head=T) )
  if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
  # Match up data
  m = match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
  m.keep = !is.na(m)
  fam = fam[m.keep,]
  pheno = pheno[m.keep,]
  m = m[m.keep]
  covar = covar[m,]
  reg = summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
  if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates\n" )
  pheno[,3] = scale(reg$resid)
  raw.pheno.file = pheno.file
  pheno.file = paste(pheno.file,".resid",sep='')
  write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)	
}

geno.file = opt$tmp
# recode to the intersection of samples and new phenotype
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --pheno ",pheno.file," --keep ",pheno.file," --make-bed --out ",geno.file,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# --- HERITABILITY ANALYSIS
if ( is.na(opt$hsq_set) ) {
  if ( opt$verbose >= 1 ) cat("### Estimating heritability\n")
  
  # 1. generate GRM
  arg = paste( opt$PATH_plink," --allow-no-sex --bfile ",opt$tmp," --make-grm-bin --out ",opt$tmp,sep='')
  system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
  
  # 2. estimate heritability
  if ( !is.na(opt$covar) ) {
    arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",raw.pheno.file," --qcovar ",opt$covar," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
  } else {
    arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",pheno.file," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
  }
  system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
  
  # 3. evaluate LRT and V(G)/Vp
  if ( !file.exists( paste(opt$tmp,".hsq",sep='') ) ) {
    cat(opt$tmp,"does not exist, likely GCTA could not converge, skipping gene\n",file=stderr())
    cleanup()
    q()
  }
  
  hsq.file = read.table(file=paste(opt$tmp,".hsq",sep=''),as.is=T,fill=T)
  hsq = as.numeric(unlist(hsq.file[hsq.file[,1] == "V(G)/Vp",2:3]))
  hsq.pv = as.numeric(unlist(hsq.file[hsq.file[,1] == "Pval",2]))
  
  if ( opt$verbose >= 1 ) cat("Heritability (se):",hsq,"LRT P-value:",hsq.pv,'\n')
  if ( opt$save_hsq ) cat( opt$out , hsq , hsq.pv , '\n' , file=paste(opt$out,".hsq",sep='') )
  
  # 4. stop if insufficient
  if ( hsq[1] < 0 || hsq.pv > opt$hsq_p ) {
    cat(opt$tmp," : heritability ",hsq[1],"; LRT P-value ",hsq.pv," : skipping gene\n",sep='',file=stderr())
    cleanup()
    q()
  }
  
} else {
  if ( opt$verbose >= 1 ) cat("### Skipping heritability estimate\n")
  hsq = opt$hsq_set
  hsq.pv = NA
}

# read in genotypes
genos = read_plink(geno.file,impute="avg")
mafs = apply(genos$bed,2,mean)/2
sds = apply(genos$bed,2,sd)
# important : genotypes are standardized and scaled here:
genos$bed = scale(genos$bed)
pheno = genos$fam[,c(1,2,6)]
pheno[,3] = scale(pheno[,3])

# check if any genotypes are NA
nasnps = apply( is.na(genos$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
  cat( "WARNING :",sum(nasnps != 0),"SNPs could not be scaled and were zeroed out, make sure all SNPs are polymorphic\n" , file=stderr())
  genos$bed[,nasnps != 0] = 0
}

N.tot = nrow(genos$bed)
if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")


# --- CL: LOAD REFERENCE WEIGHTS IF USING XPEN OR XPEN2
if ( !is.na(opt$ref_wgt) ) {
  
  ref.wgt.paths = c(unlist(strsplit(opt$ref_wgt,',')))
  W = length(ref.wgt.paths)
  
  ref.models = rep(NA, W)
  ref.model.names = rep(NA, W)
  ref.wgts = list()
  ref.cutoff = opt$ref_cutoff
  
  for (w in 1:W){
    
    load(trimws(ref.wgt.paths[w]))
    
    ### Match up and flip SNPs as needed
    m = match( genos$bim[,2] , snps[,2])
    m.keep = !is.na(m)
    snps = snps[m[m.keep],] # Now the reference snps are matched with genos
    wgt.matrix = wgt.matrix[m[m.keep], , drop = FALSE] # and so the are reference weights
    cur.bim = genos$bim[m.keep, ]
    # Flip WEIGHTS for mismatching alleles
    qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
    wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
    
    ### Get z-score for XPEN using the first set of weights
    if (w == 1){
      ref.z = rep(0, length(genos$bim[, 2]))
      names(ref.z) = genos$bim[,2]
      ref.z[m.keep] = wgt.matrix[, "top1"]
    }
    
    ### Get the rest of the XPEN2 weights
    # Choose best model according to cv
    ref.models[w] = which.max(cv.performance[1,])
    ref.model.names[w] = colnames(wgt.matrix)[ref.models[w]]
    # Load in current set of flipped/matched weights
    ref.wgts[[w]] = rep(0, length(genos$bim[, 2]))
    names(ref.wgts[[w]]) = genos$bim[,2]
    ref.wgts[[w]][m.keep] = wgt.matrix[, ref.models[w]]
    # Check if the ref.model is top1; if so, zero out everything except the best eQTL
    if (ref.model.names[w] == "top1"){
      ref.wgts[[w]][ - which.max( ref.wgts[[w]]^2 ) ] = 0
    }
    
  }
  
  # Allow user to specify ref models for XPEN2
  # if ( !is.na(opt$ref_model) ) {
  #   ref.model = which( colnames(wgt.matrix) == opt$ref_model )
  #   if ( length(ref.model) == 0 ) {
  #     cat( "WARNING : --ref_model" , ref.model ,"does not exist; choosing model based on cv metrics instead \n")
  #     ref.model = which.max(cv.performance[1,])
  #   }	
  # } else {
  #   ref.model = which.max(cv.performance[1,])
  # }
}


# --- CROSSVALIDATION ANALYSES
set.seed(1)
cv.performance = matrix(NA,nrow=2,ncol=M)
rownames(cv.performance) = c("rsq","pval")
colnames(cv.performance) = models

if ( opt$crossval <= 1 ) {
  if ( opt$verbose >= 1 ) cat("### Skipping cross-validation\n")
} else {
  if ( opt$verbose >= 1 ) cat("### Performing",opt$crossval,"fold cross-validation\n")
  cv.all = pheno
  N = nrow(cv.all)
  cv.sample = sample(N)
  cv.all = cv.all[ cv.sample , ]
  folds = cut(seq(1,N),breaks=opt$crossval,labels=FALSE)
  
  cv.calls = matrix(NA,nrow=N,ncol=M)
  
  for ( i in 1:opt$crossval ) {
    if ( opt$verbose >= 1 ) cat("- Crossval fold",i,"\n")
    indx = which(folds==i,arr.ind=TRUE)
    cv.train = cv.all[-indx,]
    # store intercept
    intercept = mean( cv.train[,3] )
    cv.train[,3] = scale(cv.train[,3])
    
    # hide current fold
    cv.file = paste(opt$tmp,".cv",sep='')
    write.table( cv.train , quote=F , row.names=F , col.names=F , file=paste(cv.file,".keep",sep=''))	
    arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$tmp," --keep ",cv.file,".keep --out ",cv.file," --make-bed",sep='')
    system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
    
    for ( mod in 1:M ) {
      if ( models[mod] == "blup" ) {
        pred.wgt = weights.bslmm( cv.file , bv_type=2 , snp=genos$bim[,2] )
      }
      else if ( models[mod] == "bslmm" ) {
        pred.wgt = weights.bslmm( cv.file , bv_type=1 , snp=genos$bim[,2] )
      }		
      else if ( models[mod] == "lasso" ) {
        pred.wgt = weights.lasso( cv.file , hsq[1] , snp=genos$bim[,2] )
      }
      else if ( models[mod] == "enet" ) {
        pred.wgt = weights.enet( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3]) , alpha=0.5 )
      }
      else if ( models[mod] == "xpen" ) {
        pred.wgt = weights.xpen( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3]) , ref = ref.z , cutoff = ref.cutoff , alpha=0.5 )
      }		
      else if ( models[mod] == "xpen2" ) {
        res = weights.xpen2( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3]) , ref = ref.wgts , alpha=0.5 )
        pred.wgt = res[[1]]
        pred.coefs = res[[2]]
      }		
      else if ( models[mod] == "top1" ) {
        pred.wgt = weights.marginal( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
        pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
      }
      
      pred.wgt[ is.na(pred.wgt) ] = 0
      
      if ( models[mod] == "xpen2" ) {
        # If using XPEN2, perform prediction in two steps
        ref.preds = matrix(nrow = length(indx), ncol = W)
        for (w in 1:W){
          ref.preds[, w] = scale(genos$bed[ cv.sample[ indx ] , ] %*% ref.wgts[[w]])
        }
        ref.preds[ is.na(ref.preds) ] = 0 # Prevent NAs from propagating if ref preds are bad
        
        cv.calls[ indx , mod ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt + ref.preds %*% pred.coefs
      } else {
        
        # Otherwise, predict from weights into sample as usual
        cv.calls[ indx , mod ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt
        
      }
      
    }
  }
  
  # compute rsq + P-value for each model
  for ( mod in 1:M ) {
    reg = summary(lm( cv.all[,3] ~ cv.calls[,mod] ))
    cv.performance[ 1, mod ] = reg$adj.r.sq
    # If predictions have zero variance, fill in 1 as the p-value; otherwise, extract p-value
    if (nrow(reg$coef) < 2){
      cv.performance[ 2, mod ] = 1
    } else {
      cv.performance[ 2, mod ] = reg$coef[2,4]
    }
    
  }
  if ( opt$verbose >= 1 ) write.table(cv.performance,quote=F,sep='\t')
}

# --- FULL ANALYSES
if ( opt$verbose >= 1 ) cat("Computing full-sample weights\n")

# call models to get weights
wgt.matrix = matrix(0,nrow=nrow(genos$bim),ncol=M)
ref.coefs = NA
colnames(wgt.matrix) = models
rownames(wgt.matrix) = genos$bim[,2]
for ( mod in 1:M ) {
  if ( models[mod] == "blup" ) {
    wgt.matrix[,mod] = weights.bslmm( geno.file , bv_type=2 , snp=genos$bim[,2] , out=opt$tmp )
  }
  else if ( models[mod] == "bslmm" ) {
    wgt.matrix[,mod] = weights.bslmm( geno.file , bv_type=1 , snp=genos$bim[,2] , out=opt$tmp )
  }		
  else if ( models[mod] == "lasso" ) {
    wgt.matrix[,mod] = weights.lasso( geno.file , hsq[1] , snp=genos$bim[,2] , out=opt$tmp )
  }
  else if ( models[mod] == "enet" ) {
    wgt.matrix[,mod] = weights.enet( genos$bed , as.matrix(pheno[,3]) , alpha=0.5 )
  }	
  else if ( models[mod] == "xpen" ) {
    wgt.matrix[,mod] = weights.xpen( genos$bed , as.matrix(pheno[,3]) , ref = ref.z , cutoff = ref.cutoff , alpha=0.5 )
  }	
  else if ( models[mod] == "xpen2" ) {
    res = weights.xpen2( genos$bed , as.matrix(pheno[,3]) , ref = ref.wgts , alpha=0.5 )
    wgt.matrix[,mod] = res[[1]]
    ref.coefs = res[[2]]
  }	
  else if ( models[mod] == "top1" ) {
    wgt.matrix[,mod] = weights.marginal( genos$bed , as.matrix(pheno[,3]) , beta=F )
  }
}

# CL: If weights from reference populations were passed in, keep them in the output
if ( !is.na(opt$ref_wgt)  ){
  for (w in 1:W){
    wgt.matrix = cbind(wgt.matrix, ref.wgts[[w]])
    colnames(wgt.matrix)[ncol(wgt.matrix)] = paste0(ref.model.names[w], ".ref.", w)
  }
}

# save weights, rsq, p-value for each model, and hsq to output
snps = genos$bim
save( wgt.matrix , ref.coefs , snps , cv.performance , hsq, hsq.pv, N.tot , file = paste( opt$out , ".wgt.RDat" , sep='' ) )
# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()
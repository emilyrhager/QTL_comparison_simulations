#!/usr/bin/env Rscript



# Given a set of effect sizes, and a marker location, simulate the resulting LOD scores in cross 2. 

# Notation: 'A', 'B': alleles from rqtl object
# 1 vs 2: 1 refers to the cross experiment where the QTL was identified
# 2 refers to the cross experiment where we are testing the detection probability.

# For now, at least, we assume (in the future these could be relaxed but would require changes to the code):
#   -both crosses are intercrosses (ie genotypes AA, AB, BB all exist in the cross object)
#   -the trait is normally distributed ( we pull phenotypes from a random normal distribution)
#   -marker names are of the format chrX:1234, chrX_1234, or similar (ie chr, chromosome #, sep, basepair)
#   -no covariates we care about
#   -the within-genotype variance is not much different from the overall variance for the phenotype of interest
#   -the cross objects are stored as .rds files and have a column 'pgm', and geno.probs are already calculated


# Inputs:
# Path to cross object (cross B)
# Path to sims object (cross A)
# nsims per effect size
# Phenotype name in cross B
# Whether 'A' vs 'B' allele is the same direction in both crosses
# Ratio of effect sizes calculated in the file to effect sizes you wish to calculate (default 1)
# Whether to save all sims. 

# ERH cross: SW is A, BK is B.
# EPK cross: NB is B, BW is A. 

require(qtl)
require(argparse)
require(dplyr)
require(tidyr)
require(stringr)

# Parse user inputs.
parser <- ArgumentParser()

# Required
parser$add_argument("-o", "--cross_object", help="Path to rqtl cross object (saved as .rds) to use for simulations", required = T)
parser$add_argument("-p","--phenotype", help = "Phenotype to test (must match name of the cross object pheno column)", required = T)
parser$add_argument("-f", "--folder", help = "Path to folder to put output files", required = T)

# Optional
parser$add_argument("-e", "--effect_ratio", default = 1, help = "Ratio of effect size in original cross to effect size in this cross. default %(default)s.")
parser$add_argument("--file_header", help = "First part of the output file names. Default is: 'sims_{phenotype}_effect_ratio_{effect_ratio}_'")
parser$add_argument("--geno_reversed", action = "store_true", help = "Use this flag if you wish to flip the effect direction between the two crosses - in other words, if allele A in cross 1 was named allele B in cross 2.")
parser$add_argument("--sim_file", help = "optional path to the file to use for the simulations of effect sizes from cross 1. Default is to take the first file in the output folder that is named sims*prob_distrib.csv")
parser$add_argument("--lose_sims", action = "store_true", help = "Do not save the output of all simulations, just save the resulting probability distribution file. Default is to save all simulations but this file can be large.")
parser$add_argument("--round", action = "store_true", help = "Round simulated data to the nearest integer. This is for simulating count data as a rounded normal distribution - use with care!")
parser$add_argument("-d", "--dom_model", help = "Dominance model to use; if not specified will assume there is only one sims file in the folder")
parser$add_argument("-v", "--total_variance_constant", action = 'store_true', help = "Hold total variance constant; default is to hold per-genotype variance constant.")
parser$add_argument("-i", "--identical_cross_object", action = 'store_true', help = "Cross object is identical to the one in step 1 - use original peak marker name.")
args <- parser$parse_args()

# Required
cross.path <- args$cross_object
pheno.name.2 <- args$phenotype
save.folder <- args$folder

# Optional
effect.ratio <- as.numeric(args$effect_ratio)
geno.reversed <- FALSE
if(args$geno_reversed == TRUE){
  geno.reversed <- TRUE
}
if(!is.null(args$dom_model)){
  dlabel = paste0('_dom_',args$dom_model,'_')
}
if(is.null(args$dom_model)){
  dlabel = ''
}
if(args$total_variance_constant == TRUE){
    vartag = 'total'
}
if(args$total_variance_constant == FALSE){
    vartag = 'pergeno'
}
if(args$identical_cross_object == TRUE){
  sameobj = T
}
if(args$identical_cross_object == FALSE){
  sameobj = F
}


file.header <- paste0('sims_',pheno.name.2,dlabel,vartag,'_var_constant_effect_ratio_',effect.ratio)
if(!is.null(args$file_header)){
  file.header <- args$file_header
}
save.all.sims <- TRUE
if(args$lose_sims == TRUE){
  save.all.sims <- FALSE
}

if(is.null(args$sim_file)){
  sim_list <- sort(Sys.glob(file.path(save.folder,paste0('sims*',dlabel,'*',vartag,'*prob_distrib.csv'))))
  
  if(length(sim_list)==0){
    print('No sim file available! Should be supplied, or else a file sims*prob_distrib.csv must exist in the project folder')
  }
  
  if(length(sim_list)>1){
    print(paste0('Possible error: more than one sim file to pick from. Taking the first one, ', sim_list[1], ', of the list ', sim_list))
  }
  
  sim.path <- sim_list[1]
}
if(!is.null(args$sim_file)){
  sim.path <- args$sim_file
}
round = FALSE
if(args$round){
  round = TRUE
}

save.path <- file.path(save.folder,file.header)

print(paste0('Beginning... Files will be saved with prefix ',save.path))

# Read data
cross.obj <- readRDS(cross.path)
sim.df <- read.csv(sim.path)

# Method, model, nsim, chr, bayes interval should be pulled from sim.df;
# note that this assumes only one value per simulation dataframe. 
model.list <- unique(sim.df$model)
method.list <- unique(sim.df$method)
nsim.list <- unique(sim.df$sims.per.effect.size)
orig.marker.list <- unique(sim.df$orig.peak.marker)
orig.marker.name.list <- unique(sim.df$orig.peak.marker.name)
chr.list <- unique(sim.df$chromosome)
lower.bayes.list <- unique(sim.df$lower.bayes)
upper.bayes.list <- unique(sim.df$upper.bayes)
var.constant.list <- unique(sim.df$constant_variance)

if(length(model.list)>1){
  print('Problem: sim.df has multiple models; expected only one. Taking first model.')
}
if(length(method.list)>1){
  print('Problem: sim.df has multiple methods; expected only one. Taking first method.')
}
if(length(nsim.list)>1){
  print('Problem: sim.df has multiple values for nsim; expected only one. Taking first nsim.')
}
if(length(orig.marker.list)>1){
  print('Problem: sim.df has multiple original markers; expected only one. Taking first marker.')
}
if(length(orig.marker.name.list)>1){
  print('Problem: sim.df has multiple original markers; expected only one. Taking first marker.')
}
if(length(chr.list)>1){
  print('Problem: sim.df has multiple chromosomes; expected only one. Taking first chr.')
}
if(length(lower.bayes.list)>1){
  print('Problem: sim.df has multiple lower bayes limits; expected only one. Taking first.')
}
if(length(upper.bayes.list)>1){
  print('Problem: sim.df has multiple upper bayes limits; expected only one. Taking first.')
}
if(length(var.constant.list)>1){
  print('Problem: sim.df has multiple constant variance settings; expected only one. Taking first.')
}

model.2 <- as.character(model.list[[1]])
method.2 <- as.character(method.list[[1]])
nsim <- as.integer(nsim.list[[1]])
orig.marker.pos.1 <- orig.marker.list[[1]]
orig.marker.name.1 <- orig.marker.name.list[[1]]
chr.2 <- chr.list[[1]]
upper.bayes.1 <- upper.bayes.list[[1]]
lower.bayes.1 <- lower.bayes.list[[1]]
if(var.constant.list[[1]]=='total_variance_constant'){
  total.var.constant = TRUE
}
if(var.constant.list[[1]]=='per_genotype_variance_constant'){
  total.var.constant = FALSE
}

n.effect <- length(sim.df$sim.a)

# Data from cross object 2. 
pheno.sd.2 = sd(cross.obj$pheno[[pheno.name.2]],na.rm=T)
ninds = dim(cross.obj$pheno)[[1]]


if(sameobj == FALSE){
  # Find the two flanking markers in the new cross object. (Or, allow for only one if both crosses happen to have the same marker!)
  # Figure out how far they are from the original marker.
  map.2 <- cross.obj$geno[[chr.2]]$map
  
  # Check for the possibility that the original marker is further to the end
  # than any marker in cross 2:
  n.right <- sum(map.2 >= orig.marker.pos.1)
  n.left <- sum(map.2 <= orig.marker.pos.1)
  if(n.right>0){
    marker.right <- min(map.2[map.2>=orig.marker.pos.1])
  }
  if(n.left>0){
    marker.left <- max(map.2[map.2<=orig.marker.pos.1])
  }
  if(n.left==0){
    # original marker is left of all cross 2 markers
    marker.left = marker.right
  }
  if(n.right==0){
    # original marker is right of all cross 2 markers
    marker.right = marker.left
  }
  # Note: this way of handling it ^ doesn't account for genetic distance between end-most marker and end
  # of the chromosome. Improve later? 
  
  marker.dists <- c(orig.marker.pos.1 - marker.left, marker.right - orig.marker.pos.1)
  if(marker.right == marker.left){
    prop.dists <- c(0.5,0.5)
  }
  if(marker.right != marker.left){
    # note that this is low if close, high if far. 
    prop.dists <- marker.dists/sum(marker.dists)
  }
  
  # Get separator from cross object.
  # Assumes names are in the format chr[number or 'X'][something that's not a digit][more digits]
  # and extracts the [something that's not a digit] from the first marker name on the first chromosome.
  sep = str_replace(names(cross.obj$geno[[1]]$map)[1],"chr[\\dX]+([^\\d]+)\\d+", "\\1")
  
  # Get genotype probs at these markers, take a weighted average by phyiscal distance.
  pAA_L = cross.obj$geno[[chr.2]]$prob[ ,paste0('chr',chr.2,sep,marker.left),'AA']
  pAB_L = cross.obj$geno[[chr.2]]$prob[ ,paste0('chr',chr.2,sep,marker.left),'AB']
  pBB_L = cross.obj$geno[[chr.2]]$prob[ ,paste0('chr',chr.2,sep,marker.left),'BB']
  
  pAA_R = cross.obj$geno[[chr.2]]$prob[ ,paste0('chr',chr.2,sep,marker.right),'AA']
  pAB_R = cross.obj$geno[[chr.2]]$prob[ ,paste0('chr',chr.2,sep,marker.right),'AB']
  pBB_R = cross.obj$geno[[chr.2]]$prob[ ,paste0('chr',chr.2,sep,marker.right),'BB']
  
  # take weighted average; if prop_dist_left is 0.1, weight_left is 1-0.1 = 0.9. 
  weight_L = 1 - prop.dists[1]
  weight_R = 1 - prop.dists[2]
  pAA = weight_L * pAA_L + weight_R * pAA_R
  pAB = weight_L * pAB_L + weight_R * pAB_R
  pBB = weight_L * pBB_L + weight_R * pBB_R
  # ^ SHOULD I DO THIS IN A MORE COMPLEX WAY? 
}
if(sameobj == TRUE){
  # Pull probs from the original peak marker.
  
  if(!is.na(str_split(orig.marker.name.1,paste0('c',chr.2,'.'))[[1]][2])){
    orig.marker.name.1 = str_split(orig.marker.name.1,paste0('c',chr.2,'.'))[[1]][2]
  }
  pAA = cross.obj$geno[[chr.2]]$prob[ ,names(cross.obj$geno[[chr.2]]$prob[1,,'AA'])==orig.marker.name.1,'AA']
  pAB = cross.obj$geno[[chr.2]]$prob[ ,names(cross.obj$geno[[chr.2]]$prob[1,,'AB'])==orig.marker.name.1,'AB']
  pBB = cross.obj$geno[[chr.2]]$prob[ ,names(cross.obj$geno[[chr.2]]$prob[1,,'BB'])==orig.marker.name.1,'BB']
  marker.left = orig.marker.name.1
  marker.right = orig.marker.name.1
}
# Run simulations
report.interval = ceiling(nsim/5)
all.sim.peaks <- data.frame()
for(row in 1:nrow(sim.df)){
  a.orig = sim.df[row,'sim.a']
  d.orig = sim.df[row,'sim.d']
  a = a.orig * effect.ratio
  d = d.orig * effect.ratio
  print(paste0('Running ',nsim,' simulations for effect size ',round(a,2)))
  
  # Calculate mean phenotypes for each genotype class
  if(!geno.reversed){
    meanAA = 0 # arbitrary baseline.
    meanAB = a + d
    meanBB = 2*a
  }
  if(geno.reversed){
    # If the corresponding genotype labels (A vs B) are flipped, then flip the 
    # effect here. This applies, for instance, if we're interested in a forest
    # allele, and SW is allele A but NB is allele B. Then all 
    # alleles "A" should be "B" for purposes of the simulation and vice versa. 
    meanBB = 0 # arbitrary baseline.
    meanAB = a + d
    meanAA = 2*a
  }
  
  # Simulate phenotypes
  simphenos = data.frame()
  simphenos = data.frame('pgm' = cross.obj$pheno$pgm)
  
  for(i in 1:nsim){
    
    # For each simulation assign genotypes to each individual according to genotype probabilities.
    pick_geno_vals = runif(ninds,0,1)
    
    # Set up phenotypic variance
    if(total.var.constant){
      # Hold total variance constant
      nAA = sum(pick_geno_vals >= (pBB + pAB))
      nBB = sum(pick_geno_vals < pBB)
      nAB = ninds - nAA - nBB
      xgbar_value = (1/(ninds**2))*(nAA * nBB * (2*a)**2 + nAA*nAB*(a+d)**2 + nAB*nBB*(a-d)**2)
      pheno.var = (pheno.sd.2**2 - xgbar_value)
      if(pheno.var >= 0){
        pheno.sd = pheno.var**(0.5)
      }
      if(pheno.var < 0){
        pheno.sd = -1 # just a flag. 
      }
      
    }
    if(!(total.var.constant)){
      # Hold per-genotype variance constant and equal to total variance in the cross object
      pheno.sd = pheno.sd.2
    }
    if(pheno.sd >= 0){
      
      # sim all from AA distrib.
      phenos_sim = rnorm(ninds, mean = meanAA, sd = pheno.sd)
      
      # replace AB and BB genos with AB distrib.
      phenos_sim[pick_geno_vals < (pBB + pAB)] = rnorm(sum(pick_geno_vals < (pBB + pAB)), mean = meanAB, sd = pheno.sd)
      
      # replace BB genos with BB distrib.
      phenos_sim[pick_geno_vals < pBB ] = rnorm(sum(pick_geno_vals < pBB), mean = meanBB, sd = pheno.sd)
      
      # if NA in original cross object, should remain NA.
      phenos_sim[is.na(cross.obj$pheno[[pheno.name.2]])] = NA
      
      # round for counts if required
      if(round){
        phenos_sim = round(phenos_sim, 0)
      }
      
    }
    
    if(pheno.sd < 0){
      phenos_sim = rep(NA, ninds)
      
    }
    
    simphenos[[paste0('sim',i)]] = phenos_sim
    
  }
  
  # Run QTL mapping for each simulated dataset.
  simB = cross.obj
  simB$pheno = simphenos
  sim.peaks.out = data.frame('sim' = 1:nsim, 'sim.a.orig'= a.orig, 'sim.d.orig' = d.orig, 
                             'sim.a' = a, 'sim.d' = d, 
                              'sim.pheno.sd' = pheno.sd.2, 'effect.ratio' = effect.ratio,
                             'trait' = pheno.name.2, 'pos.left' = marker.left, 'pos.right' = marker.right, 
                             'orig.marker.pos' = orig.marker.pos.1,
                             'orig.marker.name' = orig.marker.name.1,
                             'chr'= NA,'pos'= NA,'lod'= NA, # chromosome-wide max lod
                             'lod.left' = NA, 'lod.right' = NA, # lod at left and right flanking markers
                             'pos.bayesint' = NA, 'lod.bayesint' = NA) # max lod in bayes interval from cross 1
  for(i in 1:nsim){
    if(i%%report.interval==0){ # print progress.
      print(paste0('    sim ',i,'(',100*i/nsim,'%)'))
    }
    if(!(all(is.na(simB$pheno[[paste0('sim',i)]])))){
      
      out <- scanone(simB, pheno.col = paste0('sim',i), method = method.2, model = model.2, chr = chr.2)
    
      # save chromosome-wide max lod:
      sim.peaks.out[sim.peaks.out$sim == i,c('pos','lod')] = summary(out)[c('pos','lod')]

      # save lod at the flanking markers:
      outdf <-data.frame(out)
      if(sameobj==FALSE){
        sim.peaks.out[sim.peaks.out$sim == i, 'lod.left'] = outdf[outdf$pos==marker.left,'lod']
        sim.peaks.out[sim.peaks.out$sim == i, 'lod.right'] = outdf[outdf$pos==marker.right,'lod']
      }
      if(sameobj==TRUE){
        if(substr(orig.marker.name.1,1,3)=='loc'){
          orig.marker.name.1 = paste0('c',chr.2,'.',orig.marker.name.1)
        }
        sim.peaks.out[sim.peaks.out$sim == i, 'lod.left'] = outdf[row.names(outdf)==orig.marker.name.1,'lod']
        sim.peaks.out[sim.peaks.out$sim == i, 'lod.right'] = outdf[row.names(outdf)==orig.marker.name.1,'lod']
      }
      
      # save max lod inside the original bayes interval:
      outdf <- subset(outdf, pos >= lower.bayes.1 & pos <= upper.bayes.1)
      names(outdf)<-paste0(names(outdf),'.bayesint')
      sim.peaks.out[sim.peaks.out$sim == i, c('pos.bayesint','lod.bayesint')] = outdf[which.max(outdf$lod.bayesint), c('pos.bayesint','lod.bayesint')] 
    }
  }
  
  # Combine results with previous simulations.
  all.sim.peaks <- rbind(all.sim.peaks, sim.peaks.out)

}

# Next, obtain chromosome-wide, left- and right- flanking marker, and bayes-interval LOD scores
# from the original map. 
orig.map.2 <- scanone(cross.obj, pheno.col = pheno.name.2, chr = chr.2, method = method.2, model = model.2)

# original chromosome-wide max:
all.sim.peaks$orig.cross2.pos.chr <- summary(orig.map.2)$pos
lod.2.chr <- summary(orig.map.2)$lod
all.sim.peaks$orig.cross2.lod.chr <- lod.2.chr

# original LOD at the flanking markers:
outdf <-data.frame(orig.map.2)
if(sameobj==FALSE){
  lod.2.left = outdf[outdf$pos==marker.left,'lod']
  lod.2.right = outdf[outdf$pos==marker.right,'lod']
}
if(sameobj==TRUE){
  if(substr(orig.marker.name.1,1,3)=='loc'){
    orig.marker.name.1 = paste0('c',chr.2,'.',orig.marker.name.1)
  }
  lod.2.left = outdf[row.names(outdf)==orig.marker.name.1,'lod']
  lod.2.right = outdf[row.names(outdf)==orig.marker.name.1,'lod']
}
all.sim.peaks$orig.cross2.lod.left = lod.2.left
all.sim.peaks$orig.cross2.lod.right = lod.2.right

# biggest peak inside the original bayes interval:
outdf <- subset(outdf, pos >= lower.bayes.1 & pos <= upper.bayes.1)
all.sim.peaks$orig.cross2.pos.bayesint = outdf[which.max(outdf$lod),'pos']
lod.2.bayes = outdf[which.max(outdf$lod), 'lod']
all.sim.peaks$orig.cross2.lod.bayesint = lod.2.bayes

# less than:
all.sim.peaks$lod.chr.less.than.orig <- factor(with(all.sim.peaks, lod <= orig.cross2.lod.chr), levels = c(T, F))
all.sim.peaks$lod.left.less.than.orig <- factor(with(all.sim.peaks, lod.left <= orig.cross2.lod.left), levels = c(T, F))
all.sim.peaks$lod.right.less.than.orig <- factor(with(all.sim.peaks, lod.right <= orig.cross2.lod.right), levels = c(T, F))
all.sim.peaks$lod.bayes.less.than.orig <- factor(with(all.sim.peaks, lod.bayesint <= orig.cross2.lod.bayesint), levels = c(T, F))

            
# Summarize.
chrsum <- data.frame(all.sim.peaks
                     %>% group_by(sim.a.orig, lod.chr.less.than.orig) 
                     %>% count() %>% ungroup() 
                     %>% complete(sim.a.orig, lod.chr.less.than.orig, fill = list(n=0)))

leftsum <- data.frame(all.sim.peaks
                     %>% group_by(sim.a.orig, lod.left.less.than.orig) 
                     %>% count() %>% ungroup() 
                     %>% complete(sim.a.orig, lod.left.less.than.orig, fill = list(n=0)))

rightsum <- data.frame(all.sim.peaks
                      %>% group_by(sim.a.orig, lod.right.less.than.orig) 
                      %>% count() %>% ungroup() 
                      %>% complete(sim.a.orig, lod.right.less.than.orig, fill = list(n=0)))

bayessum <- data.frame(all.sim.peaks
                      %>% group_by(sim.a.orig, lod.bayes.less.than.orig) 
                      %>% count() %>% ungroup() 
                      %>% complete(sim.a.orig, lod.bayes.less.than.orig, fill = list(n=0)))

chrsum$lod.compared = 'max.over.chr'
leftsum$lod.compared = 'left.flanking.marker'
rightsum$lod.compared = 'right.flanking.marker'
bayessum$lod.compared = 'max.over.bayesint'

names(chrsum)[2]<-'lod.less.than.orig'
names(leftsum)[2]<-'lod.less.than.orig'
names(rightsum)[2]<-'lod.less.than.orig'
names(bayessum)[2]<-'lod.less.than.orig'

allsum <- rbind(chrsum, leftsum, rightsum, bayessum)
allsum <- subset(allsum, lod.less.than.orig==T)
allsum$prob.less.than.orig <- allsum$n / nsim

d.to.a = d.orig/a.orig # whatever the last value was is fine - all should be the same.

# add other useful info to save
allsum$orig.d = d.to.a * allsum$sim.a.orig
allsum$nsim = nsim
allsum$effect.ratio = effect.ratio
allsum$pheno2 = pheno.name.2
allsum$chr = chr.2
allsum$method = method.2
allsum$model = model.2
allsum$lower.bayes = lower.bayes.1
allsum$upper.bayes = upper.bayes.1
allsum$orig.marker.pos.1 = orig.marker.pos.1
allsum$left.pos = marker.left
allsum$right.pos = marker.right
allsum$pheno.sd.2 = pheno.sd.2
allsum$lod.2.chr = lod.2.chr
allsum$lod.2.bayes = lod.2.bayes
allsum$lod.2.left = lod.2.left
allsum$lod.2.right = lod.2.right
if(total.var.constant){
  allsum$constant_variance = 'total_variance_constant'
}
if(!(total.var.constant)){
  allsum$constant_variance = 'per_genotype_variance_constant'
}


# Save
if(save.all.sims){
  write.csv(all.sim.peaks, file = paste0(save.path,'_all_sims.csv') , row.names = F)
}           
write.csv(allsum, file = paste0(save.path, '_summary_sims.csv'), row.names = F)

# Viz:
# ggplot(all.sim.peaks)+geom_histogram(aes(x=pos.bayesint,fill='max_bayesint'),alpha = 0.2)+geom_histogram(aes(x=pos,fill='max_chr'),alpha = 0.2)+geom_vline(xintercept = c(marker.left,marker.right))+theme_classic()+geom_rug(aes(x=pos),data=out)
# ggplot(all.sim.peaks,aes(x=lod))+geom_point(aes(y=lod.right,color='right'),alpha = 0.2)+geom_point(aes(y=lod.left,color='left'),alpha = 0.2)+geom_point(aes(y=lod.bayesint,color='in_bayes'),alpha = 0.2)+geom_abline(slope=1,intercept=0) + theme_classic()


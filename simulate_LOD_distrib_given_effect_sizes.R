#!/usr/bin/env Rscript

# Author: Emily R. Hager

# Simulations Step 2:

# Use simulations to approximate the probability of obtaining as or more extreme LOD scores 
# than actual data, given an effect size e: 
#         P(LOD < L | effect size == e)


# Given a set of effect sizes, and a marker location, simulate the resulting LOD scores in cross 2. 

# Notation: 'A', 'B': alleles from rqtl object
# 1 vs 2: 1 refers to the cross experiment where the QTL was identified
# 2 refers to the cross experiment where we are testing the detection probability.

require(qtl)
require(argparse)
require(dplyr)
require(tidyr)

source('simulation_functions.R')

# Parse user inputs.
parser <- ArgumentParser()

# Required
parser$add_argument("-o", "--cross_object", help="Path to rqtl cross object (saved as .rds) to use for simulations. Genotype probabilities required; if these do not exist, run calc.genoprob first.", required = T)
parser$add_argument("-p","--phenotype", help = "Phenotype to test (must match name of the cross object pheno column)", required = T)
parser$add_argument("-f", "--folder", help = "Path to folder to put output files", required = T)

# Optional
parser$add_argument("-e", "--effect_ratio", default = 1.0, help = "Ratio of effect size in original cross to effect size in this cross. default %(default)s.")
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
effect.ratio <- args$effect_ratio
geno.reversed <- args$geno_reversed

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
  sameobj = TRUE
}
if(args$identical_cross_object == FALSE){
  sameobj = FALSE
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

round = args$round

save.path <- file.path(save.folder,file.header)

print(paste0('Beginning... Files will be saved with prefix ',save.path))

# Read data
cross.obj <- readRDS(cross.path)
sim.df <- read.csv(sim.path)

# Method, model, nsim, chr, bayes interval should be pulled from sim.df;
# note that this assumes only one value per simulation dataframe. 
model.2 <- as.character(take_first_unique_input(sim.df,'model'))
method.2 <- as.character(take_first_unique_input(sim.df,'method'))
nsim <- as.integer(take_first_unique_input(sim.df,'sims.per.effect.size'))
orig.marker.pos.1 <- take_first_unique_input(sim.df,'orig.peak.marker')
orig.marker.name.1 <- take_first_unique_input(sim.df,'orig.peak.marker.name')
chr.2 <- take_first_unique_input(sim.df,'chromosome')
lower.bayes.1 <- take_first_unique_input(sim.df,'lower.bayes')
upper.bayes.1 <- take_first_unique_input(sim.df,'upper.bayes')
var.constant <- take_first_unique_input(sim.df,'constant_variance')

if(var.constant == 'total_variance_constant'){
  total.var.constant = TRUE
}
if(var.constant == 'per_genotype_variance_constant'){
  total.var.constant = FALSE
}

n.effect <- length(sim.df$sim.a)

# Data from cross object 2. 
pheno.sd.2 = sd(cross.obj$pheno[[pheno.name.2]],na.rm=T)
ninds = dim(cross.obj$pheno)[[1]]


if(sameobj == FALSE){
  # Find the two flanking markers in the new cross object.
  # (Or, allow for only one if both crosses happen to have the same marker!)
  # Figure out how far they are from the original marker.
  map.2 <- cross.obj$geno[[chr.2]]$map
  
  # Check for the possibility that the original marker is further to the end
  # than any marker in cross 2:
  # Requires both maps be in the same units (e.g. pos in bp)
  n.right <- sum(map.2 >= orig.marker.pos.1)
  n.left <- sum(map.2 <= orig.marker.pos.1)
  if(n.right>0){
    marker.pos.right <- min(map.2[map.2>=orig.marker.pos.1])
  }
  if(n.left>0){
    marker.pos.left <- max(map.2[map.2<=orig.marker.pos.1])
  }
  
  if(n.right==0){
    # original marker is right of all cross 2 markers
    marker.pos.right = marker.pos.left
  }
  if(n.left==0){
    # original marker is left of all cross 2 markers
    marker.pos.left = marker.pos.right
  }

  marker.dists <- c(orig.marker.pos.1 - marker.pos.left, marker.pos.right - orig.marker.pos.1)
  if(marker.pos.right == marker.pos.left){
    prop.dists <- c(0.5,0.5)
  }
  if(marker.pos.right != marker.pos.left){
    # note that this is low if close, high if far. 
    prop.dists <- marker.dists/sum(marker.dists)
  }
  
  # Get labels for each position
  marker.name.right <- names(map.2[map.2==marker.pos.right])
  marker.name.left <- names(map.2[map.2==marker.pos.left])
  
  # Get genotype probs at these markers, take a weighted average by phyiscal distance.
  pAA_L = get_genotype_probabilities_at_marker(cross.obj, chr.2, marker.name.left,'AA')
  pAB_L = get_genotype_probabilities_at_marker(cross.obj, chr.2, marker.name.left,'AB')
  pBB_L = get_genotype_probabilities_at_marker(cross.obj, chr.2, marker.name.left,'BB')
  
  pAA_R = get_genotype_probabilities_at_marker(cross.obj, chr.2, marker.name.right,'AA')
  pAB_R = get_genotype_probabilities_at_marker(cross.obj, chr.2, marker.name.right,'AB')
  pBB_R = get_genotype_probabilities_at_marker(cross.obj, chr.2, marker.name.right,'BB')  

  # take weighted average; if prop_dist_left is 0.1, weight_left is 1-0.1 = 0.9. 
  weight_L = 1 - prop.dists[1]
  weight_R = 1 - prop.dists[2]
  pAA = weight_L * pAA_L + weight_R * pAA_R
  pAB = weight_L * pAB_L + weight_R * pAB_R
  pBB = weight_L * pBB_L + weight_R * pBB_R
}

if(sameobj == TRUE){
  # Pull probabilities at the original peak marker.
  pAA = get_genotype_probabilities_at_marker(cross.obj, chr.2, orig.marker.name.1, 'AA')
  pAB = get_genotype_probabilities_at_marker(cross.obj, chr.2, orig.marker.name.1, 'AB')
  pBB = get_genotype_probabilities_at_marker(cross.obj, chr.2, orig.marker.name.1, 'BB')
  
  marker.name.left = orig.marker.name.1
  marker.name.right = orig.marker.name.1
  marker.pos.left = orig.marker.pos.1
  marker.pos.right = orig.marker.pos.1
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
    # allele, and in experiment 1 forest ancestry is allele A but in experiment 2, 
    # forest ancestry is allele B. Then all 
    # alleles "A" should be "B" for purposes of the simulation and vice versa. 
    meanBB = 0 # arbitrary baseline.
    meanAB = a + d
    meanAA = 2*a
  }
  
  # Simulate phenotypes
  simphenos = data.frame('pgm' = cross.obj$pheno$pgm)
  
  for(i in 1:nsim){
    
    # For each simulation assign genotypes to each individual according to genotype probabilities.
    pick_geno_vals = runif(ninds,0,1)
    
   
    if(total.var.constant){
      # Hold total variance constant; thus, calculate per-geno variance
      # based on effect size a and d, and the number of individuals in each genotype class.
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
        # if the simulated effect size is too large to exist while maintaining the 
        # total variance equal to that from the original experiment,
        # the variance would be negative here, so we watch for it. 
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
                             'trait' = pheno.name.2, 
                             'marker.pos.left' = marker.pos.left,  'marker.name.left' = marker.name.left,
                             'marker.pos.right' = marker.pos.right,  'marker.name.right' = marker.name.right,
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
        sim.peaks.out[sim.peaks.out$sim == i, 'lod.left'] = outdf[outdf$pos==marker.pos.left,'lod']
        sim.peaks.out[sim.peaks.out$sim == i, 'lod.right'] = outdf[outdf$pos==marker.pos.right,'lod']
      }
      if(sameobj==TRUE){
        sim.peaks.out[sim.peaks.out$sim == i, 'lod.left'] = outdf[outdf$pos==orig.marker.pos.1,'lod']
        sim.peaks.out[sim.peaks.out$sim == i, 'lod.right'] = outdf[outdf$pos==orig.marker.pos.1,'lod']
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
  lod.2.left = outdf[outdf$pos==marker.pos.left,'lod']
  lod.2.right = outdf[outdf$pos==marker.pos.right,'lod']
}
if(sameobj==TRUE){
  lod.2.left = outdf[outdf$pos==orig.marker.pos.1,'lod']
  lod.2.right = outdf[outdf$pos==orig.marker.pos.1,'lod']
}
all.sim.peaks$orig.cross2.lod.left = lod.2.left
all.sim.peaks$orig.cross2.lod.right = lod.2.right

# biggest peak inside the original bayes interval:
outdf <- subset(outdf, pos >= lower.bayes.1 & pos <= upper.bayes.1)
all.sim.peaks$orig.cross2.pos.bayesint = outdf[which.max(outdf$lod),'pos']
lod.2.bayes = outdf[which.max(outdf$lod), 'lod']
all.sim.peaks$orig.cross2.lod.bayesint = lod.2.bayes

# which sims are less than the threshold value, for each condition:
all.sim.peaks$lod.chr.less.than.orig <- factor(with(all.sim.peaks, lod <= orig.cross2.lod.chr), levels = c(T, F))
all.sim.peaks$lod.left.less.than.orig <- factor(with(all.sim.peaks, lod.left <= orig.cross2.lod.left), levels = c(T, F))
all.sim.peaks$lod.right.less.than.orig <- factor(with(all.sim.peaks, lod.right <= orig.cross2.lod.right), levels = c(T, F))
all.sim.peaks$lod.bayes.less.than.orig <- factor(with(all.sim.peaks, lod.bayesint <= orig.cross2.lod.bayesint), levels = c(T, F))

# Summarize.
chrsum <- count_sims_below_lod(all.sim.peaks, 'max.over.chr')
leftsum <- count_sims_below_lod(all.sim.peaks, 'left.flanking.marker')
rightsum <- count_sims_below_lod(all.sim.peaks, 'right.flanking.marker')
bayessum <- count_sims_below_lod(all.sim.peaks, 'max.over.bayesint')

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
allsum$left.pos = marker.pos.left
allsum$right.pos = marker.pos.right
allsum$name.pos = marker.name.left
allsum$name.pos = marker.name.right
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


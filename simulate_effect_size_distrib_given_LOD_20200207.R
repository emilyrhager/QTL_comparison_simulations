#!/usr/bin/env Rscript



# Use simulations to approximate the probability distribution of 
# effect sizes that could lead to a given detected LOD score: P(E|LOD.1).

# Inputs:
# path to cross object
# chromosome of interest (at the moment, assumes a single peak per chromosome)
# trait of interest
# mapping method (default: ehk) and model (default: normal); note that the phenotypes are simulated as random normal draws so if not normal you might have a problem. 
# number of simulations to run per effect size
# range of effect sizes to run, as a function of the effect size of the QTL
# number of effect sizes to run (together these determine the resolution);  if you want to include 1 this should be odd.
# dominance assumption (options: additive-only, set ratio manually, take ratio from the map)
# window around lod score to include as 'hits' when calculating prob distribution.
# save folder: folder to save output. 

require(qtl)
require(argparse)
require(dplyr)
require(stringr)

# Parse user inputs.
parser <- ArgumentParser()

# Required
parser$add_argument("-o", "--cross_object", help="Path to rqtl cross object (saved as .rds) to use for simulations", required = T)
parser$add_argument("-c", "--chromosome", help="Chromosome with the locus of interest", required = T)
parser$add_argument("-p","--phenotype", help = "Phenotype to test (must match name of the cross object pheno column", required = T)
parser$add_argument("-f", "--folder", help = "Path to folder to put output files", required = T)

# Optional
parser$add_argument("-n", "--nsim", default = 100, help = "Number of simulations per effect size; default %(default)s.")
parser$add_argument("-e", "--neffect", default = 5, help = "Number of effect sizes to test; default %(default)s.")
parser$add_argument("--min_effect", default = 0.5, help = "Smallest effect size to test, as a function of the detected QTL effect size (i.e., --min.effect 0.5 will set the minimum effect size to 0.5 times the detected effect; default %(default)s")
parser$add_argument("--max_effect", default = 2, help = "Largest effect size to test, as a function of the detected QTL effect size (i.e., --max.effect 2 will set the minimum effect size to 2 times the detected effect; default %(default)s")
parser$add_argument("-d", "--dom_model", default = 'additive', help = "Dominance model to use for sims. Options are: (1) (default) 'additive': sims use additive effects only; (2) 'experiment': sims use the d:a ratio detected in the original map for all effect sizes; or (3) a number: the number specified here will be taken as the ratio of d to a.")
parser$add_argument("-w", "--lod_window", default = 1, help = "The half-width of the window around the detected lod score that will be used to generate the probability distribution. For example, -w 1 (the default) gives the probability distribution of effect sizes that would generate a lod score within +/- 1 of the detected score")
parser$add_argument("-m", "--method", default = 'ehk', help = "Method for QTL mapping; default ehk")
parser$add_argument("--model", default = 'normal', help = "Model to use for QTL mapping; default normal. CAREFUL: data for sims generated from random normal distribution regardless, so modify with care.")
parser$add_argument("--file_header", help = "First part of the output file names. Default is: 'sims_{phenotype}_dom_{dominance model}_window_{lod window}_'")
parser$add_argument("--lose_sims", action = "store_true", help = "Do not save the output of all simulations, just save the resulting probability distribution file. Default is to save all simulations but this file can be large.")
parser$add_argument("--round", action = "store_true", help = "Round simulated data to the nearest integer. This is for simulating count data as a rounded normal distribution - use with care!")
parser$add_argument("--total_variance_constant", "-v", action = "store_true", help = "Hold total variance constant; thus the per-genotype variance shrinks with increasing effect size. Default is to hold per-genotype variance constant.")

args <- parser$parse_args()





# Required
cross.path <- args$cross_object #'/n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/data/SWxBK_rqtl_all_object_20191118.RDS'
chr.1 <- as.numeric(as.character(args$chromosome))

# special case: X
if(is.na(chr.1)){
  print('Sims need to be adjusted for the X chromosome.')
  chr.1 = as.character(args$chromosome)
  if(!chr.1=='X'){
    print(paste0('Specified chromosome, ',args$chromosome,' is neither numeric nor X - are you sure?'))
  }
}


pheno.name.1 <- args$phenotype
save.folder <- args$folder

nsim <- as.integer(as.character(args$nsim))

n.effect <- as.integer(as.character(args$neffect))

min.effect <- as.numeric(as.character(args$min_effect))

max.effect <- as.numeric(as.character(args$max_effect))
dom.model <- as.character(args$dom_model)
if(!(dom.model %in% c('additive','experiment'))){
  dom.model <- as.numeric(dom.model)
}
if(is.na(dom.model)){
  print(paste0("Error: dom_model, ",args$dom_model," is not allowed. It must be 'additive','experiment', or numeric. Reverting to dom_model = 'additive" ))
  dom.model <- 'additive'
  # If don't revert, could result in very long runs that fail silently by being full of NAs.
  }

lod.window <- as.numeric(as.character(args$lod_window))

method.1 <- args$method

model.1 <- args$model

save.all.sims <- TRUE
if(args$lose_sims == TRUE){
  save.all.sims <- FALSE
}
total.var.constant = FALSE
var_tag = '_pergeno_var_const'
if(args$total_variance_constant == TRUE){
  total.var.constant = TRUE
  var_tag = '_total_var_const'
}

file.header <- paste0('sims_',pheno.name.1,'_dom_',dom.model,'_window_',lod.window,var_tag)
if(!is.null(args$file_header)){
  file.header <- args$file_header
}

report.interval = ceiling(nsim/5) # report progress 5x per set of sims.
save.path <- file.path(save.folder, file.header)

round = FALSE
if(args$round){
  round = TRUE
}

print('Starting up... Files will be saved in with prefix:')
print(save.path)

##### BEGIN RUNNING THE FUNCTIONS #####

# Set of effect sizes to test
effect.list <- seq( min.effect, max.effect, (max.effect - min.effect)/n.effect )

# Load cross object
cross.obj <- readRDS(cross.path)

# Generate the original map, and report approximate expected run time. 
start.time <- Sys.time()
orig.map.1  = scanone(cross.obj, pheno.col = pheno.name.1, chr = chr.1, method = method.1, model = model.1)
end.time <- Sys.time()
print(paste0('Performing ',nsim,' simulations for each of ',n.effect,' effect sizes. Each simulation will run approx. ', round(end.time-start.time,2),
      ' seconds, thus total expected run time is about ', round(n.effect * nsim * (end.time-start.time)/60, 0),' minutes.'))

# Find original LOD score, peak marker position, number of animals, overall phenotype standard deviation, genotype probabilities.
# Note the assumptions used here: 
#    (1) Overall sd in the phenotype will be used as the per-genotype sd
#    (2) Only handles one QTL peak per chromosome at this point. 
#    (3) Cross object must have genotype probabilities already calculated (use MSG or rqtl's calc.genoprob.)
orig.peaks.1 <- summary(orig.map.1)
LOD.1 = orig.peaks.1$lod
marker.pos.1 = orig.peaks.1$pos
marker.name.1 = row.names(orig.peaks.1)
# check for 'loc':
if(!is.na(str_split(marker.name.1,paste0('c',chr.1,'.'))[[1]][2])){
  marker.name.1 = str_split(marker.name.1,paste0('c',chr.1,'.'))[[1]][2]
}
pheno.sd.1 = sd(cross.obj$pheno[[pheno.name.1]],na.rm=T)
ninds = dim(cross.obj$pheno)[[1]]
BI = bayesint(orig.map.1)
lower.bayes.1 = min(BI$pos)
upper.bayes.1 = max(BI$pos)

# Get separator from cross object.
# Assumes names are in the format chr[number or 'X'][something that's not a digit][more digits]
# and extracts the [something that's not a digit] from the first marker name on the first chromosome.

pAA = cross.obj$geno[[chr.1]]$prob[ ,names(cross.obj$geno[[chr.1]]$prob[1,,'AA'])==marker.name.1,'AA']
pAB = cross.obj$geno[[chr.1]]$prob[ ,names(cross.obj$geno[[chr.1]]$prob[1,,'AB'])==marker.name.1,'AB']
pBB = cross.obj$geno[[chr.1]]$prob[ ,names(cross.obj$geno[[chr.1]]$prob[1,,'BB'])==marker.name.1,'BB']

# Get dominance ratio
mAA = mean(cross.obj$pheno[[pheno.name.1]][pAA>=0.5],na.rm=T)
mAB = mean(cross.obj$pheno[[pheno.name.1]][pAB>=0.5],na.rm=T)
mBB = mean(cross.obj$pheno[[pheno.name.1]][pBB>=0.5],na.rm=T)

orig.a = (mBB - mAA)/2
orig.d = mAB - (mAA + mBB)/2

if(dom.model=='additive'){
  d.to.a = 0
}
if(dom.model=='experiment'){
  # use the ratio from orig.map object
  d.to.a = orig.d/orig.a
}
if(dom.model!='additive' & dom.model!='experiment'){
  # User should have used d:a ratio as the argument
  d.to.a = dom.model
}


# Run simulations
all.sim.peaks <- data.frame()
for(afactor in effect.list){
  a = afactor * orig.a
  print(paste0('Running ',nsim,' simulations for effect size ',round(a,2)))
  
  # Calculate mean phenotypes for each genotype class
  d = d.to.a * a
  
  meanAA = 0 # arbitrary baseline.
  meanAB = a + d
  meanBB = 2*a
  
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
      pheno.var = (pheno.sd.1**2 - xgbar_value)
      if(pheno.var >= 0){
        pheno.sd = pheno.var**(0.5)
      }
      if(pheno.var < 0){
        pheno.sd = -1 # just a flag. 
      }
      
    }
    if(!(total.var.constant)){
      # Hold per-genotype variance constant and equal to total variance in the cross object
      pheno.sd = pheno.sd.1
    }
    if(pheno.sd >= 0){
      
      # sim all from AA distrib.
      phenos_sim = rnorm(ninds, mean = meanAA, sd = pheno.sd)
      
      # replace AB and BB genos with AB distrib.
      phenos_sim[pick_geno_vals < (pBB + pAB)] = rnorm(sum(pick_geno_vals < (pBB + pAB)), mean = meanAB, sd = pheno.sd)
      
      # replace BB genos with BB distrib.
      phenos_sim[pick_geno_vals < pBB ] = rnorm(sum(pick_geno_vals < pBB), mean = meanBB, sd = pheno.sd)
      
      # if NA in original cross object, should remain NA.
      phenos_sim[is.na(cross.obj$pheno[[pheno.name.1]])] = NA
      
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
  simA = cross.obj
  simA$pheno = simphenos
  sim.peaks.out = data.frame('sim' = 1:nsim, 'sim.a'= a, 'sim.d' = d, 'sim.pheno.sd' = pheno.sd.1,
                             'chr'= NA,'pos'= NA,'lod'= NA,
                             'trait' = pheno.name.1, 'orig.peak.marker' = marker.pos.1, 'orig.peak.marker.name' = marker.name.1)
  for(i in 1:nsim){
    if(i%%report.interval==0){
      print(paste0('    sim ',i, ' (',100*round(i/nsim,2), '%)'))
    }
    if(!(all(is.na(simA$pheno[[paste0('sim',i)]])))){
      out <- scanone(simA, pheno.col = paste0('sim',i), method = method.1, model = model.1, chr=chr.1)
      sim.peaks.out[sim.peaks.out$sim == i,c('pos','lod')] = summary(out)[c('pos','lod')]
    }
  }
  
  all.sim.peaks <- rbind(all.sim.peaks, sim.peaks.out)
}

all.sim.peaks$chr <- chr.1

if(save.all.sims){
  write.csv(all.sim.peaks,file = paste0(save.path,'_all_sims.csv'), row.names = F)
}

sims.in.window = droplevels(subset(all.sim.peaks,(lod <= LOD.1 + lod.window) & (lod >= LOD.1 - lod.window)))
counts.in.window <- data.frame(count(sims.in.window,sim.a,sim.d))
counts.in.window$prob <- counts.in.window$n/sum(counts.in.window$n)
counts.in.window$trait <- pheno.name.1
counts.in.window$sims.per.effect.size <- nsim
counts.in.window$min.sim.effect.a <- min.effect * orig.a
counts.in.window$max.sim.effect.a <- max.effect * orig.a
counts.in.window$orig.effect.a <- orig.a
counts.in.window$orig.effect.d <- orig.d
counts.in.window$dom.model <- dom.model
counts.in.window$cross.path <- cross.path
counts.in.window$orig.peak.marker <- marker.pos.1
counts.in.window$orig.peak.marker.name <- marker.name.1
counts.in.window$model <- model.1
counts.in.window$method <- method.1
counts.in.window$chromosome <- chr.1
counts.in.window$lower.bayes <- lower.bayes.1
counts.in.window$upper.bayes <- upper.bayes.1
counts.in.window$lod.1 <- LOD.1
counts.in.window$lod.window <- lod.window
if(total.var.constant){
  counts.in.window$constant_variance = 'total_variance_constant'
}
if(!(total.var.constant)){
  counts.in.window$constant_variance = 'per_genotype_variance_constant'
}

write.csv(counts.in.window,file = paste0(save.path,'_prob_distrib.csv'), row.names=F)

# flag possible issues with sims:
if(sum(counts.in.window$sim.a == counts.in.window$min.sim.effect.a)>0 & !min.effect==0){
  print('Possible issue: minimum simulated effect size is in the final distribution, but it is not zero. The simulated LOD distribution may be artificially truncated; consider expanding the range of simulated effect sizes downwards.')
}
if(sum(counts.in.window$sim.a == counts.in.window$max.sim.effect.a)>0){
  print('Possible issue: maximum simulated effect size is in the final distribution. The simulated LOD distribution may be artificially truncated; consider expanding the range of simulated effect sizes upwards.')
}


# Note: to visualize the output:
# ggplot(all.sim.peaks,aes(x=lod,fill=sim.a))+geom_histogram() + theme_classic() + geom_vline(xintercept = LOD.1) + geom_rect(xmin = LOD.1-lod.window, xmax = LOD.1 + lod.window, ymin = -Inf, ymax = Inf, alpha = 0.01, fill = 'black')+facet_grid(sim.a~.)
# ggplot(counts.in.window,aes(x=sim.a,y=prob))+geom_bar(stat='identity')+theme_classic()+geom_vline(aes(xintercept = orig.effect.a))

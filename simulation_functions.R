### Functions for QTL simulations.


get_genotype_probabilities_at_marker<- function(cross.obj, chr.name, marker.name, geno.class){
  
  " Obtain genotype probabilities for a given genotype class (AA, AB, BB) at a given marker "
  
  # Genotype probabilities for all classes, for the chromosome of interest: 
  # dimension n_individuals x n_markers x 3
  chromosome_geno_probs = cross.obj$geno[[chr.name]]$prob
  
  # Marker names in order: 
  marker_names = names(chromosome_geno_probs[1,,geno.class])
  
  # Genotype probabilities for this genotype class, for each individual:
  probs_at_marker = chromosome_geno_probs[ ,marker_names==marker.name,geno.class]
  
  return(probs_at_marker)
}

take_first_unique_input <- function(df, col){
  
  " Return the first unique value in column col of dataframe df. If more than one value present, print a warning."
  
  vals = unique(df[[col]])
  if(length(vals)>1){
    print(paste0('Problem: dataframe column ',col,' has multiple values; expected only one. Taking first value.'))
  }
  return(vals[[1]])
}


count_sims_below_lod <- function(df, grp){
  
  " Count the number of simulations that meet the given criteria "
  
  colnames <- structure(c('lod.chr.less.than.orig','lod.left.less.than.orig','lod.right.less.than.orig','lod.bayes.less.than.orig'),
                        names=c('max.over.chr','left.flanking.marker','right.flanking.marker','max.over.bayesint'))
  
  labelcol <- colnames[[grp]]
  
  dfsum <- data.frame(df
                      %>% group_by(sim.a.orig, !!rlang::sym(labelcol))
                      %>% count()
                      %>% ungroup()
                      %>% complete(sim.a.orig, !!rlang::sym(labelcol), fill = list(n=0)))
  
  dfsum$lod.compared = labelcol
  
  names(dfsum)[2]<-'lod.less.than.orig'
  
  return(dfsum)     
  
}

merge_step1_step2_output <- function(step1df, step2df){
  
  " Merge the results of step 1 and step 2 simulations and calculate the final probability value "
  
  cdf = merge(step1df,step2df, by.x = c('sim.a','constant_variance','model','method'),
              by.y=c('sim.a.orig','constant_variance','model','method'), all = T)
  cdf$pcomb = with(cdf,prob*prob.less.than.orig)
  if(dim(cdf)[1] != dim(step2df)[1]){
    print(paste0('Dimensions do not match: dim(cdf) = ',dim(cdf),' dim step2e = ',dim(step2e)))
  }
  sdf <- data.frame( cdf[c('pcomb','constant_variance','model','method',
                           'trait','dom.model','orig.effect.a','orig.effect.d',
                           'orig.peak.marker','orig.peak.marker.name',
                           'chromosome','lod.compared','effect.ratio','left.pos',
                           'right.pos','lod.2.chr','lod.2.bayes','lod.2.left',
                           'lod.2.right')] %>% group_by(constant_variance,model,method,
                                                        trait,dom.model,orig.effect.a,orig.effect.d,
                                                        orig.peak.marker,orig.peak.marker.name,
                                                        chromosome,lod.compared,effect.ratio,left.pos,
                                                        right.pos,lod.2.chr,lod.2.bayes,lod.2.left,
                                                        lod.2.right) %>% summarize(p.final = sum(pcomb)))
  return(sdf)
}


# Note: to visualize the output:

# ggplot(all.sim.peaks)+geom_histogram(aes(x=pos.bayesint,fill='max_bayesint'),alpha = 0.2)+geom_histogram(aes(x=pos,fill='max_chr'),alpha = 0.2)+geom_vline(xintercept = c(marker.left,marker.right))+theme_classic()+geom_rug(aes(x=pos),data=out)
# ggplot(all.sim.peaks,aes(x=lod))+geom_point(aes(y=lod.right,color='right'),alpha = 0.2)+geom_point(aes(y=lod.left,color='left'),alpha = 0.2)+geom_point(aes(y=lod.bayesint,color='in_bayes'),alpha = 0.2)+geom_abline(slope=1,intercept=0) + theme_classic()
# ggplot(all.sim.peaks,aes(x=lod,fill=sim.a))+geom_histogram() + theme_classic() + geom_vline(xintercept = LOD.1) + geom_rect(xmin = LOD.1-lod.window, xmax = LOD.1 + lod.window, ymin = -Inf, ymax = Inf, alpha = 0.01, fill = 'black')+facet_grid(sim.a~.)
# ggplot(output.data,aes(x=sim.a,y=prob))+geom_bar(stat='identity')+theme_classic()+geom_vline(aes(xintercept = orig.effect.a))

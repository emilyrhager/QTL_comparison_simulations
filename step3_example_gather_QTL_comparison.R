#  Example script to gather and plot QTL simulation output for a specific project.
#  The key function is merge_step1_step2_output, in simulation_functions.R.
#    This is a more detailed example gathering and plotting the results for a specific project 
#    that involved comparing multiple settings for dominance and per-genotype variance, and five effect size ratios.
#
#  Author: Emily R Hager, 2020

library(dplyr)
library(ggplot2)
library(stringr)

source('simulation_functions.R')


# Top-level folder containing simulation output
simfold = '/Volumes/hoekstra_lab/Users/hager/QTL/compare_QTL_sims_JTG/'

# Project-specific function to gather simulation output for comparing QTL for 
# two traits (QTL detected for trait w1 on chromosome chr; testing if shared with trait w2)
# const: whether total or per-genotype variance was held constant
# dom: whether locus was assumed to be additive
# w1, w2: the two phenotypes being compared

gather_sim_data <- function(chr, w1, w2, const, dom, simfold){
  " Gather simulation output for a given chromosome/trait across all values of effect size. "
  
  # Read step 1 output (same for all effect ratios)
  f = paste0(simfold,'sims_chr_',chr,'_',w1,'_',w2)
  step1df = read.csv(paste0(f,'/sims_',w1,'_res_dom_',dom,'_window_0.1_',const,'_var_const_prob_distrib.csv'))
  
  outdf = data.frame()
  
  # elist = the effect size ratios that were tested
  if(dom=='experiment' & const == 'total'){
    elist = c(0.1,0.25,0.5,0.75,1)
  }
  if(!(dom=='experiment' & const =='total')){
    elist = 1
  }
  
  for(e in elist){ 
    
    # Read step 2 output for this effect ratio
    # and combine with step 1 output to calculate final probability
    step2e = read.csv(paste0(f,'/sims_',w2,'_res_dom_',dom,'_',const,'_var_constant_effect_ratio_',
                             e,'_summary_sims.csv'))
    sdf = merge_step1_step2_output(step1df, step2e) # function from simulation_functions.R
    sdf$whisker1 = w1
    sdf$whisker2 = w2
    outdf = rbind(outdf,sdf)
  }
  return(outdf)
}

# Gather data for each set of simulation conditions
wlist = c('alpha','beta','gamma','delta')
resdf = data.frame()
for(s1 in Sys.glob(paste0(simfold,'*step1*'))){
  # Pull the chromosome and traits from the name of the output folder.
  bn = str_split(basename(s1),'_')[[1]]
  chr = bn[4]
  w1 = bn[5]
  for(const in c('total','pergeno')){
    for(dom in c('additive','experiment')){
      for(w2 in wlist[!wlist==w1]){
        print(c(chr,w1,w2,const,dom))
        resdf = rbind(resdf, gather_sim_data(chr,w1,w2,const,dom,simfold))
      }
    }
  }
}

# Save output
write.csv(resdf,'../compare_QTL_sims_JTG/data/summary_all_sim_results_20200424.csv',row.names=F)

# The primary results are from total variance constant/ dom model experiment. Here's for effect size 1:
Rte.bayes = subset(resdf, constant_variance=='total_variance_constant' & dom.model=='experiment' & lod.compared == 'max.over.bayesint' & effect.ratio==1)
ggplot(Rte.bayes,aes(x=factor(chromosome),y=p.final,color=whisker1,shape=whisker2))+
  geom_hline(yintercept = c(0.01,0.05),color='lightgray')+
  geom_point()+
  theme_classic()+
  ylim(0,0.8)

# If we look at all four types of simulation the results are nearly identical. 
ggplot(subset(resdf, lod.compared=='max.over.bayesint' & effect.ratio==1),aes(x=factor(chromosome),y=p.final,color=whisker1,shape=whisker2))+geom_point()+
  theme_classic()+geom_hline(yintercept = c(0.01,0.05),color='lightgray')+
  facet_grid(constant_variance ~ dom.model)

# ... but just to check, not quite identical:
ggplot(subset(resdf, lod.compared=='max.over.bayesint' & effect.ratio==1),aes(x=factor(chromosome),y=p.final,color=whisker1,shape=whisker2))+geom_point()+
  theme_classic()+geom_hline(yintercept = c(0.01,0.05),color='lightgray')

# How does this vary by effect ratio?
ggplot(subset(resdf, constant_variance == 'total_variance_constant' & dom.model == 'experiment' & lod.compared == 'max.over.bayesint'))+
  theme_classic()+
  geom_point(aes(x = effect.ratio, y = p.final, color = whisker2, shape=whisker1))+
  geom_line(aes(x = effect.ratio, y = p.final, color = whisker2))+
  facet_wrap(~chromosome)+scale_x_reverse()+
  geom_hline(yintercept = c(0.01,0.05,0.1,0.2),color='lightgray')+
  ggtitle('max over original bayes confidence interval')

ggplot(subset(resdf, constant_variance == 'total_variance_constant' & dom.model == 'experiment' & lod.compared == 'right.flanking.marker'))+
  theme_classic()+
  geom_point(aes(x = effect.ratio, y = p.final, color = whisker1, shape = whisker2))+
  geom_line(aes(x = effect.ratio, y = p.final, color = whisker1, linetype = whisker2))+
  facet_wrap(~chromosome)+scale_x_reverse()+
  geom_hline(yintercept = c(0.01,0.05),color='lightgray')+
  ggtitle('exact peak marker')

ggplot(subset(resdf, constant_variance == 'total_variance_constant' & dom.model == 'experiment' & lod.compared == 'max.over.chr'))+
  theme_classic()+
  geom_point(aes(x = effect.ratio, y = p.final, color = whisker1, shape = whisker2))+
  geom_line(aes(x = effect.ratio, y = p.final, color = whisker1, linetype = whisker2))+
  facet_wrap(~chromosome)+scale_x_reverse()+
  geom_hline(yintercept = c(0.01,0.05),color='lightgray')+
  ggtitle('max over whole chromosome')

ggplot(subset(resdf, constant_variance == 'total_variance_constant' & dom.model == 'experiment' & lod.compared == 'max.over.bayesint'))+
  theme_classic()+
  geom_point(aes(x = effect.ratio, y = p.final, color = whisker1))+
  geom_line(aes(x = effect.ratio, y = p.final, color = whisker1))+
  facet_grid(whisker2~chromosome)+scale_x_reverse()+
  geom_hline(yintercept = c(0.01,0.05,0.1,0.2),color='lightgray')+
  ggtitle('max over original bayes confidence interval')

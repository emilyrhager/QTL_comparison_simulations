# Jake QTL comparison sims
library(dplyr)
library(ggplot2)
library(stringr)

chr = 8
w1 = 'alpha'
w2 = 'beta'
const = 'total' # pergeno
dom = 'experiment' # additive


simfold = '/Volumes/hoekstra_lab/Users/hager/QTL/compare_QTL_sims_JTG/'


gather_sim_data <- function(chr, w1, w2, const, dom, simfold){
  f = paste0(simfold,'sims_chr_',chr,'_',w1,'_',w2)
  step1df = read.csv(paste0(f,'/sims_',w1,'_res_dom_',dom,'_window_0.1_',const,'_var_const_prob_distrib.csv'))
  
  outdf = data.frame()
  
  if(dom=='experiment' & const == 'total'){
    elist = c(0.1,0.25,0.5,0.75,1)
  }
  if(!(dom=='experiment' & const =='total')){
    elist = 1
  }
  
  for(e in elist){
    
    step2e = read.csv(paste0(f,'/sims_',w2,'_res_dom_',dom,'_',const,'_var_constant_effect_ratio_',
                             e,'_summary_sims.csv'))
    cdf = merge(step1df,step2e, by.x = c('sim.a','constant_variance','model','method'),
                by.y=c('sim.a.orig','constant_variance','model','method'), all = T)
    cdf$pcomb = with(cdf,prob*prob.less.than.orig)
    if(dim(cdf)[1] != dim(step2e)[1]){
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
    sdf$whisker1 = w1
    sdf$whisker2 = w2
    outdf = rbind(outdf,sdf)
  }
  return(outdf)
}

sdf = gather_sim_data(chr, w1, w2, const, dom, simfold)

wlist = c('alpha','beta','gamma','delta')
resdf = data.frame()
for(s1 in Sys.glob(paste0(simfold,'*step1*'))){
  bn = str_split(basename(s1),'_')[[1]]
  chr = bn[4]
  w1 = bn[5]
  for(const in c('total','pergeno')){
    for(dom in c('additive','experiment')){
      for(w2 in wlist[!wlist==w1]){
        print(c(chr,w1,w2,const,dom))
        resdf = rbind(resdf,gather_sim_data(chr,w1,w2,const,dom,simfold))
      }
    }
  }
}

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


sum(resdf$left.pos!=resdf$orig.peak.marker.name)
sum(resdf$right.pos!=resdf$orig.peak.marker.name)


library(qtl)
JTG <- readRDS('../compare_QTL_sims_JTG/data/Gable_QTL_totlresid_04_2020.rds')
beta <- scanone(JTG, pheno.col = 'beta_res', method = 'ehk', chr =c(3,4,5,8,9,10,15))
alpha <- scanone(JTG, pheno.col = 'alpha_res', method = 'ehk', chr =c(3,4,5,8,9,10,15))
gamma <- scanone(JTG, pheno.col = 'gamma_res', method = 'ehk', chr =c(3,4,5,8,9,10,15))
delta <- scanone(JTG, pheno.col = 'delta_res', method = 'ehk', chr =c(3,4,5,8,9,10,15))
plot(delta,alpha,beta)
plot(gamma,add=T,col='purple')

# 4
plot(delta,alpha,beta,chr=4)
plot(gamma,add=T,col='purple',chr=4)
abline(v=c(133,235))

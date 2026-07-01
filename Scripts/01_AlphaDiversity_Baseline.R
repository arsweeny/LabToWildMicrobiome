#ALPHA DIVERSITY

rm(list=ls())

library(MASS)
library(Matrix)
library(ade4)
library(multtest)
library(vegan)
library(phyloseq)
library(microbiome)
library(pals)
library(microViz)
library(magrittr)
library(ggregplot)
library(ggpubr)
library(tidyverse)
library(ggforce)
library(ape)
library(picante)
library(data.table)
library(DT)
library(ggplot2)
library(patchwork)
require(RColorBrewer)
library(glmmTMB)
library(TMB)

source("Scripts/0_DataPrep_Baseline.R") 
# runs data prep & sets up all the data subsets 

#-------------------------------------------------------------------------------------------------------------------------------------------------
# RAREFACTION - Need to account for different sequencing depths and sequencing plates 
# All samples ps_all (with MF1's first)
# Hewan, Penicuik & Falkirk samples all combined into one 'wild' category
# MF1 & CD1 mice collapsed into Lab 
# ash54 wm colony = lab-like & b17 wm colony = wild-like
ps_all <- ps_all %>% 
  ps_mutate(SampleOrigin=case_when(str_detect(SampleDetailed, "MF1 lab mice") ~ "Lab:Mm",
                                   str_detect(SampleDetailed, "Pen") ~ "Wild:As", 
                                   str_detect(SampleDetailed, "Hew") ~ "Wild:As",
                                   str_detect(SampleDetailed, "Fal") ~ "Wild:As",
                                   str_detect(MouseOrigin, "Ash54 col") ~ "Lab-like:As", 
                                   str_detect(MouseOrigin, "b17 col") ~ "Wild-like:As",
                                   str_detect(MouseOrigin, "charles riv") ~ "Lab:Mm")) %>% 
  ps_mutate(SampleOrigin=fct_relevel(SampleOrigin, c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>% 
  ps_filter(!SampleDetailed %in% c("Penicuik 2020", "Falkirk 2017", "MF1 lab mice")) 

#Rarefy
#subsample reads - this is rarefaction
#(subsampling reads from each sample without replacement to a constant depth)
ps_rare <- phyloseq::rarefy_even_depth(ps_all, rngseed = 123, replace = FALSE)
#367 OTUs were removed because they are no longer present in any sample after random subsampling
head(phyloseq::sample_sums(ps_rare)) 
#sequence depth set to 15,857 

#generate tree so we can do phylogentic analysis too 
random_tree = rtree(ntaxa(ps_rare), rooted=TRUE, tip.label=taxa_names(ps_rare)) 
# #save tree so don't have to repeat this step again
write.tree(random_tree, "Output/Ps_Rare_PhyloTree")

#read in tree saved previously so no need to make a new one and change analysis 
random_tree <- read.tree("Output/Ps_Rare_PhyloTree")
#merge tree with phyloseq object
ps_rare = merge_phyloseq(ps_rare, random_tree)
# otu_table()   OTU Table:         [ 3325 taxa and 278 samples ]
# sample_data() Sample Data:       [ 278 samples by 23 sample variables ]
# tax_table()   Taxonomy Table:    [ 3325 taxa by 6 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 3325 tips and 3323 internal nodes ]

adiv_R <- data.frame("Observed" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
                      "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
                      "InvSimpson" = phyloseq::estimate_richness(ps_rare, measures = "InvSimpson"),
                      "PD" = picante::pd(samp = data.frame((data.frame(phyloseq::otu_table(ps_rare)))), 
                                         tree = phyloseq::phy_tree(ps_rare))[, 1],
                      "SampleOrigin" = phyloseq::sample_data(ps_rare)$SampleOrigin)
#then add these columns to a dataframe containing all the metadata taken from the pyloseq object 
ps_meta_rare <- meta(ps_rare)
ps_meta_rare$Observed <-adiv_R$Observed
ps_meta_rare$FaithPD <- adiv_R$PD
ps_meta_rare$Shannon <- adiv_R$Shannon
ps_meta_rare$InverseSimpson <- adiv_R$InvSimpson 
colnames(ps_meta_rare)


# sanity checks / filter populations --------------------------------------
library(ggregplot)

ps_meta_rare %>% filter(SampleOrigin == "Wild:As") %>% 
  SinaGraph(., "SampleDetailed", "Observed")

ps_meta_rare %>% filter(SampleOrigin == "Lab:Mm") %>% 
  SinaGraph(., "MouseOrigin", "Observed")

### note going to edit up above now to remove these categories 


#calculate means for each alpha metric by environment -------
#Observed richness
alpha_obs_R <- ps_meta_rare %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(Observed), se=plotrix::std.error(Observed), min=min(Observed), max=max(Observed), var=var(Observed), sd=sd(Observed))

#Pylogenetic diversity
alpha_pd_R <- ps_meta_rare %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(FaithPD), se=plotrix::std.error(FaithPD), min=min(FaithPD), max=max(FaithPD), var=var(FaithPD), sd=sd(FaithPD)) 

#shannon
alpha_shan_R <- ps_meta_rare %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(Shannon), se=plotrix::std.error(Shannon), min=min(Shannon), max=max(Shannon), var=var(Shannon), sd=sd(Shannon)) 

#inverse simpson
alpha_sim_R <- ps_meta_rare %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(InverseSimpson), se=plotrix::std.error(InverseSimpson), min=min(InverseSimpson), max=max(InverseSimpson), var=var(InverseSimpson), sd=sd(InverseSimpson)) 


#plot all cohorts in one figure figure to see if any diff
adiv_R %>%
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  gather(key = metric, value = value, c("Observed", "PD", "Shannon","InvSimpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "PD", "Shannon","InvSimpson"))) %>%
  ggplot(aes(x = SampleOrigin, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = SampleOrigin), height = 0.6, width = 0.3, show.legend = F, alpha=1.0, size=2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme_bw(base_size = 20) +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> alpha_boxplot

#Observed richness PLOT 
ps_meta_rare %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=Observed, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=alpha_obs_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=alpha_obs_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Observed Richness", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> obs_all_R

#Phylogenetic Diversity plot
ps_meta_rare %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=FaithPD, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=alpha_pd_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=alpha_pd_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Phylogenetic Diversity", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> pd_all_R

#Shannon Plot
ps_meta_rare %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=Shannon, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=alpha_shan_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=alpha_shan_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Shannon Diversity", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> shan_all_R

#Inverse Simpson Plot
ps_meta_rare %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=InverseSimpson, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=alpha_sim_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=alpha_sim_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Inverse Simpson", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> sim_all_R

#Panel plot of all alpha
alpha_pan <- obs_all_R|pd_all_R|shan_all_R|sim_all_R 

#Statistical analysis

#Observed richness - plot hist first
hist(ps_meta_rare$Observed) #slight skew
hist(log(ps_meta_rare$Observed)) #looks better 
ps_meta_rare$Log_Observed = log(ps_meta_rare$Observed)
#Run GLM with log and non log data first 
obs_glm1 <-glm(Observed ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_rare)
summary(obs_glm1)
#as expected model fit terrible

#re-order factors to set Lab-Mm as intercept and reorder seq plate so that plate 8 is intercept 
ps_meta_rare$SampleOrigin <- factor(ps_meta_rare$SampleOrigin, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))
ps_meta_rare$SeqPlate <- factor(ps_meta_rare$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))
ps_meta_rare$Sex <-as.factor(ps_meta_rare$Sex)

#run as normal glm with log data 
obsl_glm1 <-glm(Log_Observed ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_rare)
summary(obsl_glm1)

#what does seqplate significance mean 
summ <- aggregate(Log_Observed ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_rare, mean)
summ
#plot the model interaction to work out what the interaction estimate is telling you 
ggplot(summ, aes(x=SampleOrigin, y=Log_Observed, colour=SeqPlate, group=SeqPlate)) +
  geom_point() +
  geom_line() +
  theme_bw()

#run mixed model, with seq plate as random effect to account for variation here 
obsl_glmm <-glmmTMB(Log_Observed ~ SampleOrigin + Sex + (1|SeqPlate), data=ps_meta_rare)
summary(obsl_glmm)

#first make a column combining Id & SampleDetailed, otherwise random affect will treat all Ax mice as the same, even if from different rse experiemnts
ps_meta_rare <- tidyr::unite(ps_meta_rare,"MouseId", Id, SampleDetailed, remove = F) 
#try running as glmm with seqplate and mouse id as random
obsl_glmm1 <-glmmTMB(Log_Observed ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data=ps_meta_rare)
summary(obsl_glmm1)
# AIC      BIC   logLik deviance df.resid 
# 117.8    146.9    -50.9    101.8      270 
# Random effects:Conditional model:
#   Groups   Name        Variance Std.Dev.
# SeqPlate (Intercept) 0.003429 0.05856 
# MouseId  (Intercept) 0.002695 0.05191 
# Residual             0.079824 0.28253 
# Number of obs: 278, groups:  SeqPlate, 8; MouseId, 171
# Dispersion estimate for gaussian family (sigma^2): 0.0798 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               5.42349    0.04408  123.04  < 2e-16 ***
#   SampleOriginLab-like:As  -0.26276    0.04812   -5.46 4.75e-08 ***
#   SampleOriginWild-like:As -0.58552    0.04838  -12.10  < 2e-16 ***
#   SampleOriginWild:As      -0.20033    0.07162   -2.80  0.00516 ** 
#   SexM                      0.02464    0.03576    0.69  0.49080    


#phylogenetic diversity 
hist(ps_meta_rare$FaithPD) # looks pretty normal 
hist(log(ps_meta_rare$FaithPD)) # looks worse
ps_meta_rare$Log_FaithPD = log(ps_meta_rare$FaithPD)
#skewed right most values on lower size of x axis - looks more normal when logged 
fpd_glm1 <-glm(FaithPD ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_rare)
summary(fpd_glm1)
#data looked good, but fit is terrible, really overdispersed 

#run glm on logged
fpdl_glm1 <-glm(Log_FaithPD ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_rare)
summary(fpdl_glm1)
  
#run as glmm with seq plate random
fpdl_glmm <-glmmTMB(Log_FaithPD ~ SampleOrigin + Sex + (1|SeqPlate), data=ps_meta_rare)
summary(fpdl_glmm)
   
#run final glmm with seq plate and mpouse id as random 
fpdl_glmm1 <-glmmTMB(Log_FaithPD ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data=ps_meta_rare)
summary(fpdl_glmm1)
# AIC      BIC   logLik deviance df.resid 
# -77.2    -48.2     46.6    -93.2      270 
# Random effects:Conditional model:
#   Groups   Name        Variance Std.Dev.
# SeqPlate (Intercept) 0.001617 0.04022 
# MouseId  (Intercept) 0.001227 0.03502 
# Residual             0.039715 0.19929 
# Number of obs: 278, groups:  SeqPlate, 8; MouseId, 171
# Dispersion estimate for gaussian family (sigma^2): 0.0397 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               6.48667    0.03083  210.39  < 2e-16 ***
#   SampleOriginLab-like:As  -0.16702    0.03384   -4.94 7.96e-07 ***
#   SampleOriginWild-like:As -0.39112    0.03401  -11.50  < 2e-16 ***
#   SampleOriginWild:As      -0.13978    0.04990   -2.80  0.00509 ** 
#   SexM                      0.01673    0.02517    0.66  0.50638  

#shannon diversity - all groups except MF1 sig lower diversity than CD1s
hist(ps_meta_rare$Shannon) # pretty skewed
hist(log(ps_meta_rare$Shannon)) # looks even worse
#run on non-transformed data
shan_glm1 <-glm(Shannon ~ SampleOrigin  + Sex + SeqPlate, data=ps_meta_rare)
summary(shan_glm1)
#run glmm with seq plate
shan_glm <-glmmTMB(Shannon ~ SampleOrigin + Sex + (1|SeqPlate), data=ps_meta_rare)
summary(shan_glm)

#run glmm with seq plate & mouse id 
shan_glm1 <-glmmTMB(Shannon ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data=ps_meta_rare)
summary(shan_glm1)
# AIC      BIC   logLik deviance df.resid 
# 540.3    569.3   -262.2    524.3      270 
# Random effects:Conditional model:
#   Groups   Name        Variance Std.Dev.
# SeqPlate (Intercept) 0.003072 0.05542 
# MouseId  (Intercept) 0.022531 0.15010 
# Residual             0.361412 0.60118 
# Number of obs: 278, groups:  SeqPlate, 8; MouseId, 171
# Dispersion estimate for gaussian family (sigma^2): 0.361 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               3.90904    0.08302   47.08  < 2e-16 ***
#   SampleOriginLab-like:As  -0.33389    0.10364   -3.22  0.00127 ** 
#   SampleOriginWild-like:As -1.25609    0.10476  -11.99  < 2e-16 ***
#   SampleOriginWild:As      -0.93477    0.13063   -7.16 8.33e-13 ***
#   SexM                      0.01871    0.07658    0.24  0.80694 

#inverse simpson - all groups except MF1 sig lower diversity than CD1s
hist(ps_meta_rare$InverseSimpson) # really skewed
hist(log(ps_meta_rare$InverseSimpson)) # better
ps_meta_rare$Log_InverseSimpson = log(ps_meta_rare$InverseSimpson)
#skewed right most values on lower size of x axis) - looks more normal when logged 
sim_glm1 <-glm(InverseSimpson ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_rare)
summary(sim_glm1)
#Not a great fit

#lookes better with log data 
siml_glm1 <-glm(Log_InverseSimpson ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_rare)
summary(siml_glm1)

#glmm with seq plate
siml_glmm <-glmmTMB(Log_InverseSimpson ~ SampleOrigin + Sex + (1|SeqPlate), data=ps_meta_rare)
summary(siml_glmm)

#glmm with seq plate and mouseid
siml_glmm1 <-glmmTMB(Log_InverseSimpson ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data=ps_meta_rare)
summary(siml_glmm1)
# AIC      BIC   logLik deviance df.resid 
# 560.4    589.4   -272.2    544.4      270 
# Random effects:Conditional model:
#   Groups   Name        Variance Std.Dev.
# SeqPlate (Intercept) 0.01571  0.1254  
# MouseId  (Intercept) 0.01268  0.1126  
# Residual             0.39313  0.6270  
# Number of obs: 278, groups:  SeqPlate, 8; MouseId, 171
# Dispersion estimate for gaussian family (sigma^2): 0.393 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               3.04235    0.09568   31.80  < 2e-16 ***
#   SampleOriginLab-like:As  -0.31458    0.10594   -2.97  0.00298 ** 
#   SampleOriginWild-like:As -1.22626    0.10680  -11.48  < 2e-16 ***
#   SampleOriginWild:As      -1.10340    0.14885   -7.41 1.23e-13 ***
#   SexM                      0.03673    0.07831    0.47  0.63906


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#NO RAREFACTION ANALYSIS
# ALPHA DIVERSITY - ANALYSIS

# All samples ps_all (with MF1's first)
# Hewan, Penicuik & Falkirk samples all combined into one 'wild' category
# MF1 & CD1 mice collapsed into Lab 
# ash54 wm colony = lab-like & b17 wm colony = wild-like
ps_all <- ps_all %>% 
  ps_mutate(SampleOrigin=case_when(str_detect(SampleDetailed, "MF1 lab mice") ~ "Lab-Mm",
                                   str_detect(SampleDetailed, "Pen") ~ "Wild-As", 
                                   str_detect(SampleDetailed, "Hew") ~ "Wild-As",
                                   str_detect(SampleDetailed, "Fal") ~ "Wild-As",
                                   str_detect(MouseOrigin, "Ash54 col") ~ "Lab-like-As", 
                                   str_detect(MouseOrigin, "b17 col") ~ "Wild-like-As",
                                   str_detect(MouseOrigin, "charles riv") ~ "Lab-Mm"))


#first generate tree so we can do phylogentic analysis too 
# random_tree = rtree(ntaxa(ps_all), rooted=TRUE, tip.label=taxa_names(ps_all)) 
#save tree so don't have to repeat this step again
# write.tree(random_tree, "DescriptiveFigures/Ps_All_PhyloTree")

#read in tree saved previously so no need to make a new one and change analysis 
random_tree <- read.tree("DescriptiveFigures/Ps_All_PhyloTree")


ps_all = merge_phyloseq(ps_all, random_tree)
# otu_table()   OTU Table:         [ 3691 taxa and 278 samples ]
# sample_data() Sample Data:       [ 278 samples by 19 sample variables ]
# tax_table()   Taxonomy Table:    [ 3691 taxa by 6 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 3691 tips and 3690 internal nodes ]
ps_all@phy_tree #check worked properly with root included etc

#make some tables 
ps_all_asvtab <- as.data.frame(ps_all@otu_table)
ps_all_metadata_tab  <- as.data.frame(ps_all@sam_data)
ps_all_pd <- pd((ps_all_asvtab), ps_all@phy_tree, include.root=T) 

#quick look how reads and observed richness correlate
ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(ps_all),
                         "observed" = phyloseq::estimate_richness(ps_all, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")
#can see from this plot, positive correlation between observed richness and total reads

#Generate a data.frame with alpha div measures
adiv <- data.frame("Observed" = phyloseq::estimate_richness(ps_all, measures = "Observed"),
                   "Shannon" = phyloseq::estimate_richness(ps_all, measures = "Shannon"),
                   "InvSimpson" = phyloseq::estimate_richness(ps_all, measures = "InvSimpson"),
                   "PD" = picante::pd(samp = data.frame((data.frame(phyloseq::otu_table(ps_all)))), tree = phyloseq::phy_tree(ps_all))[, 1],
                   "SampleOrigin" = phyloseq::sample_data(ps_all)$SampleOrigin)
#quick look at top and bottom to look for anything odd
head(adiv)
tail(adiv)
#then add these columns to a dataframe containing all the metadata taken from the pyloseq object 
ps_meta <- meta(ps_all)
ps_meta$Observed <-adiv$Observed
ps_meta$FaithPD <- adiv$PD
ps_meta$Shannon <- adiv$Shannon
ps_meta$InverseSimpson <- adiv$InvSimpson 
colnames(ps_meta)

#calculate means for each alpha metric by environment
#Observed richness
alpha_obs <- ps_meta %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  summarise(obs_mean=mean(Observed), obs_se=plotrix::std.error(Observed))
#Pylogenetic diversity
alpha_pd <- ps_meta %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  summarise(pd_mean=mean(FaithPD), pd_se=plotrix::std.error(FaithPD))
#shannon
alpha_shan <- ps_meta %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  summarise(shan_mean=mean(Shannon), shan_se=plotrix::std.error(Shannon))
#inverse simpson
alpha_sim <- ps_meta %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  summarise(sim_mean=mean(InverseSimpson), sim_se=plotrix::std.error(InverseSimpson))
#results are weird, CD1 and MF1's have highest mean for every metric 

#combine all alpha means for each metric into 1 table 
alpha_all <- cbind(alpha_obs, alpha_pd[2:3],alpha_shan[2:3],alpha_sim[2:3])


#quick look plot to see if any patterns - plot all 4 together, box and whiskers
adiv %>%
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  gather(key = metric, value = value, c("Observed", "PD", "Shannon","InvSimpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "PD", "Shannon","InvSimpson"))) %>%
  ggplot(aes(x = SampleOrigin, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = SampleOrigin), height = 0.6, width = 0.3, show.legend = F, alpha=1.0, size=2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme_bw(base_size = 20) +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#FDB462")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> p1
#CD1 mice seem to show lots of diversity which is unexpected, maybe due to moving facilities? 
#look into this more by plotting CD1 only at each timepoint, then same for wild-like and lab-like
#can figure out of certain timepoints are skewing things and need to get rid of some?
#MF1's also seem to have high alpha diversity, but put this down to there only being 10 samples

#Observed richness PLOT 
ps_meta %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=Observed, colour=SampleOrigin), 
                     alpha=1.0, size=3, show.legend = F) + 
  geom_point(data=alpha_obs, aes(x=SampleOrigin, y=obs_mean), color="black", size=3) + 
  geom_errorbar(data=alpha_obs, aes(x=SampleOrigin, ymin=obs_mean-obs_se, ymax=obs_mean+obs_se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Observed Richness") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#FDB462")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> obs_all

#Phylogenetic Diversity plot
ps_meta %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=FaithPD, colour=SampleOrigin), 
                     alpha=1.0, size=3, show.legend = F) + 
  geom_point(data=alpha_pd, aes(x=SampleOrigin, y=pd_mean), color="black", size=3) + 
  geom_errorbar(data=alpha_pd, aes(x=SampleOrigin, ymin=pd_mean-pd_se, ymax=pd_mean+pd_se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Phylogenetic Diversity") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#FDB462")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> pd_all

#Shannon Plot
ps_meta %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=Shannon, colour=SampleOrigin), 
                     alpha=1.0, size=3, show.legend = F) + 
  geom_point(data=alpha_shan, aes(x=SampleOrigin, y=shan_mean), color="black", size=3) + 
  geom_errorbar(data=alpha_shan, aes(x=SampleOrigin, ymin=shan_mean-shan_se, ymax=shan_mean+shan_se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Shannon Diversity") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#FDB462")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> shan_all

#Inverse Simpson Plot
ps_meta %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=InverseSimpson, colour=SampleOrigin), 
                     alpha=1.0, size=3, show.legend = F) + 
  geom_point(data=alpha_sim, aes(x=SampleOrigin, y=sim_mean), color="black", size=3) + 
  geom_errorbar(data=alpha_sim, aes(x=SampleOrigin, ymin=sim_mean-sim_se, ymax=sim_mean+sim_se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Inverse Simpson") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#FDB462")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> sim_all

#Panel plot of all alpha
alpha_pan <- obs_all|pd_all|shan_all|sim_all 

#Statistical analysis

#Observed richness - all groups except MF1 sig lower diversity than CD1s
hist(ps_meta$Observed) #slight skew
hist(log(ps_meta$Observed)) #looks better 
ps_meta$Log_Observed = log(ps_meta$Observed)
#Run GLM with log and non log data first 
obs_glm1 <-glm(Observed ~ SampleOrigin, data=ps_meta)
summary(obs_glm1)
#as expected model fit terrible

obsl_glm1 <-glm(Log_Observed ~ SampleOrigin, data=ps_meta)
summary(obsl_glm1) #CD1 = intercept 
#                            Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)                 5.46072    0.03543 154.121  < 2e-16 ***
# SampleOriginLab-like mice  -0.26503    0.05028  -5.271 2.76e-07 ***
# SampleOriginMF1 mice        0.09941    0.10146   0.980    0.328    
# SampleOriginWild mice      -0.12503    0.05384  -2.322    0.021 *  
# SampleOriginWild-like mice -0.57902    0.05046 -11.474  < 2e-16 ***
# (Dispersion parameter for gaussian family taken to be 0.09038743), AIC: 127.67
# Null deviance: 38.784  on 277  degrees of freedom, Residual deviance: 24.676  on 273  degrees of freedom

#phylogenetic diversity - all groups except MF1 sig lower diversity than CD1s
hist(ps_meta$FaithPD) # looks pretty normal 
hist(log(ps_meta$FaithPD)) # looks worse
ps_meta$Log_FaithPD = log(ps_meta$FaithPD)
#skewed right most values on lower size of x axis - looks more normal when logged 
fpd_glm1 <-glm(FaithPD ~ SampleOrigin, data=ps_meta)
summary(fpd_glm1)
#data looked good, but fit is terrible, really overdispersed 

fpdl_glm1 <-glm(Log_FaithPD ~ SampleOrigin, data=ps_meta)
summary(fpdl_glm1)
#                            Estimate Std. Error t value   Pr(>|t|)    
# (Intercept)                 6.51900    0.02670 244.169  < 2e-16 ***
# SampleOriginLab-like mice  -0.19049    0.03789  -5.027 9.01e-07 ***
# SampleOriginMF1 mice        0.06252    0.07645   0.818   0.4142    
# SampleOriginWild mice      -0.08782    0.04057  -2.165   0.0313 *  
# SampleOriginWild-like mice -0.43563    0.03803 -11.456  < 2e-16 ***
# (Dispersion parameter for gaussian family taken to be 0.05132336), AIC: -29.667
# Null deviance: 21.970  on 277  degrees of freedom, Residual deviance: 14.011  on 273  degrees of freedom


#shannon diversity - all groups except MF1 sig lower diversity than CD1s
hist(ps_meta$Shannon) # pretty skewed
hist(log(ps_meta$Shannon)) # looks even worse
#run on non-transformed data
shan_glm1 <-glm(Shannon ~ SampleOrigin, data=ps_meta)
summary(shan_glm1)
#                            Estimate Std. Error t value   Pr(>|t|)    
# (Intercept)                 3.89265    0.07391  52.667  < 2e-16 ***
# SampleOriginLab-like mice  -0.30586    0.10489  -2.916  0.00384 ** 
# SampleOriginMF1 mice        0.24120    0.21165   1.140  0.25543    
# SampleOriginWild mice      -0.87596    0.11231  -7.799 1.32e-13 ***
# SampleOriginWild-like mice -1.23251    0.10527 -11.708  < 2e-16 ***
# (Dispersion parameter for gaussian family taken to be 0.3933178), AIC: 536.47
# Null deviance: 178.10  on 277  degrees of freedom, Residual deviance: 107.38  on 273  degrees of freedom

#inverse simpson - all groups except MF1 sig lower diversity than CD1s
hist(ps_meta$InverseSimpson) # really skewed
hist(log(ps_meta$InverseSimpson)) # better
ps_meta$Log_InverseSimpson = log(ps_meta$InverseSimpson)
#skewed right most values on lower size of x axis) - looks more normal when logged 
sim_glm1 <-glm(InverseSimpson ~ SampleOrigin, data=ps_meta)
summary(sim_glm1)
#Not a great fit

siml_glm1 <-glm(Log_InverseSimpson ~ SampleOrigin, data=ps_meta)
summary(siml_glm1)
#                            Estimate Std. Error t value   Pr(>|t|)    
# (Intercept)                 3.01662    0.07688  39.240  < 2e-16 ***
# SampleOriginLab-like mice  -0.28265    0.10910  -2.591   0.0101 *  
# SampleOriginMF1 mice        0.30083    0.22014   1.367   0.1729    
# SampleOriginWild mice      -0.98207    0.11682  -8.407 2.36e-15 ***
# SampleOriginWild-like mice -1.20283    0.10949 -10.986  < 2e-16 ***
# (Dispersion parameter for gaussian family taken to be 0.4255088), AIC: 558.34
# Null deviance: 190.34  on 277  degrees of freedom, Residual deviance: 116.16  on 273  degrees of freedom


#------------------------------------------------------------------------------------------------------------------------------------------------------

#EXPLORING CD1 DIVERSITY 
#Make a phyloseq to remove day -7 and day -14 samples
#so only including CD1, lab-like and wild-like samples that have been at ashworth for 2 weeks (aclimatised - hopefully)
ps_d0 <- ps_all %>% 
  ps_filter(!Timept %in% c("-14", "-7"))
#make tables and calculate alpha diversity & add to df
ps_d0_asvtab <- as.data.frame(ps_d0@otu_table)
ps_d0_metadata_tab  <- as.data.frame(ps_d0@sam_data)
ps_d0_pd <- pd((ps_d0_asvtab), ps_d0@phy_tree, include.root=T) 

adiv_d0 <- data.frame("Observed" = phyloseq::estimate_richness(ps_d0, measures = "Observed"),
                      "Shannon" = phyloseq::estimate_richness(ps_d0, measures = "Shannon"),
                      "InvSimpson" = phyloseq::estimate_richness(ps_d0, measures = "InvSimpson"),
                      "PD" = picante::pd(samp = data.frame((data.frame(phyloseq::otu_table(ps_d0)))), 
                                         tree = phyloseq::phy_tree(ps_d0))[, 1],
                      "SampleOrigin" = phyloseq::sample_data(ps_d0)$SampleOrigin)

#then add these columns to a dataframe containing all the metadata taken from the pyloseq object 
ps_meta_d0 <- meta(ps_d0)
ps_meta_d0$Observed <-adiv_d0$Observed
ps_meta_d0$FaithPD <- adiv_d0$PD
ps_meta_d0$Shannon <- adiv_d0$Shannon
ps_meta_d0$InverseSimpson <- adiv_d0$InvSimpson 
colnames(ps_meta_d0)

#calculate means for each alpha metric by environment
#Observed richness
alpha_obs_d0 <- ps_meta_d0 %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  summarise(mean=mean(Observed), se=plotrix::std.error(Observed))
#Pylogenetic diversity
alpha_pd_d0 <- ps_meta_d0 %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  summarise(mean=mean(FaithPD), se=plotrix::std.error(FaithPD))
#shannon
alpha_shan_d0 <- ps_meta_d0 %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  summarise(mean=mean(Shannon), se=plotrix::std.error(Shannon))
#inverse simpson
alpha_sim_d0 <- ps_meta_d0 %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  summarise(mean=mean(InverseSimpson), se=plotrix::std.error(InverseSimpson))
#Again MF1 and CD1's still highest - what is going on? 

#plot figure to see if any diff
adiv_d0 %>%
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild mice", "Wild-like mice", 
                                                        "Lab-like mice", "CD1 mice", "MF1 mice"))) %>%
  gather(key = metric, value = value, c("Observed", "PD", "Shannon","InvSimpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "PD", "Shannon","InvSimpson"))) %>%
  ggplot(aes(x = SampleOrigin, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = SampleOrigin), height = 0.6, width = 0.3, show.legend = F, alpha=1.0, size=2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme_bw(base_size = 20) +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD",
                                 "#FB8072","#FDB462")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> p2







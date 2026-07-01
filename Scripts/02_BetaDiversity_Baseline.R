#BETA DIVERSITY

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
library(RColorBrewer)
library(emmeans)

source("/Volumes/Rowan_HD2/Chapter2-3_Analysis/VisualisationScripts/DataPrep_Baseline.R") 
# runs data prep & sets up all the data subsets 

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BETA DIVERSITY - ANALYSIS 1  

# All samples ps_all 
# Hewan, Penicuik & Falkirk samples all combined into one 'wild' category

#Add a sample origin column 
#where Hewan, Penicuik & Falkirk samples all called 'wild'
# MF1 & CD1 mice distingushied 
# ash54 wm colony = lab-like & b17 wm colony = wild-like
ps_all <- ps_all %>% 
  ps_mutate(SampleOrigin=case_when(str_detect(SampleDetailed, "MF1 lab mice") ~ "Lab:Mm",
                                   str_detect(SampleDetailed, "Pen") ~ "Wild:As", 
                                   str_detect(SampleDetailed, "Hew") ~ "Wild:As",
                                   str_detect(SampleDetailed, "Fal") ~ "Wild:As",
                                   str_detect(MouseOrigin, "Ash54 col") ~ "Lab-like:As", 
                                   str_detect(MouseOrigin, "b17 col") ~ "Wild-like:As",
                                   str_detect(MouseOrigin, "charles riv") ~ "Lab:Mm"))
                              

# CORE FOR ALL FAECAL SAMPLES
# some initial clean up - filtering out 'rare' taxa/ASV's 
# creating 'core' phyloseq objects with phylogenetic trees for each data subset
ps_core <- core(ps_all, detection =5, prevalence = 0.1)  
#drops ASVs from 3691 down to 563 - lost over 80%

core

#generate tree so we can do phylogentic analysis too 
# random_tree_core = rtree(ntaxa(ps_core), rooted=TRUE, tip.label=taxa_names(ps_core))
#save tree so don't have to repeat this step again
# write.tree(random_tree_core, "DescriptiveFigures/Ps_Core_PhyloTree")

#read in tree saved previously so no need to make a new one and change analysis 
random_tree_core <- read.tree("DescriptiveFigures/Ps_Core_PhyloTree")
#merge tree with phyloseq object

#create a phylogenetic tree for ps_core
ps_core <- merge_phyloseq(ps_core, random_tree_core)
# otu_table()   OTU Table:         [ 913 taxa and 278 samples ]
# sample_data() Sample Data:       [ 278 samples by 19 sample variables ]
# tax_table()   Taxonomy Table:    [ 913 taxa by 6 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 913 tips and 912 internal nodes ]

# BETA DIVERSITY ORDINATION PLOTS & ANALYSIS
ps_comp <- microbiome::transform(ps_core, "compositional")
ps_bray <- ordinate(ps_core, "PCoA", "bray")
ps_unif <- ordinate(ps_core, "PCoA", "unifrac")
ps_wunif <- ordinate(ps_core, "PCoA", "wunifrac")
ps_meta <- meta(ps_core)

# View(ps_bray[["vectors"]])
# View(ps_bray[["values"]])

#define what order you want the samples to be colored first 
sample_data(ps_comp)$SampleOrigin <- factor(sample_data(ps_comp)$SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))

#PLOT Bray Curtis ordination
plot_ordination(ps_comp, ps_bray, color="SampleOrigin") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) + 
  theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top", legend.title = element_blank()) +
 stat_ellipse() + 
  ggtitle('Bray-Curtis') -> All_Bray

#PLOT Weighted Unifrac ordination
plot_ordination(ps_comp, ps_wunif, color="SampleOrigin") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) + 
  theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top", legend.title = element_blank()) +
  stat_ellipse() + 
  ggtitle('Weighted Unifrac') -> All_Wunif

#PLOT Weighted Unifrac ordination
plot_ordination(ps_comp, ps_unif, color="SampleOrigin") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) + 
  theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top", legend.title = element_blank()) +
  stat_ellipse() + 
  ggtitle('Unweighted Unifrac') -> All_Unif

beta_pan <- (All_Bray|All_Wunif|All_Unif) + plot_layout(guides = 'collect') & theme(legend.position = "top")

#STATISTICAL ANALYSIS 
#re-order factors to set Lab-Mm as intercept 
ps_meta$SampleOrigin <- factor(ps_meta$SampleOrigin, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))
#reorder seq plate so that plate 8 is intercept 
ps_meta$SeqPlate <- factor(ps_meta$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))
#then convert seq plate back to integer - so doesn't show every plate 1-8 in model output as not needed 
ps_meta$SeqPlate <-as.numeric(ps_meta$SeqPlate)

#first make a column combining Id & SampleDetailed, otherwise random affect will treat all Ax mice as the same, even if from different rse experiemnts
ps_meta <- tidyr::unite(ps_meta,"MouseId", Id, SampleDetailed, remove = F) 

#BRAY CURTIS - PERMANOVA
metadf_all <- data.frame(sample_data(ps_meta))
bray_ps <- phyloseq::distance(physeq = ps_comp, method = "bray")
set.seed(999)
adonis_bray_all <- adonis2(bray_ps ~ SampleOrigin + Sex + SeqPlate, by="onedf", data = ps_meta)
adonis_bray_all

#to do pairwise analysis 
library(pairwiseAdonis)
pairwise.adonis2(bray_ps ~ SampleOrigin + Sex + SeqPlate,  data = ps_meta)

#WEIGHTED UNIFRAC - PERMANOVA
#If run with envrio:diet interaction could use for both envrio and diet analysis? 
wunif_ps <- phyloseq::distance(physeq = ps_comp, method = "wunifrac")
set.seed(999)
adonis_wunif_all <- adonis2(wunif_ps ~ SampleOrigin + Sex + SeqPlate, by="onedf", data = ps_meta)
#               Df SumOfSqs    R2      F    Pr(>F)    
# SampleOriginLab-like:As    1    4.508 0.13034 55.2349  0.001 ***
# SampleOriginWild-like:As   1    2.879 0.08323 35.2709  0.001 ***
# SampleOriginWild:As        1    4.479 0.12950 54.8816  0.001 ***
# SexM                       1    0.164 0.00473  2.0055  0.066 .  
# SeqPlate                   1    0.358 0.01036  4.3916  0.003 ** 
# Residual                 272   22.200 0.64184                   
# Total                    277   34.588 1.00000                  

pairwise.adonis2(wunif_ps ~ SampleOrigin + Sex + SeqPlate,  data = ps_meta)


#UNWEIGHTED UNIFRAC PERMANOVA
#If run with envrio:diet interaction could use for both envrio and diet analysis? 
unif_ps <- phyloseq::distance(physeq = ps_comp, method = "unifrac")
set.seed(999)
adonis_unif_all <- adonis2(unif_ps ~ SampleOrigin + Sex + SeqPlate, by="onedf", data = ps_meta)
#               Df SumOfSqs   R2      F     Pr(>F)    
# SampleOriginLab-like:As    1    7.024 0.16725 104.3971  0.001 ***
# SampleOriginWild-like:As   1    7.579 0.18047 112.6473  0.001 ***
# SampleOriginWild:As        1    8.629 0.20546 128.2504  0.001 ***
# SexM                       1    0.109 0.00259   1.6196  0.129    
# SeqPlate                   1    0.355 0.00846   5.2829  0.003 ** 
# Residual                 272   18.300 0.43576                    
# Total                    277   41.996 1.00000                 

pairwise.adonis2(wunif_ps ~ SampleOrigin + Sex + SeqPlate,  data = ps_meta)

#PLOT PC1 (variance/dissilarity scores) VS mouse cohort 
ps_comp <- microbiome::transform(ps_core, "compositional")
ps_bray <- ordinate(ps_core, "PCoA", "bray")
ps_unif <- ordinate(ps_core, "PCoA", "unifrac")
ps_wunif <- ordinate(ps_core, "PCoA", "wunifrac")
ps_meta <- meta(ps_core)

#make PC1 columns for each beta diversity  metric 
ps_meta$PC1_Bray <-ps_bray$vectors[,1]
ps_meta$PC1_Wunif <- ps_wunif$vectors[,1]
ps_meta$PC1_Unif <- ps_unif$vectors[,1]

#PLOTS 
#calculate PC1 bray mean and se so can plot error bars
bray_mean <- ps_meta %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(PC1_Bray), se=plotrix::std.error(PC1_Bray))

#PC1_Bray plot
ps_meta %>%
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=PC1_Bray, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=bray_mean, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=bray_mean, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Bray-Curtis (PC1)", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> PC1_Bray_fig 


#re-order factors to set Lab-Mm as intercept and reorder seq plate so that plate 8 is intercept 
ps_meta$SampleOrigin <- factor(ps_meta$SampleOrigin, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))
ps_meta$SeqPlate <- factor(ps_meta$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))

#glm - pc1 bray 
bray_pc1_glm <- glm(PC1_Bray ~ SampleOrigin + Sex + SeqPlate, data=ps_meta)
summary(bray_pc1_glm)

#first make a column combining Id & SampleDetailed, otherwise random affect will treat all Ax mice as the same, even if from different rse experiemnts
ps_meta <- tidyr::unite(ps_meta,"MouseId", Id, SampleDetailed, remove = F) 

library(TMB)
library(glmmTMB)
#glmm pc1 - bray curtis 
bray_pc1_glmm <- glmmTMB(PC1_Bray ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data=ps_meta)
summary(bray_pc1_glmm)
# AIC      BIC   logLik deviance df.resid 
# -593.8   -564.7    304.9   -609.8      270 
# Random effects:Conditional model:
#   Groups   Name        Variance  Std.Dev.
# SeqPlate (Intercept) 0.0005978 0.02445 
# MouseId  (Intercept) 0.0001073 0.01036 
# Residual             0.0061814 0.07862 
# Number of obs: 278, groups:  SeqPlate, 8; MouseId, 171
# Dispersion estimate for gaussian family (sigma^2): 0.00618 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              -0.016414   0.015221  -1.078  0.28087    
# SampleOriginLab:Mm       -0.266488   0.013390 -19.902  < 2e-16 ***
#   SampleOriginWild-like:As  0.257918   0.013528  19.065  < 2e-16 ***
#   SampleOriginWild:As       0.212455   0.026238   8.097 5.63e-16 ***
#   SexM                     -0.029935   0.009774  -3.063  0.00219 ** 

#pairwise
em <- emmeans(bray_pc1_glmm, "SampleOrigin")
contrast(em, "pairwise", adjust = "Tukey")

#calculate PC1 Weighted unifrac mean and se so can plot error bars
wunif_mean <- ps_meta %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(PC1_Wunif), se=plotrix::std.error(PC1_Wunif))

#PC1 weighted unifrac plot
ps_meta %>%
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=PC1_Wunif, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=wunif_mean, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=wunif_mean, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Weighted Unifrac (PC1)", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> PC1_Wunif_fig

#glm - pc1 wunif 
wunif_pc1_glm <- glm(PC1_Wunif ~ SampleOrigin + Sex + SeqPlate, data=ps_meta)
summary(wunif_pc1_glm)

#glmm - pc1 wunif
wunif_pc1_glmm <- glmmTMB(PC1_Wunif ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data=ps_meta)
summary(wunif_pc1_glmm)
# AIC      BIC   logLik deviance df.resid 
# -341.6   -312.6    178.8   -357.6      270 
# Random effects:Conditional model:
#   Groups   Name        Variance  Std.Dev. 
# SeqPlate (Intercept) 1.445e-11 3.801e-06
# MouseId  (Intercept) 1.964e-03 4.432e-02
# Residual             1.430e-02 1.196e-01
# Number of obs: 278, groups:  SeqPlate, 8; MouseId, 171
# Dispersion estimate for gaussian family (sigma^2): 0.0143 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              -0.166247   0.017870  -9.303  < 2e-16 ***
#   SampleOriginLab:Mm        0.109907   0.021808   5.040 4.66e-07 ***
#   SampleOriginWild-like:As  0.305261   0.022742  13.423  < 2e-16 ***
#   SampleOriginWild:As       0.267222   0.023611  11.318  < 2e-16 ***
#   SexM                      0.007771   0.016075   0.483    0.629  

#pairwise
em <- emmeans(wunif_pc1_glmm, "SampleOrigin")
contrast(em, "pairwise", adjust = "Tukey")

#calculate PC1 unweighted unifrac mean and se so can plot error bars
unif_mean <- ps_meta %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(PC1_Unif), se=plotrix::std.error(PC1_Unif))

#PC1 unweighted unifrac plot
ps_meta %>%
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=PC1_Unif, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=unif_mean, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=unif_mean, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Unweighted Unifrac (PC1)", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) -> PC1_Unif_fig

#panel figure
PC1_Bray_fig|PC1_Wunif_fig|PC1_Unif_fig


#glm - pc1 unif 
unif_pc1_glm <- glm(PC1_Unif ~ SampleOrigin + Sex + SeqPlate, data=ps_meta)
summary(unif_pc1_glm)

#glmm - pc1 wunif
unif_pc1_glmm <- glmmTMB(PC1_Unif ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data=ps_meta)
summary(unif_pc1_glmm)
# AIC      BIC   logLik deviance df.resid 
# -1105.0  -1076.0    560.5  -1121.0      270 
# Random effects:Conditional model:
#   Groups   Name        Variance  Std.Dev.
# SeqPlate (Intercept) 7.901e-05 0.008889
# MouseId  (Intercept) 5.240e-05 0.007239
# Residual             9.521e-04 0.030857
# Number of obs: 278, groups:  SeqPlate, 8; MouseId, 171
# Dispersion estimate for gaussian family (sigma^2): 0.000952 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               0.116419   0.006102   19.08   <2e-16 ***
#   SampleOriginLab:Mm        0.141998   0.005379   26.40   <2e-16 ***
#   SampleOriginWild-like:As -0.394849   0.005492  -71.90   <2e-16 ***
#   SampleOriginWild:As      -0.320654   0.009687  -33.10   <2e-16 ***
#   SexM                      0.005439   0.003944    1.38    0.168  

#pairwise
em <- emmeans(unif_pc1_glmm, "SampleOrigin")
contrast(em, "pairwise", adjust = "Tukey")

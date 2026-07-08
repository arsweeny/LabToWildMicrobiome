## ALPHA DIVERSITY

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
library(ggregplot)
library(emmeans)


source("Scripts/0_DataPrep_Baseline.R") 
# runs data prep & sets up all the data subsets 



# RAREFACTION - account for different sequencing depths/plates ------------

## Rarefy - subsampling reads from each sample without replacement to a constant depth)
ps_rare <- phyloseq::rarefy_even_depth(ps_all, rngseed = 123, replace = FALSE)
## 309 OTUs were removed because they are no longer present in any sample

head(phyloseq::sample_sums(ps_rare)) 
## sequence depth set to 15,857 

# create alpha diversity metrics ------------------------------------------

adiv_R <- data.frame("Observed" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
                      "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
                      "InvSimpson" = phyloseq::estimate_richness(ps_rare, measures = "InvSimpson"),
                      "PD" = picante::pd(samp = data.frame((data.frame(phyloseq::otu_table(ps_rare)))), 
                                         tree = phyloseq::phy_tree(ps_rare))[, 1],
                      "SampleOrigin" = phyloseq::sample_data(ps_rare)$SampleOrigin)

## add to dataframe containing all metadata taken from the pyloseq object 
ps_meta_rare <- meta(ps_rare)
ps_meta_rare$Observed <-adiv_R$Observed
ps_meta_rare$FaithPD <- adiv_R$PD
ps_meta_rare$Shannon <- adiv_R$Shannon
ps_meta_rare$InverseSimpson <- adiv_R$InvSimpson 
colnames(ps_meta_rare)

# # sanity checks / filter populations --------------------------------------
# 
# ps_meta_rare %>% filter(SampleOrigin == "Wild:As") %>% 
#   SinaGraph(., "SampleDetailed", "Observed")
# 
# ps_meta_rare %>% filter(SampleOrigin == "Lab:Mm") %>% 
#   SinaGraph(., "MouseOrigin", "Observed")
# 
# ### note going to edit up above now to remove these categories 



# calculate means for each alpha metric by environment --------------------

## observed richness
alpha_obs_R <- ps_meta_rare %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(Observed), se=plotrix::std.error(Observed), min=min(Observed), max=max(Observed), var=var(Observed), sd=sd(Observed))

## pylogenetic diversity
alpha_pd_R <- ps_meta_rare %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(FaithPD), se=plotrix::std.error(FaithPD), min=min(FaithPD), max=max(FaithPD), var=var(FaithPD), sd=sd(FaithPD)) 

## shannon
alpha_shan_R <- ps_meta_rare %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(Shannon), se=plotrix::std.error(Shannon), min=min(Shannon), max=max(Shannon), var=var(Shannon), sd=sd(Shannon)) 

## inverse simpson
alpha_sim_R <- ps_meta_rare %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(InverseSimpson), se=plotrix::std.error(InverseSimpson), min=min(InverseSimpson), max=max(InverseSimpson), var=var(InverseSimpson), sd=sd(InverseSimpson)) 


# assess correlations between alpha diversity metrics ---------------------

alpha_df <- adiv_R %>%
  dplyr::select(Observed, Shannon, InvSimpson, PD)
cor(alpha_df, use = "complete.obs", method = "pearson")
library(psych)
corr <- psych::corr.test(alpha_df, use = "pairwise")
corr$p


# alpha diversity graphs --------------------------------------------------

## boxplot faceted version
# alpha_boxplot <- adiv_R %>%
#   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
#   gather(key = metric, value = value, c("Observed", "PD", "Shannon","InvSimpson")) %>%
#   mutate(metric = factor(metric, levels = c("Observed", "PD", "Shannon","InvSimpson"))) %>%
#   ggplot(aes(x = SampleOrigin, y = value)) +
#   geom_boxplot(outlier.color = NA) +
#   geom_jitter(aes(color = SampleOrigin), height = 0.6, width = 0.3, show.legend = F, alpha=1.0, size=2) +
#   labs(x = "", y = "") +
#   facet_wrap(~ metric, scales = "free", nrow = 1) +
#   theme_bw(base_size = 20) +
#   scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
#   scale_x_discrete(guide = guide_axis(angle = 45))

## observed richness plot 
obs_all_R <- ps_meta_rare %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=Observed, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=alpha_obs_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=alpha_obs_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Observed Richness", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45))

## phylogenetic Diversity plot
pd_all_R <- ps_meta_rare %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=FaithPD, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=alpha_pd_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=alpha_pd_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Phylogenetic Diversity", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45))

#Shannon Plot
shan_all_R <- ps_meta_rare %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=Shannon, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=alpha_shan_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=alpha_shan_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Shannon Diversity", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45))

#Inverse Simpson Plot
sim_all_R <- ps_meta_rare %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  ggplot() +
  ggforce::geom_sina(aes(x=SampleOrigin, y=InverseSimpson, colour=SampleOrigin), 
                     alpha=1.0, size=2, show.legend = F) + 
  geom_point(data=alpha_sim_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
  geom_errorbar(data=alpha_sim_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
  theme_bw(base_size = 20) + 
  labs(y= "Inverse Simpson", x="") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45))

#Panel plot of all alpha
alpha_pan <- obs_all_R|pd_all_R|shan_all_R|sim_all_R 


# statistical analyses ----------------------------------------------------

# observed richness
hist(ps_meta_rare$Observed) # slight skew
hist(log(ps_meta_rare$Observed)) # looks better 
ps_meta_rare$Log_Observed = log(ps_meta_rare$Observed)

## re-order factors to set Lab-Mm as intercept and reorder seq plate so that plate 8 is intercept 
ps_meta_rare$SampleOrigin <- factor(ps_meta_rare$SampleOrigin, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))
ps_meta_rare$SeqPlate <- factor(ps_meta_rare$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))
ps_meta_rare$Sex <-as.factor(ps_meta_rare$Sex)

## run GLM with log and non-log data first 
obs_glm1 <- glm(Observed ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_rare)
summary(obs_glm1)

obsl_glm1 <-glm(Log_Observed ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_rare)
summary(obsl_glm1)

summ <- aggregate(Log_Observed ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_rare, mean)
ggplot(summ, aes(x=SampleOrigin, y=Log_Observed, colour=SeqPlate, group=SeqPlate)) +
  geom_point() +
  geom_line() +
  theme_bw() # alpha diversity difference by group is similar between SeqPlates

## run mixed model with seq plate as random effect
obsl_glmm <-glmmTMB(Log_Observed ~ SampleOrigin + Sex + (1|SeqPlate), data=ps_meta_rare)
summary(obsl_glmm)

## make a column combining Id & SampleDetailed so random effect doesn't treat 
## all Ax mice as the same, even if from different RSE experiments
ps_meta_rare <- tidyr::unite(ps_meta_rare,"MouseId", Id, SampleDetailed, remove = F) 

## glmm with seqplate and mouse id as random effects
obsl_glmm1 <-glmmTMB(Log_Observed ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data=ps_meta_rare)
summary(obsl_glmm1)

emm <- emmeans(obsl_glmm1, ~ SampleOrigin)
pairs(emm, adjust = "tukey")



# phylogenetic diversity 
hist(ps_meta_rare$FaithPD) # looks pretty normal 
hist(log(ps_meta_rare$FaithPD)) # looks worse
ps_meta_rare$Log_FaithPD = log(ps_meta_rare$FaithPD)
#skewed right most values on lower size of x axis - looks more normal when logged 

fpd_glm1 <- glm(FaithPD ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_rare)
summary(fpd_glm1) # overdispersed

fpdl_glm1 <- glm(Log_FaithPD ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_rare)
summary(fpdl_glm1)

fpdl_glmm <- glmmTMB(Log_FaithPD ~ SampleOrigin + Sex + (1|SeqPlate), data = ps_meta_rare)
summary(fpdl_glmm)

fpdl_glmm1 <- glmmTMB(Log_FaithPD ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta_rare)
summary(fpdl_glmm1)

emm <- emmeans(fpdl_glmm1, ~ SampleOrigin)
pairs(emm, adjust = "tukey")  



# shannon diversity
hist(ps_meta_rare$Shannon) # pretty skewed
hist(log(ps_meta_rare$Shannon)) # looks even worse

shan_glm1 <- glm(Shannon ~ SampleOrigin  + Sex + SeqPlate, data = ps_meta_rare)
summary(shan_glm1)

shan_glm <- glmmTMB(Shannon ~ SampleOrigin + Sex + (1|SeqPlate), data = ps_meta_rare)
summary(shan_glm)

shan_glm1 <- glmmTMB(Shannon ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta_rare)
summary(shan_glm1)

emm <- emmeans(shan_glm1, ~ SampleOrigin)
pairs(emm, adjust = "tukey")  


# inverse simpson
hist(ps_meta_rare$InverseSimpson) # really skewed
hist(log(ps_meta_rare$InverseSimpson)) # better
ps_meta_rare$Log_InverseSimpson = log(ps_meta_rare$InverseSimpson)
 
sim_glm1 <- glm(InverseSimpson ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_rare)
summary(sim_glm1)

siml_glm1 <- glm(Log_InverseSimpson ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_rare)
summary(siml_glm1)

siml_glmm <- glmmTMB(Log_InverseSimpson ~ SampleOrigin + Sex + (1|SeqPlate), data = ps_meta_rare)
summary(siml_glmm)

siml_glmm1 <- glmmTMB(Log_InverseSimpson ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta_rare)
summary(siml_glmm1)

emm <- emmeans(siml_glmm1, ~ SampleOrigin)
pairs(emm, adjust = "tukey") 



# # NO RAREFACTION analysis -------------------------------------------------
# 
# # generate tree for phylogentic analysis ----------------------------------
# 
# random_tree2 = rtree(ntaxa(ps_all), rooted=TRUE, tip.label=taxa_names(ps_all)) 
# 
# ## merge tree with phyloseq object
# ps_norare = merge_phyloseq(ps_all, random_tree)
# # otu_table()   OTU Table:         [ 2554 taxa and 239 samples ]
# # sample_data() Sample Data:       [ 239 samples by 23 sample variables ]
# # tax_table()   Taxonomy Table:    [ 2554 taxa by 6 taxonomic ranks ]
# # phy_tree()    Phylogenetic Tree: [ 2554 tips and 1 internal nodes ]
# 
# 
# 
# # create alpha diversity metrics ------------------------------------------
# 
# adiv_R <- data.frame("Observed" = phyloseq::estimate_richness(ps_norare, measures = "Observed"),
#                      "Shannon" = phyloseq::estimate_richness(ps_norare, measures = "Shannon"),
#                      "InvSimpson" = phyloseq::estimate_richness(ps_norare, measures = "InvSimpson"),
#                      "PD" = picante::pd(samp = data.frame((data.frame(phyloseq::otu_table(ps_norare)))), 
#                                         tree = phyloseq::phy_tree(ps_norare))[, 1],
#                      "SampleOrigin" = phyloseq::sample_data(ps_norare)$SampleOrigin)
# 
# ## add to dataframe containing all metadata taken from the pyloseq object 
# ps_meta_norare <- meta(ps_norare)
# ps_meta_norare$Observed <-adiv_R$Observed
# ps_meta_norare$FaithPD <- adiv_R$PD
# ps_meta_norare$Shannon <- adiv_R$Shannon
# ps_meta_norare$InverseSimpson <- adiv_R$InvSimpson 
# colnames(ps_meta_norare)
# 
# # calculate means for each alpha metric by environment --------------------
# 
# ## observed richness
# alpha_obs_R <- ps_meta_norare %>% group_by(SampleOrigin) %>% 
#   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
#   summarise(mean=mean(Observed), se=plotrix::std.error(Observed), min=min(Observed), max=max(Observed), var=var(Observed), sd=sd(Observed))
# 
# ## pylogenetic diversity
# alpha_pd_R <- ps_meta_norare %>% group_by(SampleOrigin) %>% 
#   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
#   summarise(mean=mean(FaithPD), se=plotrix::std.error(FaithPD), min=min(FaithPD), max=max(FaithPD), var=var(FaithPD), sd=sd(FaithPD)) 
# 
# ## shannon
# alpha_shan_R <- ps_meta_norare %>% group_by(SampleOrigin) %>% 
#   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
#   summarise(mean=mean(Shannon), se=plotrix::std.error(Shannon), min=min(Shannon), max=max(Shannon), var=var(Shannon), sd=sd(Shannon)) 
# 
# ## inverse simpson
# alpha_sim_R <- ps_meta_norare %>% group_by(SampleOrigin) %>% 
#   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
#   summarise(mean=mean(InverseSimpson), se=plotrix::std.error(InverseSimpson), min=min(InverseSimpson), max=max(InverseSimpson), var=var(InverseSimpson), sd=sd(InverseSimpson)) 
# 
# 
# # assess correlations between alpha diversity metrics ---------------------
# 
# alpha_df <- adiv_R %>%
#   dplyr::select(Observed, Shannon, InvSimpson, PD)
# cor(alpha_df, use = "complete.obs", method = "pearson")
# library(psych)
# corr <- psych::corr.test(alpha_df, use = "pairwise")
# corr$p
# 
# 
# # alpha diversity graphs --------------------------------------------------
# 
# ## boxplot faceted version
# # alpha_boxplot <- adiv_R %>%
# #   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
# #   gather(key = metric, value = value, c("Observed", "PD", "Shannon","InvSimpson")) %>%
# #   mutate(metric = factor(metric, levels = c("Observed", "PD", "Shannon","InvSimpson"))) %>%
# #   ggplot(aes(x = SampleOrigin, y = value)) +
# #   geom_boxplot(outlier.color = NA) +
# #   geom_jitter(aes(color = SampleOrigin), height = 0.6, width = 0.3, show.legend = F, alpha=1.0, size=2) +
# #   labs(x = "", y = "") +
# #   facet_wrap(~ metric, scales = "free", nrow = 1) +
# #   theme_bw(base_size = 20) +
# #   scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
# #   scale_x_discrete(guide = guide_axis(angle = 45))
# 
# ## observed richness plot 
# obs_all_R <- ps_meta_norare %>% 
#   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
#   ggplot() +
#   ggforce::geom_sina(aes(x=SampleOrigin, y=Observed, colour=SampleOrigin), 
#                      alpha=1.0, size=2, show.legend = F) + 
#   geom_point(data=alpha_obs_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
#   geom_errorbar(data=alpha_obs_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
#   theme_bw(base_size = 20) + 
#   labs(y= "Observed Richness", x="") +
#   scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
#   scale_x_discrete(guide = guide_axis(angle = 45))
# 
# ## phylogenetic Diversity plot
# pd_all_R <- ps_meta_norare %>% 
#   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
#   ggplot() +
#   ggforce::geom_sina(aes(x=SampleOrigin, y=FaithPD, colour=SampleOrigin), 
#                      alpha=1.0, size=2, show.legend = F) + 
#   geom_point(data=alpha_pd_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
#   geom_errorbar(data=alpha_pd_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
#   theme_bw(base_size = 20) + 
#   labs(y= "Phylogenetic Diversity", x="") +
#   scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
#   scale_x_discrete(guide = guide_axis(angle = 45))
# 
# #Shannon Plot
# shan_all_R <- ps_meta_norare %>% 
#   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
#   ggplot() +
#   ggforce::geom_sina(aes(x=SampleOrigin, y=Shannon, colour=SampleOrigin), 
#                      alpha=1.0, size=2, show.legend = F) + 
#   geom_point(data=alpha_shan_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
#   geom_errorbar(data=alpha_shan_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
#   theme_bw(base_size = 20) + 
#   labs(y= "Shannon Diversity", x="") +
#   scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
#   scale_x_discrete(guide = guide_axis(angle = 45))
# 
# #Inverse Simpson Plot
# sim_all_R <- ps_meta_norare %>% 
#   mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
#   ggplot() +
#   ggforce::geom_sina(aes(x=SampleOrigin, y=InverseSimpson, colour=SampleOrigin), 
#                      alpha=1.0, size=2, show.legend = F) + 
#   geom_point(data=alpha_sim_R, aes(x=SampleOrigin, y=mean), color="black", size=3) + 
#   geom_errorbar(data=alpha_sim_R, aes(x=SampleOrigin, ymin=mean-se, ymax=mean+se), width=0.3, size=0.5) + 
#   theme_bw(base_size = 20) + 
#   labs(y= "Inverse Simpson", x="") +
#   scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
#   scale_x_discrete(guide = guide_axis(angle = 45))
# 
# #Panel plot of all alpha
# alpha_pan <- obs_all_R|pd_all_R|shan_all_R|sim_all_R 
# 
# 
# # statistical analyses ----------------------------------------------------
# 
# # observed richness
# hist(ps_meta_norare$Observed) # slight skew
# hist(log(ps_meta_norare$Observed)) # looks better 
# ps_meta_norare$Log_Observed = log(ps_meta_norare$Observed)
# 
# ## re-order factors to set Lab-Mm as intercept and reorder seq plate so that plate 8 is intercept 
# ps_meta_norare$SampleOrigin <- factor(ps_meta_norare$SampleOrigin, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))
# ps_meta_norare$SeqPlate <- factor(ps_meta_norare$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))
# ps_meta_norare$Sex <-as.factor(ps_meta_norare$Sex)
# 
# ## run GLM with log and non-log data first 
# obs_glm1 <- glm(Observed ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_norare)
# summary(obs_glm1)
# 
# obsl_glm1 <-glm(Log_Observed ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_norare)
# summary(obsl_glm1)
# 
# summ <- aggregate(Log_Observed ~ SampleOrigin + Sex + SeqPlate, data=ps_meta_norare, mean)
# ggplot(summ, aes(x=SampleOrigin, y=Log_Observed, colour=SeqPlate, group=SeqPlate)) +
#   geom_point() +
#   geom_line() +
#   theme_bw() # alpha diversity difference by group is similar between SeqPlates
# 
# ## run mixed model with seq plate as random effect
# obsl_glmm <-glmmTMB(Log_Observed ~ SampleOrigin + Sex + (1|SeqPlate), data=ps_meta_norare)
# summary(obsl_glmm)
# 
# ## make a column combining Id & SampleDetailed so random effect doesn't treat 
# ## all Ax mice as the same, even if from different RSE experiments
# ps_meta_norare <- tidyr::unite(ps_meta_norare,"MouseId", Id, SampleDetailed, remove = F) 
# 
# ## glmm with seqplate and mouse id as random effects
# obsl_glmm1 <-glmmTMB(Log_Observed ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data=ps_meta_norare)
# summary(obsl_glmm1)
# 
# emm <- emmeans(obsl_glmm1, ~ SampleOrigin)
# pairs(emm, adjust = "tukey")
# 
# 
# 
# # phylogenetic diversity 
# hist(ps_meta_norare$FaithPD) # looks pretty normal 
# hist(log(ps_meta_norare$FaithPD)) # looks worse
# ps_meta_norare$Log_FaithPD = log(ps_meta_norare$FaithPD)
# #skewed right most values on lower size of x axis - looks more normal when logged 
# 
# fpd_glm1 <- glm(FaithPD ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_norare)
# summary(fpd_glm1) # overdispersed
# 
# fpdl_glm1 <- glm(Log_FaithPD ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_norare)
# summary(fpdl_glm1)
# 
# fpdl_glmm <- glmmTMB(Log_FaithPD ~ SampleOrigin + Sex + (1|SeqPlate), data = ps_meta_norare)
# summary(fpdl_glmm)
# 
# fpdl_glmm1 <- glmmTMB(Log_FaithPD ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta_norare)
# summary(fpdl_glmm1)
# 
# emm <- emmeans(fpdl_glmm1, ~ SampleOrigin)
# pairs(emm, adjust = "tukey")  
# 
# 
# 
# # shannon diversity
# hist(ps_meta_norare$Shannon) # pretty skewed
# hist(log(ps_meta_norare$Shannon)) # looks even worse
# 
# shan_glm1 <- glm(Shannon ~ SampleOrigin  + Sex + SeqPlate, data = ps_meta_norare)
# summary(shan_glm1)
# 
# shan_glm <- glmmTMB(Shannon ~ SampleOrigin + Sex + (1|SeqPlate), data = ps_meta_norare)
# summary(shan_glm)
# 
# shan_glm1 <- glmmTMB(Shannon ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta_norare)
# summary(shan_glm1)
# 
# emm <- emmeans(shan_glm1, ~ SampleOrigin)
# pairs(emm, adjust = "tukey")  
# 
# 
# # inverse simpson
# hist(ps_meta_norare$InverseSimpson) # really skewed
# hist(log(ps_meta_norare$InverseSimpson)) # better
# ps_meta_norare$Log_InverseSimpson = log(ps_meta_norare$InverseSimpson)
# 
# sim_glm1 <- glm(InverseSimpson ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_norare)
# summary(sim_glm1)
# 
# siml_glm1 <- glm(Log_InverseSimpson ~ SampleOrigin + Sex + SeqPlate, data = ps_meta_norare)
# summary(siml_glm1)
# 
# siml_glmm <- glmmTMB(Log_InverseSimpson ~ SampleOrigin + Sex + (1|SeqPlate), data = ps_meta_norare)
# summary(siml_glmm)
# 
# siml_glmm1 <- glmmTMB(Log_InverseSimpson ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta_norare)
# summary(siml_glmm1)
# 
# emm <- emmeans(siml_glmm1, ~ SampleOrigin)
# pairs(emm, adjust = "tukey") 




# removing mice with no acclimatization time ------------------------------

ps_d0 <- ps_all %>% 
  ps_filter(!Timept %in% c("-14", "-7"))

random_tree3 = rtree(ntaxa(ps_d0), rooted=TRUE, tip.label=taxa_names(ps_d0)) 

## merge tree with phyloseq object
ps_d0 = merge_phyloseq(ps_d0, random_tree3)
# otu_table()   OTU Table:         [ 2511 taxa and 133 samples ]
# sample_data() Sample Data:       [ 133 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 2511 taxa by 6 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2511 tips and 2510 internal nodes ]


# create alpha diversity metrics ------------------------------------------

adiv_R <- data.frame("Observed" = phyloseq::estimate_richness(ps_d0, measures = "Observed"),
                     "Shannon" = phyloseq::estimate_richness(ps_d0, measures = "Shannon"),
                     "InvSimpson" = phyloseq::estimate_richness(ps_d0, measures = "InvSimpson"),
                     "PD" = picante::pd(samp = data.frame((data.frame(phyloseq::otu_table(ps_d0)))), 
                                        tree = phyloseq::phy_tree(ps_d0))[, 1],
                     "SampleOrigin" = phyloseq::sample_data(ps_d0)$SampleOrigin)

## add to dataframe containing all metadata taken from the pyloseq object 
ps_d0_2 <- meta(ps_d0)
ps_d0_2$Observed <-adiv_R$Observed
ps_d0_2$FaithPD <- adiv_R$PD
ps_d0_2$Shannon <- adiv_R$Shannon
ps_d0_2$InverseSimpson <- adiv_R$InvSimpson 
colnames(ps_d0_2)

# calculate means for each alpha metric by environment --------------------

## observed richness
alpha_obs_R <- ps_d0_2 %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(Observed), se=plotrix::std.error(Observed), min=min(Observed), max=max(Observed), var=var(Observed), sd=sd(Observed))

## pylogenetic diversity
alpha_pd_R <- ps_d0_2 %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(FaithPD), se=plotrix::std.error(FaithPD), min=min(FaithPD), max=max(FaithPD), var=var(FaithPD), sd=sd(FaithPD)) 

## shannon
alpha_shan_R <- ps_d0_2 %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(Shannon), se=plotrix::std.error(Shannon), min=min(Shannon), max=max(Shannon), var=var(Shannon), sd=sd(Shannon)) 

## inverse simpson
alpha_sim_R <- ps_d0_2 %>% group_by(SampleOrigin) %>% 
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  summarise(mean=mean(InverseSimpson), se=plotrix::std.error(InverseSimpson), min=min(InverseSimpson), max=max(InverseSimpson), var=var(InverseSimpson), sd=sd(InverseSimpson)) 


# assess correlations between alpha diversity metrics ---------------------

alpha_df <- adiv_R %>%
  dplyr::select(Observed, Shannon, InvSimpson, PD)
cor(alpha_df, use = "complete.obs", method = "pearson")
library(psych)
corr <- psych::corr.test(alpha_df, use = "pairwise")
corr$p

# alpha diversity graphs --------------------------------------------------
alpha_boxplot <- adiv_R %>%
  mutate(SampleOrigin = factor(SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>%
  gather(key = metric, value = value, c("Observed", "PD", "Shannon","InvSimpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "PD", "Shannon","InvSimpson"))) %>%
  ggplot(aes(x = SampleOrigin, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = SampleOrigin), height = 0.6, width = 0.3, show.legend = F, alpha=1.0, size=2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free", nrow = 1) +
  theme_bw(base_size = 20) +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  scale_x_discrete(guide = guide_axis(angle = 45))
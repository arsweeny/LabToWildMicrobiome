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
library(psych)
library(ggpubr)
library(car)

# Load the final preprocessed baseline phyloseq object.
# Run Scripts/00_DataPrep_Baseline.R first whenever preprocessing changes.
ps_all <- readRDS("Data/Processed/baseline_microbiome_phyloseq.rds")

# Results from this script are kept together in one self-contained folder.
results_dir <- "Results/01_Alpha_Diversity"
tables_dir <- file.path(results_dir, "Tables")
figures_dir <- file.path(results_dir, "Figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
# runs data prep & sets up all the data subsets 



# RAREFACTION - account for different sequencing depths/plates ------------

## Rarefy - subsampling reads from each sample without replacement to a constant depth)
ps_rare <- phyloseq::rarefy_even_depth(ps_all, rngseed = 123, replace = FALSE)
## 315 OTUs were removed because they are no longer present in any sample

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

corr <- psych::corr.test(alpha_df, use = "pairwise")
corr$p

library(corrplot)
cor_mat <- cor(alpha_df, use = "complete.obs", method = "pearson")

corrplot(cor_mat,
         method = "color",
         type = "upper",
         addCoef.col = "black",   # show correlation coefficients
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE)

# alpha diversity graphs --------------------------------------------------

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

# panel plot of all alpha
alpha_pan <- obs_all_R|pd_all_R|shan_all_R|sim_all_R 

# panel plot of just shannon and faiths PD
alpha_pan <- ggarrange(shan_all_R, pd_all_R, labels = c("A","B"), nrow = 1)

## save final plot
pdf(file.path(figures_dir, "alpha_diversity_by_sample_origin.pdf"), width = 6, height = 5)
print(alpha_pan)
dev.off()

# statistical analyses ----------------------------------------------------

## re-order factors to set Lab-Mm as intercept and reorder seq plate so that plate 8 is intercept 
ps_meta_rare$SampleOrigin <- factor(ps_meta_rare$SampleOrigin, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))
ps_meta_rare$SeqPlate <- factor(ps_meta_rare$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))
ps_meta_rare$Sex <- as.factor(ps_meta_rare$Sex)

## make a column combining Id & SampleDetailed so random effect doesn't treat 
## all Ax mice as the same, even if from different RSE experiments
ps_meta_rare <- tidyr::unite(ps_meta_rare,"MouseId", Id, SampleDetailed, remove = F) 

# observed richness

## glmm with seqplate and mouse id as random effects
## remember to use mouse_ID_dataset (not mouse_ID) because mice are named the 
## same thing across the three experiment types
obsl_glmm1 <-glmmTMB(Observed ~ SampleOrigin + Sex + (1|SeqPlate) + (1|mouse_ID_dataset), data=ps_meta_rare)
summary(obsl_glmm1)
## not logging Observed bc model diagnostics are fine without it

emm <- emmeans(obsl_glmm1, ~ SampleOrigin)
pairs(emm, adjust = "tukey")

Anova(obsl_glmm1, type = 3)

library(DHARMa)
res <- simulateResiduals(obsl_glmm1)
plot(res)

testDispersion(res)
testUniformity(res)
testOutliers(res)


# phylogenetic diversity 
fpdl_glmm1 <- glmmTMB(FaithPD ~ SampleOrigin + Sex + (1|SeqPlate) + (1|mouse_ID_dataset), data = ps_meta_rare)
summary(fpdl_glmm1)
## not logging Observed bc model diagnostics are fine without it

emm <- emmeans(fpdl_glmm1, ~ SampleOrigin)
pairs(emm, adjust = "tukey")  

Anova(fpdl_glmm1, type = 3)

res <- simulateResiduals(fpdl_glmm1)
plot(res)

testDispersion(res)
testUniformity(res)
testOutliers(res)

# shannon diversity
hist(ps_meta_rare$Shannon) # pretty skewed
hist(log(ps_meta_rare$Shannon)) # looks even worse

shan_glm1 <- glmmTMB(Shannon ~ SampleOrigin + Sex + (1|SeqPlate) + (1|mouse_ID_dataset), data = ps_meta_rare)
summary(shan_glm1)

emm <- emmeans(shan_glm1, ~ SampleOrigin)
pairs(emm, adjust = "tukey")  

Anova(shan_glm1, type = 3)

res <- simulateResiduals(shan_glm1)
plot(res)

testDispersion(res)
testUniformity(res)
testOutliers(res)

## look at the outlier values
## outliers are evenly spread across biological and technical groups
## outliers are all showing plausible values for shannon diversity
## will not remove any outliers
extreme <- c(
  order(res$scaledResiduals)[1:5],
  order(res$scaledResiduals, decreasing = TRUE)[1:5]
)

ps_meta_rare[extreme, c("Observed",
                        "SampleOrigin",
                        "Sex",
                        "SeqPlate",
                        "mouse_ID_dataset",
                        "SeqDepth")]

# inverse simpson
ps_meta_rare$Log_InverseSimpson = log(ps_meta_rare$InverseSimpson)

siml_glmm1 <- glmmTMB(Log_InverseSimpson ~ SampleOrigin + Sex + (1|SeqPlate) + (1|mouse_ID_dataset), data = ps_meta_rare)
summary(siml_glmm1)
## needs to be logged to achieve good model diagnostics

emm <- emmeans(siml_glmm1, ~ SampleOrigin)
pairs(emm, adjust = "tukey") 

Anova(siml_glmm1, type = 3)

res <- simulateResiduals(siml_glmm1)
plot(res)

testDispersion(res)
testUniformity(res)
testOutliers(res)



# export alpha diversity results into excel table -------------------------

library(glmmTMB)
library(car)
library(emmeans)
library(broom.mixed)
library(broom)
library(dplyr)
library(openxlsx)

## Store models
alpha_models <- list(
  "Observed richness" = obsl_glmm1,
  "Faith PD" = fpdl_glmm1,
  "Shannon" = shan_glm1,
  "log(Inverse Simpson)" = siml_glmm1
)


## GLMM fixed effects
alpha_GLMM_results <- lapply(names(alpha_models), function(metric){
  
  broom.mixed::tidy(
    alpha_models[[metric]],
    effects = "fixed",
    conf.int = TRUE
  ) %>%
    mutate(Metric = metric)
  
}) %>%
  bind_rows() %>%
  select(
    Metric,
    term,
    estimate,
    std.error,
    statistic,
    p.value,
    conf.low,
    conf.high
  ) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE            ~ ""
    )
  )


## Type III Wald chi-square tests
tidy_type3_anova <- function(model, metric){
  
  x <- as.data.frame(car::Anova(model, type = 3))
  
  x$Predictor <- rownames(x)
  rownames(x) <- NULL
  
  x %>%
    filter(Predictor != "(Intercept)") %>%
    transmute(
      Metric = metric,
      Predictor,
      Chisq,
      Df,
      p.value = `Pr(>Chisq)`,
      significance = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE            ~ ""
      )
    )
  
}

anova_results <- lapply(names(alpha_models), function(metric){
  
  tidy_type3_anova(alpha_models[[metric]], metric)
  
}) %>%
  bind_rows()


## Tukey post hoc comparisons
tukey_results <- lapply(names(alpha_models), function(metric){
  
  emmeans(alpha_models[[metric]], ~ SampleOrigin) %>%
    pairs(adjust = "tukey") %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    mutate(Metric = metric)
  
}) %>%
  bind_rows() %>%
  select(
    Metric,
    contrast,
    estimate,
    SE,
    df,
    z.ratio,
    p.value,
    asymp.LCL,
    asymp.UCL
  )


## Model statistics
model_stats <- data.frame(
  
  Metric = names(alpha_models),
  
  Transformation = c(
    "None",
    "None",
    "None",
    "Log"
  ),
  
  Family = rep("Gaussian", length(alpha_models)),
  
  Random_effects =
    rep("SeqPlate + mouse_ID_dataset", length(alpha_models)),
  
  AIC = sapply(alpha_models, AIC),
  
  BIC = sapply(alpha_models, BIC),
  
  N = sapply(alpha_models, nobs)
  
)


## Random effect variance
random_effects <- lapply(names(alpha_models), function(metric){
  
  broom.mixed::tidy(
    alpha_models[[metric]],
    effects = "ran_pars"
  ) %>%
    mutate(Metric = metric)
  
}) %>%
  bind_rows()


## Export workbook
write.xlsx(
  
  list(
    
    GLMM_fixed_effects = alpha_GLMM_results,
    
    Type_III_ANOVA = anova_results,
    
    Tukey_posthoc = tukey_results,
    
    Model_statistics = model_stats,
    
    Random_effects = random_effects
    
  ),
  
  file = file.path(tables_dir, "alpha_diversity_statistical_results.xlsx"),
  
  rowNames = FALSE
  
)

# # NON RAREFACTION analysis ------------------------------------------------
# # see comparison in Output folder / supplemental material
# # create alpha diversity metrics ------------------------------------------
# 
# adiv_R <- data.frame("Observed" = phyloseq::estimate_richness(ps_all, measures = "Observed"),
#                      "Shannon" = phyloseq::estimate_richness(ps_all, measures = "Shannon"),
#                      "InvSimpson" = phyloseq::estimate_richness(ps_all, measures = "InvSimpson"),
#                      "PD" = picante::pd(samp = data.frame((data.frame(phyloseq::otu_table(ps_all)))),
#                                         tree = phyloseq::phy_tree(ps_all))[, 1],
#                      "SampleOrigin" = phyloseq::sample_data(ps_all)$SampleOrigin)
# 
# ## add to dataframe containing all metadata taken from the pyloseq object
# ps_meta_norare <- meta(ps_all)
# ps_meta_norare$Observed <-adiv_R$Observed
# ps_meta_norare$FaithPD <- adiv_R$PD
# ps_meta_norare$Shannon <- adiv_R$Shannon
# ps_meta_norare$InverseSimpson <- adiv_R$InvSimpson
# colnames(ps_meta_norare)
# 
# 
# # assess correlations between alpha diversity metrics ---------------------
# 
# alpha_df <- adiv_R %>%
#   dplyr::select(Observed, Shannon, InvSimpson, PD)
# 
# corr <- psych::corr.test(alpha_df, use = "pairwise")
# corr$p
# 
# library(corrplot)
# cor_mat <- cor(alpha_df, use = "complete.obs", method = "pearson")
# 
# corrplot(cor_mat,
#          method = "color",
#          type = "upper",
#          addCoef.col = "black",   # show correlation coefficients
#          tl.col = "black",
#          tl.srt = 45,
#          diag = FALSE)
# 
# # statistical analyses ----------------------------------------------------
# 
# ## re-order factors to set Lab-Mm as intercept and reorder seq plate so that plate 8 is intercept 
# ps_meta_norare$SampleOrigin <- factor(ps_meta_norare$SampleOrigin, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))
# ps_meta_norare$SeqPlate <- factor(ps_meta_norare$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))
# ps_meta_norare$Sex <- as.factor(ps_meta_norare$Sex)
# 
# ## make a column combining Id & SampleDetailed so random effect doesn't treat 
# ## all Ax mice as the same, even if from different RSE experiments
# ps_meta_norare <- tidyr::unite(ps_meta_norare,"MouseId", Id, SampleDetailed, remove = F) 
# 
# # observed richness
# 
# ## glmm with seqplate and mouse id as random effects
# ## remember to use mouse_ID_dataset (not mouse_ID) because mice are named the 
# ## same thing across the three experiment types
# obsl_glmm1 <-glmmTMB(Observed ~ SampleOrigin + Sex + (1|SeqPlate) + (1|mouse_ID_dataset), data=ps_meta_norare)
# summary(obsl_glmm1)
# ## not logging Observed bc model diagnostics are fine without it
# 
# emm <- emmeans(obsl_glmm1, ~ SampleOrigin)
# pairs(emm, adjust = "tukey")
# 
# library(DHARMa)
# res <- simulateResiduals(obsl_glmm1)
# plot(res)
# 
# testDispersion(res)
# testUniformity(res)
# testOutliers(res)
# 
# 
# # phylogenetic diversity 
# fpdl_glmm1 <- glmmTMB(FaithPD ~ SampleOrigin + Sex + (1|SeqPlate) + (1|mouse_ID_dataset), data = ps_meta_norare)
# summary(fpdl_glmm1)
# ## not logging Observed bc model diagnostics are fine without it
# 
# emm <- emmeans(fpdl_glmm1, ~ SampleOrigin)
# pairs(emm, adjust = "tukey")  
# 
# res <- simulateResiduals(fpdl_glmm1)
# plot(res)
# 
# testDispersion(res)
# testUniformity(res)
# testOutliers(res)
# 
# # shannon diversity
# shan_glm1 <- glmmTMB(Shannon ~ SampleOrigin + Sex + (1|SeqPlate) + (1|mouse_ID_dataset), data = ps_meta_norare)
# summary(shan_glm1)
# 
# emm <- emmeans(shan_glm1, ~ SampleOrigin)
# pairs(emm, adjust = "tukey")  
# 
# res <- simulateResiduals(shan_glm1)
# plot(res)
# 
# testDispersion(res)
# testUniformity(res)
# testOutliers(res)
# 
# ## look at the outlier values
# ## outliers are evenly spread across biological and technical groups
# ## outliers are all showing plausible values for shannon diversity
# ## will not remove any outliers
# extreme <- c(
#   order(res$scaledResiduals)[1:5],
#   order(res$scaledResiduals, decreasing = TRUE)[1:5]
# )
# 
# ps_meta_norare[extreme, c("Observed",
#                         "SampleOrigin",
#                         "Sex",
#                         "SeqPlate",
#                         "mouse_ID_dataset",
#                         "SeqDepth")]
# 
# # inverse simpson
# ps_meta_norare$Log_InverseSimpson = log(ps_meta_norare$InverseSimpson)
# 
# siml_glmm1 <- glmmTMB(Log_InverseSimpson ~ SampleOrigin + Sex + (1|SeqPlate) + (1|mouse_ID_dataset), data = ps_meta_norare)
# summary(siml_glmm1)
# ## needs to be logged to achieve good model diagnostics
# 
# emm <- emmeans(siml_glmm1, ~ SampleOrigin)
# pairs(emm, adjust = "tukey") 
# 
# res <- simulateResiduals(siml_glmm1)
# plot(res)
# 
# testDispersion(res)
# testUniformity(res)
# testOutliers(res)




# other quick checks ------------------------------------------------------

# does acclimatization time matter?
ps_d0 <- ps_all %>% 
  ps_filter(!Timept %in% c("-14", "-7"))

## create alpha diversity metrics 
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

## calculate means for each alpha metric by environment

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

## alpha diversity graph without non-acclimatized mice
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

## plotting time before start of experiment against observed alpha diversity
## not a huge effect
ggplot(ps_meta_rare[ps_meta_rare$SampleOrigin!="Wild:As",], 
       aes(x = Timept, y = Observed, fill = SampleOrigin)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw()



# does age at sampling influence alpha diversity? -------------------------

## age doesn't seem to influence alpha diversity too much 
## across the different groups
ggplot(ps_meta_rare[ps_meta_rare$SampleOrigin!="Wild:As",], 
  aes(x = age_sampling, y = Observed, fill = SampleOrigin)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw()

ggplot(ps_meta_rare, aes(x = age_sampling, y = FaithPD, fill = SampleOrigin)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ps_meta_rare$BodyMassG <- as.numeric(ps_meta_rare$BodyMassG)
ggplot(ps_meta_rare, aes(x = BodyMassG, y = Shannon, fill = SampleOrigin)) + 
  geom_point() + 
  geom_smooth(method = "lm")



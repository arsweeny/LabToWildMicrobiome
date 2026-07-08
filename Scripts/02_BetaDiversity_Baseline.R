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
library(pairwiseAdonis)

source("Scripts/0_DataPrep_Baseline.R") 
# runs data prep & sets up all the data subsets 

# filter out 'rare' ASVs and create core ----------------------------------

ps_core <- core(ps_all, detection =5, prevalence = 0.1)  
## drops ASVs from 2863 down to 546 - lost over 80%

## remove samples with missing metadata BEFORE any ordinations or statistics
keep <- !is.na(sample_data(ps_core)$Sex)
ps_core <- prune_samples(keep, ps_core)

## prune tree to only those remaining ASVs
tree_core <- ape::drop.tip(
  phy_tree(ps_all),
  setdiff(phy_tree(ps_all)$tip.label, taxa_names(ps_core))
)

phy_tree(ps_core) <- tree_core


# beta diversity ordination plots -----------------------------------------

ps_comp <- microbiome::transform(ps_core, "compositional") # from counts to relative abundances
ps_meta <- meta(ps_comp) # metadata
ps_bray <- ordinate(ps_comp, "PCoA", "bray")
ps_unif <- ordinate(ps_core, "PCoA", "unifrac")
ps_wunif <- ordinate(ps_core, "PCoA", "wunifrac")
ps_jaccard <- ordinate(ps_comp, method = "PCoA", distance = "jaccard", binary = TRUE)

## define order of sample groups 
sample_data(ps_comp)$SampleOrigin <- factor(sample_data(ps_comp)$SampleOrigin, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))

## plot Bray Curtis ordination
All_Bray <- plot_ordination(ps_comp, ps_bray, color="SampleOrigin") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) + 
  theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top", legend.title = element_blank()) +
 stat_ellipse() + 
  ggtitle('Bray-Curtis')
All_Bray

## plot Weighted Unifrac ordination
All_Wunif <- plot_ordination(ps_comp, ps_wunif, color="SampleOrigin") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) + 
  theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top", legend.title = element_blank()) +
  stat_ellipse() + 
  ggtitle('Weighted Unifrac')
All_Wunif

## plot Weighted Unifrac ordination
All_Unif <- plot_ordination(ps_comp, ps_unif, color="SampleOrigin") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) + 
  theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top", legend.title = element_blank()) +
  stat_ellipse() + 
  ggtitle('Unweighted Unifrac')
All_Unif

## plot jaccard ordination 
All_Jaccard <- plot_ordination(ps_comp, ps_jaccard, color = "SampleOrigin") +
  geom_point(size = 2) +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD", "#FB8072")) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0, size = 22),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  stat_ellipse() +
  ggtitle("Jaccard")
All_Jaccard

## plot all ordinations together
beta_pan <- (All_Bray|All_Jaccard|All_Wunif|All_Unif) + plot_layout(guides = 'collect') & theme(legend.position = "top")
beta_pan


# statistical analyses ----------------------------------------------------

## re-order factors to set Lab-Mm as intercept 
ps_meta$SampleOrigin <- factor(ps_meta$SampleOrigin, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))

## reorder seq plate so that plate 8 is intercept 
ps_meta$SeqPlate <- factor(ps_meta$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))

## convert seq plate back to integer
ps_meta$SeqPlate <- as.numeric(ps_meta$SeqPlate)

## make a column combining Id & SampleDetailed for random affect 
## otherwise will treat all Ax mice as the same, even if from different rse experiments
ps_meta <- tidyr::unite(ps_meta,"MouseId", Id, SampleDetailed, remove = F) 


# BRAY CURTIS - PERMANOVA
metadf_all <- data.frame(sample_data(ps_meta))
bray_ps <- phyloseq::distance(physeq = ps_comp, method = "bray")
set.seed(999)
adonis_bray_all <- adonis2(bray_ps ~ SampleOrigin + Sex + SeqPlate, by = "onedf", data = ps_meta)
adonis_bray_all

pairwise.adonis2(bray_ps ~ SampleOrigin + Sex + SeqPlate,  data = ps_meta)

bd <- betadisper(bray_ps, ps_meta$SampleOrigin)
anova(betadisper(bray_ps, ps_meta$SampleOrigin))
permutest(bd, pairwise = TRUE)
plot(bd)

# WEIGHTED UNIFRAC - PERMANOVA
wunif_ps <- phyloseq::distance(physeq = ps_comp, method = "wunifrac")
set.seed(999)
adonis_wunif_all <- adonis2(wunif_ps ~ SampleOrigin + Sex + SeqPlate, by = "onedf", data = ps_meta)

pairwise.adonis2(wunif_ps ~ SampleOrigin + Sex + SeqPlate,  data = ps_meta)

bd <- betadisper(wunif_ps, ps_meta$SampleOrigin)
anova(betadisper(wunif_ps, ps_meta$SampleOrigin))
permutest(bd, pairwise = TRUE)
plot(bd)

# UNWEIGHTED UNIFRAC - PERMANOVA
unif_ps <- phyloseq::distance(physeq = ps_comp, method = "unifrac")
set.seed(999)
adonis_unif_all <- adonis2(unif_ps ~ SampleOrigin + Sex + SeqPlate, by="onedf", data = ps_meta)

pairwise.adonis2(unif_ps ~ SampleOrigin + Sex + SeqPlate,  data = ps_meta)

bd <- betadisper(unif_ps, ps_meta$SampleOrigin)
anova(betadisper(unif_ps, ps_meta$SampleOrigin))
permutest(bd, pairwise = TRUE)
plot(bd)

# JACCARD - PERMANOVA
jaccard_ps <- phyloseq::distance(physeq = ps_comp, method = "jaccard", binary = TRUE)
set.seed(999)
adonis_jaccard_all <- adonis2(jaccard_ps ~ SampleOrigin + Sex + SeqPlate, data = ps_meta, by = "onedf")

pairwise.adonis2(jaccard_ps ~ SampleOrigin + Sex + SeqPlate, data = ps_meta)

bd <- betadisper(jaccard_ps, ps_meta$SampleOrigin)
anova(betadisper(jaccard_ps, ps_meta$SampleOrigin))
permutest(bd, pairwise = TRUE)
plot(bd)

# plot PC1 (variance/dissimilarity scores) across groups ------------------

## extract PC1 coordinates from previously calculated PCoAs
ps_meta$PC1_Bray  <- ps_bray$vectors[,1]
ps_meta$PC1_Wunif <- ps_wunif$vectors[,1]
ps_meta$PC1_Unif  <- ps_unif$vectors[,1]
ps_meta$PC1_Jaccard <- ps_jaccard$vectors[,1]

## set plotting order
plot_order <- c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm")
ps_meta$SampleOrigin <- factor(ps_meta$SampleOrigin, levels = plot_order)

## colour scheme
origin_colours <- c("Wild:As" = "#74C476", "Wild-like:As" = "#6BAED6", "Lab-like:As" = "#BC80BD", "Lab:Mm" = "#FB8072")

## function to calculate mean +/- SE 
pc1_summary <- function(data, variable){
  data %>%
    group_by(SampleOrigin) %>%
    summarise(
      mean = mean(.data[[variable]], na.rm = TRUE),
      se = plotrix::std.error(.data[[variable]], na.rm = TRUE),
      .groups = "drop"
    )
}

## function to make PC1 plots -----------------------------
pc1_plot <- function(data, variable, ylabel){
  
  summary_df <- pc1_summary(data, variable)
  
  ggplot(data) +
    ggforce::geom_sina(
      aes(
        x = SampleOrigin,
        y = .data[[variable]],
        colour = SampleOrigin
      ),
      size = 2,
      show.legend = FALSE
    ) +
    geom_point(
      data = summary_df,
      aes(x = SampleOrigin, y = mean),
      colour = "black",
      size = 3
    ) +
    geom_errorbar(
      data = summary_df,
      aes(
        x = SampleOrigin,
        ymin = mean - se,
        ymax = mean + se
      ),
      width = 0.3
    ) +
    theme_bw(base_size = 20) +
    scale_colour_manual(values = origin_colours) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(
      x = "",
      y = ylabel
    )
}

## generate figures -----------------------------

PC1_Bray_fig <- pc1_plot(
  ps_meta,
  "PC1_Bray",
  "Bray-Curtis (PC1)"
)

PC1_Wunif_fig <- pc1_plot(
  ps_meta,
  "PC1_Wunif",
  "Weighted UniFrac (PC1)"
)

PC1_Unif_fig <- pc1_plot(
  ps_meta,
  "PC1_Unif",
  "Unweighted UniFrac (PC1)"
)

PC1_Jaccard_fig <- pc1_plot(
  ps_meta,
  "PC1_Jaccard",
  "Jaccard (PC1)"
)


## combine plots
PC1_Bray_fig | PC1_Wunif_fig | PC1_Unif_fig | PC1_Jaccard_fig


## set reference levels for modelling
ps_meta$SampleOrigin <- factor(ps_meta$SampleOrigin, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))
ps_meta$SeqPlate <- factor(ps_meta$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))

## create mouse ID for random effect
ps_meta <- tidyr::unite(ps_meta, "MouseId", Id, SampleDetailed, remove = FALSE)

# bray curtis model
bray_pc1_glmm <- glmmTMB(PC1_Bray ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta)
summary(bray_pc1_glmm)

emmeans(bray_pc1_glmm, "SampleOrigin") %>% contrast("pairwise", adjust="Tukey")

# jaccard model
jac_pc1_glmm <- glmmTMB(PC1_Jaccard ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta)
summary(jac_pc1_glmm)

emmeans(jac_pc1_glmm, "SampleOrigin") %>% contrast("pairwise", adjust="Tukey")

# unifrac model
unif_pc1_glmm <- glmmTMB(PC1_Unif ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta)
summary(unif_pc1_glmm)

emmeans(unif_pc1_glmm, "SampleOrigin") %>% contrast("pairwise", adjust="Tukey")

# wunifrac model
wunif_pc1_glmm <- glmmTMB(PC1_Wunif ~ SampleOrigin + Sex + (1|SeqPlate) + (1|MouseId), data = ps_meta)
summary(wunif_pc1_glmm)

emmeans(wunif_pc1_glmm, "SampleOrigin") %>% contrast("pairwise", adjust="Tukey")

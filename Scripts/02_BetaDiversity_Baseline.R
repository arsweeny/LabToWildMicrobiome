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
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0, size = 16),
        legend.position = "top", legend.title = element_blank()) +
 stat_ellipse() + 
  ggtitle('Bray-Curtis')
All_Bray

## plot Weighted Unifrac ordination
All_Wunif <- plot_ordination(ps_comp, ps_wunif, color="SampleOrigin") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) + 
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0, size = 16),
        legend.position = "top", legend.title = element_blank()) +
  stat_ellipse() + 
  ggtitle('Weighted Unifrac')
All_Wunif

## plot Weighted Unifrac ordination
All_Unif <- plot_ordination(ps_comp, ps_unif, color="SampleOrigin") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) + 
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0, size = 16),
        legend.position = "top", legend.title = element_blank()) +
  stat_ellipse() + 
  ggtitle('Unweighted Unifrac')
All_Unif

## plot jaccard ordination 
All_Jaccard <- plot_ordination(ps_comp, ps_jaccard, color = "SampleOrigin") +
  geom_point(size = 2) +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD", "#FB8072")) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  stat_ellipse() +
  ggtitle("Jaccard")
All_Jaccard

## plot all ordinations together
beta_pan <- (All_Bray|All_Jaccard|All_Wunif|All_Unif) + plot_layout(guides = 'collect') & theme(legend.position = "top")

library(ggpubr)

## remove individual legends
All_Bray_clean <- All_Bray + 
  theme(legend.position = "none")
All_Jaccard_clean <- All_Jaccard + 
  theme(legend.position = "none")
All_Wunif_clean <- All_Wunif + 
  theme(legend.position = "none")
All_Unif_clean <- All_Unif + 
  theme(legend.position = "none")

## arrange plots
beta_pan <- ggarrange(
  All_Bray_clean,
  All_Jaccard_clean,
  All_Wunif_clean,
  All_Unif_clean,
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  font.label = list(size = 16, face = "bold"),
  common.legend = TRUE,
  legend = "top",
  align = "hv"
)

## save final plot
pdf("Output/Graphs/beta_diversity.pdf", width = 10, height = 10)
print(beta_pan)
dev.off()



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

adonis_bray_all <- adonis2(
  bray_ps ~ SampleOrigin + Sex + SeqPlate,
  data = ps_meta,
  strata = ps_meta$mouse_ID_dataset,
  permutations = 9999,
  by = "margin"
)
adonis_bray_all

pairwise_bray <- pairwise.adonis2(
  bray_ps ~ SampleOrigin + Sex + SeqPlate,
  data = ps_meta,
  by = "margin",
  strata = "mouse_ID_dataset"
)
pairwise_bray

bd <- betadisper(bray_ps, ps_meta$SampleOrigin)
anova(betadisper(bray_ps, ps_meta$SampleOrigin))
permutest(bd, pairwise = TRUE)
plot(bd)

# WEIGHTED UNIFRAC - PERMANOVA
wunif_ps <- phyloseq::distance(physeq = ps_comp, method = "wunifrac")
set.seed(999)

adonis_wunif_all <- adonis2(
  wunif_ps ~ SampleOrigin + Sex + SeqPlate,
  data = ps_meta,
  strata = ps_meta$mouse_ID_dataset,
  permutations = 9999,
  by = "margin"
)
adonis_wunif_all

pairwise_wunif <- pairwise.adonis2(
  wunif_ps ~ SampleOrigin + Sex + SeqPlate,
  data = ps_meta,
  by = "margin",
  strata = "mouse_ID_dataset"
)
pairwise_wunif

bd <- betadisper(wunif_ps, ps_meta$SampleOrigin)
anova(betadisper(wunif_ps, ps_meta$SampleOrigin))
permutest(bd, pairwise = TRUE)
plot(bd)

# UNWEIGHTED UNIFRAC - PERMANOVA
unif_ps <- phyloseq::distance(physeq = ps_comp, method = "unifrac")
set.seed(999)

adonis_unif_all <- adonis2(
  unif_ps ~ SampleOrigin + Sex + SeqPlate,
  data = ps_meta,
  strata = ps_meta$mouse_ID_dataset,
  permutations = 9999,
  by = "margin"
)
adonis_unif_all

pairwise_unif <- pairwise.adonis2(
  unif_ps ~ SampleOrigin + Sex + SeqPlate,
  data = ps_meta,
  by = "margin",
  strata = "mouse_ID_dataset"
)
pairwise_unif

bd <- betadisper(unif_ps, ps_meta$SampleOrigin)
anova(betadisper(unif_ps, ps_meta$SampleOrigin))
permutest(bd, pairwise = TRUE)
plot(bd)

# JACCARD - PERMANOVA
jaccard_ps <- phyloseq::distance(physeq = ps_comp, method = "jaccard", binary = TRUE)
set.seed(999)

adonis_jaccard_all <- adonis2(
  jaccard_ps ~ SampleOrigin + Sex + SeqPlate,
  data = ps_meta,
  strata = ps_meta$mouse_ID_dataset,
  permutations = 9999,
  by = "margin"
)
adonis_jaccard_all

pairwise_jaccard <- pairwise.adonis2(
  jaccard_ps ~ SampleOrigin + Sex + SeqPlate,
  data = ps_meta,
  by = "margin",
  strata = "mouse_ID_dataset"
)
pairwise_jaccard

bd <- betadisper(jaccard_ps, ps_meta$SampleOrigin)
anova(betadisper(jaccard_ps, ps_meta$SampleOrigin))
permutest(bd, pairwise = TRUE)
plot(bd)




# extract PERMANOVA and beta dispersion results ---------------------------

## convert PERMANOVA results to dataframes and add fixed effect names
add_term_column <- function(x){
  x <- as.data.frame(x)
  x$Fixed_effect <- rownames(x)
  rownames(x) <- NULL
  x <- x[, c("Fixed_effect", setdiff(names(x), "Fixed_effect"))]
  return(x)}

## global PERMANOVA results
bray_adonis_df <- add_term_column(adonis_bray_all)
wunif_adonis_df <- add_term_column(adonis_wunif_all)
unif_adonis_df <- add_term_column(adonis_unif_all)
jaccard_adonis_df <- add_term_column(adonis_jaccard_all)

## extract pairwise PERMANOVA results
extract_pairwise <- function(x){
  results <- list()
  for(i in names(x)){
  if(i != "parent_call"){
  temp <- as.data.frame(x[[i]])
  
  # Add fixed effect names
  temp$Fixed_effect <- rownames(temp)
  rownames(temp) <- NULL
  
  # Add comparison name
  temp$Comparison <- i
      
  # Put columns in logical order
  temp <- temp[, c(
  "Comparison",
  "Fixed_effect",
  setdiff(names(temp), c("Comparison", "Fixed_effect")))]
      
  results[[i]] <- temp}}

  do.call(rbind, results)}

bray_pairwise_df <- extract_pairwise(pairwise_bray)
wunif_pairwise_df <- extract_pairwise(pairwise_wunif)
unif_pairwise_df <- extract_pairwise(pairwise_unif)
jaccard_pairwise_df <- extract_pairwise(pairwise_jaccard)

## beta dispersion analysis
bd_bray <- betadisper(bray_ps, ps_meta$SampleOrigin)
bd_wunif <- betadisper(wunif_ps, ps_meta$SampleOrigin)
bd_unif <- betadisper(unif_ps, ps_meta$SampleOrigin)
bd_jaccard <- betadisper(jaccard_ps, ps_meta$SampleOrigin)

## extract dispersion ANOVA results
bray_disp_anova <- add_term_column(anova(bd_bray))
wunif_disp_anova <- add_term_column(anova(bd_wunif))
unif_disp_anova <- add_term_column(anova(bd_unif))
jaccard_disp_anova <- add_term_column(anova(bd_jaccard))

## extract dispersion permutation tests
bray_disp_perm_test <- permutest(bd_bray, permutations = 9999)
wunif_disp_perm_test <- permutest(bd_wunif, permutations = 9999)
unif_disp_perm_test <- permutest(bd_unif, permutations = 9999)
jaccard_disp_perm_test <- permutest(bd_jaccard, permutations = 9999)

## extract tables
bray_disp_perm <- add_term_column(
  as.data.frame(bray_disp_perm_test$tab)
)

wunif_disp_perm <- add_term_column(
  as.data.frame(wunif_disp_perm_test$tab)
)

unif_disp_perm <- add_term_column(
  as.data.frame(unif_disp_perm_test$tab)
)

jaccard_disp_perm <- add_term_column(
  as.data.frame(jaccard_disp_perm_test$tab)
)

## add distance metric labels
bray_adonis_df$Distance <- "Bray-Curtis"
wunif_adonis_df$Distance <- "Weighted UniFrac"
unif_adonis_df$Distance <- "Unweighted UniFrac"
jaccard_adonis_df$Distance <- "Jaccard"

## export Excel workbook
write.xlsx(
  
  list(
    
    # Bray-Curtis
    Bray_PERMANOVA = bray_adonis_df,
    Bray_pairwise = bray_pairwise_df,
    Bray_dispersion_ANOVA = bray_disp_anova,
    Bray_dispersion_perm = bray_disp_perm,
    
    
    # Weighted UniFrac
    WUniFrac_PERMANOVA = wunif_adonis_df,
    WUniFrac_pairwise = wunif_pairwise_df,
    WUniFrac_dispersion_ANOVA = wunif_disp_anova,
    WUniFrac_dispersion_perm = wunif_disp_perm,
    
    
    # Unweighted UniFrac
    UniFrac_PERMANOVA = unif_adonis_df,
    UniFrac_pairwise = unif_pairwise_df,
    UniFrac_dispersion_ANOVA = unif_disp_anova,
    UniFrac_dispersion_perm = unif_disp_perm,
    
    
    # Jaccard
    Jaccard_PERMANOVA = jaccard_adonis_df,
    Jaccard_pairwise = jaccard_pairwise_df,
    Jaccard_dispersion_ANOVA = jaccard_disp_anova,
    Jaccard_dispersion_perm = jaccard_disp_perm
    
  ),
  
  file = "Output/Beta_PERMANOVA_results.xlsx",
  
  overwrite = TRUE
)
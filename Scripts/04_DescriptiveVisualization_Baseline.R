rm(list=ls())

library(VennDiagram)
library(UpSetR)
library(ggVennDiagram)

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

# relative abundance visualization ----------------------------------------

## family level
ps_family <- ps_all %>% tax_fix(unknowns = c("Unclassified_Family"))

CFamily <- ps_family %>%
  tax_agg("Family") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>%
  comp_barplot(
    tax_level = "Family",
    facet_by = "SampleOrigin",
    n_taxa = 15
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
CFamily

## genus level
ps_genus <- ps_all %>% tax_fix(unknowns = c("Unclassified_Genus"))

CGenus <- ps_genus %>%
  tax_agg("Genus") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>%
  comp_barplot(
    tax_level = "Genus",
    facet_by = "SampleOrigin",
    n_taxa = 15
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
CGenus

## phylum level
ps_phylum <- ps_all %>% tax_fix(unknowns = c("Unclassified_Phylum"))

CPhylum <- ps_phylum %>%
  tax_agg("Phylum") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>%
  comp_barplot(
    tax_level = "Phylum",
    facet_by = "SampleOrigin",
    n_taxa = 15
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
CPhylum


# overlap of ASVs -------------------------------------------------

## make ASV abundance table
otu <- as.data.frame(otu_table(ps_all))

# transpose if ASVs are rows
if(taxa_are_rows(ps_all)){otu <- t(otu)}
otu <- as.data.frame(otu)

# define presence threshold
asv_presence <- otu >= 5

metadata <- data.frame(sample_data(ps_all))

metadata$SampleOrigin <- factor(metadata$SampleOrigin, 
    levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))

## create ASV lists for each group
asv_sets <- lapply(
  levels(metadata$SampleOrigin),
  function(group){
    samples <- rownames(metadata)[metadata$SampleOrigin == group]
    present <- colSums(asv_presence[samples,]) > 0
    names(present[present])
  }
)

names(asv_sets) <- levels(metadata$SampleOrigin)
asv_sets <- asv_sets[c("Wild:As", "Lab:Mm", "Wild-like:As", "Lab-like:As")]

## upset plot
library(grid)
library(gridGraphics)
library(ggplotify)

UpS <- upset(
  fromList(asv_sets),
  sets = c(
    "Wild:As",
    "Wild-like:As",
    "Lab-like:As",
    "Lab:Mm"
  ),
  keep.order = T,
  order.by = "freq",
  decreasing = T
)

UpS_grob <- grid.grabExpr({
  print(UpS)
}, wrap.grobs = TRUE)

UpS_gg <- as.ggplot(UpS_grob)

## venn diagram
venn <- venn.diagram(
  x = asv_sets,
  fill = c("#74C476", "#FB8072", "#6BAED6", "#BC80BD"),
  alpha = 0.4,
  filename = NULL)

## put them together in multi-panel figure
library(ggpubr)

visual <- ggarrange(
  venn,
  UpS_gg,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B"),
  font.label = list(size = 16, face = "bold"),
  widths = c(1, 1.4)
)

## save final plot
pdf("Output/Graphs/ASV_overlaps.pdf", width = 10, height = 5)
print(visual)
dev.off()

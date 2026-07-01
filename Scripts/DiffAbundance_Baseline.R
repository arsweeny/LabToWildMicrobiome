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
library(microbiomeMarker)
library(patchwork)

source("/Volumes/Rowan_HD2/Chapter2-3_Analysis/VisualisationScripts/DataPrep_Baseline.R") # runs data prep & sets up all the data subsets 


# make trees for each ps in case useful  ----------------------------------
# library(ape)
# #Phylogenetic tree for all faecal samples
# tree_faecal = rtree(ntaxa(ps_faecal), rooted=TRUE, 
#                     tip.label=taxa_names(ps_faecal))
# ps_faecal <- merge_phyloseq(ps_faecal, tree_faecal)
# 
# #Tree for colony only faecal samples
# tree_col = rtree(ntaxa(ps_col), rooted=TRUE, 
#                     tip.label=taxa_names(ps_col))
# ps_col <- merge_phyloseq(ps_col, tree_col)
# 
# #Tree for wild samples 
# tree_wild = rtree(ntaxa(ps_wild), rooted=TRUE, 
#                  tip.label=taxa_names(ps_wild))
# ps_wild <- merge_phyloseq(ps_wild, tree_wild)
# 
# #-----------------------------------------------------------------------------------------------------------------------------------------


#ALDEX 2 diff abundance 
#For now can only do differential abundance as pairwise analysis 
#Compare all groups to 'baseline' = lab-Mm

#Add a sample origin column 
#Hewan, Penicuik & Falkirk samples all called 'wild'
# MF1 & CD1 mice = lab-Mm
#ash54 wm colony = lab-like & b17 wm colony = wild-like
ps_all <- ps_all %>% 
  ps_mutate(SampleOrigin=case_when(str_detect(SampleDetailed, "MF1 lab mice") ~ "Lab-Mm",
                                   str_detect(SampleDetailed, "Pen") ~ "Wild-As", 
                                   str_detect(SampleDetailed, "Hew") ~ "Wild-As",
                                   str_detect(SampleDetailed, "Fal") ~ "Wild-As",
                                   str_detect(MouseOrigin, "Ash54 col") ~ "Lab-like-As", 
                                   str_detect(MouseOrigin, "b17 col") ~ "Wild-like-As",
                                   str_detect(MouseOrigin, "charles riv") ~ "Lab-Mm"))
#"#74C476", "#6BAED6", "#BC80BD","#FB8072"

#Lab-Mm vs Wild-As
ps_labVwild <- ps_all %>% 
  ps_filter(SampleOrigin %in% c("Lab-Mm", "Wild-As"))

microbiomeMarker::run_aldex(ps_labVwild, 
                            group="SampleOrigin",
                            taxa_rank = "Genus",
                            transform = "identity") -> aldexOut_labVwild

aldexOut_labVwild@marker_table %>% as.matrix %>% as.data.frame() %>%  
  mutate(pvalue=as.numeric(pvalue), 
         Sig=dplyr::case_when(pvalue < 0.001 ~ "***", pvalue < 0.01 ~ "**", pvalue <0.054 ~"*",
                              pvalue > 0.055 ~ "")) -> aldexResult_labVwild
#88 observations - but filter out non genus level and unknowns 

#Remove any taxa not mapped at genus level and any unknowns - 61 observations 
aldexResult_labVwild_Genus <- aldexResult_labVwild %>% filter(!str_detect(feature, "Family"),
                                                            !str_detect(feature, "Kingdom"),
                                                            !str_detect(feature, "Phylum"),
                                                            !str_detect(feature, "Order"),
                                                            !str_detect(feature, "Class"),
                                                            !str_detect(feature, "Incertae Sedis"))
#plot lab vs wild Genus 
p1<-aldexResult_labVwild_Genus %>% 
  mutate(ef_aldex=as.numeric(ef_aldex)) %>% 
  ggplot(., aes(x=reorder(feature, ef_aldex),
                y=ef_aldex, 
                fill=enrich_group)) + 
  geom_bar(stat="identity", position="dodge") + 
  coord_flip() + 
  geom_hline(yintercept=0) + 
  theme_bw(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top") +
  geom_text(aes(label=Sig, y=4.0), size=8) + 
  labs(fill="", x="Genus", y="Aldex Effect") + 
  scale_fill_manual(values = c("#FB8072","#74C476")) +
  scale_x_discrete(labels=c("[Eubacterium] xylanophilum group"=expression("Eubacterium xylanophilum group"),
                            "[Eubacterium] siraeum group"=expression("Eubacterium siraeum group"),
                            "[Eubacterium] brachy group"=expression("Eubacterium brachy group"),
                            "[Eubacterium] nodatum group"=expression("Eubacterium nodatum group"), parse=TRUE)) 

#"#74C476", "#6BAED6", "#BC80BD","#FB8072"

#Lab-Mm vs Wild-like-As
ps_labVwildlike <- ps_all %>% 
  ps_filter(SampleOrigin %in% c("Lab-Mm", "Wild-like-As"))

microbiomeMarker::run_aldex(ps_labVwildlike, 
                            group="SampleOrigin",
                            taxa_rank = "Genus",
                            transform = "identity") -> aldexOut_labVwildlike

aldexOut_labVwildlike@marker_table %>% as.matrix %>% as.data.frame() %>%  
  mutate(pvalue=as.numeric(pvalue), 
         Sig=dplyr::case_when(pvalue < 0.001 ~ "***", pvalue < 0.01 ~ "**", pvalue <0.054 ~"*",
                              pvalue > 0.055 ~ "")) -> aldexResult_labVwildlike
#87 observations - but filter out non genus level and unknowns 

#Remove any taxa not mapped at genus level and any unknowns - 63 observations 
aldexResult_labVwildlike_Genus <- aldexResult_labVwildlike %>% filter(!str_detect(feature, "Family"),
                                                              !str_detect(feature, "Kingdom"),
                                                              !str_detect(feature, "Phylum"),
                                                              !str_detect(feature, "Order"),
                                                              !str_detect(feature, "Class"),
                                                              !str_detect(feature, "Incertae Sedis"))
#plot lab vs wild Genus 
p2<-aldexResult_labVwildlike_Genus %>% 
  mutate(ef_aldex=as.numeric(ef_aldex)) %>% 
  ggplot(., aes(x=reorder(feature, ef_aldex),
                y=ef_aldex, 
                fill=enrich_group)) + 
  geom_bar(stat="identity", position="dodge") + 
  coord_flip() + 
  geom_hline(yintercept=0) + 
  theme_bw(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top") +
  geom_text(aes(label=Sig, y=4.0), size=8) + 
  labs(fill="", x="Genus", y="Aldex Effect") + 
  scale_fill_manual(values = c("#FB8072","#6BAED6")) +
  scale_x_discrete(labels=c("[Eubacterium] xylanophilum group"=expression("Eubacterium xylanophilum group"),
                            "[Eubacterium] siraeum group"=expression("Eubacterium siraeum group"),
                            "[Eubacterium] brachy group"=expression("Eubacterium brachy group"), parse=TRUE)) 

#"#74C476", "#6BAED6", "#BC80BD","#FB8072"

#Lab-Mm vs Lab-like-As
ps_labVlablike <- ps_all %>% 
  ps_filter(SampleOrigin %in% c("Lab-Mm", "Lab-like-As")) %>%   
ps_mutate(SampleOrigin = factor(SampleOrigin, levels = c("Lab-Mm", "Lab-like-As"))) 


microbiomeMarker::run_aldex(ps_labVlablike, 
                            group="SampleOrigin",
                            taxa_rank = "Genus",
                            transform = "identity") -> aldexOut_labVlablike

aldexOut_labVlablike@marker_table %>% as.matrix %>% as.data.frame() %>%  
  mutate(pvalue=as.numeric(pvalue), 
         Sig=dplyr::case_when(pvalue < 0.001 ~ "***", pvalue < 0.01 ~ "**", pvalue <0.054 ~"*",
                              pvalue > 0.055 ~ "")) -> aldexResult_labVlablike
#75 observations - but filter out non genus level and unknowns 

#Remove any taxa not mapped at genus level and any unknowns - 57 observations 
aldexResult_labVlablike_Genus <- aldexResult_labVlablike %>% filter(!str_detect(feature, "Family"),
                                                                      !str_detect(feature, "Kingdom"),
                                                                      !str_detect(feature, "Phylum"),
                                                                      !str_detect(feature, "Order"),
                                                                      !str_detect(feature, "Class"),
                                                                      !str_detect(feature, "Incertae Sedis"))
#plot lab vs wild Genus 
p3<-aldexResult_labVlablike_Genus %>% 
  mutate(ef_aldex=as.numeric(ef_aldex)) %>% 
  ggplot(., aes(x=reorder(feature, ef_aldex),
                y=ef_aldex, 
                fill=enrich_group)) + 
  geom_bar(stat="identity", position="dodge") + 
  coord_flip() + 
  geom_hline(yintercept=0) + 
  theme_bw(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top") +
  geom_text(aes(label=Sig, y=4.0), size=8) + 
  labs(fill="", x="Genus", y="Aldex Effect") + 
  scale_fill_manual(values = c("#BC80BD","#FB8072")) +
  scale_x_discrete(labels=c("[Eubacterium] xylanophilum group"=expression("Eubacterium xylanophilum group"),
                            "[Eubacterium] siraeum group"=expression("Eubacterium siraeum group"),
                            "[Eubacterium] brachy group"=expression("Eubacterium brachy group"),
                            "[Eubacterium] nodatum group"=expression("Eubacterium nodatum group"), parse=TRUE)) 

#panel
(p1|p2|p3) + plot_layout(guides = 'collect')


#Plotting the overlap across all groups 
#combine the 3 sets of results 
aldexResult_All <- bind_rows(aldexResult_labVwild_Genus, aldexResult_labVwildlike_Genus, aldexResult_labVlablike_Genus) 

aldexResult_All %>% 
  mutate(enrich_group = factor(enrich_group, levels = c("Lab-Mm", "Lab-like-As", "Wild-like-As", "Wild-As"))) %>% 
  mutate(ef_aldex=as.numeric(ef_aldex)) %>% 
  ggplot(., aes(x=reorder(feature, ef_aldex),
                y=ef_aldex, 
                fill=enrich_group)) + 
  geom_bar(stat="identity", position="dodge") + 
  coord_flip() + 
  geom_hline(yintercept=0) + 
  theme_bw(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top") +
  geom_text(aes(label=Sig, y=4.0), size=8) + 
  labs(fill="", x="Genus", y="Aldex Effect") + 
  scale_fill_manual(values = c("#FB8072","#BC80BD","#6BAED6","#74C476"))


# overlapping 
Overlap_1<- intersect(aldexResult_labVwild_Genus$feature, aldexResult_labVwildlike_Genus$feature)
Overlap_2 <- intersect(aldexResult_labVwildlike_Genus$feature, aldexResult_labVlablike_Genus$feature)
Overlap_3<- intersect(aldexResult_labVwild_Genus$feature, aldexResult_labVlablike_Genus$feature)
Overlap_all <- c(Overlap_1, Overlap_2, Overlap_3) 

#not overlapping
Unique_1 <- setdiff(aldexResult_labVwild_Genus$feature, aldexResult_labVwildlike_Genus$feature) #15 taxa
Unique_2 <- setdiff(aldexResult_labVwildlike_Genus$feature, aldexResult_labVlablike_Genus$feature) #16 taxa
Unique_3<- setdiff(aldexResult_labVwild_Genus$feature, aldexResult_labVlablike_Genus$feature) #18 taxa


aldexResult_All %>% 
  filter(feature %in% Overlap_all) %>% 
  mutate(enrich_group = factor(enrich_group, levels = c("Lab-Mm", "Lab-like-As", "Wild-like-As", "Wild-As"))) %>% 
  mutate(ef_aldex=as.numeric(ef_aldex)) %>% 
  ggplot(., aes(x=reorder(feature, ef_aldex),
                y=ef_aldex, 
                fill=enrich_group)) + 
  geom_bar(stat="identity", position="dodge") + 
  coord_flip() + 
  geom_hline(yintercept=0) + 
  theme_bw(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top") +
  geom_text(aes(label=Sig, y=4.0), size=8) + 
  labs(fill="", x="Genus", y="Aldex Effect") + 
  scale_fill_manual(values = c("#FB8072","#BC80BD","#6BAED6","#74C476")) +
  scale_x_discrete(labels=c("[Eubacterium] xylanophilum group"=expression("Eubacterium xylanophilum group"),
                            "[Eubacterium] siraeum group"=expression("Eubacterium siraeum group"),
                            "[Eubacterium] brachy group"=expression("Eubacterium brachy group"),
                            "[Eubacterium] nodatum group"=expression("Eubacterium nodatum group"), parse=TRUE)) 







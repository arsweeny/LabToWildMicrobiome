rm(list=ls())

library(tidyverse)
library(phyloseq)
library(microbiome)
library(microViz)
library(vegan)
library(mia)
library(miaViz)
library(metagMisc)
library(magrittr)
library(ape)
library(patchwork)
library(DESeq2)

detach(package:dplyr)
library(dplyr) 


source("VisualisationScripts/DataPrep_Baseline.R")

#adds sampleGroup column with name of each mouse cohort as wanted
ps_all <- ps_all %>% 
  ps_mutate(SampleGroup=
              case_when(str_detect(SampleOrigin, "MF1|CD1") ~ "Lab:Mm",
                      str_detect(SampleOrigin, "Wild mice") ~ "Wild:As",
                      str_detect(SampleOrigin, "Lab-like") ~ "Lab-like:As",
                      str_detect(SampleOrigin, "Wild-like") ~ "Wild-like:As"))

#pulls & re-orders sample groups to how we want them & assigns to obeject Cohorts
meta(ps_all) %>% pull(SampleGroup) %>% unique -> Cohorts
Cohorts <- Cohorts[c(1,3,4,2)]
# [1] "Wild:As"      "Wild-like:As" "Lab-like:As"  "Lab:Mm"  

# Comparison DF -----------------------------------------------------------

#makes a list of the phyloseq for each sample group/mouse cohort  
#makes column called Comparison, for cohort 1-3 (all except lam-Mm), where every group is compared to lab-Mm
CohortPS <- lapply(Cohorts[1:3], function(x){
  ps <- ps_all %>% 
    ps_filter(SampleGroup %in% c(x, "Lab:Mm")) %>% 
    ps_mutate(Comparison=paste(x, "Lab:Mm", sep="-")) %>% 
    #ps_filter(SampleGroup==x) %>% 
  return(ps)
})

# can access elements of the list CohortPS[[1]] 
# check the subsetting has worked: CohortPS[[1]] %>% meta %>% pull(SampleGroup)


# Sample plots ------------------------------------------------------------
#number of samples in each cohort 
ps_all %>% meta %>%  
  mutate(SampleGroup=fct_relevel(SampleGroup, Cohorts)) %>% 
  group_by(SampleGroup) %>% count %>% 
  ggplot(aes(x=SampleGroup, y=n)) + 
  geom_bar(stat="identity") + 
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  ggtitle("Samples by origin")

# make a list of metadata & then merge
MouseMeta <- lapply(CohortPS, meta) 
MouseMeta <- bind_rows(MouseMeta)


### probably not relevant but a way to do stuff across the list of ps 
# CohortPS <- lapply(CohortPS, function(x){
#   x %>% ps_mutate(MouseGroupSurv=
#               case_when(SampleId %in% NeoInclude & Id %in% Died ~ "LambsDied",
#                         SampleId %in% NeoInclude & Id %in% Survived ~ "LambsSurvived",
#                         SampleId %in% AdultInclude ~ "Adults"))
# })
# lapply(CohortPS[[c(1,4)]]) - another way to pull elements of a list 


# Ordinations -------------------------------------------------------------

# ps_clr_genus <- MouseNeo %>% 
#   ps_filter(SampleGroup %in% Cohorts[1:3]) %>% 
#   tax_transform(trans = "clr", rank = "Genus") 
# 
# ps_clr_genus %>% 
#   # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
#   ord_calc() %>% 
#   ord_plot(color = "Died", shape = "Died", size = 3) +
#   scale_colour_brewer(palette = "Dark2") + 
#   facet_grid(~Season)

ps_core <- core(ps_all, detection = 5, prevalence=0.1)
treeAll <- rtree(ntaxa(ps_all), rooted=TRUE, tip.label = taxa_names(ps_all))
coreAll <- merge_phyloseq(ps_core, treeAll)
compAll <- microbiome::transform(coreAll, "compositional")
brayAll <- ordinate(compAll, "PCoA", "bray")
unifAll <- ordinate(compAll, "PCoA", "unifrac")
wunifAll <- ordinate(compAll, "PCoA", "wunifrac")
metaAll <-meta(coreAll) %>% as.data.frame() %>%
  mutate(PC1_Bray=brayAll$vectors[,1],
         PC1_Unif=unifAll$vectors[,1],
         PC1_Wunif=wunifAll$vectors[,1])

# Doing a quick check of variance in composition within groups 
#bray
bray_ps <- phyloseq::distance(physeq = compAll, method = "bray")
all_bray_disper <- betadisper(bray_ps, metaAll$SampleGroup)
permutest(all_bray_disper, pairwise=TRUE, permutations=1000)

all_bray_disper$distances %>% as.data.frame() %>% 
  rename(DistanceToCentroid=1) %>% 
  rownames_to_column(var="SampleId") %>% 
  left_join(metaAll, .) -> BrayDispDF 

#wunif
wunif_ps <- phyloseq::distance(physeq = compAll, method = "wunifrac")
all_wunif_disper <- betadisper(wunif_ps, metaAll$SampleGroup)
permutest(all_wunif_disper, pairwise=TRUE, permutations=1000)

all_wunif_disper$distances %>% as.data.frame() %>% 
  rename(DistanceToCentroid=1) %>% 
  rownames_to_column(var="SampleId") %>% 
  left_join(metaAll, .) -> WunifDispDF

#unif
unif_ps <- phyloseq::distance(physeq = compAll, method = "unifrac")
all_unif_disper <- betadisper(unif_ps, metaAll$SampleGroup)
permutest(all_unif_disper, pairwise=TRUE, permutations=1000)

all_unif_disper$distances %>% as.data.frame() %>% 
  rename(DistanceToCentroid=1) %>% 
  rownames_to_column(var="SampleId") %>% 
  left_join(metaAll, .) -> UnifDispDF

#set levels for plotting 
BrayDispDF$SampleGroup <- factor(BrayDispDF$SampleGroup, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))
WunifDispDF$SampleGroup <- factor(WunifDispDF$SampleGroup, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))
UnifDispDF$SampleGroup <- factor(UnifDispDF$SampleGroup, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))


#Bray Plot distance to centroid for each cohort - beta dispersion plots 
ggregplot::SinaGraph(BrayDispDF, "SampleGroup", "DistanceToCentroid") + 
  scale_x_discrete(limits=Cohorts) + theme_bw(base_size=18) + theme(legend.position = "none") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  labs(x="", y= "Distance To Centroid") +
  ggtitle('Bray-Curtis') -> bray_disp


#Wunif beta dispersion plot
ggregplot::SinaGraph(WunifDispDF, "SampleGroup", "DistanceToCentroid") + 
  scale_x_discrete(limits=Cohorts) + theme_bw(base_size=18) + theme(legend.position = "none") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  labs(x="", y= "Distance To Centroid") +
  ggtitle("Weighted Unifrac") -> wunif_disp

#Unif beta dispersion plot
ggregplot::SinaGraph(UnifDispDF, "SampleGroup", "DistanceToCentroid") + 
  scale_x_discrete(limits=Cohorts) + theme_bw(base_size=18) + theme(legend.position = "none") +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  labs(x="", y= "Distance To Centroid") + 
  ggtitle("Unweighted Unifrac") -> unif_disp 

#panel fingure
bray_disp|wunif_disp|unif_disp

#Beta diversity plots on sample groups all together 
#first define what order you want the samples to be colored in plots 
sample_data(compAll)$SampleGroup <- factor(sample_data(compAll)$SampleGroup, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))

plot_ordination(compAll, brayAll, color="SampleGroup") +
  geom_point(aes(), size= 3) +
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072")) +
  theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 18),
        legend.position = "top", legend.title = element_blank())+
  stat_ellipse() +
  ggtitle('Bray-Curtis')

#permanova on bray-curtis 
#re-order factors to set Lab-Mm as intercept and plate 8 as intercept 
sample_data(compAll)$SampleGroup <- factor(sample_data(compAll)$SampleGroup, levels = c("Lab:Mm", "Lab-like:As", "Wild-like:As","Wild:As"))
#sample_data(compAll)$SeqPlate <- factor(sample_data(compAll)$SeqPlate, levels = c("8","7","6","5","4","3","2","1"))

metadf<- data.frame(sample_data(compAll))
bray_ps <- phyloseq::distance(physeq = compAll, method = "bray")
set.seed(999)
# Adonis test
adonis1 <- adonis2(bray_ps ~ SampleGroup  +Sex + SeqPlate, by="onedf",
                   data = metadf)
adonis1

#PAIRWISE GROUP GROUPS
#Each cohort compared seperately to Lab-Mm 
#Beta Diversity on pairwise group analysis - Lab-Mm baseline to compare all to 
#Takes list of all 4 phyloseqs and prunes the rare taxa for each ps individually
#Then saves as new list called cores
cores <- lapply(CohortPS, function(x){
  ps_core <- core(x, detection =5, prevalence = .1) %>% 
  return(ps_core)
})

#adds phylogenetic tree to each ps (cores)
trees <- lapply(cores, function(x){
  random_tree = rtree(ntaxa(x), rooted=TRUE, 
                      tip.label=taxa_names(x))
})

#merges trees and phyloseqs together into ps_cores
ps_cores <- lapply(1:length(cores), function(x){
  merge_phyloseq(cores[[x]], trees[[x]])
})

#transforms microbiome data from counts in ps_cores to compositional in ps_comps
ps_comps <- lapply(ps_cores, microbiome::transform, "compositional")

#calculate the different ordination distances, for bray and both unifracs 
ps_brays <- lapply(ps_comps, function(x){
  ordinate(x, "PCoA", "bray")
})

ps_unifs <- lapply(ps_comps, function(x){
  ordinate(x, "PCoA", "unifrac")
})
                   
ps_wunifs <- lapply(ps_comps, function(x){
  ordinate(x, "PCoA", "wunifrac")
})


## check out the ordination output 
#ensure columns added for each ordination distance 
ps_metas <- lapply(1:length(ps_cores), function(x){
    meta<- meta(ps_cores[[x]]) %>% as.data.frame() %>% 
    mutate(PC1_Bray=ps_brays[[x]]$vectors[,1],
           PC1_Unif=ps_unifs[[x]]$vectors[,1],
           PC1_Wunif=ps_wunifs[[x]]$vectors[,1])
    return(meta)
})

# PAIRWISE Distances plot
#plots PC1-Bray for each cohort compared to Lab-Mm
MouseGroups <- lapply(1:length(ps_metas), function(x){
  ggregplot::SinaGraph(ps_metas[[x]] %>% as.data.frame(),
            "SampleGroup", "PC1_Bray") + theme_bw(base_size = 24) + 
    scale_colour_brewer(palette = "Set2") + 
    ggtitle(Cohorts[x]) + 
    theme(legend.position = "none") + 
    ggpubr::stat_compare_means()
})
(MouseGroups[[1]] + MouseGroups[[2]] + MouseGroups[[3]]) 
#think could maybe amend each figure directly to change appearance by MouseGroups[[1]] + ggplot comands 

## PAIRWISE Ordinations 
pal <- RColorBrewer::brewer.pal(6, "Spectral")
#Plots unweighted unifrac ordinations for each cohort compared to Lab-Mm 
MouseGroupsOrd <- 
  lapply(1:length(ps_metas), function(x){
  plot_ordination(ps_comps[[x]], 
                ps_unifs[[x]], 
                color="SampleGroup") + 
  geom_point(aes(), size= 3) + 
  scale_colour_manual(values = pal) + 
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0, size = 18),
        legend.position = "bottom")+
  stat_ellipse() + 
  ggtitle(Cohorts[x])
}) 
(MouseGroupsOrd[[1]] + MouseGroupsOrd[[2]] + MouseGroupsOrd[[3]])

#PAIRWISE DISTANCES 
# Calculate and plot group distances 
LabIds <- ps_all %>% meta %>% filter(SampleGroup=="Lab:Mm") %>% pull(SampleId) 
#creates list of all Lab-Mm mice only 

#calculates bray distances for all mouse cohorts except Lab-Mm
ps_brays_comp <- lapply(ps_brays, function(x){
  df <- x$vectors %>% as.data.frame %>% 
    filter(!row.names(.) %in% LabIds) %>% as.matrix ## filter by cohort instead 
})
#calculates bray distances for Lab-Mm only 
ps_brays_lab <- lapply(ps_brays, function(x){
  df <- x$vectors %>% as.data.frame %>% 
    filter(row.names(.) %in% LabIds) %>% as.matrix 
})

#calculates the lab-Mm centroids for each cohort (so end up with list of 3x and 3y values, 1 for each mouse cohort)
LabCentroids <- # lab = lab and all others = comp, lab here is baseline 
  lapply(ps_brays_lab, function(i){
    bray <- i
    x <- mean(bray[,1])
    y <- mean(bray[,2])
    centroid <- c(x,y) 
    return(centroid)
  })

#not really sure what these functions are for/do? 
euclidean <- function(a, b) sqrt(sum((a - b)^2))
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2)) #x1 and x2 have to be matrices 

#need to go through this slowly to figure out what doing, quite complex 
MouseDistancesBray <- lapply(1:length(ps_brays_comp), function(i){
    x1 <- LabCentroids[[i]] %>% as.matrix() %>% t 
    CompPoints <- ps_brays_comp[[i]][,1:2]
    df <- CompPoints %>% as.data.frame() %>% 
      rownames_to_column(var="SampleId") %>% 
      dplyr::select(SampleId) %>% mutate(DistanceToLab="")
    Distances <- lapply(1:dim(CompPoints)[1], function(j){
      x2 <- CompPoints[j,] %>% as.matrix() %>% t 
      euc.dist(x1, x2)
      #df[j,2] <- Distance %>% round(digits=3)
    }) %>% unlist 
    df$DistanceToLabBray <- Distances
    return(df)
  })
MouseDistancesBray %<>% bind_rows
#Essentially creates a dataframe, where the "DistanceToLabBray" has been calculated for each sample


#need to go through this slowly to figure out what doing, quite complex 
MouseDiffPC1 <- lapply(1:length(ps_brays_comp), function(i){
  LabMean <- ps_metas[[i]] %>% # LabMean 
    filter(SampleGroup=="Lab:Mm") %>% 
    pull(PC1_Bray) %>% mean # mean PC1 bray for lab 
  CompPC1 <- ps_metas[[i]] %>% 
    filter(SampleGroup!="Lab:Mm") %>% 
    dplyr::select(PC1_Bray) 
  df <- CompPC1 %>% as.data.frame() %>% 
    rownames_to_column(var="SampleId") %>% 
    dplyr::select(SampleId) %>% mutate(DiffFromLabBray="")
  Differences <- lapply(1:dim(CompPC1)[1], function(j){
    abs(LabMean - CompPC1[j,])
    #df[j,2] <- Distance %>% round(digits=3)
  }) %>% unlist 
  df$DiffFromLabBray <- Differences
  return(df)
})
MouseDiffPC1 %<>% bind_rows 
#makes another df, where "DiffFromLabBray" is calculated 

#binds MouseDistancesBray and MouseDiffPC1 into a dataframe 
MouseDF <- left_join(meta(ps_all), MouseDistancesBray) %>% 
  left_join(., MouseDiffPC1) %>% 
  mutate(SampleGroup=fct_relevel(SampleGroup, Cohorts)) 

#plots Bray PC1 difference of each cohort compared to lab-Mm mice 
SinaGraph(MouseDF, "SampleGroup", "DiffFromLabBray") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_colour_manual(values=pal) + 
  labs(y='Bray PC1 Difference', x="") -> PC1Diffs 

# save figure
# PC1Diffs + ggsave("Figures/BrayPC1Diff.png", units="mm", 
#                   height=170, width=190, dpi=450)


#GROUP DIFFERECES
#Whole data set compared to Lab-Mm as baseline 
# Group differences based on all group ordinations (not just pairwise comparing to Lab-Mm)

LabIds2 <- ps_all %>% meta %>% filter(SampleGroup=="Lab:Mm") %>% pull(SampleId) 
#creates list of all Lab-Mm mice only 

CompIds2 <- lapply(Cohorts[1:3], function(x){
  ps_all %>% meta %>% filter(SampleGroup==x) %>% 
    pull(SampleId)
})
#creates list of all cohorts except Lab-Mm 


#calculates bray distances for all mouse cohorts except Lab Mm 
ps_brays_comp2 <- lapply(CompIds2, function(x){
  df <- brayAll$vectors %>% as.data.frame %>% 
    filter(row.names(.) %in% x) %>% as.matrix ## filter by cohort instead 
})
#calculates bray distances for lab-Mm only 
ps_brays_lab2 <-  brayAll$vectors %>% as.data.frame %>% 
    filter(row.names(.) %in% LabIds2) %>% as.matrix 

xx <- mean(ps_brays_lab2[,1])
yy <- mean(ps_brays_lab2[,2])
LabCentroid2 <- c(xx,yy)  
#calculates the lab-Mm centroids entire data set, so end up with 2 values (xx and yy)

euclidean2 <- function(a, b) sqrt(sum((a - b)^2))
euc.dist2 <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2)) #x1 and x2 have to be matrices 


MouseDistancesBray2 <- lapply(1:length(ps_brays_comp2), function(i){
  x1 <- LabCentroid2 %>% as.matrix() %>% t 
  CompPoints2 <- ps_brays_comp2[[i]][,1:2]
  df <- CompPoints2 %>% as.data.frame() %>% 
    rownames_to_column(var="SampleId") %>% 
    dplyr::select(SampleId) %>% mutate(DistanceToLab="")
  Distances2 <- lapply(1:dim(CompPoints2)[1], function(j){
    x2 <- CompPoints2[j,] %>% as.matrix() %>% t 
    euc.dist2(x1, x2)
    #df[j,2] <- Distance %>% round(digits=3)
  }) %>% unlist 
  df$DistanceToLabBray <- Distances2
  return(df)
})
MouseDistancesBray2 %<>% bind_rows
#Calculates distance to lab-Mm bray for whole dataset in 1 go 


MouseDiffPC1_2 <- lapply(1:length(ps_brays_comp2), function(i){
  LabMean <- metaAll %>% # LabMean 
    filter(SampleGroup=="Lab:Mm") %>% 
    pull(PC1_Bray) %>% mean # mean PC1 for lab 
  CompPC1 <- metaAll %>% 
    filter(SampleGroup==Cohorts[i]) %>% 
    dplyr::select(PC1_Bray) 
  df <- CompPC1 %>% as.data.frame() %>% 
    rownames_to_column(var="SampleId") %>% 
    dplyr::select(SampleId) %>% mutate(DiffFromLabBray="")
  Differences <- lapply(1:dim(CompPC1)[1], function(j){
    abs(LabMean - CompPC1[j,])
    #df[j,2] <- Distance %>% round(digits=3)
  }) %>% unlist 
  df$DiffFromLabBray <- Differences
  return(df)
})
MouseDiffPC1_2 %<>% bind_rows 
#Calculates distance from lab-Mm bray for whole dataset in 1 go 

#binds MouseDistancesBray2 and MouseDiffPC1_2 into a dataframe - can use this for plotting 
MouseDF2 <- left_join(meta(ps_all), MouseDistancesBray2) %>% 
  left_join(., MouseDiffPC1_2) %>% 
  mutate(SampleGroup=fct_relevel(SampleGroup, Cohorts)) 


# LambIndiv <- LambDF %>% # wouldnt worry about this so much 
#   filter(SampleGroup %in% Cohorts[1:3]) %>% 
#   group_by(Id, Died) %>% 
#   summarise(MeanDistance=mean(DistanceToAdultsBray),
#             VarDistance=plotrix::std.error(DistanceToAdultsBray),
#             MeanDiff=mean(DiffFromAdultsBray),
#             VarDiff=plotrix::std.error(DiffFromAdultsBray))

SinaGraph(MouseDF2, "SampleGroup", "DiffFromLabBray") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_colour_manual(values=pal) + 
  labs(y='Bray PC1 Difference from Lab', x="") -> PC1Diffs_2 

SinaGraph(MouseDF2, "SampleGroup", "DistanceToLabBray") + 
  theme_bw(base_size = 18) + 
  theme(legend.position = "none") + 
  scale_x_discrete(limits=Cohorts[1:3]) -> PC1Dist_2 

#for saving figs
# PC1Diffs + ggsave("Figures/BrayPC1Diff.png", units="mm", 
#                   height=170, width=190, dpi=450)


# DIFFERENTIAL ABUNDANCES 

# phy_pal <- wesanderson::wes_palette("Zissou1", 20, type = "continuous")
# phy_pal <- RColorBrewer::brewer.pal(9, "Spectral")
# ColoursRamp <- colorRampPalette(phy_pal)
# ColoursExpandPhy <- ColoursRamp(n=20)

#can't get this to run - issue seems with the phyloseq_filter_prevalence(.,) part
#drop the filter prevalence for now because seems to work without it 
detach(package:dplyr)
library(dplyr)

CohortPSPhy <- lapply(CohortPS, function(x){
  x %>% 
    #metagMisc::phyloseq_filter_prevalence(., #works if run without this
     #prev.trh = 0.1) %>% #works if run without this
    aggregate_taxa(., 'Phylum')
})
#if run without ps_filter_prev
#Creates 3 phyloseq's (1 for each cohort) with asvs at phylum level 
#[[1]] has 17 phylum(wild-As vs Lab-Mm), [[2]] has 12 (wild-like-As vs Lab-Mm) and [[3]] has 11 (lab-like-As vs Lab-Mm)


#WORKAROUND - to make sure we get the prevalence right without using ps_filter_prev  
#Re-make CohortPS but use ps_core2 instead of ps_all as original ps
#where we filter the prevalence to 10% (leave abundance@1) in ps_core2
ps_core2 <- core(ps_all, detection=1, prevalence=0.1)

#Then re-make cohorts(2) using ps_core2 
meta(ps_core2) %>% pull(SampleGroup) %>% unique -> Cohorts2
Cohorts2<- Cohorts2[c(1,3,4,2)]
# [1] "Wild:As"      "Wild-like:As" "Lab-like:As"  "Lab:Mm"   

#re-make cohortPS(2) using ps_core2
CohortPS2 <- lapply(Cohorts2[1:3], function(x){
  ps <- ps_core2 %>% 
    ps_filter(SampleGroup %in% c(x, "Lab:Mm")) %>% 
    ps_mutate(Comparison=paste(x, "Lab:Mm", sep="-")) %>% 
    #ps_filter(SampleGroup==x) %>% 
    return(ps)
})
#then try to make phylum phyloseqs (ignoring prevalence) using CohortPS2
CohortPSPhy2 <- lapply(CohortPS2, function(x){
  x %>%
    # phyloseq_filter_prevalence(., 
    # prev.trh = 0.1) %>% 
    aggregate_taxa(., 'Phylum')
})
#Creates 3 phyloseq's (1 for each pairwise cohort, lab-mm vs wild-as etc) with asvs at phylum level 
#lose some phyla per cohort compared to CohortPSPhy but maybe his is to be expected with prevalence threshold?
#[[1]] goes from 17 to 9 taxa, [[2]] goes from 12 to 10 and [[3]] goes from 11 to 10 
#check whats happening 

MousePhy <- ps_core2 %>% # again try on ps_core2 instead of ps_all to skip prevalence filter that is not working
  ps_filter(SampleGroup%in%Cohorts[1:3]) %>%  
  #phyloseq_filter_prevalence(., prev.trh = 0.1) %>% 
  aggregate_taxa(., 'Phylum')
#creates 1 phyloseq at phylumn level for all mice except Lab-Mm 
#total of 10 phylum in 196 samples (all samples except Lab-Mm)

#deseq doesn't allow zeroes, so OTU values +1 and make sample group a factor instead of chrt
#think loop means this happens for each object in CohortPSPhy2 list?
for(i in 1:length(CohortPSPhy2)){
  otu_table(CohortPSPhy2[[i]]) <- otu_table(CohortPSPhy2[[i]]) + 1  
  #sample_data(CohortPSPhy2[[i]])$SampleGroup <-as.factor(sample_data(CohortPSPhy2[[i]])$SampleGroup) 
  CohortPSPhy2[[i]] <- CohortPSPhy2[[i]] %>%
    ps_mutate(SampleGroup=fct_relevel(SampleGroup, "Lab:Mm"))
  }
#check worked -  SampleGroup at each ps now factor with 2 lavels 
sample_data(CohortPSPhy2[[3]])$SampleGroup %>% levels #double check Lab-Mm is first in each as comparing to this 
#[1] "Lab:Mm" "Wild:As", [2] "Lab:Mm" "Wild-like:As", [3] "Lab:Mm" "Lab-like:As"

ds_mouse <- list() 
dds_mouse <- list()

#calculating the differential abundance values and Wald test for each cohort 
for(i in 1:length(CohortPSPhy2)){
  ds_mouse[[i]] <- phyloseq_to_deseq2(CohortPSPhy2[[i]], ~ SampleGroup)
  dds_mouse[[i]] <- DESeq(ds_mouse[[i]], test="Wald", fitType="parametric")
}

ResultsPhy <- 
  lapply(1:length(dds_mouse), function(x){
    as.data.frame(results(dds_mouse[[x]])) %>% 
      rownames_to_column(var="Phylum") %>% 
      mutate(SampleGroup=Cohorts2[[x]]) #%>%
    #     mutate(Sig=case_when(padj<p_correctedP ~ "*",
    #                          padj>p_correctedP ~ "")) %>%
    #     mutate(Min=log2FoldChange-lfcSE,
    #            Max=log2FoldChange+lfcSE,
    #            Starloc=case_when(Sig=="*" ~ Max+(Max-Min),
    #                              Sig=="" ~ 0))
}) 

#adjust the p-value so that it is divided by the number of tests 
#then pull significant values out from this -  
#allows you to correct for multiple testing and is more conservative than just usng 0.05
0.05 / (bind_rows(ResultsPhy) %>%
         pull(Phylum) %>% unique %>% length) -> p_correctedP 

ResultsPhy <- 
  lapply(1:length(dds_mouse), function(x){
    as.data.frame(results(dds_mouse[[x]])) %>% 
    rownames_to_column(var="Phylum") %>% 
     mutate(SampleGroup=Cohorts2[[x]]) %>%
    mutate(Sig=case_when(padj<p_correctedP ~ "*",
                         padj>p_correctedP ~ "")) %>%
    mutate(Min=log2FoldChange-lfcSE,
           Max=log2FoldChange+lfcSE,
           Starloc=case_when(Sig=="*" ~ Max+(Max-Min),
                             Sig=="" ~ 0))
 }) 
ResultsPhyAll <- bind_rows(ResultsPhy)

write.csv(ResultsPhyAll, "DescriptiveTables/Phylum_DifAbd.csv")

#Set colours 
#colour palettes - need to expand as lots of taxa at genus level 
bind_rows(ResultsPhy) %>%
  filter(pvalue<p_correctedP) %>% 
  arrange(-log2FoldChange) %>%
  pull(Phylum) %>% unique -> PhyNames #10
 

phycols <-  c("#FB8072","#FDB462","#FFED6F","#B3DE69","#41AE76", 
              "#8DD3C7","#80B1D3","#BEBADA","#BC80BD","#FA9FB5")
names(phycols) <- PhyNames

phy_rnbw <- c("#FB8072","#FDB462","#FFED6F","#B3DE69","#41AE76",
              "#8DD3C7","#80B1D3","#BEBADA","#BC80BD", "#FA9FB5")

#ALL PHYLUM PLOT
#all 3 plots compared to lab-Mm at Phylum level - x axis ordered by magnitude of the log change
lapply(1:length(ResultsPhy), function(x){
  ResultsPhy[[x]] %>%
  filter(pvalue<p_correctedP) %>% 
  ggplot(aes(x=reorder(Phylum, -log2FoldChange), #just out x=Phylum to order alphabetically 
             y=log2FoldChange,
             colour= Phylum,
             label="Sig")) + 
  geom_point(size=4)+
  geom_errorbar(aes(x=Phylum, ymax=log2FoldChange + lfcSE, ymin=log2FoldChange - lfcSE), # error bar + or - the standard error
                                   width=0.25, size=0.8) +
  geom_text(aes(label=Sig,y=Starloc, vjust=0.5 ),size=8, position=position_dodge(w=1.0), 
              show.legend = FALSE, colour="black")+
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none") + 
        lims(y=c(-5.7, 10.7))  +
  geom_hline(yintercept=0, linetype="dashed", color="#525252") + 
  ggtitle(Cohorts2[x]) +
    scale_color_manual(values = phycols)
}) -> DeSeqPlots_Phy # list of plots
#plot all 3 in panel 
(DeSeqPlots_Phy[[1]] + DeSeqPlots_Phy[[2]] + DeSeqPlots_Phy[[3]]) 


#PHYLUM FACET PLOT - Each taxa is a plot, with each cohort on X compared to Lab-Mm
ResultsPhyAll %>%
  filter(pvalue<p_correctedP) %>%
  mutate(SampleGroup=fct_relevel(SampleGroup, Cohorts[1:3])) %>% 
  ggplot(aes(x=SampleGroup, 
             y=log2FoldChange,
             colour= Phylum)) + 
  geom_point(size=4)+
  geom_errorbar(aes(x=SampleGroup, 
                    ymax=log2FoldChange + lfcSE, 
                    ymin=log2FoldChange - lfcSE), # error bar + or - the standard error
                width=0.25, size=0.9) +
  geom_text(aes(label=Sig, vjust=-0.1),size=8, position=position_dodge(w=0.5), 
            show.legend = FALSE, colour="black")+
  theme_bw(base_size = 22) + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  scale_colour_manual(values=phy_rnbw) + 
  scale_y_continuous(limits=c(-10,13), breaks=c(-10,-5,0,5,10)) + 
  geom_hline(yintercept=0, linetype="dashed", color="#525252") + 
  theme(strip.background = element_rect(fill="white")) + 
  facet_wrap(~Phylum) -> DeseqFacet_Phy 


#GENUS PLOTS
#Make CohortPS but for Genus level
CohortPSGen2 <- lapply(CohortPS2, function(x){
  x %>%
    # phyloseq_filter_prevalence(., 
    # prev.trh = 0.1) %>% 
    aggregate_taxa(., 'Genus')
})
#[[1]] has 94 taxa, [[2]] has 97 and [[3]] has 96

MouseGen <- ps_core2 %>% # again try on ps_core2 instead of ps_all to skip prevalence filter that is not working
  ps_filter(SampleGroup%in%Cohorts[1:3]) %>%  
  #phyloseq_filter_prevalence(., prev.trh = 0.1) %>% 
  aggregate_taxa(., 'Genus')
#88 genus in 196 samples (all but Lab-Mm)

#deseq doesn't allow zeroes, so OTU values +1 and make sample group a factor instead of chrt
#think loop means this happens for each object in CohortPSPhy2 list?
for(i in 1:length(CohortPSGen2)){
  otu_table(CohortPSGen2[[i]]) <- otu_table(CohortPSGen2[[i]]) + 1  
  #sample_data(CohortPSGen2[[i]])$SampleGroup <-as.factor(sample_data(CohortPSGen2[[i]])$SampleGroup) 
  CohortPSGen2[[i]] <- CohortPSGen2[[i]] %>%
    ps_mutate(SampleGroup=fct_relevel(SampleGroup, "Lab:Mm"))
}
#check worked - want Lab:Mm to be first for each pairwise comparison
sample_data(CohortPSGen2[[3]])$SampleGroup %>% levels 

#calculating differential abundance and Wald test for each cohort 
ds_mouse2 <- list() 
dds_mouse2 <- list()
library(DESeq2)
#calculating the differential abundance values and Wald test for each cohort 
for(i in 1:length(CohortPSGen2)){
  ds_mouse2[[i]] <- phyloseq_to_deseq2(CohortPSGen2[[i]], ~ SampleGroup)
  dds_mouse2[[i]] <- DESeq(ds_mouse2[[i]], test="Wald", fitType="parametric")
}

ResultsGen <- 
  lapply(1:length(dds_mouse2), function(x){
    as.data.frame(results(dds_mouse2[[x]])) %>% 
      rownames_to_column(var="Genus") %>% 
      filter(!str_detect(Genus, "Family"),
             !str_detect(Genus, "Kingdom"),
             !str_detect(Genus, "Phylum"),
             !str_detect(Genus, "Order"),
             !str_detect(Genus, "Class"),
             !str_detect(Genus, "Incertae Sedis")) %>% # drops to 208 total (max = 70 Genus) 
      mutate(Genus=str_replace(Genus, "\\[|\\]", "")) %>% 
      mutate(Genus=str_replace(Genus, "\\]|\\]", "")) %>% 
      mutate(SampleGroup=Cohorts2[[x]]) #%>% 
      # mutate(Sig=case_when(padj< p_correctedG ~ "*",
      #                      padj> p_correctedG  ~ "")) %>%
      # mutate(Min=log2FoldChange-lfcSE,
      #        Max=log2FoldChange+lfcSE,
      #        Starloc=case_when(Sig=="*" ~ Max+(Max-Min),
      #                          Sig=="" ~ 0))
  }) 

#adjust the p-value so that it is divided by the number of tests 
#then pull significant values out from this -  
#allows you to correct for multiple testing and is more conservative than just usng 0.05
0.05 / (bind_rows(ResultsGen) %>%
  pull(Genus) %>% unique %>% length) -> p_correctedG 

#this includes taxa mapped at higher order than genus- 287 total 
#add in filter to drop these & mutate to remove [] round taxa names 
#prep to plot
ResultsGen <- 
  lapply(1:length(dds_mouse2), function(x){
    as.data.frame(results(dds_mouse2[[x]])) %>% 
      rownames_to_column(var="Genus") %>% 
      filter(!str_detect(Genus, "Family"),
             !str_detect(Genus, "Kingdom"),
             !str_detect(Genus, "Phylum"),
             !str_detect(Genus, "Order"),
             !str_detect(Genus, "Class"),
             !str_detect(Genus, "Incertae Sedis")) %>% # drops to 208 total (max = 70 Genus) 
      mutate(Genus=str_replace(Genus, "\\[|\\]", "")) %>% 
      mutate(Genus=str_replace(Genus, "\\]|\\]", "")) %>% 
      mutate(SampleGroup=Cohorts2[[x]]) %>% 
      mutate(Sig=case_when(padj< p_correctedG ~ "*",
                           padj> p_correctedG  ~ "")) %>%
      mutate(Min=log2FoldChange-lfcSE,
             Max=log2FoldChange+lfcSE,
             Starloc=case_when(Sig=="*" ~ Max+(Max-Min),
                               Sig=="" ~ 0))
  }) 
ResultsAllGen <- bind_rows(ResultsGen)
write.csv(ResultsAllGen, "DescriptiveTables/Genus_DifAbd.csv")

#colour palettes - need to expand as lots of taxa at genus level 
bind_rows(ResultsGen) %>%
  filter(pvalue<p_correctedG) %>% 
  arrange(-log2FoldChange) %>%
  pull(Genus) %>% unique -> GenNames #65 

gen_pal <- RColorBrewer::brewer.pal(11, "Spectral")
ColoursRamp <- colorRampPalette(gen_pal)
ColoursExpandGen <- ColoursRamp(n=65)
names(ColoursExpandGen) <- GenNames 

gen_rnbw <-c("#FEE0D2","#FCBBA1","#FB8072","#EF3B2C","#CB181D","#A50F15","#67000D","#B30000","#FC4E2A","#FC8D59",
               "#FDAE6B","#FDD0A2","#FDB462","#FD8D3C","#F16913","#D94801","#A63603","#FE9929","#FEC44F","#FED976",
               "#FFED6F","#FFEDA0","#F7FCB9","#D9F0A3","#ADDD8E","#B3DE69","#74C476","#41AB5D","#41AE76","#006D2C",
               "#00441B","#014636","#016C59","#66C2A4","#8DD3C7","#7BCCC4","#4EB3D3","#4292C6","#045A8D","#2171B5",
               "#08519C","#253494","#08306B","#0570B0","#3690C0","#6BAED6","#80B1D3","#9ECAE1","#A6BDDB","#D0D1E6",
               "#BEBADA","#8C6BB1","#88419D","#54278F","#3F007D","#4D004B","#BC80BD","#C994C7","#D4B9DA","#FDE0DD",
               "#FCC5C0","#FA9FB5","#F768A1","#DD3497","#AE017E","#CE1256","#980043","#800026","#67001F")

#ALL Genus PLOT
#all 3 plots compared to lab-Mm at Genus level - x axis ordered by magnitude of the log change
lapply(1:length(ResultsGen), function(x){
  ResultsGen[[x]] %>%
    filter(pvalue<p_correctedG) %>% 
    ggplot(aes(x=reorder(Genus, -log2FoldChange), #just put x=Genus to order alphabetically 
               y=log2FoldChange,
               colour= Genus,
               label="Sig")) + 
    geom_point(size=4)+
    geom_errorbar(aes(x=Genus, ymax=log2FoldChange + lfcSE, ymin=log2FoldChange - lfcSE), # error bar + or - the standard error
                  width=0.25, size=0.8) +
    geom_text(aes(label=Sig,y=Starloc, vjust=-0.5),size=8, position=position_dodge(w=0.5), 
              show.legend = FALSE, colour="black")+
    theme_bw(base_size = 22) + 
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.title.x = element_blank(),
          legend.position = "none") + 
    lims(y=c(-10, 15.5))  +
    geom_hline(yintercept=0, linetype="dashed", color="#525252") + 
    ggtitle(Cohorts2[x]) +
    scale_color_manual(values = ColoursExpandGen) 
}) -> DeSeqPlots2 # list of plots

#plot all 3 in panel 
(DeSeqPlots2[[1]] / DeSeqPlots2[[2]] / DeSeqPlots2[[3]])                                                               


#GENUS FACET PLOT - Each taxa is a plot, with each cohort on X compared to Lab-Mm
ResultsAllGen %>%
  filter(pvalue<p_correctedG) %>%
  mutate(SampleGroup=fct_relevel(SampleGroup, Cohorts[1:3])) %>% 
  ggplot(aes(x=SampleGroup, 
             y=log2FoldChange,
             colour= Genus)) + 
  geom_point(size=4)+
  geom_errorbar(aes(x=SampleGroup, 
                    ymax=log2FoldChange + lfcSE, 
                    ymin=log2FoldChange - lfcSE), # error bar + or - the standard error
                width=0.25, size=0.9) +
  geom_text(aes(label=Sig, vjust=-0.1),size=8, position=position_dodge(w=0.5), 
            show.legend = FALSE, colour="black")+
  theme_bw(base_size = 22) + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  scale_colour_manual(values=gen_rnbw) + 
  #lims(y=c(-9, 14)) +
  #scale_y_continuous(limits=c(12,-8), #breaks=c(-10,-5,0,5,10)) + 
  geom_hline(yintercept=0, linetype="dashed", color="#525252") + 
  theme(strip.background = element_rect(fill="white")) + 
  facet_wrap(~Genus) -> DeseqFacet_Gen 



#Thesis corrections 
rm(list=ls())

library(phyloseq)
library(microbiome)
library(colorspace)
library(pals)
library(microViz)
library(magrittr)
library(ggregplot)
library(ggpubr)
library(tidyverse)
library(ggforce)
library(ape)
library(ggrepel)
library(vegan)
library(ggplot2)
library(patchwork)
require(RColorBrewer)
library(ape)

source("/Volumes/Rowan_HD2/Chapter2-3_Analysis/VisualisationScripts/DataPrep_Baseline.R") 
# runs data prep & sets up all the data subsets 


#Then create a phyloseq that keeps negative & positive sequencing controls 
#do not prune reads <10000 as usual
ps_all1 <- ps4 %>% ps_filter(SampleDetailed!="colony control") 
#322 samples (drops the 4 colony controls not needed)
#includes mock community, golden and water samples

#To look at info on mock community, golden colonies and water controls use ps_faecal_all
ps_all1_Reads <- sort(phyloseq::sample_sums(ps_all1)) 
view(ps_all1_Reads)
#this shows the number of reads per sample 
#save as csv file - so can look at how many reads per sample, mock community, water etc
write.csv(ps_all1_Reads, "DescriptiveTables/AllSamples_ReadCounts.csv")
#Here we can see there are 4x samples with <10k reads - 
#041F00069B.2020-10-08, A14.2021-10-07, A11.2021-10-07, 15521825.2017-09-01
#All water samples >10k too, so would be dropped

#Make new phyoseq whereby:
#remove 4 x colony samples that had reads of less than >10k 
#make sure colony and wild faecal samples are distinguished
ps_all1 %>%
  ps_filter(!SampleId =="041F00069B.2020-10-08",
            !SampleId =="A14.2021-10-07",
            !SampleId =="A11.2021-10-07",
            !SampleId =="15521825.2017-09-01") %>%
  ps_mutate(SampleOrigin=case_when(str_detect(SampleDetailed, "MF1 lab mice") ~ "Lab:Mm",
                                   str_detect(SampleDetailed, "Pen") ~ "Wild:As", 
                                   str_detect(SampleDetailed, "Hew") ~ "Wild:As",
                                   str_detect(SampleDetailed, "Fal") ~ "Wild:As",
                                   str_detect(MouseOrigin, "Ash54 col") ~ "Lab-like:As", 
                                   str_detect(MouseOrigin, "b17 col") ~ "Wild-like:As",
                                   str_detect(MouseOrigin, "charles riv") ~ "Lab:Mm")) %>%
ps_mutate(SampleType2=case_when(SampleType=="water_primers" ~ "Water and Primers",
                                SampleType=="golden_colony" ~ "Golden Colony",
                                SampleType=="golden_wild" ~ "Golden Wild",
                                SampleType== "mock_community"~ "Mock Community",
                                SampleType== "water_only" ~ "Water Only",
                                SampleOrigin=="Lab:Mm"~ "Lab:Mm",
                                SampleOrigin=="Wild:As" ~ "Wild:As",
                                SampleOrigin=="Lab-like:As"~"Lab-like:As",
                                SampleOrigin=="Wild-like:As"~"Wild-like:As")) -> ps_all2


#then use ps_all2 and SampleType2 for all corrections plots/analysis 

#Get count of genus of mock community only 
ps_mock <- ps_all2 %>% ps_filter(SampleType2 == "Mock Community")
ps_mock_genus <- table(phyloseq::tax_table(ps_mock)[, "Genus"])
view(ps_mock_genus) #35 listed, but some not at genus level

sort(ps_mock_genus, decreasing = TRUE) 
#save as csv file - so can cross ref to what should be found in ATCC and see what might be contaminants/erros in taxonomy assignment
write.csv(ps_mock_genus, "DescriptiveTables/MockCommunity_RelAbd_Genus.csv")

#Remove bacteria listed that are not at genus level 
# Clostridia Class, Clostridia vadinBB60 group Order, Lachnospiraceae Family, Muribaculaceae Family

ps_mock <- subset_taxa(ps_mock, Genus!="Clostridia Class")
ps_mock <- subset_taxa(ps_mock, Genus!="Clostridia vadinBB60 group Order") 
ps_mock <- subset_taxa(ps_mock, Genus!="Lachnospiraceae Family")
ps_mock <- subset_taxa(ps_mock, Genus!="Muribaculaceae Family")
#now 31 left (should only be 20, so 11 are poss from contamination)

#Mock community - plot the relative abundance genus level
# vertical
ps_mock %>%
tax_agg("Genus") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>% # these are the defaults
  comp_barplot(tax_level = "Genus", n_taxa = 31,
               sample_order = c("mock_community_P1", "mock_community_P2", "mock_community_P3", 
                                "mock_community_P4", "mock_community_P5", "mock_community_P6",
                                "mock_community_P7", "mock_community_P8"),
               bar_width = 0.8, bar_outline_colour = "black", 
               order_with_all_taxa = TRUE) +
  ylab("Relative Abundance") +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18),
        legend.position="bottom", legend.title = element_text(size=18),
        legend.text = element_text(size=16)) +
  scale_x_discrete(guide = guide_axis(angle = 45))



#Get count of genus of water only samples 
ps_wo <- ps_all2 %>% ps_filter(SampleType2 == "Water Only")
ps_wo_genus <- table(phyloseq::tax_table(ps_wo)[, "Genus"])
view(ps_wo_genus)
# 14 bacteria total
sort(ps_wo_genus, decreasing = TRUE) 
#save as csv file - so can cross ref to what should be found in ATCC and see what might be contaminants/erros in taxonomy assignment
write.csv(ps_wo_genus, "DescriptiveTables/WaterOnly_RelAbd_Genus.csv")
#could put table into chapter - supp figs?

#remove anything not mapped at genus level 
ps_wo <-subset_taxa(ps_wo, Genus!="Lachnospiraceae Family")
ps_wo <-subset_taxa(ps_wo, Genus!="Muribaculaceae Family")
# now only 12 

#Water only  - relative abundance, genus level
# vertical
ps_wo %>%
  tax_agg("Genus") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>% # these are the defaults
  comp_barplot(tax_level = "Genus", n_taxa = 25,
               sample_order = c("water_only_P1", "water_only_P2", "water_only_P3",
                                "water_only_P4", "water_only_P5", "water_only_P6",
                                "water_only_P7","water_only_P8"),
               bar_width = 0.8, bar_outline_colour = "black", 
               order_with_all_taxa = TRUE) +
  ylab("Relative Abundance") +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18),
        legend.position="bottom", legend.title = element_text(size=18),
        legend.text = element_text(size=16)) +
  scale_x_discrete(guide = guide_axis(angle = 45))


#Get count of genus of water plus primers samples 
ps_wp <- ps_all2 %>% ps_filter(SampleType2 == "Water and Primers")
ps_wp_genus <- table(phyloseq::tax_table(ps_wp)[, "Genus"])
view(ps_wp_genus)
#22 in total
sort(ps_wp_genus, decreasing = TRUE) # 
#save as csv file - so can cross ref to what should be found in ATCC and see what might be contaminants/erros in taxonomy assignment
write.csv(ps_wp_genus, "DescriptiveTables/WaterPrimers_RelAbd_Genus.csv")
#could put table into chapter - supp figs?

ps_wp <- subset_taxa(ps_wp, Genus!="Muribaculaceae Family")
ps_wp <- subset_taxa(ps_wp, Genus!="Clostridia vadinBB60 group Order")
ps_wp <- subset_taxa(ps_wp, Genus!= "Peptococcaceae Family")
ps_wp <- subset_taxa(ps_wp, Genus!= "Staphylococcaceae Family") 
ps_wp <- subset_taxa(ps_wp, Genus!= "Lachnospiraceae Family") 

#Water and primers   - relative abundance, genus level
# vertical
ps_wp %>%
  tax_agg("Genus") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>% # these are the defaults
  comp_barplot(tax_level = "Genus", n_taxa = 25,
               sample_order = c("water_plus_primers_P1", "water_plus_primers_P2", "water_plus_primers_P3",
                                "water_plus_primers_P4", "water_plus_primers_P5", "water_plus_primers_P6",
                                "water_plus_primers_P7", "water_plus_primers_P8"),
               bar_width = 0.8, bar_outline_colour = "black", 
               order_with_all_taxa = TRUE) +
  ylab("Relative Abundance") +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18),
        legend.position="bottom", legend.title = element_text(size=18),
        legend.text = element_text(size=16)) +
  scale_x_discrete(guide = guide_axis(angle = 45))


#PCOA PLOTS - centered log ratio transformed data
#Trying to highlight that all control samples do not overlap with main faecal samples
#transoform taxanomic data to centered log ratio at genus level
ps_clr_genus <- ps_all2 %>% 
  tax_transform(trans = "clr", rank = "Genus") 
# plot clr transformation
ps_clr_genus %>% 
  ord_calc() %>% 
  ord_plot(color = "SampleType2", size = 3) 
#need to come back and make more pretty but seems to do the trick 
#this runs centered log ratio values of each sample, check with Amy S that this is enough for corrections
#would it be better to do this on bray-curtis ordinations etc?

#In case need to use distance based ordination methods used in thesis instead of CLR transformed data
#make tree for Unifrac 
random_tree_wild = rtree(ntaxa(ps_all2), rooted=TRUE, tip.label=taxa_names(ps_all2))
#save tree so don't need to make again and change analysis 
write.tree(random_tree_wild, "DescriptiveFigures/ps_all2_PhyloTree")

#read in tree saved previously so no need to make a new one and change analysis 
random_tree_wild <- read.tree("DescriptiveFigures/ps_all2_PhyloTree")

#merge tree with phyloseq 
ps_all2 <- merge_phyloseq(ps_all2, random_tree_wild)
# otu_table()   OTU Table:         [ 3819 taxa and 318 samples ]
# sample_data() Sample Data:       [ 318 samples by 20 sample variables ]
# tax_table()   Taxonomy Table:    [ 3819 taxa by 6 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 3819 tips and 3818 internal nodes ]

#see where golden samples sit with data 
ps_goldens <- ps_all2 %>% 
  ps_filter(!SampleType2 =="Mock Community",
            !SampleType2 =="Water Only",
            !SampleType2 =="Water and Primers")


#calculate beta diversity for ps_goldens 
ps_comp <- microbiome::transform(ps_goldens, "compositional")
ps_bray <- ordinate(ps_goldens, "PCoA", "bray")
ps_unif <- ordinate(ps_goldens, "PCoA", "unifrac")
ps_wunif <- ordinate(ps_goldens, "PCoA", "wunifrac")
ps_meta <- meta(ps_goldens)

#define what order you want the samples to be colored first 
sample_data(ps_goldens)$SampleType2<- factor(sample_data(ps_goldens)$SampleType2, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm", "Golden Colony", "Golden Wild"))


# BRAY-CURTIS - goldens vs samples 
plot_ordination(ps_goldens,
                ps_bray, 
                color="SampleType2") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#BDBDBD","#000000")) +
 theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top") +
  labs(color= "") +
  stat_ellipse() + 
  ggtitle('Bray-Curtis') -> gold_bray


# WUNIF - goldens vs samples 
 plot_ordination(ps_goldens,
                 ps_wunif, 
                 color="SampleType2") +
   geom_point(aes(), size= 2) + 
   scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#BDBDBD","#000000")) +
 theme_bw(base_size = 22) +
   theme(plot.title = element_text(hjust = 0, size = 22),
         legend.position = "top") +
   labs(color= "") +
   stat_ellipse() + 
   ggtitle('Weighted Unifrac') -> gold_wuif
 
 # UNIF - goldens vs samples 
 plot_ordination(ps_goldens,
                 ps_unif, 
                 color="SampleType2") +
   geom_point(aes(), size= 2) + 
   scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#BDBDBD","#000000")) +
 theme_bw(base_size = 22) +
   theme(plot.title = element_text(hjust = 0, size = 22),
         legend.position = "top") +
   labs(color= "") +
   stat_ellipse() + 
   ggtitle('Unweighted Unifrac') -> gold_unif
 
 (gold_bray|gold_wuif|gold_unif) +
  plot_layout(guides = 'collect') & theme(legend.position = "top") 

#see where waters sit with data 
ps_waters <- ps_all2 %>% 
  ps_filter(!SampleType2 =="Mock Community",
            !SampleType2 =="Golden Colony",
            !SampleType2 =="Golden Wild")

#calculate beta diversity for ps_waters
ps_comp1 <- microbiome::transform(ps_waters, "compositional")
ps_bray1 <- ordinate(ps_waters, "PCoA", "bray")
ps_unif1 <- ordinate(ps_waters, "PCoA", "unifrac")
ps_wunif1 <- ordinate(ps_waters, "PCoA", "wunifrac")
ps_meta1 <- meta(ps_waters)

#define what order you want the samples to be colored first 
sample_data(ps_waters)$SampleType2<- factor(sample_data(ps_waters)$SampleType2, levels = c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm",
                                                                                           "Water Only", "Water and Primers"))

#BRAY-CURTIS - waters vs samples 
plot_ordination(ps_waters,
                ps_bray1, 
                color="SampleType2") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#BDBDBD","#000000"))  + theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top")+
  labs(color= "") +
  stat_ellipse() + 
  ggtitle('Bray-Curtis') -> water_bray

# WUNIF- waters vs samples 
plot_ordination(ps_waters,
                ps_wunif1, 
                color="SampleType2") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#BDBDBD","#000000"))  + theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top")+
  labs(color= "") +
  stat_ellipse() + 
  ggtitle('Weighted Unifrac') -> water_wunif

# UNIF- waters vs samples 
plot_ordination(ps_waters,
                ps_unif1, 
                color="SampleType2") +
  geom_point(aes(), size= 2) + 
  scale_colour_manual(values = c("#74C476", "#6BAED6", "#BC80BD","#FB8072","#BDBDBD","#000000"))  + theme_bw(base_size = 22) +
  theme(plot.title = element_text(hjust = 0, size = 22),
        legend.position = "top") +
  labs(color= "") +
  stat_ellipse() + 
  ggtitle('Unweighted Unifrac') -> water_unif 

(water_bray|water_wunif|water_unif) +
  plot_layout(guides = 'collect') & theme(legend.position = "top") 
rm(list=ls())

library(MASS)
library(mia)
library(Matrix)
library(ade4)
library(multtest)
library(vegan)
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
library(decontam)


# read in data  -----------------------------------------------------------


ps <- readRDS("Data/Phyloseqs/phyloseqBaseline.rds")
ps 
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6557 taxa and 326 samples (controls still included)]
# sample_data() Sample Data:       [ 326 samples by 17 sample variables ]
# tax_table()   Taxonomy Table:    [ 6557 taxa by 6 taxonomic ranks ]

#### add in the age / other variables for wild mice 

meta_full <- readxl::read_xlsx("Data/Metadata/Baseline_Data_RB_full.xlsx") 
colnames(meta_full)

meta_merge <- meta_full %>% select(sample_id, age, repro, body_mass_g, treated) %>% 
  janitor::clean_names("big_camel") 

setdiff(meta_merge$SampleId, ps %>% meta %>% pull(SampleId))

meta_old <- ps %>% meta 
meta_new <- meta_old %>% 
  dplyr::left_join(., meta_merge, by="SampleId")
rownames(meta_new) <- meta_new$SampleId

# check sample names match 
setdiff(sample_names(ps), rownames(meta_new)) # all match 

meta_new <- meta_new[sample_names(ps), ] # make the order the same  

# add new data in 
sample_data(ps) <- sample_data(meta_new)

# check 
ps %>% meta %>% colnames # done 

## dealing with negatives and contaminants ----------------------------------
contams <- isContaminant(ps, method= "prevalence",
                         neg= sample_data(ps)$SampleType == c("water_primers", "water_only"),
                         threshold= 0.5)

#saveRDS(contams, "Output/contaminants.rds")
#contams <- readRDS("Outputs/contaminants.rds")

contaminants <- rownames(contams[contams$contaminant == TRUE, ])
length(contaminants) # only 6 
clean_ps <- prune_taxa(!(taxa_names(ps) %in% contaminants), ps)
clean_ps <- prune_samples(!(sample_data(clean_ps)$SampleType %in% c("water_primers", "water_only",
                                                                    "mock_community", "golden_colony", "golden_wild")), 
                          clean_ps) # keep samples only 

clean_ps # 286 samples with samples only minus controls / technical samples 

# QC  
any(taxa_sums(clean_ps) == 0) # TRUE - some have 0? 
any(taxa_sums(clean_ps) == 1) # FALSE - none have 1  
colSums(is.na(tax_table(clean_ps))) # NA table below - lots of NA's but let's keep version with them so that e.g. genus level analyses only have those that are assigned 
# Kingdom  Phylum   Class   Order  Family   Genus 
# 2     115     203     330    1399    3903

#try to fix NA's, best way to try to keep as much info as possible is to remove any data that can't be mapped to a Kingdom?
ps1 <- prune_taxa(taxa_sums(clean_ps) >1, clean_ps)
# Drops ASVs down to 3787 (lose around 2900 ASV's that were only present once or less)

ps1 <- phyloseq_validate(ps1, remove_undetected = TRUE)
# NAs detected in phyloseq tax_table:
# Consider using tax_fix() to make taxa uniquely identifiable

colSums(is.na(tax_table(ps1)))  
# Kingdom  Phylum   Class   Order  Family   Genus 
# 1      51      68     132     770    2155
# Got rid of a fair few NAs from ps, but still quite a few left need to get rid of 

#ps1 <- tax_fix(ps1)
# removes NA's
#colSums(is.na(tax_table(ps1))) 
# Kingdom  Phylum   Class   Order  Family   Genus 
# 0       0       0       0       0       0 

#However, if we then look at the taxa table and list the phylum:
ps_sample_Phylum <- table(phyloseq::tax_table(ps1)[, "Phylum"])
#View(ps_sample_Phylum)
#18 phylum

#Now need to deal with dropping some other taxa that are not useful

#First we manually change names of completely unclassified taxa  to Incertae Sedis and can then be dropped from analysis 
#1 individual ASV this is happening for and only occurs once in the dataset, can be removed as not meaningful at all - NA at kingdom 
#tax_table(ps1) %>% as.data.frame %>% filter(is.na(Kingdom)) %>% View 

# Extract taxonomy table
tax <- tax_table(ps1)
# Find taxa where all taxonomy ranks are NA
all_na_taxa <- apply(tax, 1, function(x) all(is.na(x)))
# Remove those taxa from the phyloseq object
ps1 <- prune_taxa(!all_na_taxa, ps1) # removed 1, 3786 now 

ps1@tax_table %>% as.data.frame %>% pull(Kingdom) %>% unique # no NAs 
ps1@tax_table %>% as.data.frame %>% pull(Family) %>% unique # Mitchondria can come out 
ps1@tax_table %>% as.data.frame %>% pull(Order) %>% unique # Chloroplast can come out 

#examining tax table; in Kingdom we have 1xArchaea which we can keep 
#also, we have 77 taxa under the family Mitochondria and 13 under the order Chloroplasts - none of these are bacterial or useful 
#can we try and remove them - tried all in one step, didn't work so do stepwise
ps2 <- subset_taxa(ps1, (Family!="Mitochondria") | is.na(Family)) #Drop mitochondria (will drop NAs if don't include the na bit) - 3709 
ps3 <- subset_taxa(ps2, (Order!="Chloroplast") | is.na(Order)) #Drop chloroplast - 3696 after 
#use ps3 for further downstream analysis 

# taxa_sums(ps3) %>% as.data.frame() %>% 
#   arrange(desc(.)) %>% 
#   rename(TotalAbundance=1) -> TotalAbundances 
# summary(TotalAbundances)
 

#this bit of code here changes the ASV names to unique names ASV1 - ASV3701
#instead of their long sequence name
#makes the tax and otu tables easier to manipulate


# Calculate prevalence (number of samples in which each ASV is present)
asv <- as(otu_table(ps3), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(ps3)) {
  asv <- t(asv)
}

prev <- rowSums(asv > 0)

# Order taxa by decreasing prevalence
taxa_order <- names(sort(prev, decreasing = TRUE))

# Create new ASV names
new_names <- paste0("ASV", seq_along(taxa_order))

# Create lookup vector: old name -> new name
name_map <- setNames(new_names, taxa_order)

# Rename taxa in phyloseq object
taxa_names(ps3) <- name_map[taxa_names(ps3)]

# Inspect mapping
head(data.frame(
  old_name = names(name_map),
  new_name = unname(name_map)
))
 
 # check they're in order prevalence wise 
prevalence(ps3) %>% as.data.frame() %>%
  arrange(desc(.)) %>%
  dplyr::rename(Prevalence=1) %>% 
  mutate(Prevalence=Prevalence*100)-> PrevalenceASVs
#View(PrevalenceASVs) # all good 

asv_key <- data.frame(
  Sequence = names(name_map),
  ASV = unname(name_map),
  Prevalence = prev[names(name_map)]
)

head(asv_key)
#saveRDS(asv_key, "Output/asv_key.rds")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#SUBSETTING 
#subsetting data into relevent groups for analysis 

# clean variables 
ps3 <- ps3 %>% 
  ps_mutate(SampleDetailed=str_replace(SampleDetailed, "inFection", "infection")) %>% 
  ps_mutate(SampleExpt=str_replace(SampleExpt, "inFection", "infection")) %>% 
  ps_mutate(MouseOrigin=str_replace(MouseOrigin, "wild ", "wild")) %>% 
  ps_mutate(MouseOrigin=str_replace(MouseOrigin, "charles river ", "charles river")) %>% 
  ps_mutate(MouseOrigin=str_replace(MouseOrigin, "Ash54 colony ", "Ash54 colony"))

meta(ps3) %>% pull(SampleDetailed) %>% unique # fixed 
meta(ps3) %>% pull(MouseOrigin) %>% unique 


#Make a pyloseq object from ps3 that only has the sample of interest in (so drop the 4 colony controls)

ps3@otu_table %>% as.data.frame %>% rowSums %>% summary 

#Then prune samples to remove that have less than 10,000 reads
ps_all <- ps3 %>% ps_filter(SampleDetailed!="colony control") %>% 
  prune_samples(sample_sums(.) > 10000, .) 
#leaves 278 samples with 3686 ASVs
#4 colony controls dropped 

#check low seq read samples
SeqDepth = rowSums(otu_table(ps_all))
hist(SeqDepth)
mean(SeqDepth) #72089.1
min(SeqDepth)

ps_all %<>% 
  ps_mutate(
    SeqDepth = SeqDepth, 
  )

# phyloseq technical subsets & basic summary info 

# Analysis 1 - ALL GROUPS (drop CD1 & WM colony controls from caesarean re-derivation)
#maybe drop MF1 later & collapse all wild samples into 1 group
AllIds <- meta(ps3) %>% 
  filter(SampleDetailed!="colony control") %>%
  pull(SampleId)
#282 samples (dropped 4) 

# Analysis 2 - three main experiment lab groups, wild-like, lab-like & CD1's
RseIds <- meta(ps3) %>% 
  filter(SampleDetailed %in% c("nutrition expt","infection expt","antibiotic expt")) %>%
  pull(SampleId)
# 215 samples (just RSE baseline samples)

# Analysis 3 - All wood mice cohorts - wild-like, lab-like, Hewan, Penicuik & Falkirk
WoodIds <- meta(ps3) %>% 
  filter(MouseSpecies=="Apodemus sylvaticus") %>% 
  filter(SampleDetailed!="colony control") %>%
  pull(SampleId)
# 200 samples (all woodmice except 2 colony controls from caesarean re-derivation)

#Analysis 4 - Lab groups only (no wild samples)
LabIds <- meta(ps3) %>%
  filter(MouseOrigin!= "wild") %>% #on csv file space after word wild caused some issues so watch out for
  filter(SampleDetailed!="colony control") %>% 
  pull(SampleId)
# 225 samples (wild-like, lab-like, CD1, MF1, plus pos and neg controls) 

# Analysis 1 - ps_all has all relevant groups (with col contrls & low reads dropped)
ps_all <-ps_all %>% 
  ps_mutate(SampleOrigin=case_when(str_detect(SampleDetailed, "MF1 lab mice") ~ "MF1 mice",
                                   str_detect(SampleDetailed, "Penicuik") ~ "Wild mice", 
                                   str_detect(SampleDetailed, "Hewan") ~ "Wild mice",
                                   str_detect(SampleDetailed, "Falkirk") ~ "Wild mice",
                                   str_detect(MouseOrigin, "Ash54 col") ~ "Lab-like mice", 
                                   str_detect(MouseOrigin, "b17 col") ~ "Wild-like mice",
                                   str_detect(MouseOrigin, "charles river") ~ "CD1 mice"))

ps_all %>% meta %>% group_by(SampleOrigin) %>% dplyr::count()


#basic info on abundance & prevalence for ps_all
taxa_sums(ps_all) %>% as.data.frame() %>%
  arrange(desc(.)) %>%
  dplyr::rename(TotalAbundance=1) -> TotalAbundances_All
summary(TotalAbundances_All)
# 3691 ASVs - why 0's not gone here - probably one snuck in after a sample was cut the ASV is now zero again 
# TotalAbundance
# Min.   :      0.0
# 1st Qu.:      8.0
# Median :     81.0
# Mean   :   5516.2
# 3rd Qu.:    762.5
# Max.   :2599393.0


# Analysis 2 - ps_RSE three main experiment lab groups, wild-like, lab-like & CD1's
ps_rse<- ps_all %>% 
  ps_filter(SampleId %in% RseIds) 

ps_rse
#213 samples (2 lost with low reads), 1769 taxa

#basic info on abundance & prevalence for ps_RSE 
# taxa_sums(ps_RSE) %>% as.data.frame() %>% 
#   arrange(desc(.)) %>% 
#   rename(TotalAbundance=1) -> TotalAbundances_RSE
# summary(TotalAbundances_RSE)
# TotalAbundance   
# Min.   :      2  
# 1st Qu.:      6  
# Median :    192  
# Mean   :   8801  
# 3rd Qu.:   2308  
# Max.   :2083603
# prevalence(ps_RSE) %>% as.data.frame() %>% 
#   arrange(desc(.)) %>% 
#   rename(Prevalence=1) %>% 
#   mutate(Prevalence=Prevalence*100)-> PrevalenceASVs_RSE
# summary(PrevalenceASVs_RSE)
# Prevalence      
# Min.   :  0.4695  
# 1st Qu.:  0.4695  
# Median :  2.3474  
# Mean   : 10.6685  
# 3rd Qu.: 15.0235  
# Max.   :100.0000 


# Analysis 3 - ps_wood all wood mice cohorts - wild-like, lab-like, Hewan, Penicuik & Falkirk
ps_wood <- ps_all %>% 
  ps_filter(SampleId %in% WoodIds)

ps_wood
#196 samples (4 lost with low reads) and 3012 ASVs

#basic info on abundance & prevalence for ps_wood
# taxa_sums(ps_wood) %>% as.data.frame() %>% 
#   arrange(desc(.)) %>% 
#   rename(TotalAbundance=1) -> TotalAbundances_Wood
# summary(TotalAbundances_Wood)
# TotalAbundance     
# Min.   :      2.0  
# 1st Qu.:     10.0  
# Median :     81.0  
# Mean   :   4631.0  
# 3rd Qu.:    625.5  
# Max.   :2033706.0
# 
# prevalence(ps_wood) %>% as.data.frame() %>% 
#   arrange(desc(.)) %>% 
#   rename(Prevalence=1) %>% 
#   mutate(Prevalence=Prevalence*100)-> PrevalenceASVs_Wood
# summary(PrevalenceASVs_Wood)
# Prevalence      
# Min.   :  0.5102  
# 1st Qu.:  0.5102  
# Median :  1.0204  
# Mean   :  5.9671  
# 3rd Qu.:  4.5918  
# Max.   :100.0000

#Analysis 4 - ps_lab Lab groups only (no wild samples)
ps_lab <- ps_all %>% 
  ps_filter(SampleId %in% LabIds)

ps_lab

#223 samples (2 low reads dropped) 1911 ASVs

#basic info on abundance & prevalence for ps_lab
# taxa_sums(ps_lab) %>% as.data.frame() %>% 
#   arrange(desc(.)) %>% 
#   rename(TotalAbundance=1) -> TotalAbundances_Lab
# summary(TotalAbundances_Lab)
# TotalAbundance   
# Min.   :      2  
# 1st Qu.:      6  
# Median :    172  
# Mean   :   8614  
# 3rd Qu.:   2204  
# Max.   :2175406

# prevalence(ps_lab) %>% as.data.frame() %>% 
#   arrange(desc(.)) %>% 
#   rename(Prevalence=1) %>% 
#   mutate(Prevalence=Prevalence*100)-> PrevalenceASVs_Lab
# summary(PrevalenceASVs_Lab)
# Prevalence      
# Min.   :  0.4484  
# 1st Qu.:  0.4484  
# Median :  1.7937  
# Mean   : 10.0498  
# 3rd Qu.: 13.4529  
# Max.   :100.0000



# save all 
# saveRDS(ps_all, "Data/Phyloseqs/phyloseqAllSamples.rds")
# saveRDS(ps_rse, "Data/Phyloseqs/phyloseqRSESamples.rds")
# saveRDS(ps_wood, "Data/Phyloseqs/phyloseqWoodMiceSamples.rds")
# saveRDS(ps_lab, "Data/Phyloseqs/phyloseqLabSamples.rds")
# 


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
library(readxl)
library(dplyr)

# read in data  -----------------------------------------------------------
ps <- readRDS("Data/Phyloseqs/phyloseqBaseline.rds")
ps 

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6557 taxa and 326 samples (controls still included)]
# sample_data() Sample Data:       [ 326 samples by 17 sample variables ]
# tax_table()   Taxonomy Table:    [ 6557 taxa by 6 taxonomic ranks ]



# add in the age and other variables for wild mice  -----------------------
meta_full <- readxl::read_xlsx("Data/Metadata/Baseline_Data_RB_full.xlsx") 
colnames(meta_full)

meta_merge <- meta_full %>% select(sample_id, age, repro, body_mass_g, treated) %>% 
  janitor::clean_names("big_camel") 

setdiff(meta_merge$SampleId, ps %>% meta %>% pull(SampleId))

meta_old <- ps %>% meta 
meta_new <- meta_old %>% 
  dplyr::left_join(., meta_merge, by="SampleId")

## phyloseq sample_data requires row names to correspond to sample IDs
rownames(meta_new) <- meta_new$SampleId

setdiff(sample_names(ps), rownames(meta_new)) # all match 

## ensures metadata rows appear in the same order as phyloseq sample order
meta_new <- meta_new[sample_names(ps), ]   

## inserts updated metadata into phyloseq 
sample_data(ps) <- sample_data(meta_new)

## verify new variables were added 
ps %>% meta %>% colnames



# dealing with negatives and contaminants ---------------------------------

## identify contaminating taxa as those occurring more frequently in negative 
## controls than in true biological samples. 
contams <- isContaminant(ps, method= "prevalence",
                         neg = sample_data(ps)$SampleType == c("water_primers", "water_only"),
                         threshold= 0.5)

#saveRDS(contams, "Output/contaminants.rds")
#contams <- readRDS("Outputs/contaminants.rds")

contaminants <- rownames(contams[contams$contaminant == TRUE, ])
length(contaminants) # only 7 

## remove contaminant taxa
clean_ps <- prune_taxa(!(taxa_names(ps) %in% contaminants), ps)

## remove non-biological samples
clean_ps <- prune_samples(!(sample_data(clean_ps)$SampleType %in% 
            c("water_primers", "water_only",
            "mock_community", "golden_colony", "golden_wild")), 
            clean_ps)

clean_ps # 286 biological samples



# quality control ---------------------------------------------------------

any(taxa_sums(clean_ps) == 0) # TRUE - some have 0 after samples were removed 
any(taxa_sums(clean_ps) == 1) # FALSE - none have 1  

## NA table below - lots of NA's but let's keep version with them so that 
## e.g. genus level analyses only have those that are assigned taxa
## taxonomic assignment gets worse at lower ranks (normal)
colSums(is.na(tax_table(clean_ps))) 
## Kingdom  Phylum   Class   Order  Family   Genus 
## 2     115     203     330    1399    3903

## remove taxa with total abundance </=1
## often sequencing noise and contribute sparsity and slows down analyses
ps1 <- prune_taxa(taxa_sums(clean_ps) >1, clean_ps)
## drops ASVs down to 3786 (lose ~2900 that were only present once or less)

## check (1) taxonomy dimensions match OTU table, (2) sample IDs match, 
## (3) empty taxa removed
ps1 <- phyloseq_validate(ps1, remove_undetected = TRUE)
## NAs detected in phyloseq tax_table:
## Consider using tax_fix() to make taxa uniquely identifiable

## recheck missing taxonomy
colSums(is.na(tax_table(ps1)))  
## Kingdom  Phylum   Class   Order  Family   Genus 
## 1      51      68     132     770    2155
## still quite a few NA's left we need to get rid of 

## optional: could force clean taxonomy using tax_fix 
## ps1 <- tax_fix(ps1) # removes NA's
## colSums(is.na(tax_table(ps1))) 
## Kingdom  Phylum   Class   Order  Family   Genus 
## 0       0       0       0       0       0 

## look at the taxa table and list the phylum (18 total):
length(ps_sample_Phylum <- table(phyloseq::tax_table(ps1)[, "Phylum"]))

## first we manually change names of completely unclassified taxa  to Incertae Sedis and can then be dropped from analysis 
## 1 individual ASV this is happening for and only occurs once in the dataset, can be removed as not meaningful at all - NA at kingdom 
## tax_table(ps1) %>% as.data.frame %>% filter(is.na(Kingdom)) %>% View 

## remove completely unidentified ASVs (sequencing artifacts or non-specific alignments)
tax <- tax_table(ps1)
all_na_taxa <- apply(tax, 1, function(x) all(is.na(x)))
ps1 <- prune_taxa(!all_na_taxa, ps1) # removed 1, 3785 now 

ps1@tax_table %>% as.data.frame %>% pull(Kingdom) %>% unique # no NAs 
ps1@tax_table %>% as.data.frame %>% pull(Family) %>% unique # Mitchondria can come out 
ps1@tax_table %>% as.data.frame %>% pull(Order) %>% unique # Chloroplast can come out 

## remove host-derived and plant contaminant taxa
ps2 <- subset_taxa(ps1, (Family!="Mitochondria") | is.na(Family))
ps3 <- subset_taxa(ps2, (Order!="Chloroplast") | is.na(Order))
## have 3695 ASVs now


# rename ASVs (prevalence-based ordering) ---------------------------------

## this code changes ASV names to unique names ASV1 - ASV3695
## which makes the tax and otu tables easier to manipulate

## calculate prevalence (number of samples in which each ASV is present)
asv <- as(otu_table(ps3), "matrix")

## Ensure taxa are rows
if (!taxa_are_rows(ps3)) {asv <- t(asv)}

## compute prevalence
prev <- rowSums(asv > 0)

## order taxa by decreasing prevalence
taxa_order <- names(sort(prev, decreasing = TRUE))

## rename ASVs sequentially
new_names <- paste0("ASV", seq_along(taxa_order))
name_map <- setNames(new_names, taxa_order)
taxa_names(ps3) <- name_map[taxa_names(ps3)]

## verify mapping
head(data.frame(
  old_name = names(name_map),
  new_name = unname(name_map)
))
 
## check they're in order by prevalence 
prevalence(ps3) %>% as.data.frame() %>%
  arrange(desc(.)) %>%
  dplyr::rename(Prevalence=1) %>% 
  mutate(Prevalence=Prevalence*100)-> PrevalenceASVs
#View(PrevalenceASVs) # all good 

## create ASV key (allows us to map results back to sequences)
asv_key <- data.frame(
  Sequence = names(name_map),
  ASV = unname(name_map),
  Prevalence = prev[names(name_map)]
)

head(asv_key)
#saveRDS(asv_key, "Output/asv_key.rds")


# subsetting data for analyses --------------------------------------------

# clean metadata variables 
ps3 <- ps3 %>% 
  ps_mutate(SampleDetailed=str_replace(SampleDetailed, "inFection", "infection")) %>% 
  ps_mutate(SampleExpt=str_replace(SampleExpt, "inFection", "infection")) %>% 
  ps_mutate(MouseOrigin=str_replace(MouseOrigin, "wild ", "wild")) %>% 
  ps_mutate(MouseOrigin=str_replace(MouseOrigin, "charles river ", "charles river")) %>% 
  ps_mutate(MouseOrigin=str_replace(MouseOrigin, "Ash54 colony ", "Ash54 colony"))

meta(ps3) %>% pull(SampleDetailed) %>% unique # fixed 
meta(ps3) %>% pull(MouseOrigin) %>% unique 


## drop 4 colony control samples and prune samples with <10,000 reads 
## low-depth samples cause unstable alpha diversity estimates, inflated
## beta diversity distances, and false absence of taxa
ps3@otu_table %>% as.data.frame %>% rowSums %>% summary 
ps_all <- ps3 %>% ps_filter(SampleDetailed!="colony control") %>% 
  prune_samples(sample_sums(.) > 10000, .) 
## leaves 278 samples with 3686 ASVs

## check for low sequencing depth samples
SeqDepth = rowSums(otu_table(ps_all))
summary(SeqDepth); hist(SeqDepth)

## add sequencing depth as metadata variable
ps_all %<>% 
  ps_mutate(
    SeqDepth = SeqDepth, 
  )

## combining/renaming groups in ps_all metadata
ps_all <-ps_all %>% 
  ps_mutate(SampleOrigin=case_when(str_detect(SampleDetailed, "MF1 lab mice") ~ "MF1 mice",
                                   str_detect(SampleDetailed, "Penicuik") ~ "Wild mice", 
                                   str_detect(SampleDetailed, "Hewan") ~ "Wild mice",
                                   str_detect(SampleDetailed, "Falkirk") ~ "Wild mice",
                                   str_detect(MouseOrigin, "Ash54 col") ~ "Lab-like mice", 
                                   str_detect(MouseOrigin, "b17 col") ~ "Wild-like mice",
                                   str_detect(MouseOrigin, "charles river") ~ "CD1 mice"))
ps_all %>% meta %>% group_by(SampleOrigin) %>% dplyr::count()

## basic info on abundance & prevalence for ps_all
taxa_sums(ps_all) %>% as.data.frame() %>%
  arrange(desc(.)) %>%
  dplyr::rename(TotalAbundance=1) -> TotalAbundances_All
summary(TotalAbundances_All)
## 3686 ASVs - why 0's not gone here - probably one snuck in after a sample 
## was cut the ASV is now zero again 



# create phylogenetic tree from ASV sequences -----------------------------

# library(Biostrings)   # DNA sequence handling
# library(DECIPHER)     # Multiple sequence alignment
# library(phangorn)     # Phylogenetic tree inference
# library(ape)          # Tree manipulation
# 
# # create DNA sequence object using asv_key
# seqs <- DNAStringSet(asv_key$Sequence)
# names(seqs) <- asv_key$ASV
# 
# # align ASV sequences
# alignment <- AlignSeqs(seqs)
# 
# # convert alignment for phangorn
# phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
# 
# # build starting tree, calculate distances between ASVs
# dm <- dist.ml(phang.align)
# # create a neighbor joining starting tree
# treeNJ <- NJ(dm)
# 
# # maximum likelihood tree optimization
# fit <- pml(treeNJ, phang.align)
# 
# # NNI is appropriate for 16S amplicon data because it improves the tree
# # without the very long runtime of stochastic searches
# fitGTR <- optim.pml(
#   fit,
#   model = "GTR", # general time reversible nucleotide substitution model
#   optInv = TRUE, # estimate proportion of invariant sites
#   rearrangement = "NNI" # nearest neighbor interchange topology search
# )
# 
# # extract final tree
# tree <- fitGTR$tree
# 
# # check matching ASVs
# # any in tree but not phylogeny?
# setdiff(tree$tip.label, taxa_names(ps_all))
# # any in phyloseq but missing from tree?
# setdiff(taxa_names(ps_all), tree$tip.label)
# 
# # prune tree to match phyloseq
# tree <- keep.tip(tree, taxa_names(ps_all))
# length(tree$tip.label)
# ntaxa(ps_all)
# 
# # root the tree
# # midpoint rooting is commonly used for microbiome studies
# # when no appropriate outgroup is avaiable
# tree_rooted <- midpoint(tree)
# is.rooted(tree_rooted)

# save/load the final tree and add to ps_all
# saveRDS(tree_rooted, file = "Output/phylogenetic_tree.rds")
tree <- readRDS("Output/phylogenetic_tree.rds")
phy_tree(ps_all) <- tree



# final filtering of groups -----------------------------------------------

## Hewan, Penicuik & Falkirk samples all combined into one 'wild' category
## MF1 & CD1 mice collapsed into Lab
## ash54 wm colony = lab-like; b17 wm colony = wild-like
## filter out MF1 and all wild mice except Hewan 2020
ps_all <- ps_all %>% 
  ps_mutate(SampleOrigin=case_when(str_detect(SampleDetailed, "MF1 lab mice") ~ "Lab:Mm",
                                   str_detect(SampleDetailed, "Pen") ~ "Wild:As", 
                                   str_detect(SampleDetailed, "Hew") ~ "Wild:As",
                                   str_detect(SampleDetailed, "Fal") ~ "Wild:As",
                                   str_detect(MouseOrigin, "Ash54 col") ~ "Lab-like:As", 
                                   str_detect(MouseOrigin, "b17 col") ~ "Wild-like:As",
                                   str_detect(MouseOrigin, "charles river") ~ "Lab:Mm")) %>% 
  ps_mutate(SampleOrigin=fct_relevel(SampleOrigin, c("Wild:As", "Wild-like:As", "Lab-like:As", "Lab:Mm"))) %>% 
  ps_filter(!SampleDetailed %in% c("Penicuik 2020", "Falkirk 2017", "MF1 lab mice")) 


# 
# # phyloseq technical subsets & basic summary info -------------------------
# 
# ## Analysis 1 - ALL Groups (excluding colony controls)
# ## maybe drop MF1 later & collapse all wild samples into 1 group
# AllIds <- meta(ps3) %>% 
#   filter(SampleDetailed!="colony control") %>%
#   pull(SampleId) # 282 samples
# 
# ## Analysis 2 - RSE experimental groups (wild-like, lab-like & CD1's)
# RseIds <- meta(ps3) %>% 
#   filter(SampleDetailed %in% c("nutrition expt","infection expt","antibiotic expt")) %>%
#   pull(SampleId) # 215 samples (just RSE baseline samples)
# 
# ps_rse<- ps_all %>% 
#   ps_filter(SampleId %in% RseIds) 
# ps_rse # 213 samples (2 lost with low reads), 1769 ASVs
# 
# ## abundance & prevalence
# taxa_sums(ps_rse) %>% as.data.frame() %>%
#   arrange(desc(.)) %>%
#   rename(TotalAbundance=1) -> TotalAbundances_RSE
# summary(TotalAbundances_RSE)
# 
# prevalence(ps_rse) %>% as.data.frame() %>%
#   arrange(desc(.)) %>%
#   rename(Prevalence=1) %>%
#   mutate(Prevalence=Prevalence*100)-> PrevalenceASVs_RSE
# summary(PrevalenceASVs_RSE)
# 
# ## Analysis 3 - All wood mouse cohorts (wild-like, lab-like, Hewan, Penicuik & Falkirk)
# WoodIds <- meta(ps3) %>% 
#   filter(MouseSpecies=="Apodemus sylvaticus") %>% 
#   filter(SampleDetailed!="colony control") %>%
#   pull(SampleId) # 200 samples
# 
# ps_wood <- ps_all %>% 
#   ps_filter(SampleId %in% WoodIds)
# ps_wood # 196 samples (4 lost with low reads), 3012 ASVs
# 
# ## abundance & prevalence
# taxa_sums(ps_wood) %>% as.data.frame() %>%
#   arrange(desc(.)) %>%
#   rename(TotalAbundance=1) -> TotalAbundances_Wood
# summary(TotalAbundances_Wood)
#  
# prevalence(ps_wood) %>% as.data.frame() %>%
#   arrange(desc(.)) %>%
#   rename(Prevalence=1) %>%
#   mutate(Prevalence=Prevalence*100)-> PrevalenceASVs_Wood
# summary(PrevalenceASVs_Wood)
# 
# ## Analysis 4 - Lab groups only (no wild samples)
# LabIds <- meta(ps3) %>%
#   filter(MouseOrigin!= "wild") %>% # on csv file space after word wild caused some issues so watch out for
#   filter(SampleDetailed!="colony control") %>% 
#   pull(SampleId) # 225 samples (wild-like, lab-like, CD1, MF1, plus pos and neg controls) 
# 
# ps_lab <- ps_all %>% 
#   ps_filter(SampleId %in% LabIds)
# ps_lab # 223 samples (2 low reads dropped), 1911 ASVs
# 
# ## abundance & prevalence
# taxa_sums(ps_lab) %>% as.data.frame() %>%
#   arrange(desc(.)) %>%
#   rename(TotalAbundance=1) -> TotalAbundances_Lab
# summary(TotalAbundances_Lab)
# 
# prevalence(ps_lab) %>% as.data.frame() %>%
#   arrange(desc(.)) %>%
#   rename(Prevalence=1) %>%
#   mutate(Prevalence=Prevalence*100)-> PrevalenceASVs_Lab
# summary(PrevalenceASVs_Lab)

# save all 
# saveRDS(ps_all, "Data/Phyloseqs/phyloseqAllSamples.rds")
# saveRDS(ps_rse, "Data/Phyloseqs/phyloseqRSESamples.rds")
# saveRDS(ps_wood, "Data/Phyloseqs/phyloseqWoodMiceSamples.rds")
# saveRDS(ps_lab, "Data/Phyloseqs/phyloseqLabSamples.rds")



# Add age at sampling to ps_all -----------------------------------------------

## read mouse metadata 
dob_files <- list(
  `antibiotic expt` = "Data/Metadata/ColonyFullPerturbationData/Antibiotic_Data_RB.xlsx",
  `infection expt`  = "Data/Metadata/ColonyFullPerturbationData/Infection_Data_RB.xlsx",
  `nutrition expt`  = "Data/Metadata/ColonyFullPerturbationData/Nutrition_Data_RB.xlsx"
)

mouse_metadata <- lapply(dob_files, read_xlsx)


## function: check metadata consistency and create consensus values 

process_mouse_metadata <- function(df, dataset_name) {
  
  # Keep required columns and create unique mouse ID
  df <- df %>%
    select(mouse_ID, dob, sex) %>%
    mutate(
      dataset = dataset_name,
      mouse_ID_dataset = paste(dataset_name, mouse_ID, sep = "_")
    )
  
  
  # Check whether DOB and sex are consistent within each mouse
  inconsistent <- df %>%
    group_by(mouse_ID_dataset, mouse_ID, dataset) %>%
    summarise(
      n_unique_dob = n_distinct(dob),
      dob_values = paste(unique(dob), collapse = ", "),
      n_unique_sex = n_distinct(sex),
      sex_values = paste(unique(sex), collapse = ", "),
      .groups = "drop"
    ) %>%
    filter(n_unique_dob > 1 | n_unique_sex > 1)
  
  
  # Add consensus DOB and sex (most frequent value per mouse)
  cleaned <- df %>%
    group_by(mouse_ID_dataset, mouse_ID, dataset) %>%
    mutate(
      dob_c = names(sort(table(as.character(dob)), decreasing = TRUE))[1],
      sex_c = names(sort(table(as.character(sex)), decreasing = TRUE))[1]
    ) %>%
    ungroup()
  
  
  list(
    cleaned = cleaned,
    inconsistent = inconsistent
  )
}


## Process all datasets

processed_metadata <- Map(
  process_mouse_metadata,
  mouse_metadata,
  names(mouse_metadata)
)


## view inconsistent mice 
## these should only show mice with conflicting DOB or sex within a dataset
processed_metadata$Antibiotic$inconsistent
processed_metadata$Infection$inconsistent
processed_metadata$Nutrition$inconsistent


## create one-row-per-mouse metadata table 
mouse_metadata_final <- bind_rows(
  lapply(processed_metadata, function(x) x$cleaned)
) %>%
  select(
    mouse_ID_dataset,
    mouse_ID,
    dataset,
    dob = dob_c,
    sex = sex_c
  ) %>%
  distinct()


## confirm every experimental mouse has only one DOB/sex 
mouse_metadata_final %>%
  dplyr::count(mouse_ID_dataset) %>%
  filter(n > 1)

## add matching ID to phyloseq metadata 
## extract phyloseq metadata and preserve sample names
ps_meta <- data.frame(sample_data(ps_all))
ps_meta$SampleID <- rownames(ps_meta)

## create dataset-specific ID in phyloseq metadata
ps_meta <- ps_meta %>%
  mutate(
    mouse_ID_dataset = paste(SampleDetailed, Id, sep = "_"))

## join DOB and sex information
ps_meta <- ps_meta %>%
  dplyr::left_join(
    mouse_metadata_final,
    by = "mouse_ID_dataset")


## calculate age at sampling
ps_meta <- ps_meta %>%
  mutate(
    SampleDate = as.Date(SampleDate, format = "%d/%m/%Y"),
    dob = as.Date(dob),
    age_sampling = as.numeric(SampleDate - dob))

## check for missing DOBs - should be only Wild:As mice
ps_meta %>%
  filter(is.na(dob)) %>%
  select(Id, SampleDetailed)

## check age values
summary(ps_meta$age_sampling)
hist(ps_meta$age_sampling)

ggplot(ps_meta[ps_meta$SampleOrigin!="Wild:As",], aes(x = age_sampling)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ SampleOrigin) +
  theme_classic() +
  labs(
    x = "Age at sampling (days)",
    y = "Number of samples")

## restore sample names and return metadata to phyloseq
rownames(ps_meta) <- ps_meta$SampleID
ps_meta$SampleID <- NULL

## confirm sample names match
stopifnot(all(rownames(ps_meta) == sample_names(ps_all)))

## add updated metadata back
sample_data(ps_all) <- sample_data(ps_meta)

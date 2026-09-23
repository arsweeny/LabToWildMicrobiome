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
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(ape)
library(openxlsx)
library(patchwork)

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

meta_merge <- meta_full %>% dplyr::select(sample_id, age, repro, body_mass_g, treated) %>% 
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



# dealing with sequencing controls and contaminants ------------------------

# ------------------------------------------------------------------------------
# QUALITY-CONTROL OUTPUT DIRECTORIES
# ------------------------------------------------------------------------------

qc_results_dir <- "Results/00_Quality_Control"
qc_tables_dir <- file.path(qc_results_dir, "Tables")
qc_figures_dir <- file.path(qc_results_dir, "Figures")

dir.create(qc_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_figures_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# DEFINE CONTROL TYPES
# ------------------------------------------------------------------------------

negative_types <- c(
  "water_primers",
  "water_only"
)

mock_types <- c(
  "mock_community"
)

golden_types <- c(
  "golden_colony",
  "golden_wild"
)

all_control_types <- c(
  negative_types,
  mock_types,
  golden_types
)

qc_meta <- data.frame(
  phyloseq::sample_data(ps),
  check.names = FALSE
) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::mutate(
    ControlType = dplyr::case_when(
      SampleType %in% negative_types ~ "Negative control",
      SampleType %in% mock_types ~ "Mock community",
      SampleType == "golden_colony" ~ "Golden colony",
      SampleType == "golden_wild" ~ "Golden wild",
      TRUE ~ "Biological sample"
    ),
    SeqDepth = as.numeric(
      phyloseq::sample_sums(ps)[SampleID]
    )
  )

qc_meta$ControlType <- factor(
  qc_meta$ControlType,
  levels = c(
    "Biological sample",
    "Negative control",
    "Mock community",
    "Golden colony",
    "Golden wild"
  )
)

# Per-sample sequencing-depth table.
control_sample_table <- qc_meta |>
  dplyr::select(
    SampleID,
    SampleType,
    ControlType,
    SeqDepth,
    dplyr::everything()
  )

# Compact summary for reporting.
control_summary_table <- qc_meta |>
  dplyr::group_by(
    ControlType
  ) |>
  dplyr::summarise(
    N_samples = dplyr::n(),
    Median_reads = stats::median(
      SeqDepth,
      na.rm = TRUE
    ),
    Min_reads = min(
      SeqDepth,
      na.rm = TRUE
    ),
    Max_reads = max(
      SeqDepth,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(control_summary_table)


# ------------------------------------------------------------------------------
# QC PANEL A: SEQUENCING DEPTH OF BIOLOGICAL AND CONTROL SAMPLES
# ------------------------------------------------------------------------------

plot_control_depth <- ggplot2::ggplot(
  qc_meta,
  ggplot2::aes(
    x = ControlType,
    y = SeqDepth
  )
) +
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  ggplot2::geom_jitter(
    width = 0.12,
    height = 0,
    alpha = 0.65,
    size = 1.8
  ) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = NULL,
    y = "Sequencing depth (reads; log10 scale)"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    )
  )



# ------------------------------------------------------------------------------
# IDENTIFY CONTAMINANTS USING TRUE NEGATIVE CONTROLS ONLY
# ------------------------------------------------------------------------------

biological_ids <- qc_meta |>
  dplyr::filter(
    ControlType == "Biological sample"
  ) |>
  dplyr::pull(
    SampleID
  )

negative_ids <- qc_meta |>
  dplyr::filter(
    ControlType == "Negative control"
  ) |>
  dplyr::pull(
    SampleID
  )

decontam_ids <- c(
  biological_ids,
  negative_ids
)

ps_decontam <- phyloseq::prune_samples(
  decontam_ids,
  ps
)

negative_indicator <- phyloseq::sample_names(
  ps_decontam
) %in% negative_ids

stopifnot(
  length(negative_indicator) ==
    phyloseq::nsamples(ps_decontam),
  sum(negative_indicator) ==
    length(negative_ids)
)

# NOTE: use %in% rather than == c(...), which would recycle the vector.
contams <- decontam::isContaminant(
  ps_decontam,
  method = "prevalence",
  neg = negative_indicator,
  threshold = 0.5
)

contaminants <- rownames(
  contams[
    contams$contaminant %in% TRUE,
    ,
    drop = FALSE
  ]
)

cat(
  "\nPotential contaminant ASVs identified by decontam: ",
  length(contaminants),
  "\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# DETAILED CONTAMINANT TABLE
# ------------------------------------------------------------------------------

otu_qc <- as(
  phyloseq::otu_table(ps_decontam),
  "matrix"
)

if (!phyloseq::taxa_are_rows(ps_decontam)) {
  otu_qc <- t(otu_qc)
}

bio_cols <- colnames(otu_qc) %in%
  biological_ids

neg_cols <- colnames(otu_qc) %in%
  negative_ids

safe_row_mean_presence <- function(
    x,
    keep_cols) {

  if (sum(keep_cols) == 0) {
    return(
      rep(
        NA_real_,
        nrow(x)
      )
    )
  }

  rowMeans(
    x[
      ,
      keep_cols,
      drop = FALSE
    ] > 0
  )
}

safe_row_sum <- function(
    x,
    keep_cols) {

  if (sum(keep_cols) == 0) {
    return(
      rep(
        NA_real_,
        nrow(x)
      )
    )
  }

  rowSums(
    x[
      ,
      keep_cols,
      drop = FALSE
    ]
  )
}

tax_qc <- as.data.frame(
  phyloseq::tax_table(
    ps_decontam
  ),
  stringsAsFactors = FALSE
)

tax_qc$TaxonID <- rownames(
  tax_qc
)

decontam_table <- data.frame(
  TaxonID = rownames(contams),
  Decontam_P = contams$p,
  Contaminant = contams$contaminant,
  Prevalence_negative = safe_row_mean_presence(
    otu_qc,
    neg_cols
  )[rownames(contams)],
  Prevalence_biological = safe_row_mean_presence(
    otu_qc,
    bio_cols
  )[rownames(contams)],
  Reads_negative = safe_row_sum(
    otu_qc,
    neg_cols
  )[rownames(contams)],
  Reads_biological = safe_row_sum(
    otu_qc,
    bio_cols
  )[rownames(contams)],
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(
    tax_qc,
    by = "TaxonID"
  ) |>
  dplyr::arrange(
    dplyr::desc(Contaminant),
    Decontam_P
  )

contaminant_table <- decontam_table |>
  dplyr::filter(
    Contaminant %in% TRUE
  )


# ------------------------------------------------------------------------------
# INTERNAL QC ONLY: PREVALENCE IN NEGATIVE CONTROLS VS BIOLOGICAL SAMPLES
# ------------------------------------------------------------------------------

plot_decontam_prevalence <- ggplot2::ggplot(
  decontam_table,
  ggplot2::aes(
    x = Prevalence_biological,
    y = Prevalence_negative,
    shape = Contaminant
  )
) +
  ggplot2::geom_point(
    alpha = 0.65,
    size = 2
  ) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1)
  ) +
  ggplot2::labs(
    x = "Prevalence in biological samples",
    y = "Prevalence in negative controls",
    shape = "decontam contaminant"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  )



# ------------------------------------------------------------------------------
# MOCK AND GOLDEN REFERENCE SAMPLE QC
# ------------------------------------------------------------------------------
# These controls are evaluated descriptively rather than used by decontam.
# Without an external file describing the expected mock composition, the script
# assesses observed taxonomic profiles and repeatability but does not claim
# quantitative recovery accuracy relative to the manufacturer's target mixture.

positive_control_ids <- qc_meta |>
  dplyr::filter(
    ControlType %in% c(
      "Mock community",
      "Golden colony",
      "Golden wild"
    )
  ) |>
  dplyr::pull(
    SampleID
  )

positive_control_bray <- data.frame()
positive_control_bray_summary <- data.frame()

if (length(positive_control_ids) > 0) {

  ps_positive_controls <- phyloseq::prune_samples(
    positive_control_ids,
    ps
  )

  ps_positive_controls <- phyloseq::prune_taxa(
    phyloseq::taxa_sums(
      ps_positive_controls
    ) > 0,
    ps_positive_controls
  )

  positive_meta <- data.frame(
    phyloseq::sample_data(
      ps_positive_controls
    ),
    check.names = FALSE
  ) |>
    tibble::rownames_to_column(
      "SampleID"
    ) |>
    dplyr::mutate(
      ControlType = dplyr::case_when(
        SampleType %in% mock_types ~ "Mock community",
        SampleType == "golden_colony" ~ "Golden colony",
        SampleType == "golden_wild" ~ "Golden wild",
        TRUE ~ as.character(SampleType)
      )
    )

  # --------------------------------------------------------------------------
  # QC PANEL B: GENUS-LEVEL COMPOSITION OF MOCK / GOLDEN CONTROLS
  # --------------------------------------------------------------------------

  ps_positive_genus <- phyloseq::tax_glom(
    ps_positive_controls,
    taxrank = "Genus",
    NArm = FALSE
  )

  ps_positive_genus_rel <- phyloseq::transform_sample_counts(
    ps_positive_genus,
    function(x) {
      x / sum(x)
    }
  )

  positive_genus_long <- phyloseq::psmelt(
    ps_positive_genus_rel
  ) |>
    dplyr::mutate(
      Genus = dplyr::if_else(
        is.na(Genus) |
          Genus == "",
        "Unclassified",
        as.character(Genus)
      ),
      ControlType = dplyr::case_when(
        SampleType %in% mock_types ~ "Mock community",
        SampleType == "golden_colony" ~ "Golden colony",
        SampleType == "golden_wild" ~ "Golden wild",
        TRUE ~ as.character(SampleType)
      )
    )

  top_control_genera <- positive_genus_long |>
    dplyr::group_by(
      Genus
    ) |>
    dplyr::summarise(
      MeanAbundance = mean(
        Abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(
        MeanAbundance
      )
    ) |>
    dplyr::slice_head(
      n = 15
    ) |>
    dplyr::pull(
      Genus
    )

  positive_genus_plot_data <- positive_genus_long |>
    dplyr::mutate(
      Genus_plot = dplyr::if_else(
        Genus %in% top_control_genera,
        Genus,
        "Other"
      )
    ) |>
    dplyr::group_by(
      Sample,
      ControlType,
      Genus_plot
    ) |>
    dplyr::summarise(
      RelativeAbundance = sum(
        Abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  plot_positive_controls <- ggplot2::ggplot(
    positive_genus_plot_data,
    ggplot2::aes(
      x = Sample,
      y = RelativeAbundance,
      fill = Genus_plot
    )
  ) +
    ggplot2::geom_col(
      width = 0.85
    ) +
    ggplot2::facet_grid(
      . ~ ControlType,
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(
        accuracy = 1
      ),
      expand = c(
        0,
        0
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Relative abundance",
      fill = "Genus"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      ),
      legend.position = "right"
    )



  # --------------------------------------------------------------------------
  # WITHIN-CONTROL BRAY-CURTIS DISTANCES
  # --------------------------------------------------------------------------
  # For control types represented by >=2 samples, calculate all within-type
  # pairwise Bray-Curtis distances. Lower values indicate more reproducible
  # technical/reference community profiles.

  positive_control_bray_list <- list()

  control_levels_present <- unique(
    positive_meta$ControlType
  )

  for (control_name in control_levels_present) {

    sample_ids_control <- positive_meta |>
      dplyr::filter(
        ControlType == control_name
      ) |>
      dplyr::pull(
        SampleID
      )

    if (length(sample_ids_control) >= 2) {

      ps_control_type <- phyloseq::prune_samples(
        sample_ids_control,
        ps_positive_controls
      )

      ps_control_type <- phyloseq::prune_taxa(
        phyloseq::taxa_sums(
          ps_control_type
        ) > 0,
        ps_control_type
      )

      ps_control_type_rel <- phyloseq::transform_sample_counts(
        ps_control_type,
        function(x) {
          x / sum(x)
        }
      )

      bray_dist <- phyloseq::distance(
        ps_control_type_rel,
        method = "bray"
      )

      bray_mat <- as.matrix(
        bray_dist
      )

      bray_pairs <- as.data.frame(
        as.table(
          bray_mat
        ),
        stringsAsFactors = FALSE
      ) |>
        dplyr::rename(
          Sample1 = Var1,
          Sample2 = Var2,
          BrayCurtis = Freq
        ) |>
        dplyr::filter(
          Sample1 < Sample2
        ) |>
        dplyr::mutate(
          ControlType = control_name,
          .before = 1
        )
      
      positive_control_bray_list[[control_name]] <- bray_pairs
    }
  }
  
  if (length(positive_control_bray_list) > 0) {
    
    positive_control_bray <- dplyr::bind_rows(
      positive_control_bray_list
    )
    
    positive_control_bray_summary <- positive_control_bray |>
      dplyr::group_by(
        ControlType
      ) |>
      dplyr::summarise(
        N_pairwise_comparisons = dplyr::n(),
        Mean_BrayCurtis = mean(
          BrayCurtis,
          na.rm = TRUE
        ),
        Median_BrayCurtis = stats::median(
          BrayCurtis,
          na.rm = TRUE
        ),
        SD_BrayCurtis = stats::sd(
          BrayCurtis,
          na.rm = TRUE
        ),
        Min_BrayCurtis = min(
          BrayCurtis,
          na.rm = TRUE
        ),
        Max_BrayCurtis = max(
          BrayCurtis,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
  }
}


# ------------------------------------------------------------------------------
# QC OUTPUTS ARE WRITTEN ONCE, AT THE END OF THE SCRIPT
# ------------------------------------------------------------------------------
# The important QC tables remain in memory here. A single multi-sheet workbook
# is created after final sample filtering so it can also include low-read sample
# exclusions and the final analysis-dataset summary. Detailed intermediate
# tables (all decontam calls, per-sample control metadata, and every pairwise
# control distance) are intentionally not exported.


# ------------------------------------------------------------------------------
# SAVE ONE COMBINED QC FIGURE
# ------------------------------------------------------------------------------
# Only the two reader-facing QC panels are saved:
#   A. sequencing depth across biological/control sample types
#   B. genus-level profiles of mock and golden reference controls
#
# The decontam prevalence scatter is retained as an internal plot object but is
# not exported as a separate figure.

if (exists("plot_positive_controls")) {

  combined_qc_figure <- (
    plot_control_depth
  ) /
    (
      plot_positive_controls
    ) +
    patchwork::plot_layout(
      heights = c(0.85, 1.25)
    ) +
    patchwork::plot_annotation(
      tag_levels = "A"
    )

  ggplot2::ggsave(
    filename = file.path(
      qc_figures_dir,
      "Figure_S_QC_controls.pdf"
    ),
    plot = combined_qc_figure,
    width = 11,
    height = 11
  )
}

# ------------------------------------------------------------------------------
# REMOVE CONTAMINANT TAXA AND ALL CONTROL SAMPLES
# ------------------------------------------------------------------------------

clean_ps <- phyloseq::prune_taxa(
  !(
    phyloseq::taxa_names(ps) %in%
      contaminants
  ),
  ps
)

clean_ps <- phyloseq::prune_samples(
  !(
    phyloseq::sample_data(clean_ps)$SampleType %in%
      all_control_types
  ),
  clean_ps
)

cat(
  "\nAfter contaminant and control removal:\n",
  "  biological samples = ",
  phyloseq::nsamples(clean_ps),
  "\n  ASVs = ",
  phyloseq::ntaxa(clean_ps),
  "\n",
  sep = ""
)

clean_ps


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

# ------------------------------------------------------------------------------
# CHECK BIOLOGICAL SAMPLES REMOVED FOR LOW SEQUENCING DEPTH
# ------------------------------------------------------------------------------

# Remove colony controls first, but do NOT apply the read-depth filter yet
ps_before_depth_filter <- ps3 %>%
  ps_filter(SampleDetailed != "colony control")

# Get sequencing depth for each sample
depth_check <- data.frame(
  SampleID = phyloseq::sample_names(ps_before_depth_filter),
  Reads = phyloseq::sample_sums(ps_before_depth_filter),
  SampleType = phyloseq::sample_data(ps_before_depth_filter)$SampleType,
  stringsAsFactors = FALSE
)

# Identify biological samples below 10,000 reads
low_read_biological <- depth_check %>%
  dplyr::filter(
    !(SampleType %in% c(
      "water_primers",
      "water_only",
      "mock_community",
      "golden_colony",
      "golden_wild"
    )),
    Reads < 10000
  )

# How many?
nrow(low_read_biological)

# Which samples?
low_read_biological

# Summary of their read depths
summary(low_read_biological$Reads)

# Mean and range
mean(low_read_biological$Reads)
range(low_read_biological$Reads)

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



# phylogenetic tree is constructed after all final sample filtering


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
    dplyr::select(mouse_ID, dob, sex) %>%
    dplyr::mutate(
      dataset = dataset_name,
      mouse_ID_dataset = paste(dataset_name, mouse_ID, sep = "_")
    )
  
  # Check whether DOB and sex are consistent within each mouse
  inconsistent <- df %>%
    dplyr::group_by(mouse_ID_dataset, mouse_ID, dataset) %>%
    dplyr::summarise(
      n_unique_dob = dplyr::n_distinct(dob),
      dob_values = paste(unique(dob), collapse = ", "),
      n_unique_sex = dplyr::n_distinct(sex),
      sex_values = paste(unique(sex), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_unique_dob > 1 | n_unique_sex > 1)
  
  # Add consensus DOB and sex (most frequent value per mouse)
  cleaned <- df %>%
    dplyr::group_by(mouse_ID_dataset, mouse_ID, dataset) %>%
    dplyr::mutate(
      dob_c = names(sort(table(as.character(dob)), decreasing = TRUE))[1],
      sex_c = names(sort(table(as.character(sex)), decreasing = TRUE))[1]
    ) %>%
    dplyr::ungroup()
  
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
  dplyr::select(
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
  dplyr::select(Id, SampleDetailed)

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


# ------------------------------------------------------------------------------
# BUILD, SAVE, AND ATTACH THE FINAL PHYLOGENETIC TREE
# ------------------------------------------------------------------------------
# The tree is constructed here, AFTER all final sample filtering, so the saved
# tree and the saved ps_all object correspond to the same final ASV set.
# Downstream scripts should load this saved tree rather than rebuilding it.

processed_data_dir <- "Data/Processed"
dir.create(
  processed_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Remove any ASVs that became zero-abundance after the final sample filtering.
ps_all <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_all) > 0,
  ps_all
)

# Match the ASV sequence key to the exact ASVs retained in final ps_all.
asv_key_final <- asv_key %>%
  dplyr::filter(
    ASV %in% phyloseq::taxa_names(ps_all)
  ) %>%
  dplyr::arrange(
    match(
      ASV,
      phyloseq::taxa_names(ps_all)
    )
  )

stopifnot(
  setequal(
    asv_key_final$ASV,
    phyloseq::taxa_names(ps_all)
  )
)

# Convert ASV sequences to a DNAStringSet and preserve the ASV names as tips.
seqs <- Biostrings::DNAStringSet(
  asv_key_final$Sequence
)

names(seqs) <- asv_key_final$ASV

# Multiple-sequence alignment.
alignment <- DECIPHER::AlignSeqs(
  seqs
)

# Convert the alignment to phangorn format.
phang_align <- phangorn::phyDat(
  as(
    alignment,
    "matrix"
  ),
  type = "DNA"
)

# Build a neighbour-joining starting tree from maximum-likelihood distances.
# ape exports this function as lowercase nj().
dm <- phangorn::dist.ml(
  phang_align
)

tree_nj <- ape::nj(
  dm
)

# Fit and optimize a maximum-likelihood tree using a GTR model and NNI search.
fit <- phangorn::pml(
  tree_nj,
  phang_align
)

fit_gtr <- phangorn::optim.pml(
  fit,
  model = "GTR",
  optInv = TRUE,
  rearrangement = "NNI"
)

# Extract and midpoint-root the optimized tree.
tree_rooted <- phangorn::midpoint(
  fit_gtr$tree
)

# Ensure the final tree contains exactly the ASVs in ps_all.
tree_rooted <- ape::keep.tip(
  tree_rooted,
  phyloseq::taxa_names(ps_all)
)

stopifnot(
  setequal(
    tree_rooted$tip.label,
    phyloseq::taxa_names(ps_all)
  ),
  ape::is.rooted(tree_rooted)
)

# Attach the tree to the final phyloseq object BEFORE ps_all is saved.
phyloseq::phy_tree(ps_all) <- tree_rooted

stopifnot(
  !is.null(
    phyloseq::phy_tree(
      ps_all,
      errorIfNULL = FALSE
    )
  ),
  setequal(
    phyloseq::phy_tree(ps_all)$tip.label,
    phyloseq::taxa_names(ps_all)
  )
)

# Save the tree separately so it can be versioned in the GitHub project and
# loaded directly by downstream analyses. Newick is also written as a portable,
# human-readable representation of the same tree.
tree_rds_file <- file.path(
  processed_data_dir,
  "baseline_phylogenetic_tree.rds"
)

tree_newick_file <- file.path(
  processed_data_dir,
  "baseline_phylogenetic_tree.nwk"
)

saveRDS(
  tree_rooted,
  file = tree_rds_file
)

ape::write.tree(
  tree_rooted,
  file = tree_newick_file
)

cat(
  "\nSaved final phylogenetic tree to: ",
  tree_rds_file,
  "\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# FINAL DATASET SUMMARY FOR MANUSCRIPT
# ------------------------------------------------------------------------------

## sequencing depth per retained sample
final_depths <- phyloseq::sample_sums(ps_all)

## total reads
total_reads <- sum(final_depths)

## mean reads per sample
mean_reads <- mean(final_depths)

## range
min_reads <- min(final_depths)
max_reads <- max(final_depths)

## total number of samples
n_samples <- phyloseq::nsamples(ps_all)

## sample numbers by cohort
group_counts <- table(
  phyloseq::sample_data(ps_all)$SampleOrigin
)

## print everything
cat(
  "\nFINAL DATASET SUMMARY\n",
  "Total reads: ", format(total_reads, big.mark = ","), "\n",
  "Mean reads/sample: ", format(round(mean_reads), big.mark = ","), "\n",
  "Range: ", format(min_reads, big.mark = ","),
  " - ", format(max_reads, big.mark = ","), "\n",
  "Number of samples: ", n_samples, "\n\n",
  sep = ""
)

print(group_counts)


# ------------------------------------------------------------------------------
# SAVE ONE STREAMLINED MULTI-SHEET RESULTS WORKBOOK
# ------------------------------------------------------------------------------

# 1. Control sequencing-depth summary.
control_summary_export <- control_summary_table |>
  dplyr::mutate(
    ControlType = as.character(ControlType)
  ) |>
  dplyr::rename(
    Control_type = ControlType,
    N = N_samples,
    Median_reads = Median_reads,
    Min_reads = Min_reads,
    Max_reads = Max_reads
  )

# 2. Only ASVs actually classified as potential contaminants.
#    Retain the sequence/ID, read totals, decontam score, prevalence, and taxonomy.
contaminants_export <- contaminant_table |>
  dplyr::select(
    TaxonID,
    Reads_negative,
    Reads_biological,
    Decontam_P,
    Prevalence_negative,
    Prevalence_biological,
    dplyr::any_of(
      c(
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus"
      )
    )
  )

# 3. Compact reproducibility summary for mock/golden controls.
control_repro_export <- positive_control_bray_summary |>
  dplyr::rename(
    Control_type = ControlType,
    N_pairwise = N_pairwise_comparisons,
    Mean_Bray_Curtis = Mean_BrayCurtis,
    Median_Bray_Curtis = Median_BrayCurtis,
    SD_Bray_Curtis = SD_BrayCurtis,
    Min_Bray_Curtis = Min_BrayCurtis,
    Max_Bray_Curtis = Max_BrayCurtis
  )

# 4. Biological samples excluded for low sequencing depth.
low_read_export <- low_read_biological |>
  dplyr::select(
    SampleID,
    Reads,
    SampleType
  ) |>
  dplyr::arrange(
    Reads
  )

# 5. Final analysis-dataset summary for direct manuscript reporting.
final_sample_metadata <- data.frame(
  phyloseq::sample_data(ps_all),
  check.names = FALSE
) |>
  tibble::rownames_to_column(
    "SampleID"
  ) |>
  dplyr::mutate(
    Reads = as.numeric(
      phyloseq::sample_sums(ps_all)[SampleID]
    ),
    SampleOrigin = as.character(SampleOrigin)
  )

final_cohort_summary <- final_sample_metadata |>
  dplyr::group_by(
    SampleOrigin
  ) |>
  dplyr::summarise(
    N_samples = dplyr::n(),
    Total_reads = sum(Reads),
    Mean_reads = mean(Reads),
    Median_reads = stats::median(Reads),
    Min_reads = min(Reads),
    Max_reads = max(Reads),
    .groups = "drop"
  )

final_dataset_overall <- tibble::tibble(
  SampleOrigin = "Overall",
  N_samples = phyloseq::nsamples(ps_all),
  Total_reads = sum(final_depths),
  Mean_reads = mean(final_depths),
  Median_reads = stats::median(final_depths),
  Min_reads = min(final_depths),
  Max_reads = max(final_depths)
)

final_dataset_export <- dplyr::bind_rows(
  final_cohort_summary,
  final_dataset_overall
) |>
  dplyr::mutate(
    N_ASVs = dplyr::if_else(
      SampleOrigin == "Overall",
      as.integer(phyloseq::ntaxa(ps_all)),
      NA_integer_
    )
  )

qc_workbook <- list(
  "Control Summary" = control_summary_export,
  "Contaminants" = contaminants_export,
  "Control Reproducibility" = control_repro_export,
  "Low Read Samples" = low_read_export,
  "Final Dataset" = final_dataset_export
)

openxlsx::write.xlsx(
  qc_workbook,
  file = file.path(
    qc_tables_dir,
    "Baseline_QC_results.xlsx"
  ),
  overwrite = TRUE,
  rowNames = FALSE
)

cat(
  "\nSaved streamlined QC workbook to: ",
  file.path(
    qc_tables_dir,
    "Baseline_QC_results.xlsx"
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# SAVE THE FINAL ANALYSIS OBJECT
# ------------------------------------------------------------------------------
# All downstream scripts load this object directly instead of re-running the
# complete preprocessing workflow. Re-run this script only when the raw data,
# metadata, filtering decisions, taxonomy, or phylogenetic tree change.
saveRDS(
  ps_all,
  file = file.path(processed_data_dir, "baseline_microbiome_phyloseq.rds")
)

readr::write_csv(
  asv_key_final,
  file.path(processed_data_dir, "baseline_ASV_sequence_key.csv")
)

cat(
  "\nSaved final baseline phyloseq object (with attached tree) to: ",
  file.path(processed_data_dir, "baseline_microbiome_phyloseq.rds"),
  "\n",
  sep = ""
)

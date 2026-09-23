# ==============================================================================
# PICRUSt2 KEGG PATHWAY ANALYSIS — BASELINE MOUSE MICROBIOME
# ==============================================================================
#
# PURPOSE
# -------
# This script converts PICRUSt2-predicted KEGG Ortholog (KO) abundances into
# KEGG pathway profiles and tests whether predicted functional composition
# differs among four microbiome origins:
#
#   Wild:As       = wild Apodemus sylvaticus
#   Wild-like:As  = colony A. sylvaticus with a wild-like microbiota
#   Lab-like:As   = colony A. sylvaticus with a laboratory-like microbiota
#   Lab:Mm        = laboratory Mus musculus
#
# The analysis is intentionally focused on outputs that are useful for the
# manuscript and for reproducibility:
#
#   1. weighted NSTI prediction-quality analysis
#   2. KO -> KEGG pathway conversion and QC
#   3. >=10% pathway-prevalence filtering
#   4. KEGG functional beta diversity (Bray-Curtis)
#   5. omnibus and pairwise PERMANOVA + PERMDISP
#   6. sequencing-depth sensitivity analysis
#   7. pathway-prevalence sensitivity analysis for functional PERMANOVA
#   8. differential KEGG pathway abundance with MaAsLin3
#   9. concise descriptive summaries of significant pathways
#  10. publication-oriented figures and tables
#
# IMPORTANT INTERPRETATION
# ------------------------
# PICRUSt2 predicts gene-family content from 16S phylogenetic placement and
# reference genomes. KEGG pathways here therefore represent PREDICTED FUNCTIONAL
# POTENTIAL, not directly measured gene content, transcription, metabolites, or
# pathway activity.
#
# KEGG pathway names can also be misleading when interpreted literally. Bacterial
# KOs can map to pathways named for eukaryotic signaling, human disease, or
# particular organisms. For that reason this script does NOT select or interpret
# individual pathways solely because they have large effect sizes or appealing
# names. Full annotated differential-abundance results are saved for transparent
# reporting, while the main figures emphasize overall functional composition and
# the magnitude of between-group differences.
#
# ANALYSIS DECISIONS
# ------------------
# - KEGG pathway abundance is generated from PICRUSt2 KO predictions with
#   ggpicrust2::ko2kegg_abundance().
# - Pathways must be predicted in >=10% of samples for the PRIMARY analysis.
# - Robustness of functional beta diversity to this filtering choice is checked
#   across pathway-prevalence thresholds from 0-20%.
# - MaAsLin3 is NOT rerun across prevalence thresholds.
# - Relative pathway abundance is used for Bray-Curtis beta diversity.
# - Beta diversity uses one sample per mouse because SampleOrigin is a
#   between-mouse predictor.
# - Colony groups use the sample closest to Timept = 0; wild mice use the
#   earliest SampleDate.
# - Differential abundance uses all eligible samples and accounts for repeated
#   observations with (1 | mouse_ID_dataset).
# - Sex and log10 original 16S sequencing depth are included as fixed effects.
# - MaAsLin3 tests the abundance component only.
# - Significance is defined as BH-FDR q < 0.05.
# - KEGG Level 2 summaries are DESCRIPTIVE only; no category-enrichment claim is
#   made.
#
# REQUIRED INPUTS
# ---------------
# Data/Processed/baseline_microbiome_phyloseq.rds
# Results/PICRUSt_Outputs/results/KO_metagenome_out/
#   pred_metagenome_unstrat.tsv.gz
# Results/PICRUSt_Outputs/results/EC_metagenome_out/
#   weighted_nsti.tsv.gz
#
# OUTPUTS
# -------
# Results/05_PICRUSt2_Functional_Prediction/
#   Tables/
#   Figures/
# ==============================================================================


# ==============================================================================
# 0. SETUP
# ==============================================================================

required_packages <- c(
  "phyloseq",
  "ggpicrust2",
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "vegan",
  "lmerTest",
  "emmeans",
  "maaslin3",
  "patchwork",
  "cowplot",
  "openxlsx"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

cat(
  "\nggpicrust2 version: ",
  as.character(utils::packageVersion("ggpicrust2")),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------

phyloseq_file <- "Data/Processed/baseline_microbiome_phyloseq.rds"

ko_file <- file.path(
  "Results",
  "PICRUSt_Outputs",
  "results",
  "KO_metagenome_out",
  "pred_metagenome_unstrat.tsv.gz"
)


nsti_file <- file.path(
  "Results",
  "PICRUSt_Outputs",
  "results",
  "EC_metagenome_out",
  "weighted_nsti.tsv.gz"
)

if (!file.exists(phyloseq_file)) {
  stop("Could not find processed phyloseq object: ", phyloseq_file)
}

if (!file.exists(ko_file)) {
  stop("Could not find PICRUSt2 KO table: ", ko_file)
}

if (!file.exists(nsti_file)) {
  stop("Could not find PICRUSt2 weighted NSTI table: ", nsti_file)
}

# ------------------------------------------------------------------------------
# Output directories
# ------------------------------------------------------------------------------

output_dir <- "Results/05_PICRUSt2_Functional_Prediction"

tables_dir <- file.path(output_dir, "Tables")
figures_dir <- file.path(output_dir, "Figures")

# MaAsLin3 writes several intermediate files while fitting. These are required
# during the run but are not reader-facing outputs, so keep them in the R
# session temporary directory rather than under Results/.
maaslin_dir <- file.path(
  tempdir(),
  "PICRUSt2_MaAsLin3"
)

for (x in c(tables_dir, figures_dir, maaslin_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

# Output policy:
#   - save ONE multi-sheet workbook with the key PICRUSt2 results;
#   - save only COMPLETE figures, never individual panels from multipanel figures;
#   - retain detailed intermediate tables/plot objects in memory for analysis;
#   - remove temporary MaAsLin3 output after the final workbook/figures are saved.

# ------------------------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------------------------

minimum_pathway_prevalence <- 0.10

# Sensitivity analysis for the pathway-prevalence filter used in the
# functional beta-diversity PERMANOVA. The 10% threshold remains primary.
pathway_prevalence_sweep <- c(
  0,
  0.025,
  0.05,
  0.075,
  0.10,
  0.125,
  0.15,
  0.20
)

q_threshold <- 0.05
n_permutations <- 9999
random_seed <- 123

# MaAsLin3 must be rerun to create an in-memory fit object for contrast testing.
# TRUE will replace the existing KEGG MaAsLin3 model-output directory.
rerun_maaslin3 <- TRUE

group_levels <- c(
  "Wild:As",
  "Wild-like:As",
  "Lab-like:As",
  "Lab:Mm"
)

cohort_labels <- c(
  "Wild:As" = "Wild",
  "Wild-like:As" = "Wild-like",
  "Lab-like:As" = "Lab-like",
  "Lab:Mm" = "Lab"
)

contrast_names <- c(
  "Wild-like vs Wild",
  "Lab-like vs Wild",
  "Lab vs Wild",
  "Lab-like vs Wild-like",
  "Lab vs Wild-like",
  "Lab vs Lab-like"
)


# ------------------------------------------------------------------------------
# Consistent comparison display order for all pairwise figures
# ------------------------------------------------------------------------------
# Statistical contrasts are kept in their original MaAsLin3/PERMANOVA direction
# internally. For plotting only, comparisons are re-labelled so the wilder cohort
# is always named first along the ecological spectrum:
# Wild -> Wild-like -> Lab-like -> Lab.

directional_contrast_key <- tibble::tribble(
  ~display_test,              ~source_test,                  ~multiplier,
  "Wild vs Wild-like",        "Wild-like vs Wild",          -1,
  "Wild vs Lab-like",         "Lab-like vs Wild",           -1,
  "Wild vs Lab",              "Lab vs Wild",                -1,
  "Wild-like vs Lab-like",    "Lab-like vs Wild-like",      -1,
  "Wild-like vs Lab",         "Lab vs Wild-like",           -1,
  "Lab-like vs Lab",          "Lab vs Lab-like",            -1
)

directional_contrast_order <- c(
  "Wild vs Wild-like",
  "Wild vs Lab-like",
  "Wild vs Lab",
  "Wild-like vs Lab-like",
  "Wild-like vs Lab",
  "Lab-like vs Lab"
)

comparison_display_key <- directional_contrast_key |>
  dplyr::select(
    source_test,
    display_test
  )

comparison_display_labels <- c(
  "Wild vs Wild-like" = "Wild – Wild-like",
  "Wild vs Lab-like" = "Wild – Lab-like",
  "Wild vs Lab" = "Wild – Lab",
  "Wild-like vs Lab-like" = "Wild-like – Lab-like",
  "Wild-like vs Lab" = "Wild-like – Lab",
  "Lab-like vs Lab" = "Lab-like – Lab"
)

contrast_key <- tibble::tribble(
  ~test,                    ~first_group,   ~second_group,
  "Wild-like vs Wild",     "Wild-like:As", "Wild:As",
  "Lab-like vs Wild",      "Lab-like:As",  "Wild:As",
  "Lab vs Wild",           "Lab:Mm",       "Wild:As",
  "Lab-like vs Wild-like", "Lab-like:As",  "Wild-like:As",
  "Lab vs Wild-like",      "Lab:Mm",       "Wild-like:As",
  "Lab vs Lab-like",       "Lab:Mm",       "Lab-like:As"
)

# Consistent compact P-value formatting for tables/figures.
format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}


# ==============================================================================
# 1. LOAD PHYLOSEQ DATA AND SAMPLE METADATA
# ==============================================================================

ps_all <- readRDS(phyloseq_file)

metadata_kegg <- data.frame(
  phyloseq::sample_data(ps_all),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata_kegg$SampleID <- rownames(metadata_kegg)

# Original 16S sequencing depth, used as a sensitivity/model covariate.
sequencing_depth <- phyloseq::sample_sums(ps_all)

metadata_kegg$reads <- sequencing_depth[metadata_kegg$SampleID]
metadata_kegg$log10_reads <- log10(metadata_kegg$reads)

metadata_kegg$SampleOrigin <- factor(
  metadata_kegg$SampleOrigin,
  levels = group_levels
)

metadata_kegg$mouse_ID_dataset <- as.character(
  metadata_kegg$mouse_ID_dataset
)

stopifnot(!anyDuplicated(metadata_kegg$SampleID))
stopifnot(!anyNA(metadata_kegg$SampleOrigin))
stopifnot(!anyNA(metadata_kegg$mouse_ID_dataset))
stopifnot(!anyNA(metadata_kegg$reads))
stopifnot(!anyNA(metadata_kegg$log10_reads))

cat("\nSamples in processed phyloseq object:\n")
print(
  metadata_kegg |>
    dplyr::count(SampleOrigin, name = "n_samples")
)



# ==============================================================================
# 2. WEIGHTED NSTI PREDICTION-QUALITY ANALYSIS
# ==============================================================================
#
# Weighted nearest-sequenced-taxon index (NSTI) summarizes the phylogenetic
# distance between each sample's observed taxa and reference genomes used by
# PICRUSt2. Lower values indicate closer reference-genome representation.
#
# This section is transferred from the previous functional-analysis script. It is retained
# here because NSTI is a PICRUSt2 prediction-quality diagnostic and therefore
# applies to the KEGG functional predictions as well.
# ==============================================================================

nsti <- readr::read_tsv(
  nsti_file,
  show_col_types = FALSE
)

required_nsti_columns <- c(
  "sample",
  "weighted_NSTI"
)

if (!all(required_nsti_columns %in% names(nsti))) {
  stop(
    "weighted_nsti.tsv.gz does not contain the expected columns: ",
    paste(required_nsti_columns, collapse = ", ")
  )
}

if (anyDuplicated(nsti$sample)) {
  stop("Duplicated sample IDs detected in weighted_nsti.tsv.gz.")
}

# Restrict NSTI results to the samples in the final processed phyloseq object.
nsti <- nsti |>
  dplyr::filter(
    sample %in% metadata_kegg$SampleID
  )

if (!setequal(nsti$sample, metadata_kegg$SampleID)) {
  missing_nsti <- setdiff(
    metadata_kegg$SampleID,
    nsti$sample
  )

  stop(
    "Weighted NSTI values are missing for ",
    length(missing_nsti),
    " sample(s) in the final phyloseq object: ",
    paste(missing_nsti, collapse = ", ")
  )
}

nsti_meta <- nsti |>
  dplyr::left_join(
    metadata_kegg |>
      dplyr::select(
        SampleID,
        SampleOrigin,
        mouse_ID_dataset
      ),
    by = c("sample" = "SampleID")
  ) |>
  dplyr::mutate(
    SampleOrigin = factor(
      SampleOrigin,
      levels = group_levels
    ),
    mouse_ID_dataset = factor(
      mouse_ID_dataset
    )
  )

if (anyNA(nsti_meta$SampleOrigin)) {
  stop("At least one NSTI sample failed to receive SampleOrigin metadata.")
}

if (anyNA(nsti_meta$mouse_ID_dataset)) {
  stop("At least one NSTI sample failed to receive mouse_ID_dataset metadata.")
}

nsti_summary <- nsti_meta |>
  dplyr::group_by(
    SampleOrigin
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_NSTI = mean(
      weighted_NSTI,
      na.rm = TRUE
    ),
    sd_NSTI = stats::sd(
      weighted_NSTI,
      na.rm = TRUE
    ),
    median_NSTI = stats::median(
      weighted_NSTI,
      na.rm = TRUE
    ),
    min_NSTI = min(
      weighted_NSTI,
      na.rm = TRUE
    ),
    max_NSTI = max(
      weighted_NSTI,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Account for repeated samples from the same mouse.
nsti_model <- lmerTest::lmer(
  weighted_NSTI ~
    SampleOrigin +
    (1 | mouse_ID_dataset),
  data = nsti_meta
)

nsti_anova <- stats::anova(
  nsti_model
)

nsti_emmeans <- emmeans::emmeans(
  nsti_model,
  pairwise ~ SampleOrigin,
  adjust = "tukey"
)

nsti_emmeans_table <- as.data.frame(
  nsti_emmeans$emmeans
)

nsti_pairwise_table <- as.data.frame(
  nsti_emmeans$contrasts
)

nsti_anova_table <- data.frame(
  term = rownames(nsti_anova),
  nsti_anova,
  row.names = NULL,
  check.names = FALSE
)

cat("\nOverall weighted NSTI summary:\n")
print(
  summary(nsti_meta$weighted_NSTI)
)

cat("\nWeighted NSTI by SampleOrigin:\n")
print(
  nsti_summary
)

cat("\nNSTI mixed-model ANOVA:\n")
print(
  nsti_anova
)

cat("\nNSTI estimated marginal means and Tukey contrasts:\n")
print(
  nsti_emmeans
)

# Save NSTI QC/statistical results.




# Standalone NSTI QC figure.
p_nsti <- ggplot2::ggplot(
  nsti_meta,
  ggplot2::aes(
    x = SampleOrigin,
    y = weighted_NSTI
  )
) +
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    width = 0.55
  ) +
  ggplot2::geom_jitter(
    width = 0.15,
    alpha = 0.5,
    size = 1.5
  ) +
  ggplot2::scale_x_discrete(
    labels = cohort_labels
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Weighted NSTI"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  ) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(
      colour = "black"
    )
  )



# ==============================================================================
# 3. READ AND VALIDATE PICRUSt2 KO PREDICTIONS
# ==============================================================================

ko_raw <- readr::read_tsv(
  ko_file,
  show_col_types = FALSE
)

if (ncol(ko_raw) < 2) {
  stop("PICRUSt2 KO table must contain one KO-ID column plus sample columns.")
}

# ggpicrust2 expects KO identifiers in the FIRST COLUMN, not in row names.
ko_for_kegg <- as.data.frame(
  ko_raw,
  check.names = FALSE
)

ko_ids <- as.character(ko_for_kegg[[1]])
ko_sample_ids <- names(ko_for_kegg)[-1]

if (!all(grepl("^(ko:)?K[0-9]{5}$", ko_ids))) {
  stop("Unexpected KO identifier format detected in the first KO-table column.")
}

if (anyDuplicated(ko_ids)) {
  stop("Duplicated KO identifiers detected in PICRUSt2 KO table.")
}

ko_abundance_matrix <- as.matrix(
  ko_for_kegg[, -1, drop = FALSE]
)

storage.mode(ko_abundance_matrix) <- "numeric"

if (anyNA(ko_abundance_matrix)) {
  stop("Missing values detected in PICRUSt2 KO abundances.")
}

if (!all(is.finite(ko_abundance_matrix))) {
  stop("Non-finite values detected in PICRUSt2 KO abundances.")
}

if (any(ko_abundance_matrix < 0)) {
  stop("Negative values detected in PICRUSt2 KO abundances.")
}

cat(
  "\nPICRUSt2 KO input:\n",
  "  KOs:     ", nrow(ko_for_kegg), "\n",
  "  Samples: ", length(ko_sample_ids), "\n",
  sep = ""
)


# ==============================================================================
# 4. CONVERT KO ABUNDANCES TO KEGG PATHWAY ABUNDANCES
# ==============================================================================

# ggpicrust2 aggregates predicted KOs using its bundled KO-to-KEGG mapping.
# filter_for_prokaryotes = TRUE is the package default and is retained here.
kegg_pathways <- ggpicrust2::ko2kegg_abundance(
  data = ko_for_kegg,
  filter_for_prokaryotes = TRUE
)

kegg_pathways <- as.matrix(kegg_pathways)
storage.mode(kegg_pathways) <- "numeric"

if (!setequal(colnames(kegg_pathways), ko_sample_ids)) {
  stop("Sample membership changed during KO-to-KEGG conversion.")
}

if (anyNA(kegg_pathways) || !all(is.finite(kegg_pathways))) {
  stop("Invalid abundance values detected after KO-to-KEGG conversion.")
}

if (any(kegg_pathways < 0)) {
  stop("Negative pathway abundances detected after KO-to-KEGG conversion.")
}

cat(
  "\nKEGG pathway table after KO aggregation:\n",
  "  Pathways: ", nrow(kegg_pathways), "\n",
  "  Samples:  ", ncol(kegg_pathways), "\n",
  sep = ""
)


# ==============================================================================
# 5. ANNOTATE KEGG PATHWAYS
# ==============================================================================

# ko_to_kegg_reference has one row per KO-pathway association. Collapse it to
# one annotation row per pathway.
data(
  "ko_to_kegg_reference",
  package = "ggpicrust2"
)

kegg_annotation <- ko_to_kegg_reference |>
  dplyr::select(
    pathway_id,
    pathway_number,
    pathway_name,
    level1,
    level2,
    level3
  ) |>
  dplyr::distinct() |>
  dplyr::filter(
    pathway_id %in% rownames(kegg_pathways)
  )

if (anyDuplicated(kegg_annotation$pathway_id)) {
  stop("KEGG annotation reference produced >1 annotation row per pathway.")
}

if (!setequal(rownames(kegg_pathways), kegg_annotation$pathway_id)) {
  stop("At least one predicted KEGG pathway lacks annotation.")
}


# ==============================================================================
# 6. MATCH KEGG SAMPLES TO PHYLOSEQ METADATA
# ==============================================================================

kegg_sample_ids <- colnames(kegg_pathways)

if (!setequal(kegg_sample_ids, metadata_kegg$SampleID)) {
  stop("KEGG sample IDs do not exactly match the processed phyloseq samples.")
}

metadata_kegg <- metadata_kegg[
  match(kegg_sample_ids, metadata_kegg$SampleID),
  ,
  drop = FALSE
]

rownames(metadata_kegg) <- metadata_kegg$SampleID

stopifnot(
  identical(
    colnames(kegg_pathways),
    metadata_kegg$SampleID
  )
)


# ==============================================================================
# 7. PATHWAY PREVALENCE FILTERING AND RELATIVE ABUNDANCE
# ==============================================================================

kegg_prevalence <- tibble::tibble(
  pathway_id = rownames(kegg_pathways),
  n_present = rowSums(kegg_pathways > 0),
  prevalence = rowMeans(kegg_pathways > 0),
  total_predicted_abundance = rowSums(kegg_pathways)
) |>
  dplyr::left_join(
    kegg_annotation,
    by = "pathway_id"
  )

keep_pathways <- kegg_prevalence |>
  dplyr::filter(
    prevalence >= minimum_pathway_prevalence
  ) |>
  dplyr::pull(pathway_id)

kegg_filtered <- kegg_pathways[
  rownames(kegg_pathways) %in% keep_pathways,
  ,
  drop = FALSE
]

if (any(colSums(kegg_filtered) <= 0)) {
  stop("At least one sample has zero total abundance after pathway filtering.")
}

kegg_relative <- sweep(
  kegg_filtered,
  MARGIN = 2,
  STATS = colSums(kegg_filtered),
  FUN = "/"
)

if (!all(abs(colSums(kegg_relative) - 1) < 1e-8)) {
  stop("KEGG relative abundances do not sum to 1 within each sample.")
}

cat(
  "\nKEGG pathway prevalence filtering:\n",
  "  Before: ", nrow(kegg_pathways), " pathways\n",
  "  After:  ", nrow(kegg_filtered), " pathways (>=",
  100 * minimum_pathway_prevalence, "% prevalence)\n",
  sep = ""
)

# Save the useful annotation/QC table, including whether each pathway entered
# downstream analyses.
kegg_pathway_qc <- kegg_prevalence |>
  dplyr::mutate(
    retained_for_analysis = pathway_id %in% keep_pathways
  ) |>
  dplyr::arrange(
    dplyr::desc(retained_for_analysis),
    dplyr::desc(prevalence),
    pathway_id
  )



# ==============================================================================
# 8. SELECT ONE INDEPENDENT SAMPLE PER MOUSE FOR BETA DIVERSITY
# ==============================================================================

# Colony mice: sample closest to Timept = 0.
metadata_colony <- metadata_kegg |>
  dplyr::filter(
    SampleOrigin != "Wild:As"
  ) |>
  dplyr::mutate(
    distance_from_T0 = abs(Timept)
  ) |>
  dplyr::group_by(mouse_ID_dataset) |>
  dplyr::arrange(
    distance_from_T0,
    SampleDate,
    SampleID,
    .by_group = TRUE
  ) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

# Wild mice: earliest available sample.
metadata_wild <- metadata_kegg |>
  dplyr::filter(
    SampleOrigin == "Wild:As"
  ) |>
  dplyr::mutate(
    SampleDate_order = as.Date(SampleDate)
  ) |>
  dplyr::group_by(mouse_ID_dataset) |>
  dplyr::arrange(
    SampleDate_order,
    SampleID,
    .by_group = TRUE
  ) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

metadata_kegg_single <- dplyr::bind_rows(
  metadata_wild,
  metadata_colony
) |>
  dplyr::mutate(
    SampleOrigin = factor(
      SampleOrigin,
      levels = group_levels
    )
  )

if (anyDuplicated(metadata_kegg_single$mouse_ID_dataset)) {
  stop("Independent beta-diversity data still contain repeated mice.")
}

single_sample_ids <- metadata_kegg_single$SampleID

kegg_beta_matrix <- t(
  kegg_relative[
    ,
    single_sample_ids,
    drop = FALSE
  ]
)

metadata_kegg_single <- metadata_kegg_single[
  match(rownames(kegg_beta_matrix), metadata_kegg_single$SampleID),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(kegg_beta_matrix),
    metadata_kegg_single$SampleID
  )
)

independent_sample_counts <- metadata_kegg_single |>
  dplyr::count(
    SampleOrigin,
    name = "n_mice"
  )

cat("\nIndependent mice used for KEGG beta diversity:\n")
print(independent_sample_counts)


# ==============================================================================
# 9. KEGG FUNCTIONAL BETA DIVERSITY
# ==============================================================================

kegg_bray <- vegan::vegdist(
  kegg_beta_matrix,
  method = "bray"
)

# ------------------------------------------------------------------------------
# 8A. PCoA
# ------------------------------------------------------------------------------

kegg_pcoa <- stats::cmdscale(
  kegg_bray,
  eig = TRUE,
  k = 2
)

positive_eigenvalues <- kegg_pcoa$eig[
  kegg_pcoa$eig > 0
]

pcoa1_percent <- 100 *
  kegg_pcoa$eig[1] /
  sum(positive_eigenvalues)

pcoa2_percent <- 100 *
  kegg_pcoa$eig[2] /
  sum(positive_eigenvalues)

kegg_pcoa_df <- tibble::tibble(
  SampleID = rownames(kegg_pcoa$points),
  PCoA1 = kegg_pcoa$points[, 1],
  PCoA2 = kegg_pcoa$points[, 2]
) |>
  dplyr::left_join(
    metadata_kegg_single |>
      dplyr::select(
        SampleID,
        SampleOrigin,
        mouse_ID_dataset
      ),
    by = "SampleID"
  )

# ------------------------------------------------------------------------------
# 8B. Omnibus PERMANOVA
# ------------------------------------------------------------------------------

set.seed(random_seed)

kegg_permanova <- vegan::adonis2(
  kegg_bray ~ SampleOrigin,
  data = metadata_kegg_single,
  permutations = n_permutations
)

# ------------------------------------------------------------------------------
# 8C. PERMDISP
# ------------------------------------------------------------------------------

kegg_dispersion <- vegan::betadisper(
  kegg_bray,
  group = metadata_kegg_single$SampleOrigin
)

set.seed(random_seed)

kegg_permdisp <- vegan::permutest(
  kegg_dispersion,
  permutations = n_permutations
)

kegg_dispersion_summary <- tibble::tibble(
  SampleOrigin = metadata_kegg_single$SampleOrigin,
  distance_to_centroid = kegg_dispersion$distances
) |>
  dplyr::group_by(SampleOrigin) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_distance = mean(distance_to_centroid),
    sd_distance = stats::sd(distance_to_centroid),
    median_distance = stats::median(distance_to_centroid),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 8D. Pairwise PERMANOVAs
# ------------------------------------------------------------------------------

pairwise_comparisons <- list(
  "Wild-like vs Wild" = c("Wild-like:As", "Wild:As"),
  "Lab-like vs Wild" = c("Lab-like:As", "Wild:As"),
  "Lab vs Wild" = c("Lab:Mm", "Wild:As"),
  "Lab-like vs Wild-like" = c("Lab-like:As", "Wild-like:As"),
  "Lab vs Wild-like" = c("Lab:Mm", "Wild-like:As"),
  "Lab vs Lab-like" = c("Lab:Mm", "Lab-like:As")
)

kegg_pairwise_permanova <- lapply(
  names(pairwise_comparisons),
  function(comparison_name) {

    groups <- pairwise_comparisons[[comparison_name]]

    keep <- metadata_kegg_single$SampleOrigin %in% groups

    pair_matrix <- kegg_beta_matrix[
      keep,
      ,
      drop = FALSE
    ]

    pair_metadata <- metadata_kegg_single[
      keep,
      ,
      drop = FALSE
    ]

    pair_metadata$SampleOrigin <- droplevels(
      factor(pair_metadata$SampleOrigin)
    )

    pair_bray <- vegan::vegdist(
      pair_matrix,
      method = "bray"
    )

    set.seed(random_seed)

    pair_fit <- vegan::adonis2(
      pair_bray ~ SampleOrigin,
      data = pair_metadata,
      permutations = n_permutations
    )

    tibble::tibble(
      test = comparison_name,
      first_group = groups[1],
      second_group = groups[2],
      n_first = sum(pair_metadata$SampleOrigin == groups[1]),
      n_second = sum(pair_metadata$SampleOrigin == groups[2]),
      F = pair_fit$F[1],
      R2 = pair_fit$R2[1],
      p_value = pair_fit$`Pr(>F)`[1]
    )
  }
) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    p_adjusted = stats::p.adjust(
      p_value,
      method = "BH"
    ),
    significant = p_adjusted < 0.05,
    test = factor(test, levels = contrast_names)
  ) |>
  dplyr::arrange(test)

cat("\nOmnibus KEGG PERMANOVA:\n")
print(kegg_permanova)

cat("\nKEGG PERMDISP:\n")
print(kegg_permdisp)

cat("\nPairwise KEGG PERMANOVAs:\n")
print(kegg_pairwise_permanova)


# ==============================================================================
# 10. SEQUENCING-DEPTH SENSITIVITY ANALYSIS
# ==============================================================================

# This is a diagnostic because PICRUSt2 pathway abundance derives from 16S data.
# SampleOrigin remains the biological predictor of interest.

set.seed(random_seed)

kegg_permanova_depth <- vegan::adonis2(
  kegg_bray ~ log10_reads,
  data = metadata_kegg_single,
  permutations = n_permutations
)

set.seed(random_seed)

kegg_permanova_adjusted <- vegan::adonis2(
  kegg_bray ~ SampleOrigin + log10_reads,
  data = metadata_kegg_single,
  permutations = n_permutations,
  by = "margin"
)

kegg_depth_summary <- metadata_kegg_single |>
  dplyr::group_by(SampleOrigin) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_reads = mean(reads),
    sd_reads = stats::sd(reads),
    mean_log10_reads = mean(log10_reads),
    sd_log10_reads = stats::sd(log10_reads),
    .groups = "drop"
  )

cat("\nSequencing-depth-only PERMANOVA:\n")
print(kegg_permanova_depth)

cat("\nMarginal PERMANOVA adjusted for sequencing depth:\n")
print(kegg_permanova_adjusted)



# ==============================================================================
# 10B. PATHWAY-PREVALENCE SENSITIVITY ANALYSIS — PERMANOVA ONLY
# ==============================================================================
#
# The primary analysis excludes pathways predicted in <10% of samples.
# This sensitivity analysis asks whether the overall SampleOrigin effect on
# predicted functional composition is robust to that filtering choice.
#
# The same independent one-sample-per-mouse dataset is used at every threshold.
# At each threshold:
#   1. pathways are filtered from the UNFILTERED KEGG pathway matrix,
#   2. pathway abundances are renormalized to relative abundance,
#   3. Bray-Curtis dissimilarities are recalculated,
#   4. the omnibus PERMANOVA SampleOrigin effect is rerun.
#
# MaAsLin3 is deliberately NOT rerun here. The purpose is only to establish
# robustness of the overall functional-composition result.
# ==============================================================================

kegg_prevalence_sensitivity <- lapply(
  seq_along(pathway_prevalence_sweep),
  function(i) {

    prevalence_cutoff <- pathway_prevalence_sweep[i]

    keep_ids <- rownames(kegg_pathways)[
      rowMeans(kegg_pathways > 0) >= prevalence_cutoff
    ]

    kegg_threshold <- kegg_pathways[
      keep_ids,
      ,
      drop = FALSE
    ]

    if (nrow(kegg_threshold) == 0) {
      stop(
        "No KEGG pathways retained at prevalence cutoff ",
        prevalence_cutoff
      )
    }

    if (any(colSums(kegg_threshold) <= 0)) {
      stop(
        "At least one sample has zero total predicted pathway abundance at ",
        prevalence_cutoff * 100,
        "% prevalence."
      )
    }

    # Renormalize after filtering.
    kegg_threshold_relative <- sweep(
      kegg_threshold,
      MARGIN = 2,
      STATS = colSums(kegg_threshold),
      FUN = "/"
    )

    # Use exactly the same independent mice as the primary analysis.
    threshold_beta_matrix <- t(
      kegg_threshold_relative[
        ,
        single_sample_ids,
        drop = FALSE
      ]
    )

    stopifnot(
      identical(
        rownames(threshold_beta_matrix),
        metadata_kegg_single$SampleID
      )
    )

    threshold_bray <- vegan::vegdist(
      threshold_beta_matrix,
      method = "bray"
    )

    # Use the same random seed at every cutoff so differences reflect the
    # filtering scheme rather than different permutation draws.
    set.seed(random_seed)

    threshold_permanova <- vegan::adonis2(
      threshold_bray ~ SampleOrigin,
      data = metadata_kegg_single,
      permutations = n_permutations
    )

    tibble::tibble(
      prevalence_cutoff = prevalence_cutoff,
      prevalence_percent = prevalence_cutoff * 100,
      pathways_retained = nrow(kegg_threshold),
      F = threshold_permanova$F[1],
      R2 = threshold_permanova$R2[1],
      p_value = threshold_permanova$`Pr(>F)`[1],
      primary_threshold =
        prevalence_cutoff == minimum_pathway_prevalence
    )
  }
) |>
  dplyr::bind_rows()

cat(
  "\nKEGG pathway-prevalence sensitivity — omnibus PERMANOVA:\n"
)

print(kegg_prevalence_sensitivity)



# ------------------------------------------------------------------------------
# Supplementary pathway-prevalence sensitivity figure
# ------------------------------------------------------------------------------

p_kegg_prevalence_A <- ggplot2::ggplot(
  kegg_prevalence_sensitivity,
  ggplot2::aes(
    x = prevalence_percent,
    y = pathways_retained
  )
) +
  ggplot2::geom_line(
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    size = 2.2
  ) +
  ggplot2::geom_vline(
    xintercept = 100 * minimum_pathway_prevalence,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::labs(
    x = "Global pathway prevalence cutoff (%)",
    y = "KEGG pathways retained"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  )

p_kegg_prevalence_B <- ggplot2::ggplot(
  kegg_prevalence_sensitivity,
  ggplot2::aes(
    x = prevalence_percent,
    y = R2
  )
) +
  ggplot2::geom_line(
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    size = 2.2
  ) +
  ggplot2::geom_vline(
    xintercept = 100 * minimum_pathway_prevalence,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::labs(
    x = "Global pathway prevalence cutoff (%)",
    y = expression(PERMANOVA~R^2)
  ) +
  ggplot2::theme_classic(
    base_size = 11
  )

p_kegg_prevalence_sensitivity <- patchwork::wrap_plots(
  p_kegg_prevalence_A,
  p_kegg_prevalence_B,
  nrow = 1
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "PICRUSt2_Pathway_Filtering_Sensitivity.pdf"
  ),
  plot = p_kegg_prevalence_sensitivity,
  width = 8.5,
  height = 4.2
)


# ==============================================================================
# 11. PREPARE BETA-DIVERSITY RESULT TABLES
# ==============================================================================

# Compact summary table of the principal omnibus tests.
kegg_beta_summary <- tibble::tibble(
  analysis = c(
    "SampleOrigin omnibus PERMANOVA",
    "PERMDISP",
    "Sequencing depth only",
    "SampleOrigin adjusted for sequencing depth",
    "Sequencing depth adjusted for SampleOrigin"
  ),
  statistic = c(
    kegg_permanova$F[1],
    unname(kegg_permdisp$tab[1, "F"]),
    kegg_permanova_depth$F[1],
    kegg_permanova_adjusted["SampleOrigin", "F"],
    kegg_permanova_adjusted["log10_reads", "F"]
  ),
  R2 = c(
    kegg_permanova$R2[1],
    NA_real_,
    kegg_permanova_depth$R2[1],
    kegg_permanova_adjusted["SampleOrigin", "R2"],
    kegg_permanova_adjusted["log10_reads", "R2"]
  ),
  p_value = c(
    kegg_permanova$`Pr(>F)`[1],
    unname(kegg_permdisp$tab[1, "Pr(>F)"]),
    kegg_permanova_depth$`Pr(>F)`[1],
    kegg_permanova_adjusted["SampleOrigin", "Pr(>F)"],
    kegg_permanova_adjusted["log10_reads", "Pr(>F)"]
  )
)






# ==============================================================================
# 12. PREPARE DATA FOR DIFFERENTIAL KEGG PATHWAY ABUNDANCE
# ==============================================================================

# MaAsLin3 uses all eligible samples. Repeated samples are handled using a mouse
# random intercept.

kegg_model_metadata <- metadata_kegg |>
  dplyr::mutate(
    SampleOrigin = factor(SampleOrigin, levels = group_levels),
    Sex = factor(Sex),
    mouse_ID_dataset = factor(mouse_ID_dataset)
  ) |>
  dplyr::filter(
    stats::complete.cases(
      SampleOrigin,
      Sex,
      log10_reads,
      mouse_ID_dataset
    )
  )

kegg_model_abundance <- t(
  kegg_relative[
    ,
    kegg_model_metadata$SampleID,
    drop = FALSE
  ]
)

kegg_model_metadata <- kegg_model_metadata[
  match(rownames(kegg_model_abundance), kegg_model_metadata$SampleID),
  ,
  drop = FALSE
]

rownames(kegg_model_metadata) <- kegg_model_metadata$SampleID

stopifnot(
  identical(
    rownames(kegg_model_abundance),
    rownames(kegg_model_metadata)
  )
)

cat("\nSamples entering KEGG MaAsLin3:\n")
print(
  kegg_model_metadata |>
    dplyr::count(SampleOrigin, name = "n_samples")
)


# ==============================================================================
# 13. DIFFERENTIAL KEGG PATHWAY ABUNDANCE WITH MaAsLin3
# ==============================================================================

if (rerun_maaslin3 && dir.exists(maaslin_dir)) {
  unlink(
    maaslin_dir,
    recursive = TRUE,
    force = TRUE
  )
}

dir.create(
  maaslin_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!rerun_maaslin3 && file.exists(file.path(maaslin_dir, "all_results.tsv"))) {
  stop(
    "Existing KEGG MaAsLin3 results detected. This script needs the live fit ",
    "object to run pairwise contrasts. Set rerun_maaslin3 <- TRUE."
  )
}

set.seed(random_seed)

kegg_fit <- maaslin3::maaslin3(
  input_data = kegg_model_abundance,
  input_metadata = kegg_model_metadata,
  output = maaslin_dir,
  formula = ~ SampleOrigin +
    Sex +
    log10_reads +
    (1 | mouse_ID_dataset),
  normalization = "NONE",
  transform = "LOG",
  min_abundance = 0,
  min_prevalence = 0,
  evaluate_only = "abundance",
  warn_prevalence = FALSE,
  median_comparison_abundance = TRUE,
  standardize = TRUE,
  max_significance = q_threshold,
  plot_summary_plot = FALSE,
  plot_associations = FALSE,
  max_pngs = 0,
  save_models = TRUE,
  save_plots_rds = FALSE
)


# ==============================================================================
# 14. DEFINE AND TEST ALL SIX SAMPLE-ORIGIN CONTRASTS
# ==============================================================================

# Wild:As is the model reference level. Positive coefficients from each contrast
# mean higher predicted pathway abundance in the FIRST-NAMED group.

required_origin_coefficients <- c(
  "SampleOriginWild-like:As",
  "SampleOriginLab-like:As",
  "SampleOriginLab:Mm"
)

available_coefficient_names <- unique(
  kegg_fit$fit_data_abundance$results$name
)

if (!all(required_origin_coefficients %in% available_coefficient_names)) {
  stop(
    "Could not find the expected SampleOrigin coefficients in the MaAsLin3 fit."
  )
}

kegg_contrast_matrix <- rbind(
  "Wild-like vs Wild" = c(1, 0, 0),
  "Lab-like vs Wild" = c(0, 1, 0),
  "Lab vs Wild" = c(0, 0, 1),
  "Lab-like vs Wild-like" = c(-1, 1, 0),
  "Lab vs Wild-like" = c(-1, 0, 1),
  "Lab vs Lab-like" = c(0, -1, 1)
)

colnames(kegg_contrast_matrix) <- required_origin_coefficients

kegg_contrasts <- maaslin3::maaslin_contrast_test(
  maaslin3_fit = kegg_fit,
  contrast_mat = kegg_contrast_matrix,
  median_comparison_abundance = TRUE,
  evaluate_only = "abundance",
  correction = "BH"
)

kegg_contrast_results <- tibble::as_tibble(
  kegg_contrasts$fit_data_abundance$results
) |>
  dplyr::left_join(
    kegg_annotation |>
      dplyr::select(
        pathway_id,
        pathway_name,
        level1,
        level2,
        level3
      ),
    by = c("feature" = "pathway_id")
  ) |>
  dplyr::left_join(
    contrast_key,
    by = "test"
  ) |>
  dplyr::mutate(
    test = factor(test, levels = contrast_names),
    significant = !is.na(qval_individual) &
      qval_individual < q_threshold,
    direction = dplyr::case_when(
      coef > 0 ~ paste0("Higher in ", first_group),
      coef < 0 ~ paste0("Higher in ", second_group),
      TRUE ~ "No difference"
    ),
    model_prevalence = N_not_zero / N
  ) |>
  dplyr::arrange(
    test,
    qval_individual,
    dplyr::desc(abs(coef))
  )

if (anyNA(kegg_contrast_results$pathway_name)) {
  stop("At least one MaAsLin3 pathway result failed KEGG annotation.")
}


# ==============================================================================
# 15. SUMMARIZE DIFFERENTIAL KEGG PATHWAYS
# ==============================================================================

kegg_significant_counts <- kegg_contrast_results |>
  dplyr::group_by(test) |>
  dplyr::summarise(
    n_tested = dplyr::n(),
    n_significant = sum(significant),
    percent_significant = 100 * n_significant / n_tested,
    n_higher_first = sum(significant & coef > 0),
    n_higher_second = sum(significant & coef < 0),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    test = factor(test, levels = contrast_names)
  ) |>
  dplyr::arrange(test)

kegg_unique_significant <- kegg_contrast_results |>
  dplyr::filter(significant) |>
  dplyr::distinct(
    feature,
    pathway_name,
    level1,
    level2,
    level3
  )

# Descriptive KEGG Level 2 table. This is useful for supplementary inspection,
# but is not an enrichment analysis and should not be interpreted as one.
kegg_level2_summary <- kegg_contrast_results |>
  dplyr::group_by(
    test,
    level1,
    level2
  ) |>
  dplyr::summarise(
    n_tested = dplyr::n(),
    n_significant = sum(significant),
    percent_significant = 100 * n_significant / n_tested,
    n_higher_first = sum(significant & coef > 0),
    n_higher_second = sum(significant & coef < 0),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    test = factor(test, levels = contrast_names)
  ) |>
  dplyr::arrange(
    test,
    dplyr::desc(percent_significant),
    level2
  )

cat("\nSignificant KEGG pathways by comparison:\n")
print(kegg_significant_counts)

cat(
  "\nUnique KEGG pathways significant in >=1 contrast: ",
  nrow(kegg_unique_significant),
  "\n",
  sep = ""
)


# ==============================================================================
# 16. PREPARE DIFFERENTIAL-ABUNDANCE RESULT TABLES
# ==============================================================================

# Complete annotated pairwise results: the main reproducibility table.

# Significant subset for convenient manuscript/supplementary use.




# ==============================================================================


# 17. PRINT MANUSCRIPT-RELEVANT RESULTS TO CONSOLE
# ==============================================================================

cat("\n============================================================\n")
cat("PICRUSt2 KEGG + NSTI ANALYSIS SUMMARY\n")
cat("============================================================\n")

cat("\nWeighted NSTI by SampleOrigin:\n")
print(nsti_summary)

cat("\nNSTI mixed-model ANOVA:\n")
print(nsti_anova)

cat(
  "Pathways retained: ",
  nrow(kegg_filtered),
  "/",
  nrow(kegg_pathways),
  " (>=",
  100 * minimum_pathway_prevalence,
  "% prevalence)\n",
  sep = ""
)

cat(
  "Independent mice in beta diversity: ",
  nrow(metadata_kegg_single),
  "\n",
  sep = ""
)

cat(
  "PCoA1 + PCoA2: ",
  round(pcoa1_percent + pcoa2_percent, 1),
  "% of positive-eigenvalue variation\n",
  sep = ""
)

cat("\nOmnibus PERMANOVA:\n")
print(kegg_permanova)

cat("\nPERMDISP:\n")
print(kegg_permdisp)

cat("\nPairwise PERMANOVA results:\n")
print(kegg_pairwise_permanova)

cat("\nDifferential pathway counts:\n")
print(kegg_significant_counts)

cat(
  "\nUnique significant pathways across all contrasts: ",
  nrow(kegg_unique_significant),
  "\n",
  sep = ""
)

cat("\nSaved KEGG tables to:\n  ", tables_dir, "\n", sep = "")
cat("Saved KEGG figures to:\n  ", figures_dir, "\n", sep = "")
cat("Saved pathway-prevalence PERMANOVA sensitivity table and figure.\n")
cat("Saved MaAsLin3 model output to:\n  ", maaslin_dir, "\n", sep = "")


# ==============================================================================
# 18. FIGURE SETUP
# ==============================================================================
#
# All figure construction is kept at the end of the script so the analysis,
# model fitting, and saved statistical tables are completed first.
#
# Sample-origin colours match the palette used in the baseline alpha- and
# beta-diversity scripts:
#   Wild:As      = green
#   Wild-like:As = blue
#   Lab-like:As  = purple
#   Lab:Mm       = coral
# ==============================================================================

group_colours <- c(
  "Wild:As" = "#74C476",
  "Wild-like:As" = "#6BAED6",
  "Lab-like:As" = "#BC80BD",
  "Lab:Mm" = "#FB8072"
)

# Directional colours are not cohort colours. They indicate the direction of a
# significant pathway difference along each displayed comparison.
direction_colours <- c(
  "Higher toward lab-associated" = "#F8766D",
  "Higher toward wild-associated" = "#00BFC4"
)


# ==============================================================================
# 19. FIGURE A — KEGG FUNCTIONAL COMPOSITION PCoA
# ==============================================================================

p_kegg_pcoa <- ggplot2::ggplot(
  kegg_pcoa_df,
  ggplot2::aes(
    x = PCoA1,
    y = PCoA2,
    colour = SampleOrigin
  )
) +
  ggplot2::geom_point(
    size = 2.4,
    alpha = 0.75
  ) +
  ggplot2::stat_ellipse(
    ggplot2::aes(group = SampleOrigin),
    type = "t",
    level = 0.95,
    linewidth = 0.8
  ) +
  ggplot2::scale_colour_manual(
    values = group_colours,
    breaks = group_levels,
    labels = cohort_labels
  ) +
  ggplot2::labs(
    x = paste0("PCoA1 (", round(pcoa1_percent, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa2_percent, 1), "%)"),
    colour = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    legend.position = "top",
    axis.text = ggplot2::element_text(colour = "black")
  )



# ==============================================================================
# 20. SUPPLEMENTARY FIGURE — PAIRWISE FUNCTIONAL EFFECT SIZES
# ==============================================================================
#
# Pairwise PERMANOVA R2 values are retained as a useful standalone summary,
# although this panel is not included in the final three-panel functional figure.
# ==============================================================================

pairwise_effect_plot_data <- kegg_pairwise_permanova |>
  dplyr::mutate(
    source_test = as.character(test)
  ) |>
  dplyr::left_join(
    comparison_display_key,
    by = "source_test"
  ) |>
  dplyr::mutate(
    # Reverse factor levels because this comparison is plotted on the y-axis:
    # the first level should appear at the TOP of the horizontal plot.
    test_label = factor(
      display_test,
      levels = rev(directional_contrast_order)
    ),
    p_label = paste0(
      "BH p ",
      format_p(p_adjusted)
    )
  )

p_kegg_pairwise_effects <- ggplot2::ggplot(
  pairwise_effect_plot_data,
  ggplot2::aes(
    x = R2,
    y = test_label
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = R2,
      y = test_label,
      yend = test_label
    ),
    linewidth = 0.7
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = R2,
      label = p_label
    ),
    hjust = -0.12,
    size = 3.1
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, max(pairwise_effect_plot_data$R2) * 1.35),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::scale_y_discrete(
    labels = comparison_display_labels
  ) +
  ggplot2::labs(
    x = expression("Pairwise PERMANOVA " * R^2),
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(colour = "black")
  )



# ==============================================================================
# 21. FIGURE B — NUMBER AND PROPORTION OF DIFFERENTIAL KEGG PATHWAYS
# ==============================================================================
#
# Bars show the number of FDR-significant pathways in each pairwise comparison.
# Labels report significant / total tested and the corresponding percentage.
# ==============================================================================

pathway_count_plot_data <- kegg_significant_counts |>
  dplyr::mutate(
    source_test = as.character(test)
  ) |>
  dplyr::left_join(
    comparison_display_key,
    by = "source_test"
  ) |>
  dplyr::mutate(
    # Reverse factor levels because comparisons are plotted on the y-axis.
    test_label = factor(
      display_test,
      levels = rev(directional_contrast_order)
    ),
    count_label = paste0(
      n_significant,
      "/",
      n_tested,
      " (",
      round(percent_significant, 1),
      "%)"
    )
  )

p_kegg_pathway_counts <- ggplot2::ggplot(
  pathway_count_plot_data,
  ggplot2::aes(
    x = n_significant,
    y = test_label
  )
) +
  ggplot2::geom_col(
    width = 0.65,
    fill = "grey45"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = count_label),
    hjust = -0.10,
    size = 3.1
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      max(pathway_count_plot_data$n_significant) * 1.32
    ),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::scale_y_discrete(
    labels = comparison_display_labels
  ) +
  ggplot2::labs(
    x = "FDR-significant predicted KEGG pathways",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(colour = "black")
  )



# ==============================================================================
# 22. FIGURE C — COMPLETE KEGG LEVEL-2 DIRECTIONAL SUMMARY
# ==============================================================================
#
# This figure includes ALL KEGG Level 2 categories represented among the tested
# pathways. No category is excluded according to significance, category size, or
# biological interpretation.
#
# For each Level 2 category and comparison, bars show the percentage of all
# tested pathways in that category that were individually significant at
# BH-FDR q < 0.05. Nonsignificant pathways remain in the denominator.
#
# Comparisons are reoriented for display from more wild-associated to more
# lab-associated. Teal bars extend left for pathways higher toward the
# wild-associated side; coral bars extend right for pathways higher toward the
# lab-associated side.
#
# This is a DESCRIPTIVE organization of pathway-level results, not a formal
# Level 2 enrichment test.
# ==============================================================================

# -----------------------------------------------------------------------------
# 21A. Reuse the shared wild-to-lab comparison display definitions
# -----------------------------------------------------------------------------
# directional_contrast_key and directional_contrast_order were defined near the
# analysis settings so they are available to every pairwise figure.



# -----------------------------------------------------------------------------
# 21B. Reorient pathway-level coefficients and clean KEGG hierarchy labels
# -----------------------------------------------------------------------------

kegg_directional_results <- kegg_contrast_results |>
  dplyr::inner_join(
    directional_contrast_key,
    by = c("test" = "source_test")
  ) |>
  dplyr::mutate(
    coef_display = coef * multiplier,
    significant = !is.na(qval_individual) &
      qval_individual < q_threshold,
    direction = dplyr::case_when(
      significant & coef_display > 0 ~ "Higher toward wild-associated",
      significant & coef_display < 0 ~ "Higher toward lab-associated",
      TRUE ~ "Not significant"
    ),
    level1_clean = sub(
      "^[0-9]+\\s+",
      "",
      level1
    ),
    level2_clean = sub(
      "^[0-9]+\\s+",
      "",
      level2
    )
  )

# -----------------------------------------------------------------------------
# 21C. Summarize direction within every Level 2 category
# -----------------------------------------------------------------------------

kegg_level2_direction_summary <- kegg_directional_results |>
  dplyr::group_by(
    display_test,
    level1_clean,
    level2_clean
  ) |>
  dplyr::summarise(
    n_tested = dplyr::n(),
    n_wild = sum(
      significant & coef_display > 0,
      na.rm = TRUE
    ),
    n_lab = sum(
      significant & coef_display < 0,
      na.rm = TRUE
    ),
    n_nonsignificant = sum(
      !significant,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    percent_wild = 100 * n_wild / n_tested,
    percent_lab = 100 * n_lab / n_tested,
    percent_nonsignificant = 100 * n_nonsignificant / n_tested
  )

# -----------------------------------------------------------------------------
# 21D. Add pathway denominator to each Level 2 category label
# -----------------------------------------------------------------------------

level2_label_key <- kegg_level2_direction_summary |>
  dplyr::distinct(
    level1_clean,
    level2_clean,
    n_tested
  ) |>
  dplyr::mutate(
    level2_label = paste0(
      level2_clean,
      " (n = ",
      n_tested,
      ")"
    )
  )

kegg_level2_direction_summary <- kegg_level2_direction_summary |>
  dplyr::left_join(
    level2_label_key,
    by = c(
      "level1_clean",
      "level2_clean",
      "n_tested"
    )
  )

# -----------------------------------------------------------------------------
# 21E. Preserve KEGG Level 1 hierarchy and objective within-group ordering
# -----------------------------------------------------------------------------

level1_order <- c(
  "Metabolism",
  "Genetic Information Processing",
  "Environmental Information Processing",
  "Cellular Processes",
  "Human Diseases"
)

kegg_level2_direction_summary <- kegg_level2_direction_summary |>
  dplyr::mutate(
    level1_clean = factor(
      level1_clean,
      levels = level1_order
    ),
    display_test = factor(
      display_test,
      levels = directional_contrast_order
    )
  )

level2_order_table <- level2_label_key |>
  dplyr::mutate(
    level1_clean = factor(
      level1_clean,
      levels = level1_order
    )
  ) |>
  dplyr::arrange(
    level1_clean,
    dplyr::desc(n_tested),
    level2_clean
  ) |>
  dplyr::mutate(
    level2_order = dplyr::row_number()
  )

kegg_level2_direction_summary <- kegg_level2_direction_summary |>
  dplyr::left_join(
    level2_order_table |>
      dplyr::select(
        level1_clean,
        level2_clean,
        level2_label,
        level2_order
      ),
    by = c(
      "level1_clean",
      "level2_clean",
      "level2_label"
    )
  )

# -----------------------------------------------------------------------------
# 21F. Convert to long format for diverging bars
# -----------------------------------------------------------------------------

kegg_level2_direction_long <- kegg_level2_direction_summary |>
  dplyr::select(
    display_test,
    level1_clean,
    level2_clean,
    level2_label,
    level2_order,
    n_tested,
    percent_wild,
    percent_lab
  ) |>
  tidyr::pivot_longer(
    cols = c(
      percent_wild,
      percent_lab
    ),
    names_to = "direction",
    values_to = "percent"
  ) |>
  dplyr::mutate(
    percent_plot = dplyr::if_else(
      direction == "percent_wild",
      -percent,
      percent
    ),
    direction = dplyr::recode(
      direction,
      "percent_wild" = "Higher toward wild-associated",
      "percent_lab" = "Higher toward lab-associated"
    ),
    display_test = factor(
      display_test,
      levels = directional_contrast_order
    ),
    level2_label = factor(
      level2_label,
      levels = rev(level2_order_table$level2_label)
    )
  )

# -----------------------------------------------------------------------------
# 21G. Display labels
# -----------------------------------------------------------------------------

contrast_display_labels <- c(
  "Wild vs Wild-like" = "Wild \u2013 Wild-like",
  "Wild vs Lab-like" = "Wild \u2013 Lab-like",
  "Wild vs Lab" = "Wild \u2013 Lab",
  "Wild-like vs Lab-like" = "Wild-like \u2013 Lab-like",
  "Wild-like vs Lab" = "Wild-like \u2013 Lab",
  "Lab-like vs Lab" = "Lab-like \u2013 Lab"
)

level1_display_labels <- c(
  "Metabolism" = "Metabolism",
  "Genetic Information Processing" = "Genetic information\nprocessing",
  "Environmental Information Processing" = "Environmental\ninformation processing",
  "Cellular Processes" = "Cellular processes",
  "Human Diseases" = "Human diseases"
)

# -----------------------------------------------------------------------------
# 21H. Plot all Level 2 categories
# -----------------------------------------------------------------------------

p_kegg_level2_all <- ggplot2::ggplot(
  kegg_level2_direction_long,
  ggplot2::aes(
    x = percent_plot,
    y = level2_label,
    fill = direction
  )
) +
  ggplot2::geom_col(
    width = 0.75
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.35
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(level1_clean),
    cols = ggplot2::vars(display_test),
    scales = "free_y",
    space = "free_y",
    switch = "y",
    labeller = ggplot2::labeller(
      display_test = contrast_display_labels,
      level1_clean = level1_display_labels
    )
  ) +
  ggplot2::scale_fill_manual(
    values = direction_colours
  ) +
  ggplot2::scale_x_continuous(
    labels = function(x) abs(x),
    limits = c(-100, 100),
    breaks = c(-100, -50, 0, 50, 100)
  ) +
  ggplot2::labs(
    x = "FDR-significant pathways within KEGG Level 2 category (%)",
    y = NULL,
    fill = NULL
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(
      size = 8,
      colour = "black"
    ),
    strip.text.y.left = ggplot2::element_text(
      size = 8,
      colour = "black",
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    strip.background = ggplot2::element_rect(
      fill = "white",
      colour = "black"
    ),
    strip.placement = "outside",
    axis.text.y = ggplot2::element_text(
      size = 6.5,
      colour = "black"
    ),
    axis.text.x = ggplot2::element_text(
      size = 7,
      colour = "black"
    ),
    axis.title.x = ggplot2::element_text(
      size = 9
    ),
    legend.position = "top",
    legend.text = ggplot2::element_text(
      size = 8
    ),
    panel.spacing.x = grid::unit(
      0.5,
      "lines"
    ),
    panel.spacing.y = grid::unit(
      0.7,
      "lines"
    ),
    plot.margin = ggplot2::margin(
      t = 10,
      r = 10,
      b = 10,
      l = 20
    )
  )

p_kegg_level2_all



# ==============================================================================
# FINAL MULTIPANEL FIGURE — PICRUSt2 / KEGG
# ==============================================================================
#
# All cohort/comparison displays follow the same wild-to-lab ordering:
#   Wild -> Wild-like -> Lab-like -> Lab
#
# Panel A uses that order in the cohort legend.
# Panels B-D use the same wilder-cohort-first ordering of pairwise comparisons.
# ==============================================================================

p_final_A <- p_kegg_pcoa
p_final_B <- p_kegg_pairwise_effects
p_final_C <- p_kegg_pathway_counts
p_final_D <- p_kegg_level2_all

# A
p_final_A <- p_final_A +
  ggplot2::theme(
    legend.text = ggplot2::element_text(size = 9),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    axis.title = ggplot2::element_text(size = 10)
  )

# B
p_final_B <- p_final_B +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    axis.title.x = ggplot2::element_text(size = 10)
  )

# C
p_final_C <- p_final_C +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    axis.title.x = ggplot2::element_text(size = 10)
  )

# D
p_final_D <- p_final_D +
  ggplot2::theme(
    legend.text = ggplot2::element_text(size = 9.5),
    axis.text.y = ggplot2::element_text(size = 8.5, colour = "black"),
    axis.text.x = ggplot2::element_text(size = 8.5, colour = "black"),
    axis.title.x = ggplot2::element_text(size = 10),
    strip.text.x = ggplot2::element_text(size = 9, colour = "black"),
    strip.text.y.left = ggplot2::element_text(
      size = 9,
      colour = "black",
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    )
  )

# ==============================================================================
# TOP ROW — EXACTLY 1 ROW × 3 COLUMNS
# ==============================================================================

functional_top_row <- cowplot::plot_grid(
  p_final_A,
  p_final_B,
  p_final_C,
  nrow = 1,
  ncol = 3,
  labels = c("A", "B", "C"),
  label_size = 12,
  label_fontface = "bold"
)

# ==============================================================================
# BOTTOM ROW — PANEL D
# ==============================================================================

functional_bottom_row <- cowplot::ggdraw() +
  cowplot::draw_plot(
    p_final_D,
    x = 0.03,
    y = 0,
    width = 0.94,
    height = 1
  ) +
  cowplot::draw_label(
    "D",
    x = 0.005,
    y = 0.995,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 12
  )

# ==============================================================================
# FINAL FIGURE
# ==============================================================================

functional_figure <- cowplot::plot_grid(
  functional_top_row,
  functional_bottom_row,
  ncol = 1,
  rel_heights = c(0.8, 2.0)
)

functional_figure

# Save
ggplot2::ggsave(
  file.path(figures_dir, "PICRUSt2_Functional_Summary.pdf"),
  functional_figure,
  width = 16,
  height = 10
)



# ==============================================================================
# SAVE ONE STREAMLINED PICRUSt2 RESULTS WORKBOOK
# ==============================================================================

# Compact analysis summary for manuscript reporting and reproducibility.
picrust2_analysis_summary <- tibble::tibble(
  Metric = c(
    "Samples with NSTI values",
    "Pathways before prevalence filter",
    "Pathways retained for primary analysis",
    "Primary pathway prevalence threshold",
    "Independent mice in beta diversity",
    "PCoA1 variation (%)",
    "PCoA2 variation (%)",
    "Samples in pathway differential-abundance model",
    "Unique pathways significant in >=1 contrast"
  ),
  Value = c(
    nrow(nsti_meta),
    nrow(kegg_pathways),
    nrow(kegg_filtered),
    minimum_pathway_prevalence,
    nrow(metadata_kegg_single),
    pcoa1_percent,
    pcoa2_percent,
    nrow(kegg_model_metadata),
    nrow(kegg_unique_significant)
  )
)

# Compact NSTI table: descriptive statistics plus the omnibus mixed-model test.
nsti_test_export <- nsti_anova_table |>
  dplyr::filter(
    term == "SampleOrigin"
  )

nsti_export <- nsti_summary |>
  dplyr::mutate(
    SampleOrigin = as.character(SampleOrigin)
  )

# Functional beta-diversity omnibus + sequencing-depth sensitivity.
functional_beta_export <- kegg_beta_summary |>
  dplyr::rename(
    Analysis = analysis,
    F = statistic,
    P = p_value
  )

# Pairwise PERMANOVA in the same wild-to-lab display order used in figures.
pairwise_permanova_export <- kegg_pairwise_permanova |>
  dplyr::mutate(
    source_test = as.character(test)
  ) |>
  dplyr::left_join(
    comparison_display_key,
    by = "source_test"
  ) |>
  dplyr::transmute(
    Comparison = display_test,
    n_first,
    n_second,
    F,
    R2,
    P = p_value,
    P_BH = p_adjusted,
    significant
  ) |>
  dplyr::arrange(
    factor(
      Comparison,
      levels = directional_contrast_order
    )
  )

# Keep the full prevalence-sensitivity sweep because it directly supports the
# dedicated supplementary robustness figure.
prevalence_sensitivity_export <- kegg_prevalence_sensitivity |>
  dplyr::select(
    prevalence_percent,
    pathways_retained,
    F,
    R2,
    p_value,
    primary_threshold
  ) |>
  dplyr::rename(
    Prevalence_percent = prevalence_percent,
    Pathways_retained = pathways_retained,
    P = p_value,
    Primary_threshold = primary_threshold
  )

# Counts of significant pathways by comparison.
pathway_counts_export <- kegg_significant_counts |>
  dplyr::mutate(
    source_test = as.character(test)
  ) |>
  dplyr::left_join(
    comparison_display_key,
    by = "source_test"
  ) |>
  dplyr::transmute(
    Comparison = display_test,
    n_tested,
    n_significant,
    percent_significant,
    n_higher_first,
    n_higher_second
  ) |>
  dplyr::arrange(
    factor(
      Comparison,
      levels = directional_contrast_order
    )
  )

# Full FDR-significant pathway-level results only. This is the principal
# Supplementary Data sheet and contains KEGG annotations plus effect estimates.
significant_pathways_export <- kegg_contrast_results |>
  dplyr::filter(
    significant
  ) |>
  dplyr::mutate(
    test = as.character(test)
  ) |>
  dplyr::left_join(
    comparison_display_key,
    by = c("test" = "source_test")
  ) |>
  dplyr::mutate(
    coef_display = -coef
  ) |>
  dplyr::transmute(
    Comparison = display_test,
    pathway_id = feature,
    pathway_name,
    level1,
    level2,
    level3,
    coefficient_original = coef,
    coefficient_wild_to_lab = coef_display,
    stderr,
    p_value = pval_individual,
    q_value = qval_individual,
    model_prevalence,
    N,
    N_not_zero
  ) |>
  dplyr::arrange(
    factor(
      Comparison,
      levels = directional_contrast_order
    ),
    q_value
  )

# Underlying descriptive Level-2 summary used in panel D. This is not an
# enrichment analysis; it simply summarizes significant pathway directions.
level2_export <- kegg_level2_direction_summary |>
  dplyr::mutate(
    display_test = as.character(display_test),
    level1_clean = as.character(level1_clean)
  ) |>
  dplyr::transmute(
    Comparison = display_test,
    KEGG_Level1 = level1_clean,
    KEGG_Level2 = level2_clean,
    n_tested,
    n_higher_wild = n_wild,
    n_higher_lab = n_lab,
    n_nonsignificant,
    percent_higher_wild = percent_wild,
    percent_higher_lab = percent_lab,
    percent_nonsignificant
  ) |>
  dplyr::arrange(
    factor(
      Comparison,
      levels = directional_contrast_order
    ),
    KEGG_Level1,
    KEGG_Level2
  )

picrust2_results_workbook <- list(
  "Analysis Summary" =
    picrust2_analysis_summary,

  "NSTI Summary" =
    nsti_export,

  "NSTI Test" =
    nsti_test_export,

  "Functional Beta" =
    functional_beta_export,

  "Pairwise PERMANOVA" =
    pairwise_permanova_export,

  "Filtering Sensitivity" =
    prevalence_sensitivity_export,

  "Pathway Counts" =
    pathway_counts_export,

  "Significant Pathways" =
    significant_pathways_export,

  "Level2 Summary" =
    level2_export
)

openxlsx::write.xlsx(
  x = picrust2_results_workbook,
  file = file.path(
    tables_dir,
    "PICRUSt2_Functional_Prediction_results.xlsx"
  ),
  overwrite = TRUE,
  rowNames = FALSE
)

# Remove temporary MaAsLin3 files now that all reader-facing outputs are saved.
if (dir.exists(maaslin_dir)) {
  unlink(
    maaslin_dir,
    recursive = TRUE,
    force = TRUE
  )
}

cat(
  "\n============================================================\n",
  "PICRUSt2 OUTPUTS\n",
  "============================================================\n",
  "Workbook:\n  ",
  file.path(
    tables_dir,
    "PICRUSt2_Functional_Prediction_results.xlsx"
  ),
  "\n\nComplete figures:\n  ",
  file.path(
    figures_dir,
    "PICRUSt2_Functional_Summary.pdf"
  ),
  "\n  ",
  file.path(
    figures_dir,
    "PICRUSt2_Pathway_Filtering_Sensitivity.pdf"
  ),
  "\n============================================================\n",
  sep = ""
)

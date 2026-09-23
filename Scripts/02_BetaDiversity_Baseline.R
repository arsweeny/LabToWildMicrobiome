###############################################################################
# 02. BETA DIVERSITY ANALYSIS — BASELINE SAMPLES
# ALL-IN-ONE PRIMARY + FILTERING-SENSITIVITY + FIGURE SCRIPT
#
# This script:
#   - selects one independent baseline/pre-intervention sample per mouse,
#   - applies the primary ASV filter,
#   - loads the saved phylogenetic tree created by Script 00,
#   - calculates Bray-Curtis, Jaccard, weighted and unweighted UniFrac,
#   - runs omnibus and pairwise PERMANOVA plus PERMDISP,
#   - tests sensitivity to alternative ASV prevalence filters,
#   - loads the finalized alpha-diversity panels from Script 01,
#   - saves the combined alpha + beta main-text figure,
#   - saves the supplementary UniFrac figure and reader-facing tables.
#
# Primary analysis:
#   1. Start from the final preprocessed baseline phyloseq object.
#   2. Retain one independent pre-intervention sample per mouse FIRST.
#   3. Filter ASVs on the independent-sample dataset:
#        - >= 5 reads across the dataset
#        - present in >= 10% of independent samples
#   4. Calculate Bray–Curtis, Jaccard, weighted UniFrac and unweighted UniFrac.
#   5. Run PERMANOVA, pairwise PERMANOVA and beta-dispersion tests.
#
# Filtering sensitivity analysis:
#   Compare Bray–Curtis and Jaccard results across four ASV-filtering schemes:
#     A. >=5 total reads only
#     B. >=5 total reads + >=5% global prevalence
#     C. >=5 total reads + >=10% global prevalence (PRIMARY)
#     D. >=5 total reads + >=10% prevalence within at least one cohort
#
# Important:
#   ASV prevalence filtering is performed AFTER selecting one independent
#   sample per mouse. This avoids defining prevalence from repeated observations.
#
###############################################################################


# =============================================================================
# 0. SETUP
# =============================================================================

rm(list = ls())

library(ape)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(patchwork)
library(phyloseq)
library(tibble)
library(tidyr)
library(vegan)


# =============================================================================
# 1. LOAD DATA AND DEFINE OUTPUT DIRECTORIES
# =============================================================================

# Run Scripts/00_DataPrep_Baseline.R first whenever preprocessing changes.
ps_all <- readRDS(
  "Data/Processed/baseline_microbiome_phyloseq.rds"
)

# Load the phylogenetic tree created once in 00_DataPrep_Baseline.R.
# This file should be tracked with the project/GitHub files. The beta-diversity
# script never reconstructs the tree.
tree_file <- file.path(
  "Data",
  "Processed",
  "baseline_phylogenetic_tree.rds"
)

if (!file.exists(tree_file)) {
  stop(
    "Saved phylogenetic tree not found: ",
    tree_file,
    "\nRun 00_DataPrep_Baseline_TREE_ATTACHED.R first and make sure the saved ",
    "tree file is present in Data/Processed."
  )
}

baseline_tree <- readRDS(
  tree_file
)

# Confirm that every ASV in ps_all is represented in the saved tree.
missing_tree_taxa <- setdiff(
  phyloseq::taxa_names(ps_all),
  baseline_tree$tip.label
)

if (length(missing_tree_taxa) > 0) {
  stop(
    "The saved tree is missing ",
    length(missing_tree_taxa),
    " ASVs present in ps_all. Re-run the data-preparation script so the tree ",
    "and phyloseq object are generated from the same final ASV set."
  )
}

# If the saved tree contains any extra tips, prune them to the current ps_all
# object and attach the saved tree explicitly.
baseline_tree <- ape::keep.tip(
  baseline_tree,
  phyloseq::taxa_names(ps_all)
)

phyloseq::phy_tree(ps_all) <- baseline_tree

stopifnot(
  setequal(
    phyloseq::phy_tree(ps_all)$tip.label,
    phyloseq::taxa_names(ps_all)
  )
)

results_dir <- "Results/02_Beta_Diversity"
tables_dir <- file.path(results_dir, "Tables")
figures_dir <- file.path(results_dir, "Figures")

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Output policy:
#   - save ONE multi-sheet workbook containing only manuscript-relevant results;
#   - save only complete main/supplementary figures, never individual panels;
#   - retain detailed threshold-by-threshold objects in memory for analysis/QC,
#     but do not export bulky intermediate result files.


# =============================================================================
# 2. ANALYSIS SETTINGS
# =============================================================================

cohort_order <- c(
  "Wild:As",
  "Wild-like:As",
  "Lab-like:As",
  "Lab:Mm"
)

cohort_colours <- c(
  "Wild:As"      = "#74C476",
  "Wild-like:As" = "#6BAED6",
  "Lab-like:As"  = "#BC80BD",
  "Lab:Mm"       = "#FB8072"
)

cohort_labels <- c(
  "Wild:As"      = "Wild:As",
  "Wild-like:As" = "Wild-like:As",
  "Lab-like:As"  = "Lab-like:As",
  "Lab:Mm"       = "Lab:Mm"
)

MIN_TOTAL_READS <- 5
PRIMARY_PREVALENCE <- 0.10
N_PERMUTATIONS <- 9999
PERMUTATION_SEED <- 999


# =============================================================================
# 3. RETAIN SAMPLES WITH REQUIRED METADATA
# =============================================================================

# SeqPlate is deliberately not required or modelled here.
# The beta-diversity model uses SampleOrigin + Sex.
required_metadata <- c(
  "SampleOrigin",
  "Sex",
  "mouse_ID_dataset"
)

metadata_all <- data.frame(
  phyloseq::sample_data(ps_all),
  check.names = FALSE
)

keep_samples <- stats::complete.cases(
  metadata_all[
    ,
    required_metadata,
    drop = FALSE
  ]
)

ps_available <- phyloseq::prune_samples(
  keep_samples,
  ps_all
)

ps_available <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_available) > 0,
  ps_available
)


# =============================================================================
# 4. SELECT ONE INDEPENDENT SAMPLE PER MOUSE BEFORE ASV FILTERING
# =============================================================================

sample_metadata <- data.frame(
  phyloseq::sample_data(ps_available),
  check.names = FALSE
) |>
  tibble::rownames_to_column(
    "PhyloseqSampleID"
  ) |>
  dplyr::mutate(
    Timept_numeric = suppressWarnings(
      as.numeric(
        as.character(Timept)
      )
    ),
    SampleDate_parsed = suppressWarnings(
      as.Date(
        as.character(SampleDate)
      )
    )
  )

# -------------------------------------------------------------------------
# Wild mice
# -------------------------------------------------------------------------
# Wild mice do not have the same Day-0 experimental time-point structure.
# If multiple samples are available, retain the earliest available sample.

selected_wild <- sample_metadata |>
  dplyr::filter(
    SampleOrigin == "Wild:As"
  ) |>
  dplyr::arrange(
    mouse_ID_dataset,
    is.na(SampleDate_parsed),
    SampleDate_parsed,
    PhyloseqSampleID
  ) |>
  dplyr::group_by(
    mouse_ID_dataset
  ) |>
  dplyr::slice_head(
    n = 1
  ) |>
  dplyr::ungroup()

# -------------------------------------------------------------------------
# Lab-reared cohorts
# -------------------------------------------------------------------------
# Retain the sample closest to Day 0.
# For an exact tie around zero, positive time points are preferred, matching
# the previous analysis.

selected_other <- sample_metadata |>
  dplyr::filter(
    SampleOrigin != "Wild:As"
  ) |>
  dplyr::filter(
    !is.na(Timept_numeric)
  ) |>
  dplyr::arrange(
    mouse_ID_dataset,
    abs(Timept_numeric),
    dplyr::desc(Timept_numeric),
    PhyloseqSampleID
  ) |>
  dplyr::group_by(
    mouse_ID_dataset
  ) |>
  dplyr::slice_head(
    n = 1
  ) |>
  dplyr::ungroup()

selected_samples <- dplyr::bind_rows(
  selected_wild,
  selected_other
)

# Create the independent-sample phyloseq object BEFORE any prevalence filter.
ps_single_unfiltered <- phyloseq::prune_samples(
  selected_samples$PhyloseqSampleID,
  ps_available
)

ps_single_unfiltered <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_single_unfiltered) > 0,
  ps_single_unfiltered
)

# Confirm one sample per mouse.
mouse_counts_after_selection <- data.frame(
  phyloseq::sample_data(ps_single_unfiltered),
  check.names = FALSE
) |>
  dplyr::count(
    mouse_ID_dataset,
    name = "n_samples"
  )

stopifnot(
  all(
    mouse_counts_after_selection$n_samples == 1
  )
)

# Set factor order now and retain it in all filtered objects.
phyloseq::sample_data(ps_single_unfiltered)$SampleOrigin <- factor(
  phyloseq::sample_data(ps_single_unfiltered)$SampleOrigin,
  levels = cohort_order
)

phyloseq::sample_data(ps_single_unfiltered)$Sex <- factor(
  phyloseq::sample_data(ps_single_unfiltered)$Sex
)

cat(
  "\nIndependent samples retained:\n"
)

print(
  table(
    phyloseq::sample_data(ps_single_unfiltered)$SampleOrigin,
    useNA = "ifany"
  )
)

cat(
  "\nTotal independent samples = ",
  phyloseq::nsamples(ps_single_unfiltered),
  "\n",
  sep = ""
)


# =============================================================================
# 5. ASV FILTERING FUNCTIONS
# =============================================================================

# Helper to obtain samples x taxa matrix.
get_sample_taxa_matrix <- function(ps_object) {

  otu <- as(
    phyloseq::otu_table(ps_object),
    "matrix"
  )

  if (phyloseq::taxa_are_rows(ps_object)) {
    otu <- t(otu)
  }

  otu
}


# -------------------------------------------------------------------------
# >=5 total reads only
# -------------------------------------------------------------------------

filter_min_total_reads <- function(
    ps_object,
    min_total_reads = 5) {

  keep_taxa <- phyloseq::taxa_sums(ps_object) >=
    min_total_reads

  ps_filtered <- phyloseq::prune_taxa(
    keep_taxa,
    ps_object
  )

  phyloseq::prune_taxa(
    phyloseq::taxa_sums(ps_filtered) > 0,
    ps_filtered
  )
}


# -------------------------------------------------------------------------
# Global prevalence threshold
# -------------------------------------------------------------------------

filter_global_prevalence <- function(
    ps_object,
    min_total_reads = 5,
    prevalence_threshold = 0.10) {

  # First apply the total-read threshold.
  ps_filtered <- filter_min_total_reads(
    ps_object,
    min_total_reads = min_total_reads
  )

  otu <- get_sample_taxa_matrix(
    ps_filtered
  )

  prevalence <- colMeans(
    otu > 0
  )

  keep_taxa <- prevalence >=
    prevalence_threshold

  ps_filtered <- phyloseq::prune_taxa(
    keep_taxa,
    ps_filtered
  )

  phyloseq::prune_taxa(
    phyloseq::taxa_sums(ps_filtered) > 0,
    ps_filtered
  )
}


# -------------------------------------------------------------------------
# Group-specific prevalence threshold
# -------------------------------------------------------------------------
# Retain an ASV if it is present in >= threshold of samples within AT LEAST
# ONE microbiota-origin group.

filter_within_group_prevalence <- function(
    ps_object,
    min_total_reads = 5,
    prevalence_threshold = 0.10) {

  ps_filtered <- filter_min_total_reads(
    ps_object,
    min_total_reads = min_total_reads
  )

  otu <- get_sample_taxa_matrix(
    ps_filtered
  )

  metadata <- data.frame(
    phyloseq::sample_data(ps_filtered),
    check.names = FALSE
  )

  stopifnot(
    identical(
      rownames(metadata),
      rownames(otu)
    )
  )

  group_levels <- levels(
    droplevels(
      factor(metadata$SampleOrigin)
    )
  )

  prevalence_by_group <- vapply(
    group_levels,
    function(group_name) {

      group_rows <- metadata$SampleOrigin ==
        group_name

      colMeans(
        otu[
          group_rows,
          ,
          drop = FALSE
        ] > 0
      )
    },
    FUN.VALUE = numeric(
      ncol(otu)
    )
  )

  if (is.null(dim(prevalence_by_group))) {
    prevalence_by_group <- matrix(
      prevalence_by_group,
      ncol = 1
    )
  }

  keep_taxa <- apply(
    prevalence_by_group,
    1,
    max,
    na.rm = TRUE
  ) >= prevalence_threshold

  ps_filtered <- phyloseq::prune_taxa(
    keep_taxa,
    ps_filtered
  )

  phyloseq::prune_taxa(
    phyloseq::taxa_sums(ps_filtered) > 0,
    ps_filtered
  )
}


# =============================================================================
# 6. PRIMARY FILTER: >=5 TOTAL READS + >=10% GLOBAL PREVALENCE
# =============================================================================

ps_single <- filter_global_prevalence(
  ps_single_unfiltered,
  min_total_reads = MIN_TOTAL_READS,
  prevalence_threshold = PRIMARY_PREVALENCE
)

cat(
  "\nPrimary beta-diversity filter:\n",
  "  >= ",
  MIN_TOTAL_READS,
  " total reads\n",
  "  >= ",
  PRIMARY_PREVALENCE * 100,
  "% global prevalence\n",
  "  ASVs retained = ",
  phyloseq::ntaxa(ps_single),
  "\n",
  sep = ""
)

# Confirm tree contains all retained ASVs.
if (!is.null(
  phyloseq::phy_tree(
    ps_single,
    errorIfNULL = FALSE
  )
)) {

  retained_taxa <- phyloseq::taxa_names(
    ps_single
  )

  tree_taxa <- phyloseq::phy_tree(
    ps_single
  )$tip.label

  stopifnot(
    all(
      retained_taxa %in% tree_taxa
    )
  )
}


# =============================================================================
# 7. PRIMARY BETA-DIVERSITY OBJECTS
# =============================================================================

# Relative abundance for Bray–Curtis.
ps_comp <- phyloseq::transform_sample_counts(
  ps_single,
  function(x) {
    x / sum(x)
  }
)

ps_meta <- data.frame(
  phyloseq::sample_data(ps_single),
  check.names = FALSE
)

ps_meta$SampleOrigin <- factor(
  ps_meta$SampleOrigin,
  levels = cohort_order
)

ps_meta$Sex <- factor(
  ps_meta$Sex
)

stopifnot(
  identical(
    rownames(ps_meta),
    phyloseq::sample_names(ps_single)
  )
)

# Bray–Curtis: abundance weighted.
bray_ps <- phyloseq::distance(
  ps_comp,
  method = "bray"
)

# Jaccard: presence/absence.
jaccard_ps <- phyloseq::distance(
  ps_single,
  method = "jaccard",
  binary = TRUE
)

# Weighted UniFrac.
wunif_ps <- phyloseq::distance(
  ps_single,
  method = "wunifrac"
)

# Unweighted UniFrac.
unif_ps <- phyloseq::distance(
  ps_single,
  method = "unifrac"
)


# =============================================================================
# 8. PRIMARY PCoA FIGURES
# =============================================================================

make_pcoa_plot <- function(
    ps_object,
    ordination,
    title) {

  phyloseq::plot_ordination(
    ps_object,
    ordination,
    color = "SampleOrigin"
  ) +
    ggplot2::geom_point(
      size = 2
    ) +
    ggplot2::stat_ellipse() +
    ggplot2::scale_colour_manual(
      values = cohort_colours,
      drop = FALSE
    ) +
    ggplot2::theme_bw(
      base_size = 14
    ) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle(
      title
    )
}

ord_bray <- phyloseq::ordinate(
  ps_comp,
  method = "PCoA",
  distance = "bray"
)

ord_jaccard <- phyloseq::ordinate(
  ps_single,
  method = "PCoA",
  distance = "jaccard",
  binary = TRUE
)

ord_wunif <- phyloseq::ordinate(
  ps_single,
  method = "PCoA",
  distance = "wunifrac"
)

ord_unif <- phyloseq::ordinate(
  ps_single,
  method = "PCoA",
  distance = "unifrac"
)

plot_bray <- make_pcoa_plot(
  ps_comp,
  ord_bray,
  "Bray–Curtis"
)

plot_jaccard <- make_pcoa_plot(
  ps_single,
  ord_jaccard,
  "Jaccard"
)

plot_wunif <- make_pcoa_plot(
  ps_single,
  ord_wunif,
  "Weighted UniFrac"
)

plot_unif <- make_pcoa_plot(
  ps_single,
  ord_unif,
  "Unweighted UniFrac"
)

# -------------------------------------------------------------------------
# Main-text combined alpha + beta diversity figure
# -------------------------------------------------------------------------
#
# The alpha-diversity script saves the finalized Shannon and Faith's PD
# plotting objects after repeated rarefaction. We load those objects here so
# the main manuscript figure uses exactly the same alpha-diversity estimates
# as the statistical analysis rather than recalculating alpha diversity.

alpha_panel_file <- file.path(
  "Results",
  "01_Alpha_Diversity",
  "Objects",
  "alpha_main_panels_for_combined_figure.rds"
)

if (!file.exists(alpha_panel_file)) {
  stop(
    "Run 01_AlphaDiversity_Baseline_COMBINED_READY_FIXED.R first so that ",
    "Results/01_Alpha_Diversity/Objects/",
    "alpha_main_panels_for_combined_figure.rds exists."
  )
}

alpha_main_panels <- readRDS(
  alpha_panel_file
)

plot_shannon <- alpha_main_panels$shannon_plot
plot_faith <- alpha_main_panels$faith_pd_plot

# Apply a consistent manuscript theme to all four panels.
combined_theme <- ggplot2::theme_classic(
  base_size = 12
) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      colour = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(
      size = 12,
      colour = "black"
    ),
    axis.text = ggplot2::element_text(
      size = 10,
      colour = "black"
    ),
    plot.title = ggplot2::element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = ggplot2::margin(
        b = 6
      )
    ),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(
      size = 10
    ),
    plot.margin = ggplot2::margin(
      8,
      8,
      8,
      8
    )
  )

plot_shannon_final <- plot_shannon +
  ggplot2::scale_x_discrete(
    labels = cohort_labels
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Shannon diversity"
  ) +
  combined_theme +
  ggplot2::guides(
    colour = "none",
    fill = "none"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1,
      vjust = 1
    )
  )

plot_faith_final <- plot_faith +
  ggplot2::scale_x_discrete(
    labels = cohort_labels
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Faith's phylogenetic diversity"
  ) +
  combined_theme +
  ggplot2::guides(
    colour = "none",
    fill = "none"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1,
      vjust = 1
    )
  )

plot_bray_final <- plot_bray +
  combined_theme +
  ggplot2::theme(
    legend.position = "bottom"
  )

plot_jaccard_final <- plot_jaccard +
  combined_theme +
  ggplot2::guides(
    colour = "none"
  )

combined_diversity_figure <- patchwork::wrap_plots(
  A = plot_shannon_final,
  B = plot_faith_final,
  C = plot_bray_final,
  D = plot_jaccard_final,
  design = "AB\nCD"
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        size = 16,
        face = "bold"
      )
    )
  )

print(
  combined_diversity_figure
)

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Main_Alpha_Beta_Diversity.pdf"
  ),
  plot = combined_diversity_figure,
  width = 8,
  height = 8,
  units = "in",
  device = grDevices::cairo_pdf
)

# Supplementary UniFrac figure.
beta_unifrac <- ggpubr::ggarrange(
  plot_wunif + ggplot2::theme(legend.position = "none"),
  plot_unif + ggplot2::theme(legend.position = "none"),
  labels = c("A", "B"),
  nrow = 1,
  common.legend = TRUE,
  legend = "top"
)

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Beta_Diversity_UniFrac.pdf"
  ),
  plot = beta_unifrac,
  width = 10,
  height = 5
)


# =============================================================================
# 9. PRIMARY PERMANOVA / DISPERSION HELPERS
# =============================================================================

extract_adonis_term <- function(
    adonis_object,
    term = "SampleOrigin") {

  x <- as.data.frame(
    adonis_object
  )

  x$Term <- rownames(x)
  rownames(x) <- NULL

  x |>
    dplyr::filter(
      Term == term
    )
}


run_pairwise_permanova <- function(
    distance_object,
    metadata,
    permutations = 9999,
    seed = 999) {

  groups <- levels(
    droplevels(
      metadata$SampleOrigin
    )
  )

  comparisons <- utils::combn(
    groups,
    2,
    simplify = FALSE
  )

  results <- lapply(
    seq_along(comparisons),
    function(i) {

      comparison <- comparisons[[i]]

      keep <- metadata$SampleOrigin %in%
        comparison

      metadata_sub <- metadata[
        keep,
        ,
        drop = FALSE
      ]

      metadata_sub$SampleOrigin <- droplevels(
        metadata_sub$SampleOrigin
      )

      # Subset the distance matrix by sample IDs.
      distance_matrix <- as.matrix(
        distance_object
      )

      sample_ids <- rownames(
        metadata_sub
      )

      distance_sub <- stats::as.dist(
        distance_matrix[
          sample_ids,
          sample_ids,
          drop = FALSE
        ]
      )

      set.seed(
        seed + i
      )

      fit <- vegan::adonis2(
        distance_sub ~ SampleOrigin + Sex,
        data = metadata_sub,
        permutations = permutations,
        by = "margin"
      )

      term_result <- extract_adonis_term(
        fit,
        term = "SampleOrigin"
      )

      data.frame(
        Comparison = paste(
          comparison[1],
          "vs",
          comparison[2]
        ),
        F = term_result$F,
        R2 = term_result$R2,
        P = term_result$`Pr(>F)`,
        stringsAsFactors = FALSE
      )
    }
  )

  out <- dplyr::bind_rows(
    results
  )

  out$P_BH <- stats::p.adjust(
    out$P,
    method = "BH"
  )

  out
}


run_distance_analysis <- function(
    distance_object,
    metadata,
    distance_name,
    permutations = 9999,
    seed = 999) {

  set.seed(seed)

  overall <- vegan::adonis2(
    distance_object ~ SampleOrigin + Sex,
    data = metadata,
    permutations = permutations,
    by = "margin"
  )

  overall_origin <- extract_adonis_term(
    overall,
    term = "SampleOrigin"
  )

  pairwise <- run_pairwise_permanova(
    distance_object = distance_object,
    metadata = metadata,
    permutations = permutations,
    seed = seed
  )

  dispersion <- vegan::betadisper(
    distance_object,
    group = metadata$SampleOrigin
  )

  set.seed(seed)

  dispersion_test <- vegan::permutest(
    dispersion,
    permutations = permutations
  )

  dispersion_table <- as.data.frame(
    dispersion_test$tab
  )

  list(
    overall_full = overall,
    overall_summary = data.frame(
      Distance = distance_name,
      Term = "SampleOrigin",
      F = overall_origin$F,
      R2 = overall_origin$R2,
      P = overall_origin$`Pr(>F)`,
      stringsAsFactors = FALSE
    ),
    pairwise = dplyr::mutate(
      pairwise,
      Distance = distance_name,
      .before = 1
    ),
    dispersion = data.frame(
      Distance = distance_name,
      F = dispersion_table$F[1],
      P = dispersion_table$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  )
}


# =============================================================================
# 10. PRIMARY STATISTICAL ANALYSES
# =============================================================================

bray_results <- run_distance_analysis(
  distance_object = bray_ps,
  metadata = ps_meta,
  distance_name = "Bray-Curtis",
  permutations = N_PERMUTATIONS,
  seed = PERMUTATION_SEED
)

jaccard_results <- run_distance_analysis(
  distance_object = jaccard_ps,
  metadata = ps_meta,
  distance_name = "Jaccard",
  permutations = N_PERMUTATIONS,
  seed = PERMUTATION_SEED
)

wunif_results <- run_distance_analysis(
  distance_object = wunif_ps,
  metadata = ps_meta,
  distance_name = "Weighted UniFrac",
  permutations = N_PERMUTATIONS,
  seed = PERMUTATION_SEED
)

unif_results <- run_distance_analysis(
  distance_object = unif_ps,
  metadata = ps_meta,
  distance_name = "Unweighted UniFrac",
  permutations = N_PERMUTATIONS,
  seed = PERMUTATION_SEED
)

primary_omnibus <- dplyr::bind_rows(
  bray_results$overall_summary,
  jaccard_results$overall_summary,
  wunif_results$overall_summary,
  unif_results$overall_summary
)

primary_pairwise <- dplyr::bind_rows(
  bray_results$pairwise,
  jaccard_results$pairwise,
  wunif_results$pairwise,
  unif_results$pairwise
)

primary_dispersion <- dplyr::bind_rows(
  bray_results$dispersion,
  jaccard_results$dispersion,
  wunif_results$dispersion,
  unif_results$dispersion
)

print(primary_omnibus)
print(primary_pairwise)
print(primary_dispersion)


# =============================================================================
# 11. FILTERING SENSITIVITY ANALYSIS
# =============================================================================

# The goal is NOT to identify an "optimal" prevalence cutoff from the data.
# Instead, we ask whether the beta-diversity conclusions are robust across a
# reasonable range of global prevalence thresholds while holding the minimum
# total-read requirement constant at >=5 reads.
#
# Global prevalence thresholds tested:
#   0%, 2.5%, 5%, 7.5%, 10%, 12.5%, and 15%
#
# Sensitivity is evaluated for all four beta-diversity measures used in the
# manuscript:
#   - Bray-Curtis
#   - Jaccard
#   - Weighted UniFrac
#   - Unweighted UniFrac
#
# The 10% global-prevalence filter remains the PRIMARY analysis.
#
# We additionally retain one biologically targeted alternative:
#   >=10% prevalence within at least one cohort.
# This directly addresses the concern that a taxon can be common within one
# cohort but uncommon across the complete four-cohort dataset.

GLOBAL_PREVALENCE_SWEEP <- c(
  0,
  0.025,
  0.05,
  0.075,
  0.10,
  0.125,
  0.15
)


# -------------------------------------------------------------------------
# 11A. Build globally filtered phyloseq objects
# -------------------------------------------------------------------------

global_filtering_schemes <- lapply(
  GLOBAL_PREVALENCE_SWEEP,
  function(prevalence_cutoff) {

    if (prevalence_cutoff == 0) {

      filter_min_total_reads(
        ps_single_unfiltered,
        min_total_reads = MIN_TOTAL_READS
      )

    } else {

      filter_global_prevalence(
        ps_single_unfiltered,
        min_total_reads = MIN_TOTAL_READS,
        prevalence_threshold = prevalence_cutoff
      )
    }
  }
)

names(global_filtering_schemes) <- paste0(
  "Global prevalence >= ",
  GLOBAL_PREVALENCE_SWEEP * 100,
  "%"
)

# Targeted group-specific alternative.
ps_within_group_10 <- filter_within_group_prevalence(
  ps_single_unfiltered,
  min_total_reads = MIN_TOTAL_READS,
  prevalence_threshold = 0.10
)


# -------------------------------------------------------------------------
# 11B. Run Bray-Curtis and Jaccard across the GLOBAL threshold sweep
# -------------------------------------------------------------------------

global_sensitivity_results <- list()

for (i in seq_along(GLOBAL_PREVALENCE_SWEEP)) {

  prevalence_cutoff <- GLOBAL_PREVALENCE_SWEEP[i]
  scheme_name <- names(global_filtering_schemes)[i]
  ps_scheme <- global_filtering_schemes[[i]]

  message(
    "Global filtering sensitivity: ",
    prevalence_cutoff * 100,
    "%"
  )

  # Preserve factor structure.
  phyloseq::sample_data(ps_scheme)$SampleOrigin <- factor(
    phyloseq::sample_data(ps_scheme)$SampleOrigin,
    levels = cohort_order
  )

  phyloseq::sample_data(ps_scheme)$Sex <- factor(
    phyloseq::sample_data(ps_scheme)$Sex
  )

  metadata_scheme <- data.frame(
    phyloseq::sample_data(ps_scheme),
    check.names = FALSE
  )

  metadata_scheme$SampleOrigin <- factor(
    metadata_scheme$SampleOrigin,
    levels = cohort_order
  )

  metadata_scheme$Sex <- factor(
    metadata_scheme$Sex
  )

  # Relative abundance for Bray-Curtis.
  ps_scheme_comp <- phyloseq::transform_sample_counts(
    ps_scheme,
    function(x) {
      x / sum(x)
    }
  )

  distance_objects <- list(
    "Bray-Curtis" = phyloseq::distance(
      ps_scheme_comp,
      method = "bray"
    ),

    "Jaccard" = phyloseq::distance(
      ps_scheme,
      method = "jaccard",
      binary = TRUE
    ),

    "Weighted UniFrac" = phyloseq::distance(
      ps_scheme,
      method = "wunifrac"
    ),

    "Unweighted UniFrac" = phyloseq::distance(
      ps_scheme,
      method = "unifrac"
    )
  )

  for (distance_name in names(distance_objects)) {

    result <- run_distance_analysis(
      distance_object = distance_objects[[distance_name]],
      metadata = metadata_scheme,
      distance_name = distance_name,
      permutations = N_PERMUTATIONS,
      seed = PERMUTATION_SEED
    )

    overall_row <- result$overall_summary |>
      dplyr::mutate(
        Filtering_scheme = scheme_name,
        Prevalence_cutoff = prevalence_cutoff,
        Prevalence_percent = prevalence_cutoff * 100,
        ASVs_retained = phyloseq::ntaxa(ps_scheme),
        Test = "Overall",
        Comparison = NA_character_,
        P_BH = NA_real_,
        Dispersion_F = result$dispersion$F,
        Dispersion_P = result$dispersion$P,
        .before = 1
      ) |>
      dplyr::select(
        Filtering_scheme,
        Prevalence_cutoff,
        Prevalence_percent,
        ASVs_retained,
        Distance,
        Test,
        Comparison,
        F,
        R2,
        P,
        P_BH,
        Dispersion_F,
        Dispersion_P
      )

    pairwise_rows <- result$pairwise |>
      dplyr::mutate(
        Filtering_scheme = scheme_name,
        Prevalence_cutoff = prevalence_cutoff,
        Prevalence_percent = prevalence_cutoff * 100,
        ASVs_retained = phyloseq::ntaxa(ps_scheme),
        Test = "Pairwise",
        Dispersion_F = NA_real_,
        Dispersion_P = NA_real_,
        .before = 1
      ) |>
      dplyr::select(
        Filtering_scheme,
        Prevalence_cutoff,
        Prevalence_percent,
        ASVs_retained,
        Distance,
        Test,
        Comparison,
        F,
        R2,
        P,
        P_BH,
        Dispersion_F,
        Dispersion_P
      )

    global_sensitivity_results[[
      paste(
        prevalence_cutoff,
        distance_name,
        sep = "__"
      )
    ]] <- dplyr::bind_rows(
      overall_row,
      pairwise_rows
    )
  }
}

global_filtering_sensitivity <- dplyr::bind_rows(
  global_sensitivity_results
) |>
  dplyr::arrange(
    Distance,
    Prevalence_cutoff,
    Test,
    Comparison
  )


# -------------------------------------------------------------------------
# 11C. Run the targeted >=10% WITHIN-COHORT alternative
# -------------------------------------------------------------------------

phyloseq::sample_data(ps_within_group_10)$SampleOrigin <- factor(
  phyloseq::sample_data(ps_within_group_10)$SampleOrigin,
  levels = cohort_order
)

phyloseq::sample_data(ps_within_group_10)$Sex <- factor(
  phyloseq::sample_data(ps_within_group_10)$Sex
)

metadata_within_group <- data.frame(
  phyloseq::sample_data(ps_within_group_10),
  check.names = FALSE
)

metadata_within_group$SampleOrigin <- factor(
  metadata_within_group$SampleOrigin,
  levels = cohort_order
)

metadata_within_group$Sex <- factor(
  metadata_within_group$Sex
)

ps_within_group_comp <- phyloseq::transform_sample_counts(
  ps_within_group_10,
  function(x) {
    x / sum(x)
  }
)

within_group_distances <- list(
  "Bray-Curtis" = phyloseq::distance(
    ps_within_group_comp,
    method = "bray"
  ),

  "Jaccard" = phyloseq::distance(
    ps_within_group_10,
    method = "jaccard",
    binary = TRUE
  ),

  "Weighted UniFrac" = phyloseq::distance(
    ps_within_group_10,
    method = "wunifrac"
  ),

  "Unweighted UniFrac" = phyloseq::distance(
    ps_within_group_10,
    method = "unifrac"
  )
)

within_group_sensitivity <- list()

for (distance_name in names(within_group_distances)) {

  result <- run_distance_analysis(
    distance_object = within_group_distances[[distance_name]],
    metadata = metadata_within_group,
    distance_name = distance_name,
    permutations = N_PERMUTATIONS,
    seed = PERMUTATION_SEED
  )

  overall_row <- result$overall_summary |>
    dplyr::mutate(
      Filtering_scheme = "Within-group prevalence >=10%",
      Prevalence_cutoff = 0.10,
      Prevalence_percent = 10,
      ASVs_retained = phyloseq::ntaxa(ps_within_group_10),
      Test = "Overall",
      Comparison = NA_character_,
      P_BH = NA_real_,
      Dispersion_F = result$dispersion$F,
      Dispersion_P = result$dispersion$P,
      .before = 1
    ) |>
    dplyr::select(
      Filtering_scheme,
      Prevalence_cutoff,
      Prevalence_percent,
      ASVs_retained,
      Distance,
      Test,
      Comparison,
      F,
      R2,
      P,
      P_BH,
      Dispersion_F,
      Dispersion_P
    )

  pairwise_rows <- result$pairwise |>
    dplyr::mutate(
      Filtering_scheme = "Within-group prevalence >=10%",
      Prevalence_cutoff = 0.10,
      Prevalence_percent = 10,
      ASVs_retained = phyloseq::ntaxa(ps_within_group_10),
      Test = "Pairwise",
      Dispersion_F = NA_real_,
      Dispersion_P = NA_real_,
      .before = 1
    ) |>
    dplyr::select(
      Filtering_scheme,
      Prevalence_cutoff,
      Prevalence_percent,
      ASVs_retained,
      Distance,
      Test,
      Comparison,
      F,
      R2,
      P,
      P_BH,
      Dispersion_F,
      Dispersion_P
    )

  within_group_sensitivity[[distance_name]] <- dplyr::bind_rows(
    overall_row,
    pairwise_rows
  )
}

within_group_filtering_sensitivity <- dplyr::bind_rows(
  within_group_sensitivity
)


# -------------------------------------------------------------------------
# 11D. Compact summaries for interpretation
# -------------------------------------------------------------------------

# Overall SampleOrigin effects across the global prevalence sweep.
global_sensitivity_overall <- global_filtering_sensitivity |>
  dplyr::filter(
    Test == "Overall"
  ) |>
  dplyr::select(
    Prevalence_percent,
    ASVs_retained,
    Distance,
    F,
    R2,
    P,
    Dispersion_F,
    Dispersion_P
  )

# Pairwise robustness summary:
# Rather than forcing readers to inspect every threshold-specific pairwise test,
# summarize the range of pairwise R2 and the largest BH-adjusted P value across
# the complete global threshold sweep.
global_pairwise_robustness <- global_filtering_sensitivity |>
  dplyr::filter(
    Test == "Pairwise"
  ) |>
  dplyr::group_by(
    Distance,
    Comparison
  ) |>
  dplyr::summarise(
    R2_min = min(
      R2,
      na.rm = TRUE
    ),
    R2_max = max(
      R2,
      na.rm = TRUE
    ),
    max_P_BH = max(
      P_BH,
      na.rm = TRUE
    ),
    significant_at_all_thresholds = all(
      P_BH < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

within_group_overall <- within_group_filtering_sensitivity |>
  dplyr::filter(
    Test == "Overall"
  ) |>
  dplyr::select(
    Filtering_scheme,
    ASVs_retained,
    Distance,
    F,
    R2,
    P,
    Dispersion_F,
    Dispersion_P
  )


# -------------------------------------------------------------------------
# 11E. Compact supplementary sensitivity figure
# -------------------------------------------------------------------------

# Panel A: number of ASVs retained as the GLOBAL prevalence cutoff increases.
asv_retention_plot_data <- global_sensitivity_overall |>
  dplyr::distinct(
    Prevalence_percent,
    ASVs_retained
  )

asv_retention_plot <- ggplot2::ggplot(
  asv_retention_plot_data,
  ggplot2::aes(
    x = Prevalence_percent,
    y = ASVs_retained
  )
) +
  ggplot2::geom_line(
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    size = 2.2
  ) +
  ggplot2::geom_vline(
    xintercept = 10,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::labs(
    x = "Global prevalence cutoff (%)",
    y = "ASVs retained"
  ) +
  ggplot2::theme_bw(
    base_size = 12
  )

# Panel B: overall PERMANOVA effect size across the same thresholds.
effect_size_sensitivity_plot <- ggplot2::ggplot(
  global_sensitivity_overall,
  ggplot2::aes(
    x = Prevalence_percent,
    y = R2,
    linetype = Distance,
    shape = Distance,
    group = Distance
  )
) +
  ggplot2::geom_line(
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    size = 2.2
  ) +
  ggplot2::geom_vline(
    xintercept = 10,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Bray-Curtis" = 16,
      "Jaccard" = 17,
      "Weighted UniFrac" = 15,
      "Unweighted UniFrac" = 18
    )
  ) +
  ggplot2::labs(
    x = "Global prevalence cutoff (%)",
    y = expression(PERMANOVA~R^2),
    linetype = "Distance",
    shape = "Distance"
  ) +
  ggplot2::theme_bw(
    base_size = 12
  ) +
  ggplot2::guides(
    linetype = ggplot2::guide_legend(
      title = "Distance",
      nrow = 2,
      byrow = TRUE
    ),
    shape = ggplot2::guide_legend(
      title = "Distance",
      nrow = 2,
      byrow = TRUE
    )
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = ggplot2::element_text(
      size = 10
    ),
    legend.text = ggplot2::element_text(
      size = 9
    )
  )

filtering_sensitivity_figure <- ggpubr::ggarrange(
  asv_retention_plot,
  effect_size_sensitivity_plot,
  labels = c(
    "A",
    "B"
  ),
  nrow = 1,
  widths = c(
    1,
    1.2
  )
)

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Beta_Diversity_Filtering_Sensitivity.pdf"
  ),
  plot = filtering_sensitivity_figure,
  width = 9,
  height = 5
)

print(
  global_sensitivity_overall
)

print(
  global_pairwise_robustness
)

print(
  within_group_overall
)


# =============================================================================
# 12. SAVE ONE STREAMLINED BETA-DIVERSITY RESULTS WORKBOOK
# =============================================================================

# Keep only the sensitivity information needed to document robustness of the
# primary prevalence-filtering choice. Full threshold-specific pairwise results
# remain available in memory during the analysis but are not exported.

filtering_sensitivity_export <- dplyr::bind_rows(

  global_sensitivity_overall |>
    dplyr::transmute(
      Filtering_scheme = paste0(
        "Global prevalence >= ",
        Prevalence_percent,
        "%"
      ),
      Prevalence_percent = Prevalence_percent,
      ASVs_retained = ASVs_retained,
      Distance = Distance,
      F = F,
      R2 = R2,
      P = P,
      Dispersion_F = Dispersion_F,
      Dispersion_P = Dispersion_P
    ),

  within_group_overall |>
    dplyr::transmute(
      Filtering_scheme = Filtering_scheme,
      Prevalence_percent = 10,
      ASVs_retained = ASVs_retained,
      Distance = Distance,
      F = F,
      R2 = R2,
      P = P,
      Dispersion_F = Dispersion_F,
      Dispersion_P = Dispersion_P
    )
) |>
  dplyr::arrange(
    Distance,
    Prevalence_percent,
    Filtering_scheme
  )

# One reader-facing workbook for all key beta-diversity results.
#
# Sheets:
#   Omnibus PERMANOVA
#     Overall SampleOrigin effect for Bray-Curtis, Jaccard, weighted UniFrac,
#     and unweighted UniFrac.
#
#   Pairwise PERMANOVA
#     All six SampleOrigin comparisons for each distance metric with BH-adjusted
#     P values.
#
#   Dispersion Tests
#     PERMDISP results for each distance metric.
#
#   Filtering Sensitivity
#     Overall PERMANOVA and dispersion results across the global prevalence
#     sweep plus the targeted 10% within-cohort prevalence alternative.

beta_results_workbook <- list(
  "Omnibus PERMANOVA" =
    primary_omnibus,

  "Pairwise PERMANOVA" =
    primary_pairwise,

  "Dispersion Tests" =
    primary_dispersion,

  "Filtering Sensitivity" =
    filtering_sensitivity_export
)

openxlsx::write.xlsx(
  x = beta_results_workbook,
  file = file.path(
    tables_dir,
    "Beta_Diversity_results.xlsx"
  ),
  overwrite = TRUE,
  rowNames = FALSE
)


# =============================================================================
# 13. OUTPUT SUMMARY
# =============================================================================

cat(
  "\n============================================================\n",
  "BETA-DIVERSITY OUTPUTS\n",
  "============================================================\n",
  "Workbook:\n  ",
  file.path(
    tables_dir,
    "Beta_Diversity_results.xlsx"
  ),
  "\n\nComplete figures:\n  ",
  file.path(
    figures_dir,
    "Main_Alpha_Beta_Diversity.pdf"
  ),
  "\n  ",
  file.path(
    figures_dir,
    "Beta_Diversity_UniFrac.pdf"
  ),
  "\n  ",
  file.path(
    figures_dir,
    "Beta_Diversity_Filtering_Sensitivity.pdf"
  ),
  "\n============================================================\n",
  sep = ""
)

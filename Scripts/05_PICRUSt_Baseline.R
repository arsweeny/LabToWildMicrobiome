# ==============================================================================
# PICRUSt2 FUNCTIONAL ANALYSIS — BASELINE MOUSE MICROBIOME
# ==============================================================================
#
# PURPOSE
# -------
# This script analyzes PICRUSt2-predicted functional profiles from the baseline
# mouse microbiome dataset. It performs:
#
#   1. weighted NSTI quality control;
#   2. MetaCyc pathway preprocessing and descriptive QC;
#   3. functional beta-diversity analysis;
#   4. differential MetaCyc pathway analysis with MaAsLin3;
#   5. annotation of MetaCyc IDs with human-readable pathway names;
#   6. publication-oriented visualizations of functional differences.
#
# IMPORTANT INTERPRETATION
# ------------------------
# PICRUSt2 predicts functional potential from 16S phylogenetic placement and
# reference genomes. These are not directly measured metagenomic functions.
# Use wording such as "PICRUSt2-predicted functional potential" or
# "predicted MetaCyc pathway composition".
#
# ANALYSIS DECISIONS
# ------------------
# - MetaCyc pathways are the primary functional unit.
# - Pathways must occur in >=10% of samples.
# - Relative pathway abundances are used for Bray-Curtis beta diversity.
# - Beta diversity uses one sample per mouse because SampleOrigin is a
#   between-mouse predictor; within-mouse restricted permutations cannot test it.
# - Differential pathway analysis uses all eligible samples and accounts for
#   repeated sampling with (1 | mouse_ID_dataset).
# - Sex and log10 sequencing depth are fixed effects.
# - Samples missing Sex are excluded from MaAsLin3 only.
# - MaAsLin3 fits the abundance component only.
#
# REQUIRED INPUTS
# ---------------
# - ps_all, created by Scripts/00_DataPrep_Baseline.R
# - PICRUSt2/results/EC_metagenome_out/weighted_nsti.tsv.gz
# - PICRUSt2/results/pathways_out/path_abun_unstrat.tsv.gz
# - PICRUSt2/results/path_abun_unstrat_descriptions.tsv.gz
#   (created on Eddie with PICRUSt2 add_descriptions.py -m METACYC)
# ==============================================================================


# ==============================================================================
# 0. SETUP
# ==============================================================================

# Load the final preprocessed baseline phyloseq object.
# Run Scripts/00_DataPrep_Baseline.R first whenever preprocessing changes.
ps_all <- readRDS("Data/Processed/baseline_microbiome_phyloseq.rds")

required_packages <- c(
  "phyloseq", "dplyr", "tidyr", "tibble", "readr", "ggplot2",
  "vegan", "lmerTest", "emmeans", "maaslin3", "mirai"
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

picrust_dir <- "Results/PICRUSt_Outputs"

output_dir <- "Results/05_PICRUSt2_Functional_Prediction"
tables_dir <- file.path(output_dir, "Tables")
# Only manuscript/reproducibility-relevant tables are saved; QC/intermediate
# summaries are retained as R objects and printed to the console.
figures_dir <- file.path(output_dir, "Figures")
maaslin_dir <- file.path(output_dir, "Model_Output", "MetaCyc_MaAsLin3")

for (x in c(output_dir, tables_dir, figures_dir, maaslin_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

minimum_pathway_prevalence <- 0.10
q_threshold <- 0.05
force_rerun <- TRUE

group_levels <- c(
  "Wild:As",
  "Wild-like:As",
  "Lab-like:As",
  "Lab:Mm"
)

contrast_names <- c(
  "Wild-like:As vs Wild:As",
  "Lab-like:As vs Wild:As",
  "Lab:Mm vs Wild:As",
  "Lab-like:As vs Wild-like:As",
  "Lab:Mm vs Wild-like:As",
  "Lab:Mm vs Lab-like:As"
)

contrast_key <- tibble::tribble(
  ~contrast,                            ~first_group,   ~second_group,
  "Wild-like:As vs Wild:As",           "Wild-like:As", "Wild:As",
  "Lab-like:As vs Wild:As",            "Lab-like:As",  "Wild:As",
  "Lab:Mm vs Wild:As",                 "Lab:Mm",       "Wild:As",
  "Lab-like:As vs Wild-like:As",       "Lab-like:As",  "Wild-like:As",
  "Lab:Mm vs Wild-like:As",            "Lab:Mm",       "Wild-like:As",
  "Lab:Mm vs Lab-like:As",             "Lab:Mm",       "Lab-like:As"
)


# ==============================================================================
# 1. PREPARE SAMPLE METADATA
# ==============================================================================

meta <- data.frame(
  phyloseq::sample_data(ps_all),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta$sample <- phyloseq::sample_names(ps_all)

meta <- tibble::as_tibble(meta) |>
  dplyr::mutate(
    SampleOrigin = factor(SampleOrigin, levels = group_levels),
    mouse_ID_dataset = as.character(mouse_ID_dataset)
  )

stopifnot(!anyDuplicated(meta$sample))
stopifnot(!anyNA(meta$SampleOrigin))
stopifnot(!anyNA(meta$mouse_ID_dataset))

sample_counts <- meta |>
  dplyr::count(SampleOrigin, name = "n_samples")


# ==============================================================================
# 2. WEIGHTED NSTI QUALITY CONTROL
# ==============================================================================

nsti_file <- file.path(
  picrust_dir,
  "results",
  "EC_metagenome_out",
  "weighted_nsti.tsv.gz"
)

nsti <- readr::read_tsv(
  nsti_file,
  show_col_types = FALSE
)

required_nsti_columns <- c("sample", "weighted_NSTI")

if (!all(required_nsti_columns %in% names(nsti))) {
  stop(
    "weighted_nsti.tsv.gz does not contain the expected columns: ",
    paste(required_nsti_columns, collapse = ", ")
  )
}

if (!all(nsti$sample %in% meta$sample)) {
  stop("At least one NSTI sample does not match ps_all metadata.")
}

nsti_meta <- nsti |>
  dplyr::left_join(meta, by = "sample")

if (anyNA(nsti_meta$SampleOrigin)) {
  stop("At least one NSTI sample failed to receive SampleOrigin metadata.")
}

nsti_summary <- nsti_meta |>
  dplyr::group_by(SampleOrigin) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_NSTI = mean(weighted_NSTI, na.rm = TRUE),
    sd_NSTI = stats::sd(weighted_NSTI, na.rm = TRUE),
    median_NSTI = stats::median(weighted_NSTI, na.rm = TRUE),
    min_NSTI = min(weighted_NSTI, na.rm = TRUE),
    max_NSTI = max(weighted_NSTI, na.rm = TRUE),
    .groups = "drop"
  )

nsti_model <- lmerTest::lmer(
  weighted_NSTI ~ SampleOrigin + (1 | mouse_ID_dataset),
  data = nsti_meta
)

nsti_anova <- stats::anova(nsti_model)

nsti_emmeans <- emmeans::emmeans(
  nsti_model,
  pairwise ~ SampleOrigin,
  adjust = "tukey"
)

cat("\nOverall weighted NSTI summary:\n")
print(summary(nsti_meta$weighted_NSTI))

cat("\nWeighted NSTI by SampleOrigin:\n")
print(nsti_summary)

cat("\nNSTI mixed-model ANOVA:\n")
print(nsti_anova)

cat("\nNSTI estimated marginal means and Tukey contrasts:\n")
print(nsti_emmeans)

p_nsti <- ggplot2::ggplot(
  nsti_meta,
  ggplot2::aes(x = SampleOrigin, y = weighted_NSTI)
) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  ggplot2::labs(x = NULL, y = "Weighted NSTI") +
  ggplot2::theme_classic()

ggplot2::ggsave(
  file.path(figures_dir, "PICRUSt2_prediction_quality_NSTI.pdf"),
  p_nsti,
  width = 7,
  height = 5
)


# ==============================================================================
# 3. READ AND INSPECT METACYC PATHWAY OUTPUT
# ==============================================================================

pathway_file <- file.path(
  picrust_dir,
  "results",
  "pathways_out",
  "path_abun_unstrat.tsv.gz"
)

pathways <- readr::read_tsv(
  pathway_file,
  show_col_types = FALSE
)

if (!"pathway" %in% names(pathways)) {
  stop("MetaCyc pathway table does not contain a 'pathway' column.")
}

path_mat <- as.matrix(pathways[, -1, drop = FALSE])
storage.mode(path_mat) <- "numeric"
rownames(path_mat) <- pathways$pathway

if (!all(colnames(path_mat) %in% meta$sample)) {
  stop("At least one MetaCyc sample does not match ps_all metadata.")
}

cat("\nRaw MetaCyc pathway matrix dimensions:\n")
print(dim(path_mat))


# ------------------------------------------------------------------------------
# 3A. READ HUMAN-READABLE METACYC PATHWAY DESCRIPTIONS
# ------------------------------------------------------------------------------
#
# This file is the same PICRUSt2 MetaCyc abundance table, but with an additional
# 'description' column added using:
#
#   add_descriptions.py -m METACYC
#
# We use it only as an annotation dictionary. The numerical pathway abundance
# values used below still come from path_abun_unstrat.tsv.gz.

pathway_description_file <- file.path(
  picrust_dir,
  "results",
  "path_abun_unstrat_descriptions.tsv.gz"
)

if (!file.exists(pathway_description_file)) {
  stop(
    "Could not find the annotated MetaCyc pathway file:\n  ",
    pathway_description_file,
    "\nCopy path_abun_unstrat_descriptions.tsv.gz into PICRUSt2/results/."
  )
}

pathway_descriptions_raw <- readr::read_tsv(
  pathway_description_file,
  show_col_types = FALSE
)

required_description_columns <- c("pathway", "description")

if (!all(required_description_columns %in% names(pathway_descriptions_raw))) {
  stop(
    "Annotated MetaCyc file must contain columns named 'pathway' and ",
    "'description'. Found: ",
    paste(names(pathway_descriptions_raw)[1:min(5, ncol(pathway_descriptions_raw))],
          collapse = ", ")
  )
}

pathway_descriptions <- pathway_descriptions_raw |>
  dplyr::select(pathway, description) |>
  dplyr::distinct(pathway, .keep_all = TRUE) |>
  dplyr::mutate(
    pathway = as.character(pathway),
    description = dplyr::na_if(trimws(as.character(description)), "")
  )

description_coverage <- tibble::tibble(
  pathway = rownames(path_mat)
) |>
  dplyr::left_join(pathway_descriptions, by = "pathway") |>
  dplyr::summarise(
    n_pathways = dplyr::n(),
    n_with_description = sum(!is.na(description)),
    percent_with_description = 100 * mean(!is.na(description))
  )

cat("\nMetaCyc pathway-description coverage:\n")
print(description_coverage)


# ==============================================================================
# 4. PATHWAY PREVALENCE AND SAMPLE-TOTAL QC
# ==============================================================================

path_prev <- tibble::tibble(
  pathway = rownames(path_mat),
  n_present = rowSums(path_mat > 0),
  prevalence = rowMeans(path_mat > 0),
  total_abundance = rowSums(path_mat)
)

prevalence_threshold_summary <- tibble::tibble(
  threshold = c(">=1%", ">=5%", ">=10%", ">=25%", ">=50%", ">=75%", "100%"),
  n_pathways = c(
    sum(path_prev$prevalence >= 0.01),
    sum(path_prev$prevalence >= 0.05),
    sum(path_prev$prevalence >= 0.10),
    sum(path_prev$prevalence >= 0.25),
    sum(path_prev$prevalence >= 0.50),
    sum(path_prev$prevalence >= 0.75),
    sum(path_prev$prevalence == 1)
  )
)

cat("\nPathway prevalence distribution:\n")
print(summary(path_prev$prevalence))

cat("\nNumber of pathways at prevalence thresholds:\n")
print(prevalence_threshold_summary)

sample_path_totals <- colSums(path_mat)

path_total_df <- tibble::tibble(
  sample = names(sample_path_totals),
  pathway_total = as.numeric(sample_path_totals)
) |>
  dplyr::left_join(
    meta |>
      dplyr::select(sample, SampleOrigin, SeqDepth),
    by = "sample"
  )

depth_total_cor <- stats::cor.test(
  log10(path_total_df$pathway_total),
  log10(path_total_df$SeqDepth),
  method = "spearman",
  exact = FALSE
)

cat("\nRaw predicted pathway total vs sequencing depth:\n")
print(depth_total_cor)


# ==============================================================================
# 5. FILTER PATHWAYS AND CONVERT TO RELATIVE ABUNDANCE
# ==============================================================================

keep_pathways <- path_prev$pathway[
  path_prev$prevalence >= minimum_pathway_prevalence
]

path_mat_filt <- path_mat[
  rownames(path_mat) %in% keep_pathways,
  ,
  drop = FALSE
]

path_rel <- sweep(
  path_mat_filt,
  MARGIN = 2,
  STATS = colSums(path_mat_filt),
  FUN = "/"
)

if (!all(abs(colSums(path_rel) - 1) < 1e-10)) {
  stop("Relative pathway abundances do not sum to 1 for every sample.")
}

if (any(!is.finite(path_rel))) {
  stop("Non-finite values detected in relative pathway abundance matrix.")
}

fraction_retained <- colSums(path_mat_filt) / colSums(path_mat)

path_mean_rel <- tibble::tibble(
  pathway = rownames(path_rel),
  mean_relative_abundance = rowMeans(path_rel),
  median_relative_abundance = apply(
    path_rel,
    MARGIN = 1,
    FUN = stats::median
  ),
  prevalence = rowMeans(path_rel > 0)
)

cat("\nFiltered MetaCyc matrix dimensions:\n")
print(dim(path_mat_filt))

cat("\nPercent of predicted abundance retained after filtering:\n")
print(summary(fraction_retained * 100))

cat("\nMean relative abundance of retained pathways:\n")
print(summary(path_mean_rel$mean_relative_abundance))


# ==============================================================================
# 6. PREPARE ONE-SAMPLE-PER-MOUSE FUNCTIONAL BETA DIVERSITY
# ==============================================================================

path_rel_samples <- t(path_rel)

meta_path <- meta[
  match(rownames(path_rel_samples), meta$sample),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(path_rel_samples),
    meta_path$sample
  )
)

# SampleOrigin is constant within mouse, so do NOT use permutations restricted
# within mouse to test SampleOrigin. Instead select one independent sample/mouse.
#
# Colony/lab groups: sample closest to Timept = 0.
# Wild mice: earliest SampleDate because Timept is unavailable.
# Remaining ties: SampleDate then sample ID.
meta_path_ind <- meta_path |>
  dplyr::group_by(mouse_ID_dataset) |>
  dplyr::arrange(
    dplyr::if_else(
      SampleOrigin == "Wild:As",
      as.numeric(SampleDate - min(SampleDate, na.rm = TRUE)),
      abs(Timept)
    ),
    SampleDate,
    sample,
    .by_group = TRUE
  ) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

if (max(table(meta_path_ind$mouse_ID_dataset)) != 1) {
  stop("Independent beta-diversity dataset still contains repeated mice.")
}

path_rel_ind <- path_rel_samples[
  meta_path_ind$sample,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(path_rel_ind),
    meta_path_ind$sample
  )
)

independent_sample_counts <- meta_path_ind |>
  dplyr::count(SampleOrigin, name = "n_mice")

cat("\nIndependent-mouse sample counts:\n")
print(independent_sample_counts)


# ==============================================================================
# 7. FUNCTIONAL BETA DIVERSITY: BRAY-CURTIS, PERMANOVA, PERMDISP
# ==============================================================================

path_bray_ind <- vegan::vegdist(
  path_rel_ind,
  method = "bray"
)

set.seed(123)

path_permanova_ind <- vegan::adonis2(
  path_bray_ind ~ SampleOrigin,
  data = meta_path_ind,
  permutations = 9999
)

cat("\nOmnibus functional PERMANOVA:\n")
print(path_permanova_ind)

path_disp_ind <- vegan::betadisper(
  path_bray_ind,
  group = meta_path_ind$SampleOrigin
)

path_disp_anova <- stats::anova(path_disp_ind)

set.seed(123)

path_disp_perm_ind <- vegan::permutest(
  path_disp_ind,
  permutations = 9999
)

path_disp_tukey <- stats::TukeyHSD(path_disp_ind)

path_disp_df <- tibble::tibble(
  distance = path_disp_ind$distances,
  SampleOrigin = meta_path_ind$SampleOrigin
)

path_disp_summary <- path_disp_df |>
  dplyr::group_by(SampleOrigin) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_distance = mean(distance),
    sd_distance = stats::sd(distance),
    median_distance = stats::median(distance),
    .groups = "drop"
  )

cat("\nFunctional PERMDISP ANOVA:\n")
print(path_disp_anova)

cat("\nFunctional PERMDISP permutation test:\n")
print(path_disp_perm_ind)

cat("\nPairwise Tukey comparisons of dispersion:\n")
print(path_disp_tukey)


# ==============================================================================
# 8. PAIRWISE FUNCTIONAL PERMANOVA
# ==============================================================================

origin_pairs <- combn(
  group_levels,
  2,
  simplify = FALSE
)

pairwise_path_permanova <- lapply(
  origin_pairs,
  function(groups) {

    keep <- meta_path_ind$SampleOrigin %in% groups

    mat_sub <- path_rel_ind[
      keep,
      ,
      drop = FALSE
    ]

    meta_sub <- meta_path_ind[
      keep,
      ,
      drop = FALSE
    ]

    meta_sub$SampleOrigin <- droplevels(meta_sub$SampleOrigin)

    bray_sub <- vegan::vegdist(
      mat_sub,
      method = "bray"
    )

    set.seed(123)

    fit <- vegan::adonis2(
      bray_sub ~ SampleOrigin,
      data = meta_sub,
      permutations = 9999
    )

    tibble::tibble(
      group1 = groups[1],
      group2 = groups[2],
      n1 = sum(meta_sub$SampleOrigin == groups[1]),
      n2 = sum(meta_sub$SampleOrigin == groups[2]),
      F = fit$F[1],
      R2 = fit$R2[1],
      p = fit$`Pr(>F)`[1]
    )
  }
) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    p_adj = stats::p.adjust(
      p,
      method = "BH"
    )
  )

readr::write_csv(
  pairwise_path_permanova,
  file.path(tables_dir, "MetaCyc_pairwise_PERMANOVA_results.csv")
)

cat("\nPairwise functional PERMANOVAs:\n")
print(pairwise_path_permanova)


# ==============================================================================
# 9. PCoA AND DISPERSION FIGURES
# ==============================================================================

path_pcoa_ind <- stats::cmdscale(
  path_bray_ind,
  eig = TRUE,
  k = 2
)

positive_eig_ind <- path_pcoa_ind$eig[
  path_pcoa_ind$eig > 0
]

pcoa1_pct_ind <- 100 *
  path_pcoa_ind$eig[1] /
  sum(positive_eig_ind)

pcoa2_pct_ind <- 100 *
  path_pcoa_ind$eig[2] /
  sum(positive_eig_ind)

path_pcoa_ind_df <- tibble::tibble(
  sample = rownames(path_pcoa_ind$points),
  PCoA1 = path_pcoa_ind$points[, 1],
  PCoA2 = path_pcoa_ind$points[, 2]
) |>
  dplyr::left_join(
    meta_path_ind |>
      dplyr::select(sample, SampleOrigin, mouse_ID_dataset),
    by = "sample"
  )

p_functional_pcoa <- ggplot2::ggplot(
  path_pcoa_ind_df,
  ggplot2::aes(
    x = PCoA1,
    y = PCoA2,
    colour = SampleOrigin
  )
) +
  ggplot2::geom_point(size = 2.2, alpha = 0.75) +
  ggplot2::stat_ellipse(linewidth = 0.8) +
  ggplot2::labs(
    x = paste0("PCoA1 (", round(pcoa1_pct_ind, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa2_pct_ind, 1), "%)"),
    colour = "Sample origin"
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  file.path(figures_dir, "MetaCyc_pathway_composition_PCoA.pdf"),
  p_functional_pcoa,
  width = 7,
  height = 5.5
)

p_functional_dispersion <- ggplot2::ggplot(
  path_disp_df,
  ggplot2::aes(
    x = SampleOrigin,
    y = distance,
    colour = SampleOrigin
  )
) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  ggplot2::geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.2) +
  ggplot2::labs(
    x = NULL,
    y = "Distance to group centroid",
    colour = "Sample origin"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  file.path(figures_dir, "MetaCyc_pathway_dispersion_by_sample_origin.pdf"),
  p_functional_dispersion,
  width = 7,
  height = 5.5
)


# ==============================================================================
# 10. SEQUENCING-DEPTH DIAGNOSTIC FOR FUNCTIONAL COMPOSITION
# ==============================================================================

meta_path_ind <- meta_path_ind |>
  dplyr::mutate(
    log10_reads = log10(SeqDepth)
  )

set.seed(123)

path_depth_test_ind <- vegan::adonis2(
  path_bray_ind ~ log10_reads,
  data = meta_path_ind,
  permutations = 9999
)

set.seed(123)

path_depth_adjusted_ind <- vegan::adonis2(
  path_bray_ind ~ SampleOrigin + log10_reads,
  data = meta_path_ind,
  permutations = 9999,
  by = "margin"
)

cat("\nSequencing-depth diagnostic — unadjusted:\n")
print(path_depth_test_ind)

cat("\nSequencing-depth diagnostic — adjusted for SampleOrigin:\n")
print(path_depth_adjusted_ind)


# ==============================================================================
# 11. PREPARE METACYC DATA FOR MAASLIN3
# ==============================================================================

# Differential pathway analysis uses all eligible samples.
# Repeated samples are handled by the mouse random intercept.
#
# MaAsLin3 expects samples in rows and features in columns.
# We pass the already filtered, relative-abundance MetaCyc table.

pathway_feature_data <- as.data.frame(
  t(path_rel),
  check.names = FALSE
)

pathway_metadata <- meta[
  match(
    rownames(pathway_feature_data),
    meta$sample
  ),
  ,
  drop = FALSE
] |>
  dplyr::mutate(
    SampleOrigin = factor(SampleOrigin, levels = group_levels),
    Sex = factor(Sex),
    mouse_ID_dataset = factor(mouse_ID_dataset),
    log10_reads = log10(SeqDepth)
  )

stopifnot(
  identical(
    rownames(pathway_feature_data),
    pathway_metadata$sample
  )
)

# Exclude samples missing any model covariate.
keep_model_samples <- stats::complete.cases(
  pathway_metadata[
    c(
      "SampleOrigin",
      "Sex",
      "log10_reads",
      "mouse_ID_dataset"
    )
  ]
)

pathway_metadata <- pathway_metadata[
  keep_model_samples,
  ,
  drop = FALSE
]

pathway_feature_data <- pathway_feature_data[
  pathway_metadata$sample,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(pathway_feature_data),
    pathway_metadata$sample
  )
)

# MaAsLin3 requires sample IDs as metadata row names.
metadata <- pathway_metadata |>
  dplyr::select(-sample) |>
  as.data.frame(check.names = FALSE)

rownames(metadata) <- pathway_metadata$sample

metadata <- metadata[
  rownames(pathway_feature_data),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(pathway_feature_data),
    rownames(metadata)
  )
)

cat("\nMaAsLin3 feature-table dimensions:\n")
print(dim(pathway_feature_data))

cat("\nMaAsLin3 metadata dimensions:\n")
print(dim(metadata))

cat("\nSamples entering MaAsLin3 by SampleOrigin:\n")
print(
  metadata |>
    tibble::rownames_to_column("sample") |>
    dplyr::count(SampleOrigin)
)


# ==============================================================================
# 12. MAASLIN3 DIFFERENTIAL METACYC PATHWAY ANALYSIS
# ==============================================================================

primary_formula <- stats::as.formula(
  "~ SampleOrigin + Sex + log10_reads + (1 | mouse_ID_dataset)"
)

run_abundance_analysis <- function(
  feature_data,
  metadata,
  model_dir,
  feature_level
) {

  results_file <- file.path(model_dir, "all_results.tsv")

  if (force_rerun && dir.exists(model_dir)) {
    unlink(model_dir, recursive = TRUE, force = TRUE)
  }

  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(results_file) && !force_rerun) {
    stop(
      "Existing ", feature_level, " results found in ", model_dir,
      ". Set force_rerun <- TRUE to replace them."
    )
  }

  stopifnot(
    identical(
      rownames(feature_data),
      rownames(metadata)
    )
  )

  detected_cores <- parallel::detectCores(logical = TRUE)
  n_workers <- max(1L, min(6L, detected_cores - 2L))

  try(mirai::daemons(0), silent = TRUE)
  mirai::daemons(n_workers)

  fit <- tryCatch(
    {
      maaslin3::maaslin3(
        input_data = feature_data,
        input_metadata = metadata,
        output = model_dir,
        formula = primary_formula,

        # Relative abundances and prevalence filtering were done upstream.
        normalization = "NONE",
        transform = "LOG",
        min_abundance = 0,
        min_prevalence = 0,

        median_comparison_abundance = TRUE,
        standardize = TRUE,

        evaluate_only = "abundance",
        warn_prevalence = FALSE,

        max_significance = q_threshold,
        plot_summary_plot = FALSE,
        plot_associations = FALSE,
        max_pngs = 0,
        save_models = TRUE,
        save_plots_rds = FALSE
      )
    },
    finally = {
      mirai::daemons(0)
    }
  )

  model_results <- readr::read_tsv(
    results_file,
    show_col_types = FALSE
  )

  origin_key <- model_results |>
    dplyr::filter(
      metadata == "SampleOrigin",
      !is.na(value),
      !is.na(name)
    ) |>
    dplyr::distinct(value, name)

  coef_names <- c(
    origin_key |>
      dplyr::filter(value == "Wild-like:As") |>
      dplyr::pull(name) |>
      dplyr::first(),

    origin_key |>
      dplyr::filter(value == "Lab-like:As") |>
      dplyr::pull(name) |>
      dplyr::first(),

    origin_key |>
      dplyr::filter(value == "Lab:Mm") |>
      dplyr::pull(name) |>
      dplyr::first()
  )

  if (anyNA(coef_names)) {
    stop(
      "Could not identify all SampleOrigin coefficients for ",
      feature_level
    )
  }

  contrast_matrix <- rbind(
    "Wild-like:As vs Wild:As"      = c( 1,  0,  0),
    "Lab-like:As vs Wild:As"       = c( 0,  1,  0),
    "Lab:Mm vs Wild:As"            = c( 0,  0,  1),
    "Lab-like:As vs Wild-like:As"  = c(-1,  1,  0),
    "Lab:Mm vs Wild-like:As"       = c(-1,  0,  1),
    "Lab:Mm vs Lab-like:As"        = c( 0, -1,  1)
  )

  colnames(contrast_matrix) <- coef_names

  contrast_fit <- maaslin3::maaslin_contrast_test(
    maaslin3_fit = fit,
    contrast_mat = contrast_matrix,
    evaluate_only = "abundance"
  )

  if (
    is.list(contrast_fit) &&
    !is.null(contrast_fit$fit_data_abundance) &&
    !is.null(contrast_fit$fit_data_abundance$results)
  ) {
    pairwise <- tibble::as_tibble(
      contrast_fit$fit_data_abundance$results
    )
  } else if (is.data.frame(contrast_fit)) {
    pairwise <- tibble::as_tibble(contrast_fit)
  } else {
    stop(
      "Unexpected maaslin_contrast_test() result for ",
      feature_level
    )
  }

  if ("test" %in% names(pairwise)) {
    pairwise <- pairwise |>
      dplyr::rename(contrast = test)
  }

  pairwise <- pairwise |>
    dplyr::mutate(
      feature = as.character(feature),
      feature_level = feature_level,
      contrast = factor(contrast, levels = contrast_names),
      significant =
        !is.na(qval_individual) &
        qval_individual < q_threshold,
      lower_ci = coef - 1.96 * stderr,
      upper_ci = coef + 1.96 * stderr,
      fold_difference = 2^coef
    ) |>
    dplyr::left_join(contrast_key, by = "contrast") |>
    dplyr::mutate(
      direction = dplyr::case_when(
        coef > 0 ~ paste0("Higher in ", first_group),
        coef < 0 ~ paste0("Higher in ", second_group),
        TRUE ~ "No difference"
      )
    ) |>
    dplyr::arrange(contrast, qval_individual)

  list(
    fit = fit,
    pairwise = pairwise,
    significant = pairwise |>
      dplyr::filter(significant)
  )
}

message("\nRunning MetaCyc pathway abundance analysis...")

pathway_results <- run_abundance_analysis(
  feature_data = pathway_feature_data,
  metadata = metadata,
  model_dir = maaslin_dir,
  feature_level = "MetaCyc pathway"
)


# ==============================================================================
# 13. ANNOTATE, SUMMARIZE, AND SAVE DIFFERENTIAL PATHWAY RESULTS
# ==============================================================================

pathway_pairwise <- pathway_results$pairwise |>
  dplyr::left_join(
    pathway_descriptions,
    by = c("feature" = "pathway")
  ) |>
  dplyr::mutate(
    # Use the readable MetaCyc description wherever possible.
    pathway_name = dplyr::if_else(
      !is.na(description),
      description,
      feature
    )
  )

# If two MetaCyc IDs share the same text description, append the pathway ID so
# that the heatmap row labels remain unique.
duplicate_names <- pathway_pairwise |>
  dplyr::distinct(feature, pathway_name) |>
  dplyr::count(pathway_name, name = "n_ids") |>
  dplyr::filter(n_ids > 1) |>
  dplyr::pull(pathway_name)

pathway_pairwise <- pathway_pairwise |>
  dplyr::mutate(
    pathway_label = dplyr::if_else(
      pathway_name %in% duplicate_names,
      paste0(pathway_name, " [", feature, "]"),
      pathway_name
    )
  )

pathway_significant <- pathway_pairwise |>
  dplyr::filter(significant)

pathway_totals <- pathway_significant |>
  dplyr::count(
    contrast,
    name = "n_significant_pathways"
  ) |>
  tidyr::complete(
    contrast = factor(
      contrast_names,
      levels = contrast_names
    ),
    fill = list(
      n_significant_pathways = 0
    )
  )

pathway_directions <- pathway_significant |>
  dplyr::count(
    contrast,
    direction,
    name = "n_pathways"
  )

pathway_recurrence <- pathway_significant |>
  dplyr::group_by(feature, description, pathway_name, pathway_label) |>
  dplyr::summarise(
    n_significant_contrasts = dplyr::n(),
    min_q = min(qval_individual, na.rm = TRUE),
    max_abs_coef = max(abs(coef), na.rm = TRUE),
    significant_contrasts = paste(
      as.character(contrast),
      collapse = "; "
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(n_significant_contrasts),
    min_q,
    dplyr::desc(max_abs_coef)
  )

cat("\nSignificant predicted MetaCyc pathways by comparison:\n")
print(pathway_totals)

cat(
  "\nUnique significant predicted pathways across all contrasts: ",
  dplyr::n_distinct(pathway_significant$feature),
  "\n",
  sep = ""
)

cat("\nMost recurrent significant pathways:\n")
print(pathway_recurrence, n = min(50, nrow(pathway_recurrence)))


# ==============================================================================
# 14. FIGURE 1 — NUMBER OF DIFFERENTIAL PATHWAYS BY COMPARISON
# ==============================================================================

p_pathway_counts <- ggplot2::ggplot(
  pathway_totals,
  ggplot2::aes(
    x = contrast,
    y = n_significant_pathways
  )
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_significant_pathways),
    hjust = -0.2,
    size = 3.5
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    title = "Differentially abundant predicted MetaCyc pathways",
    subtitle = "MaAsLin3 abundance model; individual FDR q < 0.05",
    x = NULL,
    y = "Number of significant pathways"
  ) +
  ggplot2::theme_classic(base_size = 11)

ggplot2::ggsave(
  file.path(figures_dir, "MetaCyc_significant_pathway_counts.pdf"),
  p_pathway_counts,
  width = 8,
  height = 5.5
)



# ==============================================================================
# 15. ASSIGN SIGNIFICANT METACYC PATHWAYS TO BROAD FUNCTIONAL CATEGORIES
# ==============================================================================
#
# PURPOSE
# -------
# The pathway-level MaAsLin3 analysis can identify many significant MetaCyc
# pathways. MetaCyc pathways overlap biologically and are not statistically
# independent, so these categories are used for ORGANIZATION and INTERPRETATION,
# not as a replacement for the pathway-level statistical tests.
#
# Classification is based on human-readable MetaCyc descriptions and broad
# biological themes relevant to gut microbial metabolism.
#
# Some unusual MetaCyc labels (for example, photosynthesis-specific pathway
# names) can arise because PICRUSt2 predicts gene-family content and reconstructs
# pathways from shared enzyme sets. These labels should not be interpreted
# literally as evidence that the named organism-specific process is occurring
# in the mouse gut.
# ==============================================================================

pathway_category_table <- pathway_recurrence |>
  dplyr::mutate(
    description_lower = stringr::str_to_lower(description),

    functional_category = dplyr::case_when(

      # Fermentation and short-chain-fatty-acid-related metabolism
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "fermentation", "butanoate", "butyrate", "propanoate",
            "propionate", "acetate and butanoate", "acetate and lactate",
            "butanediol", "propanediol"
          ),
          collapse = "|"
        )
      ) ~ "Fermentation / SCFA metabolism",

      # Carbohydrate use and central sugar pathways
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "glucose", "fructose", "galactose", "galactitol", "galactarate",
            "glucarate", "mannan", "mannitol", "mannose", "xylose",
            "arabinose", "rhamnose", "fucose", "starch", "glycogen",
            "sucrose", "lactose", "hexitol", "inositol", "stachyose",
            "fructan", "chitin", "erythronate", "threonate", "ascorbate",
            "glycolysis", "gluconeogenesis", "pentose phosphate",
            "entner-doudoroff"
          ),
          collapse = "|"
        )
      ) ~ "Carbohydrate utilization",

      # Amino-acid biosynthesis and degradation
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "lysine", "arginine", "ornithine", "citrulline", "glutamate",
            "glutamine", "aspartate", "asparagine", "methionine", "cysteine",
            "serine", "glycine", "histidine", "tryptophan", "tyrosine",
            "phenylalanine", "threonine", "isoleucine", "leucine", "valine",
            "alanine", "amino acid", "chorismate", "cyanophycin"
          ),
          collapse = "|"
        )
      ) ~ "Amino-acid metabolism",

      # Nucleotide / nucleoside metabolism
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "purine", "pyrimidine", "adenosine", "adenine", "guanosine",
            "guanine", "inosine", "uridine", "ump", "nucleotide",
            "nucleoside", "deoxyribonucleotide", "queuosine", "preq0"
          ),
          collapse = "|"
        )
      ) ~ "Nucleotide metabolism",

      # Vitamins and cofactors
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "thiamine", "pyridox", "folate", "tetrahydrofolate",
            "riboflavin", "flavin", "biotin", "cobalamin", "coenzyme",
            "nad ", "nad salvage", "heme", "menaquin", "ubiquin",
            "molybdopterin", "pantothen", "tetrapyrrole"
          ),
          collapse = "|"
        )
      ) ~ "Vitamin / cofactor metabolism",

      # Lipids, membrane lipids and isoprenoids
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "fatty acid", "stearate", "palmit", "phospholipid", "phosphatid",
            "lipid", "vaccenate", "dodecenoate", "gondoate", "plasmalogen",
            "carnitine", "phytol", "farnesol", "isoprene", "isoprenoid",
            "polyisoprenoid", "methylerythritol", "mevalonate"
          ),
          collapse = "|"
        )
      ) ~ "Lipid / isoprenoid metabolism",

      # Bacterial cell envelope and surface structures
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "peptidoglycan", "lipid a", "o-antigen", "cell wall",
            "muropeptide", "lipopolysaccharide", "kdo", "colanic acid",
            "common antigen", "n-acetylmuram", "udp-n-acetyl",
            "cmp-pseudaminate", "cmp-legionaminate", "glycosylation",
            "heptose"
          ),
          collapse = "|"
        )
      ) ~ "Cell envelope / surface structures",

      # Central carbon and energy metabolism
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "tca cycle", "tricarboxylic", "glyoxylate", "acetyl coenzyme a",
            "acetyl-coa", "electron transport", "respiration", "rubisco",
            "calvin", "photosynthetic carbon assimilation",
            "formaldehyde assimilation", "2-methylcitrate"
          ),
          collapse = "|"
        )
      ) ~ "Energy / central carbon metabolism",

      # Nitrogen, sulfur and phosphorus metabolism
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "nitrate", "denitrification", "sulfate", "sulfur", "phosphonate",
            "phytate", "urea cycle", "ureide", "allantoin"
          ),
          collapse = "|"
        )
      ) ~ "Nitrogen / sulfur / phosphorus metabolism",

      # Polyamines
      stringr::str_detect(
        description_lower,
        paste(
          c("polyamine", "putrescine", "spermidine", "norspermidine"),
          collapse = "|"
        )
      ) ~ "Polyamine metabolism",

      # Aromatic / xenobiotic degradation
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "catechol", "protocatechuate", "toluene", "benzoate",
            "hydroxybenzoate", "octane oxidation", "glyphosate",
            "pesticide", "androstenedione", "aromatic degradation"
          ),
          collapse = "|"
        )
      ) ~ "Aromatic / xenobiotic degradation",

      # Amines and neuroactive compounds
      stringr::str_detect(
        description_lower,
        paste(
          c(
            "dopamine", "serotonin", "gaba", "aminobutanoate",
            "ethanolamine", "creatinine"
          ),
          collapse = "|"
        )
      ) ~ "Amine / neuroactive compound metabolism",

      # RNA / translation / stringent-response-related pathways
      stringr::str_detect(
        description_lower,
        paste(c("trna", "ppgpp", "stringent response"), collapse = "|")
      ) ~ "Translation / regulatory metabolism",

      TRUE ~ "Other / unclear"
    )
  ) |>
  dplyr::select(-description_lower)


# Summarize the number of UNIQUE significant pathways assigned to each category.
category_summary <- pathway_category_table |>
  dplyr::count(
    functional_category,
    name = "n_unique_significant_pathways",
    sort = TRUE
  )

cat("\nSignificant pathways by broad functional category:\n")
print(category_summary, n = Inf)

# Keep a separate table of pathways that remain difficult to classify.
unclear_pathways <- pathway_category_table |>
  dplyr::filter(functional_category == "Other / unclear") |>
  dplyr::select(
    feature,
    description,
    n_significant_contrasts,
    min_q,
    max_abs_coef
  )

# ==============================================================================
# 16. ADD FUNCTIONAL CATEGORIES TO ALL PAIRWISE MAASLIN3 RESULTS
# ==============================================================================
#
# The category annotation above was built from pathways significant in at least
# one contrast. We now attach those categories to the complete pairwise result
# table. Pathways never significant in any contrast are labelled "Not selected /
# not significant" because they were not part of the biological categorization.
# ==============================================================================

pathway_pairwise_categorized <- pathway_pairwise |>
  dplyr::select(
    -dplyr::any_of("functional_category")
  ) |>
  dplyr::left_join(
    pathway_category_table |>
      dplyr::select(feature, functional_category),
    by = "feature"
  ) |>
  dplyr::mutate(
    functional_category = dplyr::coalesce(
      functional_category,
      "Not selected / not significant"
    )
  )

readr::write_csv(
  pathway_pairwise_categorized,
  file.path(tables_dir, "MetaCyc_differential_abundance_all_pairwise_results.csv")
)


# Category-by-contrast summary among FDR-significant pathways.
#
# IMPORTANT: This is descriptive. A pathway can overlap biologically with other
# pathways, so counts within categories are not independent statistical tests.

category_contrast_summary <- pathway_pairwise_categorized |>
  dplyr::filter(significant) |>
  dplyr::count(
    contrast,
    functional_category,
    name = "n_significant_pathways"
  ) |>
  dplyr::arrange(contrast, dplyr::desc(n_significant_pathways))

readr::write_csv(
  category_contrast_summary,
  file.path(tables_dir, "MetaCyc_functional_category_counts_by_comparison.csv")
)


# ==============================================================================
# 17. FIGURE 2 — DESCRIPTIVE FUNCTIONAL-CATEGORY HEATMAP
# ==============================================================================
#
# This figure answers:
#   "Which broad kinds of predicted functions contain the most significant
#    pathways in each comparison?"
#
# It is a descriptive organization of pathway-level results, NOT a category-
# level enrichment test.
# ==============================================================================

category_plot_order <- c(
  "Carbohydrate utilization",
  "Fermentation / SCFA metabolism",
  "Amino-acid metabolism",
  "Nucleotide metabolism",
  "Vitamin / cofactor metabolism",
  "Lipid / isoprenoid metabolism",
  "Cell envelope / surface structures",
  "Energy / central carbon metabolism",
  "Nitrogen / sulfur / phosphorus metabolism",
  "Polyamine metabolism",
  "Aromatic / xenobiotic degradation",
  "Amine / neuroactive compound metabolism",
  "Translation / regulatory metabolism",
  "Other / unclear"
)

category_contrast_plot <- category_contrast_summary |>
  tidyr::complete(
    contrast = factor(contrast_names, levels = contrast_names),
    functional_category = category_plot_order,
    fill = list(n_significant_pathways = 0)
  ) |>
  dplyr::mutate(
    contrast = factor(contrast, levels = contrast_names),
    functional_category = factor(
      functional_category,
      levels = rev(category_plot_order)
    )
  )

p_category_counts <- ggplot2::ggplot(
  category_contrast_plot,
  ggplot2::aes(
    x = contrast,
    y = functional_category,
    fill = n_significant_pathways
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_significant_pathways),
    size = 2.8
  ) +
  ggplot2::labs(
    title = "Functional categories represented among differential pathways",
    subtitle = paste0(
      "Number of FDR-significant MetaCyc pathways per category; ",
      "descriptive summary only"
    ),
    x = NULL,
    y = NULL,
    fill = "Significant\npathways"
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

ggplot2::ggsave(
  file.path(figures_dir, "MetaCyc_functional_category_summary_heatmap.pdf"),
  p_category_counts,
  width = 10.5,
  height = 7.5
)


# ==============================================================================
# 18. MANUALLY CURATED REPRESENTATIVE PATHWAYS FOR MAIN FUNCTIONAL HEATMAP
# ==============================================================================
#
# WHY MANUAL CURATION?
# --------------------
# Ranking pathways only by q-value or effect size tends to prioritize unusual
# database-specific pathway labels and can produce a biologically incoherent
# figure. Here we instead show four interpretable representative pathways from
# each of eight major functional themes.
#
# This selection is ONLY for visualization. Statistical inference remains based
# on the complete pathway-level MaAsLin3 analysis above.
# ==============================================================================

selected_functional_pathways <- tibble::tribble(
  ~functional_category, ~feature,

  # Carbohydrate utilization
  "Carbohydrate utilization", "P562-PWY",
  "Carbohydrate utilization", "PWY-7456",
  "Carbohydrate utilization", "PWY-6531",
  "Carbohydrate utilization", "RHAMCAT-PWY",

  # Fermentation / SCFA metabolism
  "Fermentation / SCFA metabolism", "PWY-5677",
  "Fermentation / SCFA metabolism", "PWY-5484",
  "Fermentation / SCFA metabolism", "GOLPDLCAT-PWY",
  "Fermentation / SCFA metabolism", "P163-PWY",

  # Amino-acid metabolism
  "Amino-acid metabolism", "CITRULBIO-PWY",
  "Amino-acid metabolism", "PWY-2942",
  "Amino-acid metabolism", "HISDEG-PWY",
  "Amino-acid metabolism", "HOMOSER-METSYN-PWY",

  # Nucleotide metabolism
  "Nucleotide metabolism", "PWY-7220",
  "Nucleotide metabolism", "PWY-7221",
  "Nucleotide metabolism", "DENOVOPURINE2-PWY",
  "Nucleotide metabolism", "PWY-6126",

  # Vitamin / cofactor metabolism
  "Vitamin / cofactor metabolism", "NAD-BIOSYNTHESIS-II",
  "Vitamin / cofactor metabolism", "PYRIDOXSYN-PWY",
  "Vitamin / cofactor metabolism", "HEMESYN2-PWY",
  "Vitamin / cofactor metabolism", "PWY-7357",

  # Lipid / isoprenoid metabolism
  "Lipid / isoprenoid metabolism", "PWY-5989",
  "Lipid / isoprenoid metabolism", "POLYISOPRENSYN-PWY",
  "Lipid / isoprenoid metabolism", "PWY-7783",
  "Lipid / isoprenoid metabolism", "NAGLIPASYN-PWY",

  # Cell envelope / surface structures
  "Cell envelope / surface structures", "PWY-5265",
  "Cell envelope / surface structures", "UDPNAGSYN-PWY",
  "Cell envelope / surface structures", "PWY0-1241",
  "Cell envelope / surface structures", "PWY-6143",

  # Energy / central carbon metabolism
  "Energy / central carbon metabolism", "PWY-5913",
  "Energy / central carbon metabolism", "P23-PWY",
  "Energy / central carbon metabolism", "PWY-6969",
  "Energy / central carbon metabolism", "REDCITCYC"
)


# Confirm every manually selected pathway exists in the PICRUSt2 result set.
missing_selected_pathways <- selected_functional_pathways |>
  dplyr::anti_join(
    pathway_category_table,
    by = "feature"
  )

if (nrow(missing_selected_pathways) > 0) {
  print(missing_selected_pathways)
  stop(
    "At least one manually selected MetaCyc pathway ID was not found. ",
    "See missing_selected_pathways above."
  )
}


selected_functional_pathways_named <- selected_functional_pathways |>
  dplyr::left_join(
    pathway_category_table |>
      dplyr::select(
        feature,
        description,
        pathway_label,
        n_significant_contrasts,
        min_q,
        max_abs_coef
      ),
    by = "feature"
  )

readr::write_csv(
  selected_functional_pathways_named,
  file.path(tables_dir, "MetaCyc_curated_pathways_shown_in_heatmap.csv")
)


# ==============================================================================
# 19. BUILD DATA FOR THE CURATED FUNCTIONAL HEATMAP
# ==============================================================================

# Remove annotation columns already present in pathway_pairwise before joining,
# preventing duplicate description.x / description.y columns.

pathway_pairwise_for_heatmap <- pathway_pairwise |>
  dplyr::select(
    -dplyr::any_of(
      c(
        "description",
        "functional_category",
        "pathway_label",
        "pathway_name"
      )
    )
  )

functional_heatmap_data <- pathway_pairwise_for_heatmap |>
  dplyr::inner_join(
    selected_functional_pathways_named |>
      dplyr::select(
        feature,
        functional_category,
        description,
        pathway_label
      ),
    by = "feature"
  ) |>
  dplyr::mutate(
    significant_heatmap =
      !is.na(qval_individual) &
      qval_individual < q_threshold,

    # Non-significant cells are displayed in grey.
    plot_coef = dplyr::if_else(
      significant_heatmap,
      as.numeric(coef),
      NA_real_
    ),

    significance_symbol = dplyr::case_when(
      !significant_heatmap ~ "",
      qval_individual < 0.001 ~ "***",
      qval_individual < 0.01  ~ "**",
      qval_individual < 0.05  ~ "*",
      TRUE ~ ""
    )
  )


functional_category_order <- c(
  "Carbohydrate utilization",
  "Fermentation / SCFA metabolism",
  "Amino-acid metabolism",
  "Nucleotide metabolism",
  "Vitamin / cofactor metabolism",
  "Lipid / isoprenoid metabolism",
  "Cell envelope / surface structures",
  "Energy / central carbon metabolism"
)

functional_heatmap_data <- functional_heatmap_data |>
  dplyr::mutate(
    functional_category = factor(
      functional_category,
      levels = functional_category_order
    ),
    contrast = factor(
      contrast,
      levels = contrast_names
    )
  )


# QC: there should be 32 unique pathways: 4 in each of 8 categories.
heatmap_pathway_qc <- functional_heatmap_data |>
  dplyr::distinct(
    functional_category,
    feature,
    description
  ) |>
  dplyr::count(
    functional_category,
    name = "n_pathways"
  )

cat("\nCurated pathways per functional category:\n")
print(heatmap_pathway_qc, n = Inf)

if (
  dplyr::n_distinct(functional_heatmap_data$feature) != 32L ||
  any(heatmap_pathway_qc$n_pathways != 4L)
) {
  stop(
    "Curated heatmap QC failed: expected 32 pathways total and 4 per category."
  )
}


# ==============================================================================
# 20. FIGURE 3 — CURATED FUNCTIONAL EFFECT-SIZE HEATMAP
# ==============================================================================
#
# Positive MaAsLin3 coefficient:
#   higher predicted pathway abundance in the FIRST-NAMED group.
#
# Negative MaAsLin3 coefficient:
#   higher predicted pathway abundance in the SECOND-NAMED group.
#
# Grey:
#   pathway contrast is not significant at q < 0.05.
# ==============================================================================

coef_limit <- max(
  abs(functional_heatmap_data$plot_coef),
  na.rm = TRUE
)

if (!is.finite(coef_limit) || coef_limit == 0) {
  coef_limit <- 1
}

p_curated_functional_heatmap <- ggplot2::ggplot(
  functional_heatmap_data,
  ggplot2::aes(
    x = contrast,
    y = description,
    fill = plot_coef
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.30
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance_symbol),
    size = 2.7
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#B7433F",
    midpoint = 0,
    limits = c(-coef_limit, coef_limit),
    na.value = "grey92",
    name = "MaAsLin3\ncoefficient"
  ) +
  ggplot2::facet_grid(
    functional_category ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::labs(
    title = "Predicted functional pathways associated with microbiome origin",
    subtitle = paste0(
      "Representative MetaCyc pathways grouped by biological function; ",
      "grey cells are not significant at FDR q < 0.05"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Positive coefficients = higher predicted abundance in the first-named ",
      "group; negative coefficients = higher in the second-named group. ",
      "* q < 0.05, ** q < 0.01, *** q < 0.001."
    )
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.y = ggplot2::element_text(size = 7.5),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.y.left = ggplot2::element_text(
      angle = 0,
      face = "bold",
      hjust = 1,
      size = 8
    ),
    panel.spacing.y = grid::unit(0.55, "lines"),
    legend.position = "right",
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 7.5
    )
  )

ggplot2::ggsave(
  file.path(figures_dir, "MetaCyc_curated_differential_pathway_heatmap.pdf"),
  p_curated_functional_heatmap,
  width = 12,
  height = 15
)


# ==============================================================================
# 21. CATEGORY-LEVEL DIRECTION SUMMARY
# ==============================================================================
#
# For each contrast and functional category, summarize how many significant
# pathways are higher in the first-named versus second-named group.
#
# Again, this is DESCRIPTIVE. It should not be interpreted as a formal
# enrichment test because MetaCyc pathways are overlapping/non-independent.
# ==============================================================================

category_direction_summary <- pathway_pairwise_categorized |>
  dplyr::filter(
    significant,
    functional_category != "Not selected / not significant"
  ) |>
  dplyr::mutate(
    direction_simple = dplyr::case_when(
      coef > 0 ~ "Higher in first group",
      coef < 0 ~ "Higher in second group",
      TRUE ~ "No direction"
    )
  ) |>
  dplyr::count(
    contrast,
    functional_category,
    direction_simple,
    name = "n_pathways"
  ) |>
  tidyr::pivot_wider(
    names_from = direction_simple,
    values_from = n_pathways,
    values_fill = 0
  )

# Make sure both expected columns exist even if a future dataset has zero
# pathways in one direction.
if (!"Higher in first group" %in% names(category_direction_summary)) {
  category_direction_summary[["Higher in first group"]] <- 0L
}

if (!"Higher in second group" %in% names(category_direction_summary)) {
  category_direction_summary[["Higher in second group"]] <- 0L
}

category_direction_summary <- category_direction_summary |>
  dplyr::mutate(
    n_total =
      `Higher in first group` +
      `Higher in second group`,

    proportion_higher_first = dplyr::if_else(
      n_total > 0,
      `Higher in first group` / n_total,
      NA_real_
    ),

    direction_balance = dplyr::if_else(
      n_total > 0,
      (
        `Higher in first group` -
        `Higher in second group`
      ) / n_total,
      NA_real_
    )
  )

# ==============================================================================
# 22. SAVE A COMPACT BIOLOGICAL-INTERPRETATION TABLE
# ==============================================================================
#
# This is the most convenient table to open when writing the Results section.
# It contains only significant pathway-level contrasts, human-readable names,
# functional categories, effect sizes, confidence intervals, FDR q-values and
# direction.
# ==============================================================================

interpretation_table <- pathway_pairwise_categorized |>
  dplyr::filter(significant) |>
  dplyr::select(
    feature,
    description,
    pathway_label,
    functional_category,
    contrast,
    first_group,
    second_group,
    coef,
    lower_ci,
    upper_ci,
    fold_difference,
    qval_individual,
    direction
  ) |>
  dplyr::arrange(
    functional_category,
    contrast,
    qval_individual
  )

readr::write_csv(
  interpretation_table,
  file.path(
    tables_dir,
    "MetaCyc_differential_abundance_significant_results.csv"
  )
)


# ==============================================================================
# 23. FINAL CONSOLE SUMMARY
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("PICRUSt2 FUNCTIONAL ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat(
  "\nMetaCyc pathways before prevalence filtering: ",
  nrow(path_mat),
  "\n",
  sep = ""
)

cat(
  "MetaCyc pathways retained at >= ",
  minimum_pathway_prevalence * 100,
  "% prevalence: ",
  nrow(path_rel),
  "\n",
  sep = ""
)

cat(
  "Samples entering functional beta diversity: ",
  nrow(meta_path_ind),
  " (one sample per mouse)\n",
  sep = ""
)

cat(
  "Samples entering MaAsLin3: ",
  nrow(metadata),
  "\n",
  sep = ""
)

cat(
  "Unique significant MetaCyc pathways: ",
  dplyr::n_distinct(pathway_significant$feature),
  "\n",
  sep = ""
)

cat(
  "Significant pathways assigned to 'Other / unclear': ",
  nrow(unclear_pathways),
  "\n",
  sep = ""
)

cat(
  "Curated pathways shown in main functional heatmap: ",
  dplyr::n_distinct(functional_heatmap_data$feature),
  "\n",
  sep = ""
)

cat("\nMain outputs:\n")
cat("  Tables:  ", tables_dir, "\n", sep = "")
cat("  Figures: ", figures_dir, "\n", sep = "")
cat("  MaAsLin3 model: ", maaslin_dir, "\n", sep = "")

cat("\nKey figures generated:\n")
cat("  - PICRUSt2_prediction_quality_NSTI.pdf\n")
cat("  - MetaCyc_pathway_composition_PCoA.pdf\n")
cat("  - MetaCyc_pathway_dispersion_by_sample_origin.pdf\n")
cat("  - MetaCyc_significant_pathway_counts.pdf\n")
cat("  - MetaCyc_functional_category_summary_heatmap.pdf\n")
cat("  - MetaCyc_curated_differential_pathway_heatmap.pdf\n")

cat("\nKey interpretation table:\n")
cat("  - MetaCyc_differential_abundance_significant_results.csv\n")

cat("\nInterpretation reminder:\n")
cat(
  "PICRUSt2 estimates predicted functional potential from 16S data; ",
  "these are not directly measured metagenomic functions or pathway activities.\n"
)

cat(
  "Functional categories are descriptive organizational labels; ",
  "the inferential statistics remain the pathway-level MaAsLin3 tests.\n"
)

cat("============================================================\n")

# ==============================================================================
# MAASLIN 3 DIFFERENTIAL ABUNDANCE — ASV + GENUS + FAMILY ANALYSES
# ==============================================================================
# Biological question:
#   How do microbial abundances differ among the four SampleOrigin groups?
#
# Analysis strategy:
#   PRIMARY STATISTICAL RESOLUTION = ASV level
#     - Preserves the full sequence-level resolution of the 16S dataset.
#
#   MAIN-TEXT TAXONOMIC PRESENTATION = genus level
#     - Raw ASV counts are summed within named genera before modelling.
#     - Only genuinely genus-assigned ASVs enter the genus analysis.
#     - The main heatmap shows the 20 strongest/recurrent significant genera.
#
#   TAXONOMIC-SCALE SENSITIVITY = family level
#     - Raw ASV counts are summed within named families before modelling.
#     - Provides broader taxonomic coverage where genus assignment is incomplete.
#
#   The script also quantifies the fraction of ASVs AND total reads represented
#   by valid genus and family assignments.
#
# For ALL taxonomic levels:
#   - all eligible samples are used;
#   - repeated samples are handled with (1 | mouse_ID_dataset);
#   - Sex and log10 sequencing depth are included as fixed effects;
#   - SeqPlate is not included because it is confounded with SampleOrigin;
#   - only MaAsLin3's abundance component is fitted;
#   - six pairwise SampleOrigin contrasts are obtained from one fitted model;
#   - significance is individual FDR q < 0.05.
#
# Required object:
#   ps_all : phyloseq object created by 00_DataPrep_Baseline.R
# ==============================================================================

# Load the final preprocessed baseline phyloseq object.
# Run Scripts/00_DataPrep_Baseline.R first whenever preprocessing changes.
ps_all <- readRDS("Data/Processed/baseline_microbiome_phyloseq.rds")

# Packages ---------------------------------------------------------------------
required_packages <- c(
  "phyloseq", "maaslin3", "dplyr", "tidyr", "tibble",
  "readr", "ggplot2", "mirai", "openxlsx", "cowplot"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE,
          FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop("Install the following packages before running: ",
       paste(missing_packages, collapse = ", "))
}

# User settings ----------------------------------------------------------------
output_dir <- "Results/03_ASV_Differential_Abundance"
tables_dir <- file.path(output_dir, "Tables")
figures_dir <- file.path(output_dir, "Figures")

# MaAsLin3 writes several intermediate files to its output directory. These are
# required during model fitting/contrast testing but are not manuscript outputs,
# so place them in the session temporary directory rather than Results/.
model_dir <- file.path(
  tempdir(),
  "Baseline_Differential_Abundance_MaAsLin3"
)
asv_model_dir <- file.path(model_dir, "ASV")
genus_model_dir <- file.path(model_dir, "Genus")
family_model_dir <- file.path(model_dir, "Family")

q_threshold <- 0.05
minimum_prevalence <- 0.10
force_rerun <- FALSE

# Main-text genus heatmap: intentionally compact for publication.
# Ranking = number of significant contrasts, then minimum q-value, then effect size.
max_genera_heatmap <- 20

for (x in c(output_dir, tables_dir, figures_dir, model_dir,
            asv_model_dir, genus_model_dir, family_model_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}


# Output policy:
#   - save ONE multi-sheet workbook with the key differential-abundance results;
#   - save ONE complete integrated figure;
#   - do not save component panels or large intermediate MaAsLin3 tables;
#   - detailed model objects/results remain available in memory during the run.

# Group order ------------------------------------------------------------------
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
  ~contrast, ~first_group, ~second_group,
  "Wild-like:As vs Wild:As",     "Wild-like:As", "Wild:As",
  "Lab-like:As vs Wild:As",      "Lab-like:As",  "Wild:As",
  "Lab:Mm vs Wild:As",           "Lab:Mm",       "Wild:As",
  "Lab-like:As vs Wild-like:As", "Lab-like:As",  "Wild-like:As",
  "Lab:Mm vs Wild-like:As",      "Lab:Mm",       "Wild-like:As",
  "Lab:Mm vs Lab-like:As",       "Lab:Mm",       "Lab-like:As"
)


# Display-only comparison orientation -------------------------------------------
# Statistical contrasts above are retained exactly as fitted. For figures, all
# comparisons are reoriented along the same ecological gradient:
#
#   Wild -> Wild-like -> Lab-like -> Lab
#
# Every fitted contrast is currently written with the more laboratory-associated
# cohort first. Therefore all displayed heatmap coefficients are multiplied by
# -1. This is a DISPLAY-ONLY transformation:
#
#   negative displayed coefficient = higher in the first (wilder) cohort
#   positive displayed coefficient = higher in the second (more lab-associated)
#
# Model coefficients, q-values, significance calls, and exported statistical
# result tables remain unchanged.

display_contrast_key <- tibble::tribble(
  ~contrast,                           ~display_comparison,       ~coef_multiplier,
  "Wild-like:As vs Wild:As",           "Wild vs Wild-like",       -1,
  "Lab-like:As vs Wild:As",            "Wild vs Lab-like",        -1,
  "Lab:Mm vs Wild:As",                 "Wild vs Lab",             -1,
  "Lab-like:As vs Wild-like:As",       "Wild-like vs Lab-like",   -1,
  "Lab:Mm vs Wild-like:As",            "Wild-like vs Lab",        -1,
  "Lab:Mm vs Lab-like:As",             "Lab-like vs Lab",         -1
)

display_comparison_order <- c(
  "Wild vs Wild-like",
  "Wild vs Lab-like",
  "Wild vs Lab",
  "Wild-like vs Lab-like",
  "Wild-like vs Lab",
  "Lab-like vs Lab"
)

display_comparison_labels <- c(
  "Wild vs Wild-like"     = "Wild vs\nWild-like",
  "Wild vs Lab-like"      = "Wild vs\nLab-like",
  "Wild vs Lab"           = "Wild vs\nLab",
  "Wild-like vs Lab-like" = "Wild-like vs\nLab-like",
  "Wild-like vs Lab"      = "Wild-like vs\nLab",
  "Lab-like vs Lab"       = "Lab-like vs\nLab"
)

# Extract counts and taxonomy --------------------------------------------------
otu_matrix <- as(phyloseq::otu_table(ps_all), "matrix")
if (phyloseq::taxa_are_rows(ps_all)) otu_matrix <- t(otu_matrix)
storage.mode(otu_matrix) <- "numeric"

asv_counts <- as.data.frame(otu_matrix, check.names = FALSE)

if (anyNA(asv_counts) || any(as.matrix(asv_counts) < 0)) {
  stop("ASV table contains NA or negative values.")
}

taxonomy <- data.frame(
  phyloseq::tax_table(ps_all),
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  tibble::rownames_to_column("ASV")

if (!all(c("Genus", "Family") %in% names(taxonomy))) {
  stop("The taxonomy table must contain both Genus and Family columns.")
}

# Prepare metadata -------------------------------------------------------------
metadata <- data.frame(
  phyloseq::sample_data(ps_all),
  check.names = FALSE
) |>
  tibble::rownames_to_column("sample_id") |>
  dplyr::mutate(
    SampleOrigin = factor(SampleOrigin, levels = group_levels),
    Sex = factor(Sex),
    mouse_ID_dataset = factor(mouse_ID_dataset)
  ) |>
  dplyr::filter(
    !is.na(SampleOrigin),
    !is.na(Sex),
    !is.na(mouse_ID_dataset)
  ) |>
  tibble::column_to_rownames("sample_id")

shared_samples <- intersect(rownames(asv_counts), rownames(metadata))
asv_counts <- asv_counts[shared_samples, , drop = FALSE]
metadata <- metadata[shared_samples, , drop = FALSE]
stopifnot(identical(rownames(asv_counts), rownames(metadata)))

# Remove empty samples and add sequencing depth --------------------------------
# Sequencing depth is always calculated from the COMPLETE ASV count table.
sample_reads <- rowSums(asv_counts)
keep_nonempty <- sample_reads > 0

asv_counts <- asv_counts[keep_nonempty, , drop = FALSE]
metadata <- metadata[keep_nonempty, , drop = FALSE]
sample_reads <- sample_reads[keep_nonempty]
metadata$log10_reads <- log10(sample_reads)

# Compact design summaries -----------------------------------------------------
sample_summary <- metadata |>
  tibble::rownames_to_column("sample_id") |>
  dplyr::count(SampleOrigin, name = "n_samples")

repeat_summary <- metadata |>
  tibble::rownames_to_column("sample_id") |>
  dplyr::count(mouse_ID_dataset, name = "n_samples") |>
  dplyr::summarise(
    n_mice = dplyr::n(),
    n_mice_repeated = sum(n_samples > 1),
    mean_samples_per_mouse = mean(n_samples),
    median_samples_per_mouse = median(n_samples),
    maximum_samples_per_mouse = max(n_samples)
  )

# TSS denominator from complete community -------------------------------------
# We normalize BEFORE feature filtering/aggregation so every feature abundance
# remains relative to the full observed community, not to the subset retained
# for a particular taxonomic analysis.
asv_tss_all <- sweep(
  as.matrix(asv_counts),
  MARGIN = 1,
  STATS = sample_reads,
  FUN = "/"
)

# ASV input --------------------------------------------------------------------
asv_prevalence <- colMeans(asv_counts > 0)
asv_keep <- asv_prevalence >= minimum_prevalence

asv_filter_summary <- tibble::tibble(
  ASV = colnames(asv_counts),
  prevalence = as.numeric(asv_prevalence),
  retained = as.logical(asv_keep)
) |>
  dplyr::left_join(taxonomy, by = "ASV")

asv_input <- as.data.frame(
  asv_tss_all[, asv_keep, drop = FALSE],
  check.names = FALSE
)

# Genus and family inputs -------------------------------------------------------
# Helper defining an informative taxonomic assignment. Generic placeholders are
# excluded because pooling unrelated ASVs under labels such as "uncultured" or
# "unclassified" would create biologically invalid aggregate features.
is_valid_taxon <- function(x) {
  x <- trimws(as.character(x))
  !(is.na(x) |
      x == "" |
      toupper(x) == "NA" |
      tolower(x) %in% c(
        "uncultured",
        "unclassified",
        "unknown",
        "uncultured bacterium",
        "unidentified"
      ))
}

# Aggregate raw ASV counts to one taxonomic rank, then divide by sequencing
# depth from the COMPLETE ASV community. This means genus/family abundances are
# still expressed relative to the whole observed community, not merely the
# classified subset.
make_taxon_input <- function(rank_name) {

  taxon_raw <- trimws(as.character(taxonomy[[rank_name]]))
  valid_taxon <- is_valid_taxon(taxon_raw)

  taxon_for_asv <- taxon_raw[
    match(colnames(asv_counts), taxonomy$ASV)
  ]
  valid_for_asv <- valid_taxon[
    match(colnames(asv_counts), taxonomy$ASV)
  ]

  if (anyNA(valid_for_asv)) {
    stop("Some ASVs in the count table are missing from the taxonomy table.")
  }

  classified_counts <- as.matrix(
    asv_counts[, valid_for_asv, drop = FALSE]
  )
  classified_taxa <- taxon_for_asv[valid_for_asv]

  aggregated_counts <- t(rowsum(
    t(classified_counts),
    group = classified_taxa,
    reorder = FALSE
  ))

  taxon_prevalence <- colMeans(aggregated_counts > 0)
  taxon_keep <- taxon_prevalence >= minimum_prevalence

  taxon_tss <- sweep(
    aggregated_counts,
    MARGIN = 1,
    STATS = sample_reads,
    FUN = "/"
  )

  taxon_input <- as.data.frame(
    taxon_tss[, taxon_keep, drop = FALSE],
    check.names = FALSE
  )

  filter_summary <- tibble::tibble(
    taxon = colnames(aggregated_counts),
    prevalence = as.numeric(taxon_prevalence),
    retained = as.logical(taxon_keep)
  )
  names(filter_summary)[1] <- rank_name

  # Coverage metrics. Read coverage is calculated from raw reads across all
  # samples, so it indicates what proportion of the observed community is
  # represented at this taxonomic rank.
  read_coverage <- sum(as.matrix(asv_counts[, valid_for_asv, drop = FALSE])) /
    sum(as.matrix(asv_counts))

  retained_asv_rank_coverage <- mean(valid_for_asv[asv_keep])

  list(
    input = taxon_input,
    filter_summary = filter_summary,
    valid_for_asv = valid_for_asv,
    n_taxa_before_filter = ncol(aggregated_counts),
    n_taxa_after_filter = ncol(taxon_input),
    read_coverage = read_coverage,
    retained_asv_rank_coverage = retained_asv_rank_coverage
  )
}

genus_data <- make_taxon_input("Genus")
family_data <- make_taxon_input("Family")

genus_input <- genus_data$input
family_input <- family_data$input

# ASV-to-taxonomy mapping retained for traceability.
taxonomy_mapping <- taxonomy |>
  dplyr::transmute(
    ASV,
    Genus = as.character(Genus),
    Family = as.character(Family),
    valid_genus = is_valid_taxon(Genus),
    valid_family = is_valid_taxon(Family)
  )

# Taxonomic coverage summary ---------------------------------------------------
# This table is especially important because the proportion of ASVs assigned to
# a rank can look low when many rare ASVs are unresolved. Read coverage shows
# how much of the actual sequenced community those classified ASVs represent.
taxonomic_resolution_summary <- tibble::tibble(
  taxonomic_level = c("Genus", "Family"),
  n_total_ASVs = ncol(asv_counts),
  n_ASVs_with_valid_assignment = c(
    sum(genus_data$valid_for_asv),
    sum(family_data$valid_for_asv)
  ),
  proportion_all_ASVs_with_valid_assignment = c(
    mean(genus_data$valid_for_asv),
    mean(family_data$valid_for_asv)
  ),
  proportion_retained_ASVs_with_valid_assignment = c(
    genus_data$retained_asv_rank_coverage,
    family_data$retained_asv_rank_coverage
  ),
  proportion_total_reads_with_valid_assignment = c(
    genus_data$read_coverage,
    family_data$read_coverage
  ),
  n_taxa_before_prevalence_filter = c(
    genus_data$n_taxa_before_filter,
    family_data$n_taxa_before_filter
  ),
  n_taxa_after_prevalence_filter = c(
    genus_data$n_taxa_after_filter,
    family_data$n_taxa_after_filter
  )
)


# Shared mixed-model formula ---------------------------------------------------
primary_formula <- stats::as.formula(
  "~ SampleOrigin + Sex + log10_reads + (1 | mouse_ID_dataset)"
)

# Helper: run one abundance-only MaAsLin3 analysis -----------------------------
run_abundance_analysis <- function(feature_data, model_dir, feature_level) {

  results_file <- file.path(model_dir, "all_results.tsv")

  if (force_rerun && dir.exists(model_dir)) {
    unlink(model_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  # Always refit unless a completed output is deliberately being reused.
  # Contrast tests need the in-memory fit object, so a fresh run is the safest
  # and clearest workflow for this compact script.
  if (file.exists(results_file) && !force_rerun) {
    stop(
      "Existing ", feature_level, " results found in ", model_dir,
      ". Set force_rerun <- TRUE to replace them."
    )
  }

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

        # Relative abundances and prevalence filtering were done explicitly.
        normalization = "NONE",
        transform = "LOG",
        min_abundance = 0,
        min_prevalence = 0,

        median_comparison_abundance = TRUE,
        standardize = TRUE,

        # Abundance inference only: no logistic prevalence component.
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

  model_results <- readr::read_tsv(results_file, show_col_types = FALSE)

  # Identify treatment-coded SampleOrigin coefficients relative to Wild:As.
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
    stop("Could not identify all SampleOrigin coefficients for ", feature_level)
  }

  # One fitted model -> six formal linear contrasts.
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

  if (is.list(contrast_fit) &&
      !is.null(contrast_fit$fit_data_abundance) &&
      !is.null(contrast_fit$fit_data_abundance$results)) {
    pairwise <- tibble::as_tibble(contrast_fit$fit_data_abundance$results)
  } else if (is.data.frame(contrast_fit)) {
    pairwise <- tibble::as_tibble(contrast_fit)
  } else {
    stop("Unexpected maaslin_contrast_test() result for ", feature_level)
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
      significant = !is.na(qval_individual) &
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
    significant = pairwise |> dplyr::filter(significant)
  )
}

# Run ASV analysis -------------------------------------------------------------
message("\nRunning ASV-level abundance analysis...")
asv_analysis <- run_abundance_analysis(
  feature_data = asv_input,
  model_dir = asv_model_dir,
  feature_level = "ASV"
)

# Add taxonomy to ASV result tables.
asv_pairwise <- asv_analysis$pairwise |>
  dplyr::rename(ASV = feature) |>
  dplyr::left_join(taxonomy, by = "ASV")

asv_significant <- asv_pairwise |>
  dplyr::filter(significant)

# Run genus analysis -----------------------------------------------------------
message("\nRunning genus-level abundance analysis...")
genus_analysis <- run_abundance_analysis(
  feature_data = genus_input,
  model_dir = genus_model_dir,
  feature_level = "Genus"
)

genus_pairwise <- genus_analysis$pairwise |>
  dplyr::rename(Genus = feature)

genus_significant <- genus_pairwise |>
  dplyr::filter(significant)


# Run family analysis ----------------------------------------------------------
message("\nRunning family-level abundance analysis...")
family_analysis <- run_abundance_analysis(
  feature_data = family_input,
  model_dir = family_model_dir,
  feature_level = "Family"
)

family_pairwise <- family_analysis$pairwise |>
  dplyr::rename(Family = feature)

family_significant <- family_pairwise |>
  dplyr::filter(significant)

# Compact summaries ------------------------------------------------------------
make_contrast_totals <- function(results) {
  results |>
    dplyr::filter(significant) |>
    dplyr::count(contrast, name = "n_significant_features") |>
    tidyr::complete(
      contrast = factor(contrast_names, levels = contrast_names),
      fill = list(n_significant_features = 0)
    )
}

make_direction_summary <- function(results) {
  results |>
    dplyr::filter(significant) |>
    dplyr::count(contrast, direction, name = "n_features")
}

make_recurrence <- function(results, feature_column) {
  results |>
    dplyr::filter(significant) |>
    dplyr::group_by(.data[[feature_column]]) |>
    dplyr::summarise(
      n_significant_contrasts = dplyr::n(),
      min_q = min(qval_individual, na.rm = TRUE),
      max_abs_coef = max(abs(coef), na.rm = TRUE),
      significant_contrasts = paste(as.character(contrast), collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_significant_contrasts), min_q)
}

asv_totals <- make_contrast_totals(asv_pairwise)
asv_directions <- make_direction_summary(asv_pairwise)
asv_recurrence <- make_recurrence(asv_pairwise, "ASV")

genus_totals <- make_contrast_totals(genus_pairwise)
genus_directions <- make_direction_summary(genus_pairwise)
genus_recurrence <- make_recurrence(genus_pairwise, "Genus")

family_totals <- make_contrast_totals(family_pairwise)
family_directions <- make_direction_summary(family_pairwise)
family_recurrence <- make_recurrence(family_pairwise, "Family")

# Reader-facing result tables are written once, at the end of the script, into
# a single multi-sheet workbook. Full model results remain in memory here.

# INTERNAL PLOT OBJECT: number of significant genera ---------------------------
p_genus_counts <- ggplot2::ggplot(
  genus_totals,
  ggplot2::aes(x = contrast, y = n_significant_features)
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_significant_features),
    hjust = -0.2,
    size = 3.5
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    title = "Differentially abundant genera among sample origins",
    subtitle = "MaAsLin3 abundance model; individual FDR q < 0.05",
    x = NULL,
    y = "Number of significant genera"
  ) +
  ggplot2::theme_classic(base_size = 11)


# MAIN-FIGURE COMPONENT: top 20 genus coefficient heatmap ----------------------
# Every row is the same taxonomic rank (Genus). Genera are ranked objectively by
# recurrence across significant contrasts, then minimum q-value, then effect size.
selected_genera <- genus_recurrence |>
  # Exclude taxonomically uninformative placeholder groups from the main figure.
  # The next highest-ranked informative genus is then included automatically so
  # the heatmap still contains max_genera_heatmap genera.
  dplyr::filter(
    !is.na(Genus),
    trimws(Genus) != "",
    tolower(trimws(Genus)) != "incertae sedis"
  ) |>
  dplyr::slice_head(
    n = max_genera_heatmap
  ) |>
  dplyr::pull(Genus)

cat("\nGenera selected for main heatmap after excluding 'Incertae Sedis':\n")
print(selected_genera)

genus_heatmap_data <- genus_pairwise |>
  dplyr::filter(Genus %in% selected_genera) |>
  dplyr::mutate(
    contrast = as.character(contrast)
  ) |>
  dplyr::left_join(
    display_contrast_key,
    by = "contrast"
  ) |>
  dplyr::mutate(
    significance_marker = dplyr::case_when(
      qval_individual < 0.001 ~ "***",
      qval_individual < 0.01  ~ "**",
      qval_individual < 0.05  ~ "*",
      TRUE ~ ""
    ),
    coef_display = coef * coef_multiplier,
    Genus = factor(Genus, levels = rev(selected_genera)),
    display_comparison = factor(
      display_comparison,
      levels = display_comparison_order
    )
  )

p_genus_heatmap <- ggplot2::ggplot(
  genus_heatmap_data,
  ggplot2::aes(
    x = display_comparison,
    y = Genus,
    fill = coef_display
  )
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance_marker),
    size = 3
  ) +
  ggplot2::scale_fill_gradient2(
    midpoint = 0,
    name = "Coefficient"
  ) +
  ggplot2::scale_x_discrete(
    limits = display_comparison_order,
    labels = display_comparison_labels,
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Pairwise genus-level abundance coefficients",
    subtitle = "Top 20 recurrent/strong genera; * q<0.05, ** q<0.01, *** q<0.001",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 8),
    panel.grid = ggplot2::element_blank()
  )

genus_heatmap_height <- max(6, 0.28 * length(selected_genera) + 2)


# INTERNAL QC PLOT: complete significant-genus heatmap -------------------------
all_significant_genera <- genus_recurrence$Genus

full_genus_heatmap_data <- genus_pairwise |>
  dplyr::filter(Genus %in% all_significant_genera) |>
  dplyr::mutate(
    significance_marker = dplyr::case_when(
      qval_individual < 0.001 ~ "***",
      qval_individual < 0.01  ~ "**",
      qval_individual < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Genus = factor(Genus, levels = rev(all_significant_genera)),
    contrast = factor(contrast, levels = contrast_names)
  )

p_full_genus_heatmap <- ggplot2::ggplot(
  full_genus_heatmap_data,
  ggplot2::aes(x = contrast, y = Genus, fill = coef)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
  ggplot2::geom_text(ggplot2::aes(label = significance_marker), size = 2.7) +
  ggplot2::scale_fill_gradient2(midpoint = 0, name = "Coefficient") +
  ggplot2::labs(
    title = "Supplementary genus-level abundance coefficients",
    subtitle = "All genera significant in at least one comparison",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 7),
    panel.grid = ggplot2::element_blank()
  )


# INTERNAL QC PLOT: family coefficient heatmap ---------------------------------
# Family aggregation gives broader taxonomic coverage and tests whether the same
# origin effects persist at a coarser taxonomic scale.
all_significant_families <- family_recurrence$Family

if (length(all_significant_families) > 0) {
  family_heatmap_data <- family_pairwise |>
    dplyr::filter(Family %in% all_significant_families) |>
    dplyr::mutate(
      significance_marker = dplyr::case_when(
        qval_individual < 0.001 ~ "***",
        qval_individual < 0.01  ~ "**",
        qval_individual < 0.05  ~ "*",
        TRUE ~ ""
      ),
      Family = factor(Family, levels = rev(all_significant_families)),
      contrast = factor(contrast, levels = contrast_names)
    )

  p_family_heatmap <- ggplot2::ggplot(
    family_heatmap_data,
    ggplot2::aes(x = contrast, y = Family, fill = coef)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = significance_marker), size = 2.8) +
    ggplot2::scale_fill_gradient2(midpoint = 0, name = "Coefficient") +
    ggplot2::labs(
      title = "Supplementary family-level abundance coefficients",
      subtitle = "All families significant in at least one comparison",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 7.5),
      panel.grid = ggplot2::element_blank()
    )

}

# INTERNAL QC PLOT: complete significant-ASV heatmap ---------------------------
# ASVs are still tested individually here. Taxonomy is retained in the table,
# but the figure label uses Genus where available and ASV ID for traceability.
# ASVs lacking a genus assignment remain in this supplementary figure as ASV ID.
asv_recurrence_order <- asv_recurrence$ASV

asv_label_key <- taxonomy |>
  dplyr::transmute(
    ASV,
    Genus = trimws(as.character(Genus)),
    display_label = dplyr::if_else(
      !is.na(Genus) & Genus != "" & toupper(Genus) != "NA",
      paste0(Genus, " [", ASV, "]"),
      ASV
    )
  )

asv_heatmap_data <- asv_pairwise |>
  dplyr::filter(ASV %in% asv_recurrence_order) |>
  dplyr::left_join(
    asv_label_key |> dplyr::select(ASV, display_label),
    by = "ASV"
  ) |>
  dplyr::mutate(
    significance_marker = dplyr::case_when(
      qval_individual < 0.001 ~ "***",
      qval_individual < 0.01  ~ "**",
      qval_individual < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

ordered_asv_labels <- asv_label_key$display_label[
  match(asv_recurrence_order, asv_label_key$ASV)
]

asv_heatmap_data <- asv_heatmap_data |>
  dplyr::mutate(
    display_label = factor(display_label, levels = rev(ordered_asv_labels)),
    contrast = factor(contrast, levels = contrast_names)
  )

p_asv_heatmap <- ggplot2::ggplot(
  asv_heatmap_data,
  ggplot2::aes(x = contrast, y = display_label, fill = coef)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance_marker),
    size = 2.6
  ) +
  ggplot2::scale_fill_gradient2(midpoint = 0, name = "Coefficient") +
  ggplot2::labs(
    title = "Supplementary ASV-level abundance coefficients",
    subtitle = "All ASVs significant in at least one comparison",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 6.5),
    panel.grid = ggplot2::element_blank()
  )

asv_heatmap_height <- max(7, 0.22 * length(asv_recurrence_order) + 2)


# Console summary --------------------------------------------------------------
cat("\n============================================================\n")
cat("MAASLIN3 ASV + GENUS + FAMILY ABUNDANCE SUMMARY\n")
cat("============================================================\n")
cat("Samples analysed: ", nrow(asv_input), "\n", sep = "")
cat("ASVs retained for ASV model: ", ncol(asv_input), "\n", sep = "")
cat("Named genera retained for genus model: ", ncol(genus_input), "\n", sep = "")
cat("Named families retained for family model: ", ncol(family_input), "\n", sep = "")
cat("\nTaxonomic coverage:\n")
print(taxonomic_resolution_summary)
cat("\nASV-level results (supplementary):\n")
cat("  Unique significant ASVs: ",
    dplyr::n_distinct(asv_significant$ASV), "\n", sep = "")
cat("  Significant ASV-contrast associations: ",
    nrow(asv_significant), "\n", sep = "")
cat("\nGenus-level results (main presentation):\n")
cat("  Unique significant genera: ",
    dplyr::n_distinct(genus_significant$Genus), "\n", sep = "")
cat("  Significant genus-contrast associations: ",
    nrow(genus_significant), "\n", sep = "")
cat("\nGenus counts by contrast:\n")
print(genus_totals)
cat("\nFamily-level results (taxonomic-scale sensitivity):\n")
cat("  Unique significant families: ",
    dplyr::n_distinct(family_significant$Family), "\n", sep = "")
cat("  Significant family-contrast associations: ",
    nrow(family_significant), "\n", sep = "")
cat("\nFamily counts by contrast:\n")
print(family_totals)
cat("\nOutputs written to: ", output_dir, "\n", sep = "")
cat("============================================================\n")


# ============================================================
# INTEGRATED DIFFERENTIAL ABUNDANCE FIGURE
#
# Top:
#   Top-20 genus-level differential abundance heatmap
#
# Bottom:
#   Number of significant ASVs / genera / families
#
# The six comparison columns are panel-aligned exactly between plots and displayed from
# wilder to more laboratory-associated cohorts. Heatmap coefficients are
# reoriented for display so negative values indicate higher abundance in the
# first (wilder) cohort and positive values indicate higher abundance in the
# second (more lab-associated) cohort.
# Comparison labels appear only once, below the barplot.
# Legends are placed horizontally below the full figure.
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(cowplot)
library(grid)


# ------------------------------------------------------------
# 1. Define comparison display order
# ------------------------------------------------------------
# Use the same wild-to-lab gradient as the other manuscript figures.

comparison_order <- display_comparison_order
comparison_labels <- display_comparison_labels


# ------------------------------------------------------------
# 2. Number of significant taxa
# ------------------------------------------------------------
# Generate counts directly from the current model results rather than hard-coding
# values, so this panel always reflects the latest ASV/genus/family analyses.

count_data <- dplyr::bind_rows(
  asv_totals |>
    dplyr::mutate(taxonomic_level = "ASV"),
  genus_totals |>
    dplyr::mutate(taxonomic_level = "Genus"),
  family_totals |>
    dplyr::mutate(taxonomic_level = "Family")
) |>
  dplyr::mutate(
    contrast = as.character(contrast)
  ) |>
  dplyr::left_join(
    display_contrast_key |>
      dplyr::select(
        contrast,
        display_comparison
      ),
    by = "contrast"
  ) |>
  dplyr::transmute(
    comparison = factor(
      display_comparison,
      levels = comparison_order
    ),
    taxonomic_level = factor(
      taxonomic_level,
      levels = c(
        "ASV",
        "Genus",
        "Family"
      )
    ),
    n_significant = n_significant_features
  )


# ------------------------------------------------------------
# 3. Modify existing top-20 genus heatmap
# ------------------------------------------------------------
#
# p_genus_heatmap is your existing top-20 genus heatmap.
# These changes apply ONLY to the integrated figure and do
# not alter the original standalone p_genus_heatmap object.
# ------------------------------------------------------------

heatmap_integrated <- p_genus_heatmap +
  
  # Force columns into exactly the desired comparison order
  ggplot2::scale_x_discrete(
    limits = comparison_order,
    labels = comparison_labels,
    drop = FALSE
  ) +
  
  # Remove title/subtitle and x-axis title from this version
  ggplot2::labs(
    title = NULL,
    subtitle = NULL,
    x = NULL
  ) +
  
  ggplot2::theme(
    
    # Comparison labels appear only below barplot
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    
    # Make genus names more readable
    axis.text.y = ggplot2::element_text(
      size = 11
    ),
    
    # Thin black separator immediately beneath heatmap
    axis.line.x.bottom = ggplot2::element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    # Small space between heatmap separator and barplot
    plot.margin = ggplot2::margin(
      t = 3,
      r = 3,
      b = 4,
      l = 3
    )
  )


# ------------------------------------------------------------
# 4. Create compact barplot
# ------------------------------------------------------------

barplot_integrated <- ggplot2::ggplot(
  count_data,
  ggplot2::aes(
    x = comparison,
    y = n_significant,
    fill = taxonomic_level
  )
) +
  
  # Bars
  ggplot2::geom_col(
    position = ggplot2::position_dodge(
      width = 0.78
    ),
    width = 0.70
  ) +
  
  # Number of significant taxa above each bar
  ggplot2::geom_text(
    ggplot2::aes(
      label = n_significant
    ),
    position = ggplot2::position_dodge(
      width = 0.78
    ),
    vjust = -0.3,
    size = 3.5
  ) +
  
  # Same comparison order as heatmap
  ggplot2::scale_x_discrete(
    limits = comparison_order,
    labels = comparison_labels,
    drop = FALSE
  ) +
  
  ggplot2::scale_y_continuous(
    limits = c(0, 36),
    breaks = c(
      0,
      10,
      20,
      30
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0)
    )
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Significant taxa",
    fill = "Taxonomic level"
  ) +
  
  ggplot2::theme_classic(
    base_size = 11
  ) +
  
  # Add right-hand vertical border to match plotting area
  ggplot2::annotate(
    "segment",
    x = Inf,
    xend = Inf,
    y = -Inf,
    yend = Inf,
    colour = "black",
    linewidth = 0.4
  ) +
  
  # Allow border at edge of plotting area to render
  ggplot2::coord_cartesian(
    clip = "off"
  ) +
  
  ggplot2::theme(
    
    # Larger comparison labels
    axis.text.x = ggplot2::element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 1,
      size = 11,
      lineheight = 0.9
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 10
    ),
    
    axis.title.y = ggplot2::element_text(
      size = 11
    ),
    
    # Keep normal left and bottom axes
    axis.line.y.left = ggplot2::element_line(
      colour = "black",
      linewidth = 0.4
    ),
    
    axis.line.x.bottom = ggplot2::element_line(
      colour = "black",
      linewidth = 0.4
    ),
    
    # No top border:
    # this keeps the two plots visually distinct
    axis.line.x.top = ggplot2::element_blank(),
    
    # Legend formatting
    legend.position = "bottom",
    legend.direction = "horizontal",
    
    legend.title = ggplot2::element_text(
      size = 10
    ),
    
    legend.text = ggplot2::element_text(
      size = 10
    ),
    
    # Small top margin prevents barplot y-axis from
    # visually connecting to heatmap y-axis
    plot.margin = ggplot2::margin(
      t = 5,
      r = 3,
      b = 3,
      l = 3
    )
  )


# ------------------------------------------------------------
# 5. Extract legends for placement LEFT of barplot
# ------------------------------------------------------------

# Horizontal coefficient legend
heatmap_legend <- cowplot::get_legend(
  heatmap_integrated +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title = "Coefficient\n(- = wilder,\n+ = lab-associated)",
        title.position = "top",
        title.hjust = 0,
        barwidth = grid::unit(2.2, "cm"),
        barheight = grid::unit(0.35, "cm")
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.box.just = "left",
      legend.title = ggplot2::element_text(
        size = 9,
        lineheight = 0.95
      ),
      legend.text = ggplot2::element_text(
        size = 8.5
      ),
      legend.margin = ggplot2::margin(
        t = 0,
        r = 0,
        b = 0,
        l = 4
      ),
      legend.box.margin = ggplot2::margin(
        t = 0,
        r = 0,
        b = 0,
        l = 4
      )
    )
)


# Compact taxonomic-level legend
bar_legend <- cowplot::get_legend(
  barplot_integrated +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        title.position = "top",
        title.hjust = 0,
        ncol = 1
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.title = ggplot2::element_text(
        size = 10
      ),
      legend.text = ggplot2::element_text(
        size = 9
      ),
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    )
)


# ------------------------------------------------------------
# 6. Remove legends from plots
# ------------------------------------------------------------

heatmap_no_legend <- heatmap_integrated +
  ggplot2::theme(
    legend.position = "none"
  )

barplot_no_legend <- barplot_integrated +
  ggplot2::theme(
    legend.position = "none"
  )


# ------------------------------------------------------------
# 7. Stack the two legends vertically
# ------------------------------------------------------------

legend_block <- cowplot::plot_grid(
  heatmap_legend,
  bar_legend,
  ncol = 1,

  rel_heights = c(
    1.0,
    1.2
  ),

  align = "v"
)


# ------------------------------------------------------------
# 8. ALIGN HEATMAP AND BARPLOT PLOTTING PANELS
# ------------------------------------------------------------
#
# IMPORTANT:
# The legend block must NOT take up horizontal layout space beside the barplot.
# If it does, the barplot panel becomes narrower and its six comparison groups
# no longer sit directly below the six heatmap columns.
#
# cowplot::align_plots() aligns the actual plotting panels of the heatmap and
# barplot using their left and right axes. This ensures that the centre of each
# comparison in the barplot is directly beneath the corresponding heatmap
# column, even though the heatmap has long genus labels on its y-axis.
# ------------------------------------------------------------

aligned_DA_plots <- cowplot::align_plots(
  heatmap_no_legend,
  barplot_no_legend,
  align = "v",
  axis = "lr"
)

heatmap_aligned <- aligned_DA_plots[[1]]
barplot_aligned <- aligned_DA_plots[[2]]


# ------------------------------------------------------------
# 9. STACK THE ALIGNED HEATMAP AND BARPLOT
# ------------------------------------------------------------

combined_core <- cowplot::plot_grid(
  heatmap_aligned,
  barplot_aligned,

  ncol = 1,

  rel_heights = c(
    3.6,
    1.6
  ),

  align = "v",
  axis = "lr"
)


# ------------------------------------------------------------
# 10. ADD LEGENDS INTO THE LOWER-LEFT GUTTER
# ------------------------------------------------------------
#
# Once the two plotting panels are aligned, the long genus labels in the
# heatmap create a natural empty gutter to the left of the lower barplot.
# Place the legends in that gutter as an overlay rather than as a separate
# column. This preserves exact heatmap/barplot column alignment.
#
# If the journal/export changes font metrics slightly, only the four values
# below should need small visual adjustment; they do not affect the data panels.
# ------------------------------------------------------------

combined_DA_figure <- cowplot::ggdraw() +

  cowplot::draw_plot(
    combined_core,
    x = 0.01,
    y = 0.01,
    width = 0.98,
    height = 0.98
  ) +

  cowplot::draw_plot(
    legend_block,
    x = 0.005,
    y = 0.015,
    width = 0.205,
    height = 0.255
  )


# ------------------------------------------------------------
# 11. Save
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Differential_Abundance_Summary.pdf"
  ),
  plot = combined_DA_figure,
  width = 8.4,
  height = 8,
  units = "in",
  device = cairo_pdf
)


# ==============================================================================
# SAVE ONE STREAMLINED DIFFERENTIAL-ABUNDANCE RESULTS WORKBOOK
# ==============================================================================

# Compact summary of the design and the number of significant results.
analysis_summary <- tibble::tibble(
  Metric = c(
    "Samples analysed",
    "Independent mice",
    "Repeated-measure mice",
    "Minimum prevalence",
    "ASVs retained for modelling",
    "Named genera retained for modelling",
    "Named families retained for modelling",
    "Unique significant ASVs",
    "Unique significant genera",
    "Unique significant families",
    "Significant ASV-contrast associations",
    "Significant genus-contrast associations",
    "Significant family-contrast associations"
  ),
  Value = c(
    nrow(asv_input),
    repeat_summary$n_mice,
    repeat_summary$n_mice_repeated,
    minimum_prevalence,
    ncol(asv_input),
    ncol(genus_input),
    ncol(family_input),
    dplyr::n_distinct(asv_significant$ASV),
    dplyr::n_distinct(genus_significant$Genus),
    dplyr::n_distinct(family_significant$Family),
    nrow(asv_significant),
    nrow(genus_significant),
    nrow(family_significant)
  )
)

# Counts shown in the lower panel of the main figure, using the same wild-to-lab
# display order as the figure.
significant_counts_export <- count_data |>
  dplyr::mutate(
    comparison = as.character(comparison),
    taxonomic_level = as.character(taxonomic_level)
  ) |>
  tidyr::pivot_wider(
    names_from = taxonomic_level,
    values_from = n_significant
  ) |>
  dplyr::arrange(
    factor(
      comparison,
      levels = display_comparison_order
    )
  )

# Keep only FDR-significant associations in the reader-facing workbook.
# Complete model results are available in memory while this script runs.
asv_significant_export <- asv_significant |>
  dplyr::mutate(
    contrast = as.character(contrast)
  ) |>
  dplyr::arrange(
    contrast,
    qval_individual
  )

genus_significant_export <- genus_significant |>
  dplyr::mutate(
    contrast = as.character(contrast)
  ) |>
  dplyr::arrange(
    contrast,
    qval_individual
  )

family_significant_export <- family_significant |>
  dplyr::mutate(
    contrast = as.character(contrast)
  ) |>
  dplyr::arrange(
    contrast,
    qval_individual
  )

# One reader-facing workbook for the entire differential-abundance analysis.
#
# Sheets:
#   Analysis Summary
#     Sample/model counts and total numbers of significant results.
#
#   Taxonomic Coverage
#     Proportion of ASVs and reads represented by valid genus/family assignments.
#
#   Significant Counts
#     Numbers of significant ASVs, genera, and families for each pairwise contrast.
#
#   ASV Significant
#     All FDR-significant ASV-level associations with taxonomy.
#
#   Genus Significant
#     All FDR-significant genus-level associations.
#
#   Family Significant
#     All FDR-significant family-level associations.

da_results_workbook <- list(
  "Analysis Summary" =
    analysis_summary,

  "Taxonomic Coverage" =
    taxonomic_resolution_summary,

  "Significant Counts" =
    significant_counts_export,

  "ASV Significant" =
    asv_significant_export,

  "Genus Significant" =
    genus_significant_export,

  "Family Significant" =
    family_significant_export
)

openxlsx::write.xlsx(
  x = da_results_workbook,
  file = file.path(
    tables_dir,
    "Differential_Abundance_results.xlsx"
  ),
  overwrite = TRUE,
  rowNames = FALSE
)

cat(
  "\n============================================================\n",
  "DIFFERENTIAL-ABUNDANCE OUTPUTS\n",
  "============================================================\n",
  "Workbook:\n  ",
  file.path(
    tables_dir,
    "Differential_Abundance_results.xlsx"
  ),
  "\n\nComplete figure:\n  ",
  file.path(
    figures_dir,
    "Differential_Abundance_Summary.pdf"
  ),
  "\n============================================================\n",
  sep = ""
)

# Remove temporary MaAsLin3 output files created during this run.
if (dir.exists(model_dir)) {
  unlink(
    model_dir,
    recursive = TRUE,
    force = TRUE
  )
}

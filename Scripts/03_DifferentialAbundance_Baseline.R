
# ==============================================================================
# ANCOM-BC2 RESULTS: TABLES, SUMMARIES, AND FIGURES
#
# Required existing objects:
#   da_fixed : completed ANCOM-BC2 result
#   ps_ancom : filtered phyloseq object used in ANCOM-BC2
#
# Assumptions:
#   SampleOrigin is an unordered factor with levels:
#     Wild_As, Wildlike_As, Lablike_As, Lab_Mm
#
#   da_fixed was run with:
#     group = "SampleOrigin"
#     global = TRUE
#     pairwise = TRUE
# ==============================================================================


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(forcats)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(patchwork)

# ==============================================================================
# 3. LOAD DATA
# ==============================================================================

source("Scripts/0_DataPrep_Baseline.R")

# The sourced script should create ps_all.
stopifnot(exists("ps_all"))
stopifnot(inherits(ps_all, "phyloseq"))

ancom <- readRDS("C:/Users/wolfs/OneDrive/Documents/1. Ongoing Projects/LabToWildMicrobiome/Output/ANCOMBC2_SampleOrigin_Sex_fixed.rds")


# # prepare for ANCOMBC2 analysis -------------------------------------------
#  
# ## convert sample metadata to a regular data.frame.
# meta_ancom <- data.frame(
#   phyloseq::sample_data(ps_all),
#   check.names = FALSE
# )
# 
# ## define SampleOrigin as an UNORDERED factor.
# meta_ancom$SampleOrigin <- factor(
#   as.character(meta_ancom$SampleOrigin),
#   levels = c(
#     "Wild:As",
#     "Wild-like:As",
#     "Lab-like:As",
#     "Lab:Mm"
#   ),
#   labels = c(
#     "Wild_As",
#     "Wildlike_As",
#     "Lablike_As",
#     "Lab_Mm"
#   ),
#   ordered = FALSE
# )
# 
# 
# ## prepare other model variables.
# meta_ancom$Sex <- droplevels(
#   factor(meta_ancom$Sex)
# )
# 
# meta_ancom$mouse_ID_dataset <- droplevels(
#   factor(meta_ancom$mouse_ID_dataset)
# )
# 
# ## SeqPlate is retained in the metadata for descriptive checks, but it is not
# ## included in the final model because it is completely confounded with the
# ## Wild group in these data.
# meta_ancom$SeqPlate <- droplevels(
#   factor(meta_ancom$SeqPlate)
# )
# 
# ## remove samples missing variables required for the final fixed-effects model.
# #
# # mouse_ID_dataset is not required by this fixed model, but retaining only
# # samples with a valid mouse identifier keeps ps_ancom compatible with other
# # repeated-measures analyses.
# keep_samples <- rownames(meta_ancom)[
#   !is.na(meta_ancom$SampleOrigin) &
#     !is.na(meta_ancom$Sex) &
#     !is.na(meta_ancom$mouse_ID_dataset)
# ]
# 
# 
# ps_ancom <- phyloseq::prune_samples(
#   keep_samples,
#   ps_all
# )
# 
# 
# ## replace the sample metadata in ps_ancom with the cleaned metadata.
# meta_ancom <- meta_ancom[
#   phyloseq::sample_names(ps_ancom),
#   ,
#   drop = FALSE
# ]
# 
# phyloseq::sample_data(ps_ancom) <-
#   phyloseq::sample_data(meta_ancom)
# 
# ## remove zero-count taxa after sample pruning.
# ps_ancom <- phyloseq::prune_taxa(
#   phyloseq::taxa_sums(ps_ancom) > 0,
#   ps_ancom
# )
# 
# ## drop any unused levels again after pruning.
# phyloseq::sample_data(ps_ancom)$SampleOrigin <-
#   droplevels(
#     factor(
#       phyloseq::sample_data(ps_ancom)$SampleOrigin,
#       levels = c(
#         "Wild_As",
#         "Wildlike_As",
#         "Lablike_As",
#         "Lab_Mm"
#       ),
#       ordered = FALSE
#     )
#   )
# 
# phyloseq::sample_data(ps_ancom)$Sex <-
#   droplevels(
#     factor(
#       phyloseq::sample_data(ps_ancom)$Sex
#     )
#   )
# 
# phyloseq::sample_data(ps_ancom)$mouse_ID_dataset <-
#   droplevels(
#     factor(
#       phyloseq::sample_data(ps_ancom)$mouse_ID_dataset
#     )
#   )
# 
# ## VALIDATE MODEL INPUT
# meta_ancom <- data.frame(
#   phyloseq::sample_data(ps_ancom),
#   check.names = FALSE
# )
# 
# expected_origin_levels <- c(
#   "Wild_As",
#   "Wildlike_As",
#   "Lablike_As",
#   "Lab_Mm"
# )
# 
# if (!identical(
#   levels(meta_ancom$SampleOrigin),
#   expected_origin_levels
# )) {
#   stop(
#     paste0(
#       "SampleOrigin levels are incorrect.\n",
#       "Current levels: ",
#       paste(
#         levels(meta_ancom$SampleOrigin),
#         collapse = ", "
#       ),
#       "\nExpected levels: ",
#       paste(
#         expected_origin_levels,
#         collapse = ", "
#       )
#     )
#   )
# }
# 
# 
# if (is.ordered(meta_ancom$SampleOrigin)) {
#   stop(
#     "SampleOrigin is ordered. It must be an unordered factor."
#   )
# }
# 
# 
# if (nlevels(meta_ancom$Sex) < 2) {
#   stop(
#     "Sex has fewer than two observed levels after sample filtering."
#   )
# }
# 
# 
# ## check sample counts by group.
# sample_counts_check <- meta_ancom %>%
#   dplyr::count(
#     SampleOrigin,
#     name = "n_samples"
#   )
# 
# print(sample_counts_check)
# 
# 
# ## check sex distribution within groups.
# sex_counts_check <- meta_ancom %>%
#   dplyr::count(
#     SampleOrigin,
#     Sex,
#     name = "n_samples"
#   )
# 
# print(sex_counts_check)
# 
# 
# ## confirm that treatment contrasts—not polynomial contrasts—will be used.
# design_columns <- colnames(
#   stats::model.matrix(
#     ~ SampleOrigin + Sex,
#     data = meta_ancom,
#     contrasts.arg = list(
#       SampleOrigin = stats::contr.treatment(
#         nlevels(meta_ancom$SampleOrigin)
#       ),
#       Sex = stats::contr.treatment(
#         nlevels(meta_ancom$Sex)
#       )
#     )
#   )
# )
# 
# print(design_columns)
# 
# 
# # Expected SampleOrigin coefficient names include:
# # 
# #   SampleOriginWildlike_As
# #   SampleOriginLablike_As
# #   SampleOriginLab_Mm
# # 
# # They should not include:
# # 
# #   SampleOrigin.L
# #   SampleOrigin.Q
# #   SampleOrigin.C
# 
# if (any(grepl(
#   "SampleOrigin\\.[LQC]",
#   design_columns
# ))) {
#   stop(
#     "Polynomial SampleOrigin contrasts were detected."
#   )
# }
# 
# ## RUN FINAL FIXED-EFFECTS ANCOM-BC2 MODEL
# set.seed(999)
# 
# da_fixed <- ANCOMBC::ancombc2(
#   data = ps_ancom,
#   
#   # phyloseq count table
#   assay_name = "counts",
#   
#   # ASV-level analysis
#   tax_level = NULL,
#   
#   # Final fixed-effects model
#   fix_formula = "SampleOrigin + Sex",
#   
#   # No random effect in the final model
#   rand_formula = NULL,
#   
#   # Multiple-testing correction
#   p_adj_method = "BH",
#   
#   # Retain taxa present in at least 10% of samples
#   prv_cut = 0.10,
#   
#   # Variable used for global and pairwise comparisons
#   group = "SampleOrigin",
#   
#   # Run all six pairwise group comparisons
#   pairwise = TRUE,
#   
#   # Run the overall four-group test
#   global = TRUE,
#   
#   # Do not restrict comparisons to the reference group
#   dunnet = FALSE,
#   
#   # No ordered-trend test
#   trend = FALSE,
#   
#   alpha = 0.05,
#   
#   # Use one core for reproducibility and easier debugging
#   n_cl = 1
# )
# 
# stopifnot(
#   is.list(da_fixed),
#   is.data.frame(da_fixed$res),
#   is.data.frame(da_fixed$res_global),
#   is.data.frame(da_fixed$res_pair)
# )
# 
# cat(
#   "\nANCOM-BC2 MODEL COMPLETE\n",
#   "========================\n",
#   "Samples analysed: ",
#   phyloseq::nsamples(ps_ancom),
#   "\nASVs retained by ANCOM-BC2: ",
#   nrow(da_fixed$res),
#   "\nGlobally robust differential ASVs: ",
#   sum(
#     da_fixed$res_global$diff_robust_abn,
#     na.rm = TRUE
#   ),
#   "\n\n",
#   sep = ""
# )
# 
# print(
#   summary(
#     da_fixed$res_global
#   )
# )
# 
# ## confirm that pairwise columns have treatment-coded names.
# pairwise_lfc_columns <- grep(
#   "^lfc_",
#   names(da_fixed$res_pair),
#   value = TRUE
# )
# 
# print(pairwise_lfc_columns)
# 
# expected_pairwise_lfc_columns <- c(
#   "lfc_SampleOriginWildlike_As",
#   "lfc_SampleOriginLablike_As",
#   "lfc_SampleOriginLab_Mm",
#   paste0(
#     "lfc_SampleOriginLablike_As_",
#     "SampleOriginWildlike_As"
#   ),
#   paste0(
#     "lfc_SampleOriginLab_Mm_",
#     "SampleOriginWildlike_As"
#   ),
#   paste0(
#     "lfc_SampleOriginLab_Mm_",
#     "SampleOriginLablike_As"
#   )
# )
# 
# missing_pairwise_lfc_columns <- setdiff(
#   expected_pairwise_lfc_columns,
#   names(da_fixed$res_pair)
# )
# 
# if (length(missing_pairwise_lfc_columns) > 0) {
#   stop(
#     paste0(
#       "The expected pairwise output columns are missing:\n",
#       paste(
#         missing_pairwise_lfc_columns,
#         collapse = "\n"
#       ),
#       "\n\nActual lfc columns:\n",
#       paste(
#         pairwise_lfc_columns,
#         collapse = "\n"
#       )
#     )
#   )
# }
# 
# 
# ## SAVE MODEL AND ANALYSIS OBJECT
# ancom <- da_fixed
# saveRDS(
#   ancom, file = "Output/ANCOMBC2_SampleOrigin_Sex_fixed.rds")

# ==============================================================================
# 3. VALIDATE INPUT OBJECTS
# ==============================================================================

stopifnot(
  exists("ancom"),
  exists("ps_ancom"),
  inherits(ps_ancom, "phyloseq"),
  is.data.frame(ancom$res),
  is.data.frame(ancom$res_global),
  is.data.frame(ancom$res_pair)
)


expected_origin_levels <- c(
  "Wild_As",
  "Wildlike_As",
  "Lablike_As",
  "Lab_Mm"
)

origin_levels <- levels(
  phyloseq::sample_data(ps_ancom)$SampleOrigin
)

if (!identical(origin_levels, expected_origin_levels)) {
  stop(
    paste0(
      "SampleOrigin levels are not as expected.\n",
      "Current levels: ",
      paste(origin_levels, collapse = ", "),
      "\nExpected levels: ",
      paste(expected_origin_levels, collapse = ", ")
    )
  )
}

if (is.ordered(
  phyloseq::sample_data(ps_ancom)$SampleOrigin
)) {
  stop(
    "SampleOrigin is ordered. It must be an unordered factor ",
    "for named treatment and pairwise contrasts."
  )
}


# ==============================================================================
# 4. PLOTTING LABELS AND ORDER
# ==============================================================================

origin_labels <- c(
  Wild_As = "Wild",
  Wildlike_As = "Wild-like",
  Lablike_As = "Lab-like",
  Lab_Mm = "Lab"
)

origin_plot_levels <- unname(origin_labels)


# ==============================================================================
# 5. SAMPLE METADATA SUMMARY
# ==============================================================================

meta_ancom <- data.frame(
  phyloseq::sample_data(ps_ancom),
  check.names = FALSE
)


sample_counts <- meta_ancom %>%
  dplyr::count(
    SampleOrigin,
    name = "n_samples"
  ) %>%
  dplyr::mutate(
    SampleOrigin_label =
      origin_labels[as.character(SampleOrigin)]
  )


sample_sex_counts <- meta_ancom %>%
  dplyr::count(
    SampleOrigin,
    Sex,
    name = "n_samples"
  ) %>%
  dplyr::mutate(
    SampleOrigin_label =
      origin_labels[as.character(SampleOrigin)]
  )


plate_origin_table <- stats::xtabs(
  ~ SeqPlate + SampleOrigin,
  data = meta_ancom
)


# ==============================================================================
# 6. PREPARE TAXONOMY
# ==============================================================================

taxonomy <- as.data.frame(
  phyloseq::tax_table(ps_ancom),
  stringsAsFactors = FALSE
) %>%
  tibble::rownames_to_column("taxon")


taxonomy_ranks <- intersect(
  c(
    "Genus",
    "Family",
    "Order",
    "Class",
    "Phylum"
  ),
  names(taxonomy)
)


clean_taxonomy <- function(x) {
  
  x <- as.character(x)
  
  invalid <- is.na(x) |
    x == "" |
    tolower(x) %in% c(
      "na",
      "unknown",
      "unclassified",
      "uncultured",
      "uncultured bacterium",
      "unidentified"
    )
  
  x[invalid] <- NA_character_
  x
}


taxonomy[taxonomy_ranks] <- lapply(
  taxonomy[taxonomy_ranks],
  clean_taxonomy
)


choose_taxon_name <- function(x) {
  
  x <- as.character(x)
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return("Unclassified")
  }
  
  x[[1]]
}


taxonomy$taxon_name <- apply(
  taxonomy[, taxonomy_ranks, drop = FALSE],
  1,
  choose_taxon_name
)


taxonomy <- taxonomy %>%
  dplyr::mutate(
    taxon_label = paste0(
      taxon_name,
      " [",
      taxon,
      "]"
    )
  )


# ==============================================================================
# 7. GLOBAL RESULTS
# ==============================================================================

global_results <- ancom$res_global %>%
  dplyr::left_join(
    taxonomy,
    by = "taxon"
  ) %>%
  dplyr::arrange(q_val)


global_significant <- global_results %>%
  dplyr::filter(
    diff_robust_abn %in% TRUE
  )


global_summary <- tibble::tibble(
  n_samples = phyloseq::nsamples(ps_ancom),
  n_ASVs_in_phyloseq = phyloseq::ntaxa(ps_ancom),
  n_ASVs_tested = nrow(global_results),
  n_global_significant = sum(
    global_results$diff_robust_abn,
    na.rm = TRUE
  ),
  percent_global_significant = round(
    100 *
      n_global_significant /
      n_ASVs_tested,
    1
  ),
  prevalence_cutoff = 0.10,
  adjustment_method = "Benjamini-Hochberg",
  model = "SampleOrigin + Sex"
)

print(global_summary)


# ==============================================================================
# 8. DEFINE PAIRWISE CONTRASTS
# ==============================================================================

contrast_map <- tibble::tribble(
  ~suffix,
  ~group_1,
  ~group_2,
  ~comparison,
  ~comparison_short,
  ~comparison_order,
  ~comparison_type,
  
  "SampleOriginWildlike_As",
  "Wild-like",
  "Wild",
  "Wild-like vs Wild",
  "Wild-like / Wild",
  1,
  "Adjacent",
  
  "SampleOriginLablike_As_SampleOriginWildlike_As",
  "Lab-like",
  "Wild-like",
  "Lab-like vs Wild-like",
  "Lab-like / Wild-like",
  2,
  "Adjacent",
  
  "SampleOriginLab_Mm_SampleOriginLablike_As",
  "Lab",
  "Lab-like",
  "Lab vs Lab-like",
  "Lab / Lab-like",
  3,
  "Adjacent",
  
  "SampleOriginLablike_As",
  "Lab-like",
  "Wild",
  "Lab-like vs Wild",
  "Lab-like / Wild",
  4,
  "Non-adjacent",
  
  "SampleOriginLab_Mm_SampleOriginWildlike_As",
  "Lab",
  "Wild-like",
  "Lab vs Wild-like",
  "Lab / Wild-like",
  5,
  "Non-adjacent",
  
  "SampleOriginLab_Mm",
  "Lab",
  "Wild",
  "Lab vs Wild",
  "Lab / Wild",
  6,
  "Full gradient"
)


# ==============================================================================
# 9. VALIDATE PAIRWISE COLUMNS
# ==============================================================================

metric_prefixes <- c(
  "lfc_",
  "se_",
  "W_",
  "p_",
  "q_",
  "diff_",
  "passed_ss_",
  "diff_robust_"
)


required_pairwise_columns <- unlist(
  lapply(
    contrast_map$suffix,
    function(suffix_i) {
      paste0(metric_prefixes, suffix_i)
    }
  )
)


missing_pairwise_columns <- setdiff(
  required_pairwise_columns,
  names(da_fixed$res_pair)
)


if (length(missing_pairwise_columns) > 0) {
  stop(
    paste0(
      "Missing expected pairwise columns:\n",
      paste(
        missing_pairwise_columns,
        collapse = "\n"
      )
    )
  )
}


# ==============================================================================
# 10. PAIRWISE RESULTS IN LONG FORMAT
# ==============================================================================

pairwise_results <- purrr::map_dfr(
  seq_len(nrow(contrast_map)),
  function(i) {
    
    suffix_i <- contrast_map$suffix[[i]]
    
    da_fixed$res_pair %>%
      dplyr::transmute(
        taxon,
        
        suffix = suffix_i,
        group_1 = contrast_map$group_1[[i]],
        group_2 = contrast_map$group_2[[i]],
        comparison = contrast_map$comparison[[i]],
        comparison_short =
          contrast_map$comparison_short[[i]],
        comparison_order =
          contrast_map$comparison_order[[i]],
        comparison_type =
          contrast_map$comparison_type[[i]],
        
        lfc = .data[[paste0("lfc_", suffix_i)]],
        se = .data[[paste0("se_", suffix_i)]],
        W = .data[[paste0("W_", suffix_i)]],
        p_value = .data[[paste0("p_", suffix_i)]],
        q_value = .data[[paste0("q_", suffix_i)]],
        significant =
          .data[[paste0("diff_", suffix_i)]],
        passed_sensitivity =
          .data[[paste0("passed_ss_", suffix_i)]],
        robust_significant =
          .data[[paste0("diff_robust_", suffix_i)]]
      )
  }
) %>%
  dplyr::left_join(
    taxonomy,
    by = "taxon"
  ) %>%
  dplyr::mutate(
    comparison = factor(
      comparison,
      levels = contrast_map$comparison[
        order(contrast_map$comparison_order)
      ]
    ),
    
    comparison_short = factor(
      comparison_short,
      levels = contrast_map$comparison_short[
        order(contrast_map$comparison_order)
      ]
    ),
    
    enriched_in = dplyr::case_when(
      robust_significant %in% TRUE &
        lfc > 0 ~ group_1,
      
      robust_significant %in% TRUE &
        lfc < 0 ~ group_2,
      
      TRUE ~ NA_character_
    ),
    
    depleted_in = dplyr::case_when(
      robust_significant %in% TRUE &
        lfc > 0 ~ group_2,
      
      robust_significant %in% TRUE &
        lfc < 0 ~ group_1,
      
      TRUE ~ NA_character_
    ),
    
    direction = dplyr::case_when(
      robust_significant %in% TRUE &
        lfc > 0 ~ "Positive",
      
      robust_significant %in% TRUE &
        lfc < 0 ~ "Negative",
      
      TRUE ~ "Not significant"
    ),
    
    ci_low = lfc - 1.96 * se,
    ci_high = lfc + 1.96 * se
  ) %>%
  dplyr::arrange(
    comparison_order,
    q_value
  )


significant_pairwise_results <- pairwise_results %>%
  dplyr::filter(
    robust_significant %in% TRUE
  )


# ==============================================================================
# 11. PAIRWISE SUMMARY TABLES
# ==============================================================================

pairwise_summary <- pairwise_results %>%
  dplyr::group_by(
    comparison,
    comparison_short,
    comparison_order,
    comparison_type,
    group_1,
    group_2
  ) %>%
  dplyr::summarise(
    n_tested = dplyr::n(),
    
    n_significant = sum(
      robust_significant,
      na.rm = TRUE
    ),
    
    n_higher_in_group_1 = sum(
      robust_significant %in% TRUE &
        lfc > 0,
      na.rm = TRUE
    ),
    
    n_higher_in_group_2 = sum(
      robust_significant %in% TRUE &
        lfc < 0,
      na.rm = TRUE
    ),
    
    median_abs_lfc_significant = stats::median(
      abs(lfc[robust_significant %in% TRUE]),
      na.rm = TRUE
    ),
    
    max_abs_lfc_significant = max(
      abs(lfc[robust_significant %in% TRUE]),
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(comparison_order)


pairwise_direction_summary <- significant_pairwise_results %>%
  dplyr::count(
    comparison,
    comparison_short,
    comparison_order,
    enriched_in,
    name = "n_significant"
  ) %>%
  dplyr::arrange(comparison_order)


# Top taxa by q-value in each contrast.
top_pairwise_by_q <- significant_pairwise_results %>%
  dplyr::group_by(comparison) %>%
  dplyr::arrange(q_value, .by_group = TRUE) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::ungroup()


# Largest positive effects.
top_pairwise_increases <- significant_pairwise_results %>%
  dplyr::filter(lfc > 0) %>%
  dplyr::group_by(comparison) %>%
  dplyr::arrange(
    dplyr::desc(lfc),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::ungroup()


# Largest negative effects.
top_pairwise_decreases <- significant_pairwise_results %>%
  dplyr::filter(lfc < 0) %>%
  dplyr::group_by(comparison) %>%
  dplyr::arrange(
    lfc,
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::ungroup()


# ==============================================================================
# 12. RELATIVE ABUNDANCE FOR DESCRIPTIVE FIGURES
# ==============================================================================

ps_relative <- phyloseq::transform_sample_counts(
  ps_ancom,
  function(x) {
    if (sum(x) == 0) {
      return(x)
    }
    
    x / sum(x)
  }
)


abundance_long <- phyloseq::psmelt(
  ps_relative
) %>%
  tibble::as_tibble() %>%
  dplyr::rename(
    taxon = OTU,
    relative_abundance = Abundance
  ) %>%
  dplyr::left_join(
    taxonomy,
    by = "taxon"
  ) %>%
  dplyr::mutate(
    SampleOrigin_label = factor(
      origin_labels[as.character(SampleOrigin)],
      levels = origin_plot_levels
    )
  )


group_abundance_summary <- abundance_long %>%
  dplyr::group_by(
    taxon,
    taxon_name,
    taxon_label,
    SampleOrigin_label
  ) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    mean_relative_abundance = mean(
      relative_abundance,
      na.rm = TRUE
    ),
    median_relative_abundance = stats::median(
      relative_abundance,
      na.rm = TRUE
    ),
    prevalence = mean(
      relative_abundance > 0,
      na.rm = TRUE
    ),
    sd_relative_abundance = stats::sd(
      relative_abundance,
      na.rm = TRUE
    ),
    se_relative_abundance =
      sd_relative_abundance /
      sqrt(n_samples),
    .groups = "drop"
  )


group_abundance_wide <- group_abundance_summary %>%
  dplyr::select(
    taxon,
    taxon_name,
    taxon_label,
    SampleOrigin_label,
    mean_relative_abundance,
    prevalence
  ) %>%
  tidyr::pivot_wider(
    names_from = SampleOrigin_label,
    values_from = c(
      mean_relative_abundance,
      prevalence
    ),
    names_sep = "__"
  )


# Add abundance information to global results.
global_results_with_abundance <- global_results %>%
  dplyr::left_join(
    group_abundance_wide,
    by = c(
      "taxon",
      "taxon_name",
      "taxon_label"
    )
  )


# ==============================================================================
# 13. TRANSITION-PATTERN CLASSIFICATION
# ==============================================================================

adjacent_results <- pairwise_results %>%
  dplyr::filter(
    comparison %in% c(
      "Wild-like vs Wild",
      "Lab-like vs Wild-like",
      "Lab vs Lab-like"
    )
  ) %>%
  dplyr::select(
    taxon,
    comparison,
    lfc,
    robust_significant
  ) %>%
  tidyr::pivot_wider(
    names_from = comparison,
    values_from = c(
      lfc,
      robust_significant
    ),
    names_sep = "__"
  )


transition_patterns <- group_abundance_wide %>%
  dplyr::left_join(
    adjacent_results,
    by = "taxon"
  ) %>%
  dplyr::left_join(
    global_results %>%
      dplyr::select(
        taxon,
        global_q = q_val,
        global_significant =
          diff_robust_abn
      ),
    by = "taxon"
  ) %>%
  dplyr::mutate(
    wild =
      `mean_relative_abundance__Wild`,
    
    wildlike =
      `mean_relative_abundance__Wild-like`,
    
    lablike =
      `mean_relative_abundance__Lab-like`,
    
    lab =
      `mean_relative_abundance__Lab`,
    
    monotonic_increase =
      wild <= wildlike &
      wildlike <= lablike &
      lablike <= lab,
    
    monotonic_decrease =
      wild >= wildlike &
      wildlike >= lablike &
      lablike >= lab,
    
    sig_early =
      `robust_significant__Wild-like vs Wild`
    %in% TRUE,
    
    sig_middle =
      `robust_significant__Lab-like vs Wild-like`
    %in% TRUE,
    
    sig_late =
      `robust_significant__Lab vs Lab-like`
    %in% TRUE,
    
    n_significant_adjacent_steps =
      rowSums(
        cbind(
          sig_early,
          sig_middle,
          sig_late
        ),
        na.rm = TRUE
      ),
    
    first_significant_transition =
      dplyr::case_when(
        sig_early ~
          "Wild to Wild-like",
        
        !sig_early & sig_middle ~
          "Wild-like to Lab-like",
        
        !sig_early &
          !sig_middle &
          sig_late ~
          "Lab-like to Lab",
        
        TRUE ~
          "No significant adjacent step"
      ),
    
    transition_pattern =
      dplyr::case_when(
        monotonic_increase &
          global_significant %in% TRUE ~
          "Progressive increase",
        
        monotonic_decrease &
          global_significant %in% TRUE ~
          "Progressive decrease",
        
        sig_early &
          !sig_middle &
          !sig_late ~
          "Early shift only",
        
        !sig_early &
          sig_middle &
          !sig_late ~
          "Middle shift only",
        
        !sig_early &
          !sig_middle &
          sig_late ~
          "Late shift only",
        
        sig_early &
          sig_middle &
          sig_late ~
          "Persistent stepwise change",
        
        (sig_early + sig_middle + sig_late) >= 2 ~
          "Multiple transition stages",
        
        global_significant %in% TRUE ~
          "Non-monotonic or non-adjacent",
        
        TRUE ~
          "Not globally significant"
      ),
    
    lab_vs_wild_log2_ratio = log2(
      (lab + 1e-6) /
        (wild + 1e-6)
    )
  )


transition_summary <- transition_patterns %>%
  dplyr::filter(
    global_significant %in% TRUE
  ) %>%
  dplyr::count(
    transition_pattern,
    name = "n_ASVs",
    sort = TRUE
  )


first_transition_summary <- transition_patterns %>%
  dplyr::filter(
    global_significant %in% TRUE
  ) %>%
  dplyr::count(
    first_significant_transition,
    name = "n_ASVs",
    sort = TRUE
  )


persistence_summary <- transition_patterns %>%
  dplyr::filter(
    global_significant %in% TRUE
  ) %>%
  dplyr::count(
    n_significant_adjacent_steps,
    name = "n_ASVs"
  )


# ==============================================================================
# 14. TAXONOMIC SUMMARY TABLES
# ==============================================================================

make_taxonomic_summary <- function(
    results,
    rank_name
) {
  
  if (!rank_name %in% names(results)) {
    return(NULL)
  }
  
  results %>%
    dplyr::filter(
      !is.na(.data[[rank_name]])
    ) %>%
    dplyr::count(
      .data[[rank_name]],
      name = "n_significant_ASVs",
      sort = TRUE
    )
}


global_family_summary <- make_taxonomic_summary(
  global_significant,
  "Family"
)

global_genus_summary <- make_taxonomic_summary(
  global_significant,
  "Genus"
)

global_phylum_summary <- make_taxonomic_summary(
  global_significant,
  "Phylum"
)


pairwise_family_summary <- significant_pairwise_results %>%
  dplyr::filter(!is.na(Family)) %>%
  dplyr::count(
    comparison,
    Family,
    name = "n_significant_ASVs"
  ) %>%
  dplyr::arrange(
    comparison,
    dplyr::desc(n_significant_ASVs)
  )


# ==============================================================================
# 15. FIGURE 1: SIGNIFICANT ASV COUNTS AND DIRECTION
# ==============================================================================

p_pairwise_counts <- ggplot2::ggplot(
  pairwise_direction_summary,
  ggplot2::aes(
    x = comparison,
    y = n_significant,
    fill = enriched_in
  )
) +
  ggplot2::geom_col(
    width = 0.75
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "Number of robust differential ASVs",
    fill = "Higher abundance in",
    title = "Differential ASVs across the wild-to-lab gradient",
    subtitle = paste0(
      global_summary$n_global_significant,
      " of ",
      global_summary$n_ASVs_tested,
      " ASVs were globally significant"
    )
  ) +
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    plot.title.position = "plot",
    legend.position = "right"
  )

print(p_pairwise_counts)


# ==============================================================================
# 16. VOLCANO PLOT FUNCTION
# ==============================================================================
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

label_n <- 8

comparison_order <- c(
  "Wild-like vs Wild",
  "Lab-like vs Wild-like",
  "Lab vs Lab-like",
  "Lab-like vs Wild",
  "Lab vs Wild-like",
  "Lab vs Wild"
)

volcano_colours <- c(
  "Positive LFC" = "#F8766D",
  "Negative LFC" = "#00BA38",
  "Not significant" = "#619CFF"
)

missing_comparisons <- setdiff(
  comparison_order,
  unique(pairwise_results$comparison)
)

if (length(missing_comparisons) > 0) {
  stop(
    paste0(
      "Missing comparisons:\n",
      paste(missing_comparisons, collapse = "\n")
    )
  )
}

## volcano function
make_volcano <- function(
    comparison_name,
    label_n = 8
) {
  
  plot_data <- pairwise_results %>%
    dplyr::filter(
      comparison == comparison_name
    ) %>%
    dplyr::mutate(
      q_plot = pmax(
        q_value,
        1e-300
      ),
      
      plot_group = dplyr::case_when(
        robust_significant %in% TRUE & lfc > 0 ~
          "Positive LFC",
        
        robust_significant %in% TRUE & lfc < 0 ~
          "Negative LFC",
        
        TRUE ~
          "Not significant"
      ),
      
      plot_group = factor(
        plot_group,
        levels = c(
          "Positive LFC",
          "Negative LFC",
          "Not significant"
        )
      )
    )
  
  
  group_1_name <- unique(
    as.character(plot_data$group_1)
  )
  
  group_2_name <- unique(
    as.character(plot_data$group_2)
  )
  
  
  if (
    length(group_1_name) != 1 ||
    length(group_2_name) != 1
  ) {
    stop(
      paste0(
        "Comparison '",
        comparison_name,
        "' does not have exactly one group_1 and group_2."
      )
    )
  }
  
  
  label_data <- plot_data %>%
    dplyr::filter(
      robust_significant %in% TRUE,
      !is.na(taxon_name),
      taxon_name != ""
    ) %>%
    dplyr::arrange(
      q_value,
      dplyr::desc(abs(lfc))
    ) %>%
    dplyr::slice_head(
      n = label_n
    )
  
  
  n_robust <- sum(
    plot_data$robust_significant %in% TRUE,
    na.rm = TRUE
  )
  
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = lfc,
      y = -log10(q_plot),
      colour = plot_group
    )
  ) +
    
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = 2,
      linewidth = 0.35,
      colour = "grey35"
    ) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.35,
      colour = "grey35"
    ) +
    
    ggplot2::geom_point(
      size = 1.8,
      alpha = 0.75,
      na.rm = TRUE
    ) +
    
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(
        label = taxon_name
      ),
      size = 2.6,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.20,
      segment.linewidth = 0.25,
      segment.colour = "grey55",
      show.legend = FALSE,
      seed = 999
    ) +
    
    # This scale must be exactly identical in every panel.
    ggplot2::scale_colour_manual(
      values = volcano_colours,
      breaks = c(
        "Positive LFC",
        "Negative LFC"
      ),
      labels = c(
        "Higher in first-named group",
        "Higher in second-named group"
      ),
      drop = FALSE
    ) + 
    
    ggplot2::labs(
      x = paste0(
        "Log fold change\n"
      ),
      
      y = expression(
        -log[10]("q-value")
      ),
      
      colour = NULL,
      
      title = comparison_name,
      
    ) +
    
    ggplot2::theme_classic(
      base_size = 10
    ) +
    
    ggplot2::theme(
      plot.title.position = "plot",
      
      plot.title = ggplot2::element_text(
        size = 11,
        margin = ggplot2::margin(b = 1)
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 9,
        margin = ggplot2::margin(b = 4)
      ),
      
      axis.title = ggplot2::element_text(
        size = 9
      ),
      
      axis.text = ggplot2::element_text(
        size = 8
      ),
      
      legend.position = "bottom",
      
      plot.margin = ggplot2::margin(
        t = 4,
        r = 6,
        b = 4,
        l = 4
      )
    )
}

## create plates
volcano_plots <- stats::setNames(
  lapply(
    comparison_order,
    make_volcano,
    label_n = label_n
  ),
  comparison_order
)

## combine plots
p_all_volcanoes <- patchwork::wrap_plots(
  volcano_plots[comparison_order],
  ncol = 3,
  guides = "collect"
) +
  
  patchwork::plot_annotation(
    tag_levels = "A",
    
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold",
        hjust = 0
      ),
      
      plot.tag = ggplot2::element_text(
        size = 12,
        face = "plain"
      )
    )
  ) &
  
  ggplot2::theme(
    legend.position = "bottom",
    
    legend.box = "horizontal",
    
    legend.justification = "center",
    
    legend.text = ggplot2::element_text(
      size = 8
    ),
    
    legend.key.width = grid::unit(
      0.6,
      "cm"
    ),
    
    legend.key.height = grid::unit(
      0.35,
      "cm"
    ),
    
    legend.spacing.x = grid::unit(
      0.08,
      "cm"
    ),
    
    legend.margin = ggplot2::margin(
      t = 2,
      r = 0,
      b = 0,
      l = 0
    )
  )

print(p_all_volcanoes)


# ==============================================================================
# 17. FIGURE 3: GLOBAL GROUP-MEAN HEATMAP
# ==============================================================================

n_heatmap_taxa <- min(
  50,
  nrow(global_significant)
)


top_global_taxa <- global_significant %>%
  dplyr::slice_min(
    q_val,
    n = n_heatmap_taxa,
    with_ties = FALSE
  ) %>%
  dplyr::pull(taxon)


heatmap_data <- group_abundance_summary %>%
  dplyr::filter(
    taxon %in% top_global_taxa
  ) %>%
  dplyr::mutate(
    log_mean_abundance = log10(
      mean_relative_abundance + 1e-6
    )
  ) %>%
  dplyr::group_by(taxon) %>%
  dplyr::mutate(
    row_z_score = as.numeric(
      scale(log_mean_abundance)
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    global_results %>%
      dplyr::select(
        taxon,
        global_q = q_val
      ),
    by = "taxon"
  ) %>%
  dplyr::mutate(
    taxon_label = forcats::fct_reorder(
      taxon_label,
      global_q,
      .desc = TRUE
    )
  )


p_global_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(
    x = SampleOrigin_label,
    y = taxon_label,
    fill = row_z_score
  )
) +
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.2
  ) +
  ggplot2::scale_fill_gradient2(
    midpoint = 0,
    name = "Row\nz-score"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    title = paste0(
      "Top ",
      n_heatmap_taxa,
      " globally differential ASVs"
    ),
    subtitle = "Row-scaled log10 mean relative abundance"
  ) +
  ggplot2::theme_minimal(
    base_size = 10
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 30,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_text(
      size = 6.5
    ),
    plot.title.position = "plot"
  )

print(p_global_heatmap)


# ==============================================================================
# 18. FIGURE 4: PAIRWISE SIGNIFICANCE TILE PLOT
# ==============================================================================

n_tile_taxa <- 50


tile_taxa <- pairwise_results %>%
  dplyr::group_by(taxon) %>%
  dplyr::summarise(
    n_significant = sum(
      robust_significant,
      na.rm = TRUE
    ),
    min_q = min(
      q_value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_significant > 0) %>%
  dplyr::arrange(
    dplyr::desc(n_significant),
    min_q
  ) %>%
  dplyr::slice_head(n = n_tile_taxa) %>%
  dplyr::pull(taxon)


tile_data <- pairwise_results %>%
  dplyr::filter(
    taxon %in% tile_taxa
  ) %>%
  dplyr::mutate(
    tile_value = dplyr::case_when(
      robust_significant %in% TRUE &
        lfc > 0 ~ 1,
      
      robust_significant %in% TRUE &
        lfc < 0 ~ -1,
      
      TRUE ~ 0
    )
  ) %>%
  dplyr::group_by(taxon) %>%
  dplyr::mutate(
    n_significant = sum(
      robust_significant,
      na.rm = TRUE
    ),
    strongest_effect = max(
      abs(lfc),
      na.rm = TRUE
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    taxon_label = forcats::fct_reorder(
      taxon_label,
      n_significant * 100 +
        strongest_effect
    )
  )


p_significance_tile <- ggplot2::ggplot(
  tile_data,
  ggplot2::aes(
    x = comparison_short,
    y = taxon_label,
    fill = tile_value
  )
) +
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.15
  ) +
  ggplot2::scale_fill_gradient2(
    low = "steelblue",
    high = "firebrick",
    midpoint = 0,
    breaks = c(-1, 1),
    labels = c(
      "Higher in group 2",
      "Higher in group 1"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    title = "Direction and timing of differential abundance",
    subtitle = paste0(
      "Top ",
      n_tile_taxa,
      " ASVs ranked by number of significant contrasts"
    )
  ) +
  ggplot2::theme_minimal(
    base_size = 9
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 40,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_text(
      size = 5.5
    ),
    legend.position = "bottom"
  )

print(p_significance_tile)

# ==============================================================================
# 20. FIGURE 6: REPRESENTATIVE TAXON TRAJECTORIES
# ==============================================================================

trajectory_taxa <- transition_patterns %>%
  dplyr::filter(
    global_significant %in% TRUE,
    transition_pattern %in% c(
      "Progressive increase",
      "Progressive decrease",
      "Early shift only",
      "Middle shift only",
      "Late shift only",
      "Persistent stepwise change"
    )
  ) %>%
  dplyr::group_by(transition_pattern) %>%
  dplyr::arrange(
    global_q,
    dplyr::desc(
      abs(lab_vs_wild_log2_ratio)
    ),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 2) %>%
  dplyr::ungroup() %>%
  dplyr::pull(taxon)


trajectory_data <- group_abundance_summary %>%
  dplyr::filter(
    taxon %in% trajectory_taxa
  ) %>%
  dplyr::left_join(
    transition_patterns %>%
      dplyr::select(
        taxon,
        transition_pattern
      ),
    by = "taxon"
  ) %>%
  dplyr::mutate(
    mean_percent =
      100 * mean_relative_abundance,
    
    se_percent =
      100 * se_relative_abundance
  )


p_taxon_trajectories <- ggplot2::ggplot(
  trajectory_data,
  ggplot2::aes(
    x = SampleOrigin_label,
    y = mean_percent,
    group = taxon
  )
) +
  ggplot2::geom_line(
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    size = 2
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = pmax(
        mean_percent - se_percent,
        0
      ),
      ymax =
        mean_percent + se_percent
    ),
    width = 0.1
  ) +
  ggplot2::facet_wrap(
    ~ taxon_label + transition_pattern,
    scales = "free_y",
    ncol = 4
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Mean relative abundance (%)",
    title = "Representative microbial trajectories",
    subtitle = "Group means ± standard error"
  ) +
  ggplot2::theme_classic(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    ),
    strip.text = ggplot2::element_text(
      size = 7
    )
  )

print(p_taxon_trajectories)


# ==============================================================================
# 24. FIGURE 10: TRANSITION-PATTERN COUNTS
# ==============================================================================

p_transition_counts <- transition_summary %>%
  dplyr::mutate(
    transition_pattern = forcats::fct_reorder(
      transition_pattern,
      n_ASVs
    )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = transition_pattern,
      y = n_ASVs
    )
  ) +
  ggplot2::geom_col(
    width = 0.75
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "Number of globally significant ASVs",
    title = "Timing and shape of microbial responses"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  )

print(p_transition_counts)


# ==============================================================================
# 28. SAVE ALL FIGURES
# ==============================================================================

pdf("Output/Graphs/ANCOM_pairwise_counts.pdf", width = 5, height = 4)
print(p_pairwise_counts)
dev.off()

pdf("Output/Graphs/ANCOM_volcano.pdf", width = 10, height = 6)
print(p_all_volcanoes)
dev.off()

pdf("Output/Graphs/ANCOM_global_heatmap.pdf", width = 7, height = 7)
print(p_global_heatmap)
dev.off()

pdf("Output/Graphs/ANCOM_significance_tiles.pdf", width = 7, height = 7)
print(p_significance_tile)
dev.off()

# ==============================================================================
# 29. EXCEL EXPORT
# ==============================================================================

wb <- openxlsx::createWorkbook()


add_sheet <- function(
    workbook,
    sheet_name,
    data
) {
  
  if (is.null(data)) {
    return(invisible(NULL))
  }
  
  openxlsx::addWorksheet(
    workbook,
    sheet_name
  )
  
  openxlsx::writeData(
    workbook,
    sheet_name,
    data
  )
  
  openxlsx::freezePane(
    workbook,
    sheet_name,
    firstRow = TRUE
  )
  
  openxlsx::setColWidths(
    workbook,
    sheet_name,
    cols = seq_len(ncol(data)),
    widths = "auto"
  )
}


add_sheet(
  wb,
  "Analysis summary",
  global_summary
)

add_sheet(
  wb,
  "Sample counts",
  sample_counts
)

add_sheet(
  wb,
  "Sample by sex",
  sample_sex_counts
)

add_sheet(
  wb,
  "Global all",
  global_results_with_abundance
)

add_sheet(
  wb,
  "Global significant",
  global_results_with_abundance %>%
    dplyr::filter(
      diff_robust_abn %in% TRUE
    )
)

add_sheet(
  wb,
  "Pairwise summary",
  pairwise_summary
)

add_sheet(
  wb,
  "Pairwise all",
  pairwise_results
)

add_sheet(
  wb,
  "Pairwise significant",
  significant_pairwise_results
)

add_sheet(
  wb,
  "Top pairwise by q",
  top_pairwise_by_q
)

add_sheet(
  wb,
  "Largest increases",
  top_pairwise_increases
)

add_sheet(
  wb,
  "Largest decreases",
  top_pairwise_decreases
)

add_sheet(
  wb,
  "Transition patterns",
  transition_patterns
)

add_sheet(
  wb,
  "Transition summary",
  transition_summary
)

add_sheet(
  wb,
  "First transition",
  first_transition_summary
)

add_sheet(
  wb,
  "Persistence summary",
  persistence_summary
)

add_sheet(
  wb,
  "Family global",
  global_family_summary
)

add_sheet(
  wb,
  "Genus global",
  global_genus_summary
)

add_sheet(
  wb,
  "Phylum global",
  global_phylum_summary
)

add_sheet(
  wb,
  "Family pairwise",
  pairwise_family_summary
)

add_sheet(
  wb,
  "Group abundances",
  group_abundance_summary
)


openxlsx::saveWorkbook(
  wb,
  file = paste0(
    "Output/ANCOMBC2/",
    "ANCOMBC2_complete_results.xlsx"
  ),
  overwrite = TRUE
)


# ==============================================================================
# 30. EXPORT ANCOM-BC2 RESULTS TO ONE EXCEL WORKBOOK
# ==============================================================================

library(openxlsx)
library(dplyr)

# Prepare the tables to export
excel_tables <- list(
  "Global_all" =
    global_results_with_abundance,
  
  "Global_significant" =
    global_results_with_abundance %>%
    dplyr::filter(
      diff_robust_abn %in% TRUE
    ),
  
  "Pairwise_all" =
    pairwise_results,
  
  "Pairwise_significant" =
    significant_pairwise_results,
  
  "Pairwise_summary" =
    pairwise_summary,
  
  "Transition_patterns" =
    transition_patterns
)

# Create workbook
wb <- openxlsx::createWorkbook()

# Header style
header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom",
  wrapText = TRUE
)

# Add each table as a separate worksheet
for (sheet_name in names(excel_tables)) {
  
  table_i <- excel_tables[[sheet_name]]
  
  openxlsx::addWorksheet(
    wb,
    sheetName = sheet_name
  )
  
  openxlsx::writeData(
    wb,
    sheet = sheet_name,
    x = table_i,
    headerStyle = header_style,
    withFilter = TRUE
  )
  
  # Freeze the header row
  openxlsx::freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
  
  # Automatically size columns
  openxlsx::setColWidths(
    wb,
    sheet = sheet_name,
    cols = seq_len(ncol(table_i)),
    widths = "auto"
  )
}

# Save the workbook
openxlsx::saveWorkbook(
  wb,
  file = "Output/ANCOMBC2_results.xlsx",
  overwrite = TRUE
)

# ==============================================================================
# 31. PRINT RESULTS SUMMARY
# ==============================================================================

cat(
  "\nANCOM-BC2 RESULTS SUMMARY\n",
  "========================\n",
  "Samples: ",
  global_summary$n_samples,
  "\nASVs tested: ",
  global_summary$n_ASVs_tested,
  "\nGlobally significant ASVs: ",
  global_summary$n_global_significant,
  " (",
  global_summary$percent_global_significant,
  "%)\n\n",
  sep = ""
)


print(
  pairwise_summary %>%
    dplyr::select(
      comparison,
      n_significant,
      n_higher_in_group_1,
      n_higher_in_group_2,
      median_abs_lfc_significant
    )
)


cat(
  "\nTransition-pattern summary\n",
  "--------------------------\n"
)

print(transition_summary)


cat(
  "\nEarliest significant adjacent transition\n",
  "----------------------------------------\n"
)

print(first_transition_summary)

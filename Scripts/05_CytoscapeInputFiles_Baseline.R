# ==============================================================================
# CYTOSCAPE GROUP-FAMILY NETWORK
# THRESHOLD-CONSISTENT PIE CHARTS
#
# Required object:
#   ps_ancom
#
# Optional object:
#   da_fixed
#
# This script creates:
#
#   1. Four group nodes:
#        Wild
#        Wild-like
#        Lab-like
#        Lab
#
#   2. Family nodes connected to every group in which they meet the defined
#      prevalence and abundance thresholds.
#
#   3. Exact overlap classes such as:
#        Wild
#        Wild__Wildlike
#        Wildlike__Lablike__Lab
#        Wild__Wildlike__Lablike__Lab
#
#   4. Threshold-consistent pie-chart values:
#        - a group colour appears only when the family passes the presence
#          definition in that group;
#        - therefore, a "Wild only" family has a fully Wild-coloured pie.
#
#   5. Family node size based on overall mean relative abundance.
#
# No X/Y positions are exported. Layout is done interactively in Cytoscape.
#
# Main Cytoscape files:
#   cytoscape_group_family_nodes.csv
#   cytoscape_group_family_edges.csv
#
# Supporting files:
#   family_group_summary.csv
#   family_node_attributes.csv
#   overlap_summary.csv
#   pie_membership_check.csv
#   network_settings.csv
# ==============================================================================


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(readr)


# ==============================================================================
# 2. OUTPUT DIRECTORY
# ==============================================================================

output_dir <- "Output/Cytoscape_group_family_final"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. USER SETTINGS
# ==============================================================================

# Taxonomic rank represented by the network nodes.
taxonomic_rank <- "Family"


# ------------------------------------------------------------------------------
# Family-retention thresholds
#
# A family is retained if it meets BOTH thresholds in at least one group.
#
# Relative abundance is represented as a proportion:
#   0.001  = 0.1%
#   0.0001 = 0.01%
# ------------------------------------------------------------------------------

min_group_prevalence <- 0.10
min_group_mean_abundance <- 0.001


# ------------------------------------------------------------------------------
# Group-presence thresholds
#
# These determine:
#   - whether a group-to-family edge is created;
#   - the family's exact overlap class;
#   - which colours appear in the family pie chart.
#
# Keeping these equal to the retention thresholds is usually clearest.
# ------------------------------------------------------------------------------

presence_prevalence_cutoff <- 0.10
presence_abundance_cutoff <- 0.001


# ------------------------------------------------------------------------------
# Unclassified-family handling
#
# FALSE:
#   Keep each unresolved ASV as a separate short node, e.g. UC_ASV17.
#
# TRUE:
#   Pool all unresolved family assignments into one node called Unclassified.
#
# Keeping them separate avoids combining unrelated ASVs.
# ------------------------------------------------------------------------------

collapse_unclassified <- FALSE


# ------------------------------------------------------------------------------
# Remove unresolved family nodes entirely
#
# FALSE:
#   Keep UC_ASV nodes.
#
# TRUE:
#   Remove all unresolved family nodes and their edges.
# ------------------------------------------------------------------------------

remove_unclassified_nodes <- FALSE


# ------------------------------------------------------------------------------
# Node-size range exported for Cytoscape
#
# Node diameter is scaled from sqrt(overall mean relative abundance), so node
# AREA is approximately proportional to abundance.
# ------------------------------------------------------------------------------

minimum_family_size <- 22
maximum_family_size <- 82

# Fixed size for the four central group nodes.
group_node_size <- 100


# ------------------------------------------------------------------------------
# Edge-width range exported for Cytoscape
# ------------------------------------------------------------------------------

minimum_edge_width <- 0.5
maximum_edge_width <- 6


# Small constant used only for logarithms.
pseudocount <- 1e-6


# ==============================================================================
# 4. VALIDATE INPUT
# ==============================================================================

stopifnot(
  exists("ps_ancom"),
  inherits(ps_ancom, "phyloseq")
)


available_ranks <- colnames(
  phyloseq::tax_table(ps_ancom)
)


if (!taxonomic_rank %in% available_ranks) {
  
  stop(
    paste0(
      "Taxonomic rank '",
      taxonomic_rank,
      "' is unavailable.\nAvailable ranks: ",
      paste(
        available_ranks,
        collapse = ", "
      )
    )
  )
}


expected_groups <- c(
  "Wild_As",
  "Wildlike_As",
  "Lablike_As",
  "Lab_Mm"
)


observed_groups <- unique(
  as.character(
    phyloseq::sample_data(ps_ancom)$SampleOrigin
  )
)


missing_groups <- setdiff(
  expected_groups,
  observed_groups
)


if (length(missing_groups) > 0) {
  
  stop(
    paste0(
      "Missing SampleOrigin groups: ",
      paste(
        missing_groups,
        collapse = ", "
      )
    )
  )
}


# ==============================================================================
# 5. GROUP DEFINITIONS
# ==============================================================================

group_lookup <- tibble::tribble(
  ~SampleOrigin, ~group_id,         ~group_label, ~group_order,
  
  "Wild_As",     "GROUP__Wild",     "Wild",       1,
  "Wildlike_As", "GROUP__Wildlike", "Wild-like",  2,
  "Lablike_As",  "GROUP__Lablike",  "Lab-like",   3,
  "Lab_Mm",      "GROUP__Lab",      "Lab",        4
)


# ==============================================================================
# 6. EXTRACT SAMPLE METADATA
# ==============================================================================

metadata <- data.frame(
  phyloseq::sample_data(ps_ancom),
  check.names = FALSE
) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::mutate(
    SampleOrigin = factor(
      SampleOrigin,
      levels = expected_groups,
      ordered = FALSE
    )
  ) %>%
  dplyr::left_join(
    group_lookup,
    by = "SampleOrigin"
  )


if (any(is.na(metadata$group_id))) {
  stop("Some samples could not be assigned to a SampleOrigin group.")
}


# ==============================================================================
# 7. EXTRACT ASV COUNT TABLE
# ==============================================================================

asv_counts <- as(
  phyloseq::otu_table(ps_ancom),
  "matrix"
)


# Ensure rows are taxa and columns are samples.
if (!phyloseq::taxa_are_rows(ps_ancom)) {
  asv_counts <- t(asv_counts)
}


missing_count_samples <- setdiff(
  metadata$sample_id,
  colnames(asv_counts)
)


if (length(missing_count_samples) > 0) {
  
  stop(
    paste0(
      "The following metadata samples are absent from the count table:\n",
      paste(
        missing_count_samples,
        collapse = "\n"
      )
    )
  )
}


asv_counts_long <- as.data.frame(
  asv_counts,
  check.names = FALSE
) %>%
  tibble::rownames_to_column("taxon") %>%
  tidyr::pivot_longer(
    cols = -taxon,
    names_to = "sample_id",
    values_to = "count"
  )


# ==============================================================================
# 8. EXTRACT AND CLEAN TAXONOMY
# ==============================================================================

taxonomy <- as.data.frame(
  phyloseq::tax_table(ps_ancom),
  stringsAsFactors = FALSE
) %>%
  tibble::rownames_to_column("taxon")


# Extract a numeric part from an ASV ID.
#
# ASV30 becomes 30.
# If no number is present, the original taxon identifier is retained in a
# shortened Cytoscape-safe form.

make_short_asv_label <- function(taxon_id) {
  
  taxon_id <- as.character(taxon_id)
  
  asv_number <- stringr::str_extract(
    taxon_id,
    "[0-9]+"
  )
  
  fallback <- stringr::str_replace_all(
    taxon_id,
    "[^A-Za-z0-9]",
    ""
  )
  
  dplyr::if_else(
    !is.na(asv_number),
    paste0("UC_ASV", asv_number),
    paste0(
      "UC_",
      stringr::str_sub(
        fallback,
        1,
        15
      )
    )
  )
}


make_network_taxon <- function(
    rank_value,
    taxon_id,
    collapse = FALSE
) {
  
  rank_value <- as.character(rank_value)
  
  invalid <- is.na(rank_value) |
    rank_value == "" |
    tolower(rank_value) %in% c(
      "na",
      "unknown",
      "unclassified",
      "uncultured",
      "uncultured bacterium",
      "unidentified",
      "unassigned"
    )
  
  
  if (collapse) {
    
    rank_value[invalid] <- "Unclassified"
    
  } else {
    
    rank_value[invalid] <- make_short_asv_label(
      taxon_id[invalid]
    )
  }
  
  
  # Make node labels safe for Cytoscape identifiers.
  stringr::str_replace_all(
    rank_value,
    "[^A-Za-z0-9_.-]",
    "_"
  )
}


taxonomy <- taxonomy %>%
  dplyr::mutate(
    network_taxon = make_network_taxon(
      .data[[taxonomic_rank]],
      taxon,
      collapse = collapse_unclassified
    ),
    
    unclassified_taxon = stringr::str_detect(
      network_taxon,
      "^UC_ASV|^UC_|^Unclassified$"
    ),
    
    family_id = paste0(
      "FAMILY__",
      network_taxon
    )
  )


# ==============================================================================
# 9. AGGREGATE ASVs TO FAMILY WITHIN EACH SAMPLE
# ==============================================================================

family_sample_counts <- asv_counts_long %>%
  dplyr::left_join(
    taxonomy %>%
      dplyr::select(
        taxon,
        family_id,
        network_taxon,
        unclassified_taxon
      ),
    by = "taxon"
  ) %>%
  dplyr::group_by(
    sample_id,
    family_id,
    network_taxon,
    unclassified_taxon
  ) %>%
  dplyr::summarise(
    count = sum(
      count,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ==============================================================================
# 10. SAMPLE-LEVEL RELATIVE ABUNDANCE
# ==============================================================================

family_sample_abundance <- family_sample_counts %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(
    library_size = sum(
      count,
      na.rm = TRUE
    ),
    
    relative_abundance = dplyr::if_else(
      library_size > 0,
      count / library_size,
      0
    )
  ) %>%
  dplyr::ungroup()


# ==============================================================================
# 11. COMPLETE SAMPLE-BY-FAMILY TABLE
#
# Explicit zeros are required for correct prevalence and group means.
# ==============================================================================

all_families <- family_sample_abundance %>%
  dplyr::distinct(
    family_id,
    network_taxon,
    unclassified_taxon
  )


family_sample_complete <- tidyr::crossing(
  sample_id = metadata$sample_id,
  family_id = all_families$family_id
) %>%
  dplyr::left_join(
    all_families,
    by = "family_id"
  ) %>%
  dplyr::left_join(
    family_sample_abundance %>%
      dplyr::select(
        sample_id,
        family_id,
        count,
        relative_abundance
      ),
    by = c(
      "sample_id",
      "family_id"
    )
  ) %>%
  dplyr::mutate(
    count = tidyr::replace_na(
      count,
      0
    ),
    
    relative_abundance = tidyr::replace_na(
      relative_abundance,
      0
    )
  ) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(
        sample_id,
        SampleOrigin,
        group_id,
        group_label,
        group_order
      ),
    by = "sample_id"
  )


# ==============================================================================
# 12. GROUP-LEVEL FAMILY SUMMARY
# ==============================================================================

family_group_summary <- family_sample_complete %>%
  dplyr::group_by(
    family_id,
    network_taxon,
    unclassified_taxon,
    SampleOrigin,
    group_id,
    group_label,
    group_order
  ) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    
    total_count = sum(
      count,
      na.rm = TRUE
    ),
    
    mean_relative_abundance = mean(
      relative_abundance,
      na.rm = TRUE
    ),
    
    median_relative_abundance = stats::median(
      relative_abundance,
      na.rm = TRUE
    ),
    
    prevalence = mean(
      count > 0,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# ==============================================================================
# 13. OVERALL FAMILY SUMMARY
#
# Each sample receives equal weight.
# ==============================================================================

overall_family_summary <- family_sample_complete %>%
  dplyr::group_by(
    family_id,
    network_taxon,
    unclassified_taxon
  ) %>%
  dplyr::summarise(
    overall_mean_relative_abundance = mean(
      relative_abundance,
      na.rm = TRUE
    ),
    
    overall_median_relative_abundance = stats::median(
      relative_abundance,
      na.rm = TRUE
    ),
    
    overall_prevalence = mean(
      count > 0,
      na.rm = TRUE
    ),
    
    total_count_all_samples = sum(
      count,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# ==============================================================================
# 14. RETAIN FAMILIES
# ==============================================================================

family_retention <- family_group_summary %>%
  dplyr::group_by(
    family_id,
    network_taxon,
    unclassified_taxon
  ) %>%
  dplyr::summarise(
    max_group_prevalence = max(
      prevalence,
      na.rm = TRUE
    ),
    
    max_group_mean_abundance = max(
      mean_relative_abundance,
      na.rm = TRUE
    ),
    
    retain_family = any(
      prevalence >= min_group_prevalence &
        mean_relative_abundance >=
        min_group_mean_abundance
    ),
    
    .groups = "drop"
  )


retained_family_ids <- family_retention %>%
  dplyr::filter(
    retain_family
  ) %>%
  dplyr::pull(
    family_id
  )


# ==============================================================================
# 15. CREATE GROUP-WIDE FAMILY TABLES
# ==============================================================================

family_group_retained <- family_group_summary %>%
  dplyr::filter(
    family_id %in% retained_family_ids
  )


make_group_wide <- function(
    data,
    value_column,
    prefix
) {
  
  data %>%
    dplyr::select(
      family_id,
      network_taxon,
      unclassified_taxon,
      group_label,
      value = {{ value_column }}
    ) %>%
    tidyr::pivot_wider(
      names_from = group_label,
      values_from = value,
      values_fill = 0,
      names_prefix = prefix
    )
}


family_mean_wide <- make_group_wide(
  family_group_retained,
  mean_relative_abundance,
  "mean_abundance_"
)


family_prevalence_wide <- make_group_wide(
  family_group_retained,
  prevalence,
  "prevalence_"
)


family_count_wide <- make_group_wide(
  family_group_retained,
  total_count,
  "count_"
)


family_attributes <- family_mean_wide %>%
  dplyr::left_join(
    family_prevalence_wide,
    by = c(
      "family_id",
      "network_taxon",
      "unclassified_taxon"
    )
  ) %>%
  dplyr::left_join(
    family_count_wide,
    by = c(
      "family_id",
      "network_taxon",
      "unclassified_taxon"
    )
  ) %>%
  dplyr::left_join(
    family_retention,
    by = c(
      "family_id",
      "network_taxon",
      "unclassified_taxon"
    )
  ) %>%
  dplyr::left_join(
    overall_family_summary,
    by = c(
      "family_id",
      "network_taxon",
      "unclassified_taxon"
    )
  ) %>%
  dplyr::rename(
    mean_abundance_Wildlike =
      `mean_abundance_Wild-like`,
    
    mean_abundance_Lablike =
      `mean_abundance_Lab-like`,
    
    prevalence_Wildlike =
      `prevalence_Wild-like`,
    
    prevalence_Lablike =
      `prevalence_Lab-like`,
    
    count_Wildlike =
      `count_Wild-like`,
    
    count_Lablike =
      `count_Lab-like`
  )


# ==============================================================================
# 16. DETERMINE GROUP PRESENCE
#
# Presence is defined using BOTH prevalence and mean abundance.
# These variables drive:
#
#   - overlap_class
#   - visible edges
#   - pie-chart colours
# ==============================================================================

family_attributes <- family_attributes %>%
  dplyr::mutate(
    present_Wild =
      prevalence_Wild >=
      presence_prevalence_cutoff &
      mean_abundance_Wild >=
      presence_abundance_cutoff,
    
    present_Wildlike =
      prevalence_Wildlike >=
      presence_prevalence_cutoff &
      mean_abundance_Wildlike >=
      presence_abundance_cutoff,
    
    present_Lablike =
      prevalence_Lablike >=
      presence_prevalence_cutoff &
      mean_abundance_Lablike >=
      presence_abundance_cutoff,
    
    present_Lab =
      prevalence_Lab >=
      presence_prevalence_cutoff &
      mean_abundance_Lab >=
      presence_abundance_cutoff,
    
    n_groups_present = as.integer(
      present_Wild +
        present_Wildlike +
        present_Lablike +
        present_Lab
    )
  )


# ==============================================================================
# 17. CREATE THRESHOLD-CONSISTENT PIE VALUES
#
# A group contributes to a family pie only when the family is classified as
# present in that group.
#
# This guarantees:
#   Wild overlap class          -> Wild-only pie
#   Wild__Wildlike              -> Wild and Wild-like slices only
#   Shared by all four groups   -> up to four slices
# ==============================================================================

family_attributes <- family_attributes %>%
  dplyr::mutate(
    pie_abundance_Wild = dplyr::if_else(
      present_Wild,
      mean_abundance_Wild,
      0
    ),
    
    pie_abundance_Wildlike = dplyr::if_else(
      present_Wildlike,
      mean_abundance_Wildlike,
      0
    ),
    
    pie_abundance_Lablike = dplyr::if_else(
      present_Lablike,
      mean_abundance_Lablike,
      0
    ),
    
    pie_abundance_Lab = dplyr::if_else(
      present_Lab,
      mean_abundance_Lab,
      0
    ),
    
    pie_abundance_total =
      pie_abundance_Wild +
      pie_abundance_Wildlike +
      pie_abundance_Lablike +
      pie_abundance_Lab,
    
    pie_Wild = dplyr::if_else(
      pie_abundance_total > 0,
      pie_abundance_Wild /
        pie_abundance_total,
      0
    ),
    
    pie_Wildlike = dplyr::if_else(
      pie_abundance_total > 0,
      pie_abundance_Wildlike /
        pie_abundance_total,
      0
    ),
    
    pie_Lablike = dplyr::if_else(
      pie_abundance_total > 0,
      pie_abundance_Lablike /
        pie_abundance_total,
      0
    ),
    
    pie_Lab = dplyr::if_else(
      pie_abundance_total > 0,
      pie_abundance_Lab /
        pie_abundance_total,
      0
    )
  )


# ==============================================================================
# 18. DEFINE EXACT OVERLAP CLASSES
# ==============================================================================

family_attributes$overlap_class <- purrr::pmap_chr(
  list(
    family_attributes$present_Wild,
    family_attributes$present_Wildlike,
    family_attributes$present_Lablike,
    family_attributes$present_Lab
  ),
  
  function(
    wild,
    wildlike,
    lablike,
    lab
  ) {
    
    groups <- c(
      if (wild) "Wild",
      if (wildlike) "Wildlike",
      if (lablike) "Lablike",
      if (lab) "Lab"
    )
    
    
    if (length(groups) == 0) {
      return("Below_threshold")
    }
    
    
    paste(
      groups,
      collapse = "__"
    )
  }
)


family_attributes <- family_attributes %>%
  dplyr::mutate(
    overlap_display_class = dplyr::case_when(
      overlap_class == "Wild" ~
        "Wild only",
      
      overlap_class == "Wildlike" ~
        "Wild-like only",
      
      overlap_class == "Lablike" ~
        "Lab-like only",
      
      overlap_class == "Lab" ~
        "Lab only",
      
      overlap_class == "Wild__Wildlike" ~
        "Wild + Wild-like",
      
      overlap_class == "Wildlike__Lablike" ~
        "Wild-like + Lab-like",
      
      overlap_class == "Lablike__Lab" ~
        "Lab-like + Lab",
      
      overlap_class ==
        "Wild__Wildlike__Lablike" ~
        "Wild through Lab-like",
      
      overlap_class ==
        "Wildlike__Lablike__Lab" ~
        "Wild-like through Lab",
      
      overlap_class ==
        "Wild__Wildlike__Lablike__Lab" ~
        "Shared by all four",
      
      overlap_class == "Wild__Lablike" ~
        "Wild + Lab-like",
      
      overlap_class == "Wild__Lab" ~
        "Wild + Lab",
      
      overlap_class == "Wildlike__Lab" ~
        "Wild-like + Lab",
      
      overlap_class ==
        "Wild__Wildlike__Lab" ~
        "Wild + Wild-like + Lab",
      
      overlap_class ==
        "Wild__Lablike__Lab" ~
        "Wild + Lab-like + Lab",
      
      TRUE ~
        overlap_class
    )
  ) %>%
  dplyr::filter(
    overlap_class != "Below_threshold"
  )


# ==============================================================================
# 19. OPTIONALLY REMOVE UNCLASSIFIED NODES
# ==============================================================================

if (remove_unclassified_nodes) {
  
  family_attributes <- family_attributes %>%
    dplyr::filter(
      !unclassified_taxon
    )
}


# ==============================================================================
# 20. SCALE FAMILY NODE SIZE
#
# Cytoscape maps node_size_scaled directly to Node Size.
#
# We scale from sqrt(overall mean relative abundance), so node area is
# approximately proportional to abundance.
# ==============================================================================

size_input <- sqrt(
  family_attributes$overall_mean_relative_abundance
)


if (
  length(
    unique(size_input)
  ) <= 1
) {
  
  family_attributes$node_size_scaled <-
    mean(
      c(
        minimum_family_size,
        maximum_family_size
      )
    )
  
} else {
  
  family_attributes$node_size_scaled <-
    scales::rescale(
      size_input,
      to = c(
        minimum_family_size,
        maximum_family_size
      )
    )
}


family_attributes <- family_attributes %>%
  dplyr::mutate(
    log10_overall_mean_abundance =
      log10(
        overall_mean_relative_abundance +
          pseudocount
      )
  )


# ==============================================================================
# 21. OPTIONAL ANCOM-BC2 FAMILY ANNOTATIONS
#
# These summarize ASV-level ANCOM-BC2 results within each family.
# They are not formal family-level differential-abundance tests.
# ==============================================================================

safe_min <- function(x) {
  
  x <- x[
    is.finite(x) &
      !is.na(x)
  ]
  
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  
  min(x)
}


if (
  exists("da_fixed") &&
  is.data.frame(da_fixed$res_global)
) {
  
  family_ancom_summary <- da_fixed$res_global %>%
    dplyr::select(
      taxon,
      global_q_value = q_val,
      global_significant =
        diff_robust_abn
    ) %>%
    dplyr::left_join(
      taxonomy %>%
        dplyr::select(
          taxon,
          family_id
        ),
      by = "taxon"
    ) %>%
    dplyr::group_by(
      family_id
    ) %>%
    dplyr::summarise(
      n_ASVs_tested = dplyr::n(),
      
      n_global_significant_ASVs = sum(
        global_significant,
        na.rm = TRUE
      ),
      
      proportion_significant_ASVs =
        n_global_significant_ASVs /
        n_ASVs_tested,
      
      min_ASV_global_q =
        safe_min(global_q_value),
      
      .groups = "drop"
    )
  
  
  family_attributes <- family_attributes %>%
    dplyr::left_join(
      family_ancom_summary,
      by = "family_id"
    )
}


# Create optional annotation columns when da_fixed is unavailable.
if (!"n_ASVs_tested" %in% names(family_attributes)) {
  family_attributes$n_ASVs_tested <- NA_integer_
}

if (!"n_global_significant_ASVs" %in%
    names(family_attributes)) {
  family_attributes$n_global_significant_ASVs <-
    NA_integer_
}

if (!"proportion_significant_ASVs" %in%
    names(family_attributes)) {
  family_attributes$proportion_significant_ASVs <-
    NA_real_
}

if (!"min_ASV_global_q" %in%
    names(family_attributes)) {
  family_attributes$min_ASV_global_q <-
    NA_real_
}


# ==============================================================================
# 22. CREATE FAMILY NODE TABLE
# ==============================================================================

family_nodes <- family_attributes %>%
  dplyr::transmute(
    id = family_id,
    name = network_taxon,
    node_type = "Family",
    
    # Family nodes do not belong to one single group.
    group = "",
    
    unclassified_taxon,
    
    overlap_class,
    overlap_display_class,
    n_groups_present,
    
    node_size_scaled,
    
    overall_mean_relative_abundance,
    overall_median_relative_abundance,
    overall_prevalence,
    total_count_all_samples,
    log10_overall_mean_abundance,
    
    present_Wild,
    present_Wildlike,
    present_Lablike,
    present_Lab,
    
    # Threshold-consistent pie-chart values.
    pie_Wild,
    pie_Wildlike,
    pie_Lablike,
    pie_Lab,
    
    # Thresholded abundance values used to create the pies.
    pie_abundance_Wild,
    pie_abundance_Wildlike,
    pie_abundance_Lablike,
    pie_abundance_Lab,
    
    # Original unthresholded group summaries.
    prevalence_Wild,
    prevalence_Wildlike,
    prevalence_Lablike,
    prevalence_Lab,
    
    mean_abundance_Wild,
    mean_abundance_Wildlike,
    mean_abundance_Lablike,
    mean_abundance_Lab,
    
    n_ASVs_tested,
    n_global_significant_ASVs,
    proportion_significant_ASVs,
    min_ASV_global_q
  )


# ==============================================================================
# 23. CREATE FOUR GROUP NODES
# ==============================================================================

group_sample_counts <- metadata %>%
  dplyr::count(
    group_id,
    group_label,
    group_order,
    name = "n_samples"
  )


group_nodes <- group_sample_counts %>%
  dplyr::transmute(
    id = group_id,
    name = group_label,
    node_type = "Group",
    
    group = group_label,
    
    unclassified_taxon = FALSE,
    
    overlap_class = "",
    overlap_display_class = "",
    n_groups_present = 0L,
    
    node_size_scaled =
      as.numeric(group_node_size),
    
    overall_mean_relative_abundance = 0,
    overall_median_relative_abundance = 0,
    overall_prevalence = 0,
    total_count_all_samples = 0,
    log10_overall_mean_abundance = 0,
    
    present_Wild = FALSE,
    present_Wildlike = FALSE,
    present_Lablike = FALSE,
    present_Lab = FALSE,
    
    pie_Wild = 0,
    pie_Wildlike = 0,
    pie_Lablike = 0,
    pie_Lab = 0,
    
    pie_abundance_Wild = 0,
    pie_abundance_Wildlike = 0,
    pie_abundance_Lablike = 0,
    pie_abundance_Lab = 0,
    
    prevalence_Wild = 0,
    prevalence_Wildlike = 0,
    prevalence_Lablike = 0,
    prevalence_Lab = 0,
    
    mean_abundance_Wild = 0,
    mean_abundance_Wildlike = 0,
    mean_abundance_Lablike = 0,
    mean_abundance_Lab = 0,
    
    n_ASVs_tested = 0L,
    n_global_significant_ASVs = 0L,
    proportion_significant_ASVs = 0,
    min_ASV_global_q = 1
  )


nodes <- dplyr::bind_rows(
  group_nodes,
  family_nodes
)


# ==============================================================================
# 24. CREATE GROUP-TO-FAMILY EDGES
#
# An edge exists only when the family passes the same thresholds used for the
# overlap class and pie colours.
# ==============================================================================

edges <- family_group_retained %>%
  dplyr::filter(
    prevalence >=
      presence_prevalence_cutoff,
    
    mean_relative_abundance >=
      presence_abundance_cutoff,
    
    family_id %in%
      family_nodes$id
  ) %>%
  dplyr::transmute(
    source = group_id,
    target = family_id,
    
    interaction =
      "group_contains_family",
    
    group = group_label,
    group_order,
    
    mean_relative_abundance,
    median_relative_abundance,
    prevalence,
    total_count,
    
    sqrt_mean_relative_abundance =
      sqrt(mean_relative_abundance),
    
    log10_mean_relative_abundance =
      log10(
        mean_relative_abundance +
          pseudocount
      )
  )


# Scale edge width for direct passthrough mapping in Cytoscape.
if (
  nrow(edges) == 0
) {
  
  stop(
    "No group-to-family edges passed the selected thresholds."
  )
}


if (
  length(
    unique(
      edges$sqrt_mean_relative_abundance
    )
  ) <= 1
) {
  
  edges$edge_width_scaled <-
    mean(
      c(
        minimum_edge_width,
        maximum_edge_width
      )
    )
  
} else {
  
  edges$edge_width_scaled <-
    scales::rescale(
      edges$sqrt_mean_relative_abundance,
      to = c(
        minimum_edge_width,
        maximum_edge_width
      )
    )
}


# Cytoscape transparency uses approximately:
#   0   = transparent
#   255 = opaque
#
# Map prevalence onto a useful opacity range.
edges$edge_transparency_scaled <- scales::rescale(
  edges$prevalence,
  to = c(
    35,
    220
  )
)


# ==============================================================================
# 25. PIE AND MEMBERSHIP CONSISTENCY CHECKS
# ==============================================================================

pie_tolerance <- 1e-10


pie_membership_check <- family_nodes %>%
  dplyr::mutate(
    invalid_Wild =
      !present_Wild &
      pie_Wild > pie_tolerance,
    
    invalid_Wildlike =
      !present_Wildlike &
      pie_Wildlike > pie_tolerance,
    
    invalid_Lablike =
      !present_Lablike &
      pie_Lablike > pie_tolerance,
    
    invalid_Lab =
      !present_Lab &
      pie_Lab > pie_tolerance,
    
    pie_sum =
      pie_Wild +
      pie_Wildlike +
      pie_Lablike +
      pie_Lab,
    
    invalid_pie_sum =
      abs(pie_sum - 1) >
      pie_tolerance,
    
    any_problem =
      invalid_Wild |
      invalid_Wildlike |
      invalid_Lablike |
      invalid_Lab |
      invalid_pie_sum
  )


pie_membership_problems <- pie_membership_check %>%
  dplyr::filter(
    any_problem
  )


if (nrow(pie_membership_problems) > 0) {
  
  warning(
    nrow(pie_membership_problems),
    " family nodes failed the pie-membership consistency checks."
  )
  
} else {
  
  message(
    "PASS: Pie colours agree with overlap-class membership."
  )
}


# ==============================================================================
# 26. EDGE-MEMBERSHIP CONSISTENCY CHECK
# ==============================================================================

edge_membership_check <- edges %>%
  dplyr::left_join(
    family_nodes %>%
      dplyr::select(
        target = id,
        present_Wild,
        present_Wildlike,
        present_Lablike,
        present_Lab
      ),
    by = "target"
  ) %>%
  dplyr::mutate(
    edge_matches_membership = dplyr::case_when(
      group == "Wild" ~
        present_Wild,
      
      group == "Wild-like" ~
        present_Wildlike,
      
      group == "Lab-like" ~
        present_Lablike,
      
      group == "Lab" ~
        present_Lab,
      
      TRUE ~
        FALSE
    )
  )


if (any(!edge_membership_check$edge_matches_membership)) {
  
  stop(
    "At least one exported edge disagrees with the family membership columns."
  )
  
} else {
  
  message(
    "PASS: All edges agree with group-presence membership."
  )
}


# ==============================================================================
# 27. OVERLAP SUMMARY
# ==============================================================================

overlap_summary <- family_nodes %>%
  dplyr::count(
    overlap_class,
    overlap_display_class,
    n_groups_present,
    name = "n_families",
    sort = TRUE
  )


# ==============================================================================
# 28. EXPORT CYTOSCAPE FILES
# ==============================================================================

readr::write_csv(
  nodes,
  file.path(
    output_dir,
    "cytoscape_group_family_nodes.csv"
  ),
  na = ""
)


readr::write_csv(
  edges,
  file.path(
    output_dir,
    "cytoscape_group_family_edges.csv"
  ),
  na = ""
)


# ==============================================================================
# 29. EXPORT SUPPORTING FILES
# ==============================================================================

readr::write_csv(
  family_group_summary,
  file.path(
    output_dir,
    "family_group_summary.csv"
  ),
  na = ""
)


readr::write_csv(
  family_attributes,
  file.path(
    output_dir,
    "family_node_attributes.csv"
  ),
  na = ""
)


readr::write_csv(
  overlap_summary,
  file.path(
    output_dir,
    "overlap_summary.csv"
  ),
  na = ""
)


readr::write_csv(
  pie_membership_check,
  file.path(
    output_dir,
    "pie_membership_check.csv"
  ),
  na = ""
)


network_settings <- tibble::tibble(
  setting = c(
    "taxonomic_rank",
    "min_group_prevalence",
    "min_group_mean_abundance",
    "presence_prevalence_cutoff",
    "presence_abundance_cutoff",
    "collapse_unclassified",
    "remove_unclassified_nodes",
    "minimum_family_size",
    "maximum_family_size",
    "group_node_size",
    "pie_definition",
    "node_size_definition",
    "edge_definition"
  ),
  
  value = c(
    taxonomic_rank,
    min_group_prevalence,
    min_group_mean_abundance,
    presence_prevalence_cutoff,
    presence_abundance_cutoff,
    collapse_unclassified,
    remove_unclassified_nodes,
    minimum_family_size,
    maximum_family_size,
    group_node_size,
    paste(
      "Group mean abundance restricted to groups",
      "passing the presence thresholds"
    ),
    paste(
      "sqrt(overall mean relative abundance),",
      "rescaled to Cytoscape node diameter"
    ),
    paste(
      "Family passes prevalence and abundance",
      "thresholds in connected group"
    )
  )
)


readr::write_csv(
  network_settings,
  file.path(
    output_dir,
    "network_settings.csv"
  ),
  na = ""
)


# ==============================================================================
# 30. FINAL VALIDATION
# ==============================================================================

stopifnot(
  nrow(group_nodes) == 4,
  !anyDuplicated(nodes$id),
  all(edges$source %in% nodes$id),
  all(edges$target %in% nodes$id),
  all(is.finite(family_nodes$node_size_scaled)),
  all(is.finite(edges$edge_width_scaled)),
  all(is.finite(edges$edge_transparency_scaled))
)


# ==============================================================================
# 31. PRINT SUMMARY
# ==============================================================================

cat(
  "\nCYTOSCAPE GROUP-FAMILY EXPORT COMPLETE\n",
  "======================================\n",
  "Family nodes retained: ",
  nrow(family_nodes),
  "\n",
  "Group nodes: ",
  nrow(group_nodes),
  "\n",
  "Group-to-family edges: ",
  nrow(edges),
  "\n",
  "Overlap classes observed: ",
  nrow(overlap_summary),
  "\n",
  "Unclassified family nodes retained: ",
  sum(
    family_nodes$unclassified_taxon
  ),
  "\n",
  "Pie consistency problems: ",
  nrow(pie_membership_problems),
  "\n",
  "Output directory: ",
  output_dir,
  "\n\n",
  sep = ""
)


print(
  overlap_summary
)
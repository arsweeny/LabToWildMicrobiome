# ==============================================================================
# CORE MICROBIOME ANALYSIS ACROSS SAMPLE-ORIGIN GROUPS
# ==============================================================================
#
# Main questions:
#
#   1. Which ASVs occur in at least 90% of individuals within each group?
#
#   2. Which core ASVs are:
#        - shared among all four groups;
#        - shared among subsets of groups;
#        - unique to a particular group?
#
#   3. Does the degree of overlap change when taxa are examined at broader
#      taxonomic levels:
#        ASV -> Genus -> Family -> Order -> Class -> Phylum?
#
#   4. What is the genus-, family-, and higher-level composition of the
#      ASV-defined core microbiome of each group?
#
#   5. How abundant are core ASVs in each group?
#
#   6. How sensitive are conclusions to the chosen core-prevalence threshold?
#
#
# IMPORTANT DISTINCTION
# ---------------------
#
# This script performs two related analyses:
#
#   A. Independently defined cores at each taxonomic level
#
#      Counts are aggregated to each level first, and prevalence is then
#      calculated at that level.
#
#      Example:
#        A genus can be core even when no single ASV belonging to that genus
#        is core, because different ASVs from that genus occur in different
#        individuals.
#
#   B. Taxonomic composition of the ASV-defined cores
#
#      Core ASVs are identified first, and those ASVs are subsequently
#      classified by genus, family, order, class, and phylum.
#
#      This asks which broader bacterial groups contain the core ASVs.
#
#
# INPUT
# -----
#
# Required source script:
#
#   Scripts/00_DataPrep_Baseline.R
#
# Required phyloseq object:
#
#   ps_all
#
#
# OUTPUT
# ------
#
# The script writes:
#
#   Output/Core_microbiome/Tables/
#   Output/Core_microbiome/Graphs/
#
# including:
#
#   - prevalence and abundance tables for every taxonomic level;
#   - core-membership tables;
#   - core-size summaries;
#   - shared and group-specific core taxa;
#   - ASV prevalence heatmap;
#   - ASV abundance heatmap;
#   - UpSet plots at each taxonomic level;
#   - taxonomic-composition plots of ASV cores;
#   - core-size comparison across taxonomic levels;
#   - threshold-sensitivity plots.
#
# ==============================================================================


# ==============================================================================
# 1. CLEAN WORKSPACE AND LOAD PACKAGES
# ==============================================================================

rm(list = ls())

library(phyloseq)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(purrr)

library(ggplot2)
library(scales)
library(ComplexUpset)
library(patchwork)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

ps_all <- readRDS("Data/Processed/baseline_microbiome_phyloseq.rds")


stopifnot(
  exists("ps_all"),
  inherits(ps_all, "phyloseq")
)


# ==============================================================================
# 3. USER SETTINGS
# ==============================================================================

# ------------------------------------------------------------------------------
# Core definition
# ------------------------------------------------------------------------------

# Minimum number of reads required for a taxon to be considered detected in
# one sample.
detection_threshold <- 5


# Minimum proportion of individuals in which a taxon must be detected for it
# to be considered part of that group's core microbiome.
core_prevalence_threshold <- 0.90


# Additional thresholds used for the sensitivity analysis.
sensitivity_thresholds <- c(
  0.50,
  0.75,
  0.90,
  1.00
)


# ------------------------------------------------------------------------------
# Unit used to calculate prevalence
# ------------------------------------------------------------------------------

# Ideally this should identify the biological individual rather than the
# sequencing sample.
#
# The script will use mouse_ID_dataset if that column exists. Otherwise, each
# sample will be treated as an independent individual.
preferred_individual_column <- "mouse_ID_dataset"


# When an individual has more than one sample, the individual is considered
# positive when the taxon is detected in at least one of its samples.
#
# This is appropriate for the baseline dataset if each individual appears only
# once. If multiple biological timepoints are included, consider subsetting to
# one timepoint before this analysis.
individual_presence_rule <- "any"


# ------------------------------------------------------------------------------
# Groups and plotting order
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# Taxonomic levels
# ------------------------------------------------------------------------------

requested_taxonomic_levels <- c(
  "ASV",
  "Genus",
  "Family",
  "Order",
  "Class",
  "Phylum"
)


# ------------------------------------------------------------------------------
# Figure settings
# ------------------------------------------------------------------------------

# Maximum number of categories shown separately in taxonomic composition plots.
maximum_taxa_in_composition_plot <- 15


# Maximum number of ASV labels shown on heatmap axes.
#
# If there are more ASVs, labels are omitted while the complete heatmap is
# retained.
maximum_ASV_labels <- 80


# Small value added before log-transforming mean relative abundance for plotting.
abundance_pseudocount_percent <- 0.0001


# ------------------------------------------------------------------------------
# Output folders
# ------------------------------------------------------------------------------

output_directory <- "Results/04_Descriptive_Visualization/Core_Microbiome"

table_directory <- file.path(
  output_directory,
  "Tables"
)

graph_directory <- file.path(
  output_directory,
  "Figures"
)


dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  graph_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 4. EXTRACT AND VALIDATE METADATA
# ==============================================================================

metadata <- data.frame(
  phyloseq::sample_data(ps_all),
  check.names = FALSE
) %>%
  
  tibble::rownames_to_column(
    "sample_id"
  )


if (!"SampleOrigin" %in% names(metadata)) {
  
  stop(
    "The sample metadata does not contain a SampleOrigin column."
  )
  
}


metadata <- metadata %>%
  
  dplyr::mutate(
    
    SampleOrigin = factor(
      as.character(SampleOrigin),
      levels = cohort_order
    )
  )


if (any(is.na(metadata$SampleOrigin))) {
  
  invalid_groups <- unique(
    as.character(
      data.frame(
        phyloseq::sample_data(ps_all),
        check.names = FALSE
      )$SampleOrigin[
        is.na(metadata$SampleOrigin)
      ]
    )
  )
  
  stop(
    paste0(
      "Some SampleOrigin values do not match cohort_order: ",
      paste(invalid_groups, collapse = ", ")
    )
  )
  
}


missing_groups <- setdiff(
  cohort_order,
  unique(as.character(metadata$SampleOrigin))
)


if (length(missing_groups) > 0) {
  
  stop(
    paste0(
      "The following expected SampleOrigin groups are absent: ",
      paste(missing_groups, collapse = ", ")
    )
  )
  
}


# Use the preferred biological-individual column when available.
if (preferred_individual_column %in% names(metadata)) {
  
  individual_column <- preferred_individual_column
  
  message(
    "Calculating prevalence at the individual level using: ",
    individual_column
  )
  
} else {
  
  individual_column <- "sample_id"
  
  warning(
    paste0(
      preferred_individual_column,
      " was not found in the metadata. Each sample will therefore be treated ",
      "as a separate individual."
    )
  )
  
}


metadata <- metadata %>%
  
  dplyr::mutate(
    
    individual_id = as.character(
      .data[[individual_column]]
    )
  )


if (any(is.na(metadata$individual_id))) {
  
  stop(
    "At least one sample has a missing individual identifier."
  )
  
}


# Confirm that each individual belongs to only one SampleOrigin group.
individual_group_check <- metadata %>%
  
  dplyr::distinct(
    individual_id,
    SampleOrigin
  ) %>%
  
  dplyr::count(
    individual_id,
    name = "n_groups"
  ) %>%
  
  dplyr::filter(
    n_groups > 1
  )


if (nrow(individual_group_check) > 0) {
  
  stop(
    paste0(
      nrow(individual_group_check),
      " individual identifiers occur in more than one SampleOrigin group. ",
      "Confirm that individual_id uniquely identifies individuals across ",
      "datasets and groups."
    )
  )
  
}


metadata_core <- metadata %>%
  
  dplyr::select(
    sample_id,
    individual_id,
    SampleOrigin
  ) %>%
  
  dplyr::distinct()


# Summarize sample and individual numbers.
group_sample_sizes <- metadata_core %>%
  
  dplyr::group_by(
    SampleOrigin
  ) %>%
  
  dplyr::summarise(
    
    n_samples = dplyr::n_distinct(sample_id),
    
    n_individuals = dplyr::n_distinct(individual_id),
    
    .groups = "drop"
  )




print(group_sample_sizes)


# ==============================================================================
# 5. EXTRACT THE ASV COUNT MATRIX
#
# Rows must represent samples and columns must represent ASVs.
# ==============================================================================

otu_matrix <- as(
  phyloseq::otu_table(ps_all),
  "matrix"
)


if (phyloseq::taxa_are_rows(ps_all)) {
  
  otu_matrix <- t(otu_matrix)
  
}


storage.mode(otu_matrix) <- "numeric"


if (any(otu_matrix < 0, na.rm = TRUE)) {
  
  stop(
    "The ASV count matrix contains negative values."
  )
  
}


if (any(!is.finite(otu_matrix))) {
  
  stop(
    "The ASV count matrix contains non-finite values."
  )
  
}


samples_missing_metadata <- setdiff(
  rownames(otu_matrix),
  metadata_core$sample_id
)


metadata_missing_samples <- setdiff(
  metadata_core$sample_id,
  rownames(otu_matrix)
)


if (length(samples_missing_metadata) > 0) {
  
  stop(
    paste0(
      "Samples occur in the ASV table but not in metadata: ",
      paste(samples_missing_metadata, collapse = ", ")
    )
  )
  
}


if (length(metadata_missing_samples) > 0) {
  
  stop(
    paste0(
      "Samples occur in metadata but not in the ASV table: ",
      paste(metadata_missing_samples, collapse = ", ")
    )
  )
  
}


# Reorder metadata to match the ASV matrix.
metadata_core <- metadata_core %>%
  
  dplyr::slice(
    match(
      rownames(otu_matrix),
      sample_id
    )
  )


stopifnot(
  identical(
    metadata_core$sample_id,
    rownames(otu_matrix)
  )
)


# ==============================================================================
# 6. EXTRACT AND CLEAN TAXONOMY
# ==============================================================================

taxonomy <- as.data.frame(
  phyloseq::tax_table(ps_all),
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  
  tibble::rownames_to_column(
    "ASV"
  )


if (!setequal(taxonomy$ASV, colnames(otu_matrix))) {
  
  stop(
    "The taxonomy table and ASV count table do not contain the same ASV IDs."
  )
  
}


# Reorder taxonomy to match the count matrix.
taxonomy <- taxonomy %>%
  
  dplyr::slice(
    match(
      colnames(otu_matrix),
      ASV
    )
  )


standard_taxonomy_ranks <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)


available_taxonomy_ranks <- intersect(
  standard_taxonomy_ranks,
  names(taxonomy)
)


if (length(available_taxonomy_ranks) == 0) {
  
  stop(
    "No recognised taxonomic rank columns were found in the taxonomy table."
  )
  
}


# Standardize empty and uninformative taxonomy labels.
unknown_taxonomy_terms <- c(
  "",
  "na",
  "n/a",
  "unknown",
  "unclassified",
  "uncultured",
  "uncultured bacterium",
  "unidentified",
  "unassigned",
  "none"
)


for (rank_name in available_taxonomy_ranks) {
  
  raw_values <- as.character(
    taxonomy[[rank_name]]
  )
  
  cleaned_values <- stringr::str_trim(
    raw_values
  )
  
  unknown_values <- is.na(cleaned_values) |
    stringr::str_to_lower(cleaned_values) %in% unknown_taxonomy_terms
  
  cleaned_values[unknown_values] <- paste0(
    "Unclassified_",
    rank_name
  )
  
  taxonomy[[rank_name]] <- cleaned_values
  
}


# Determine which requested levels can be analysed.
taxonomic_levels <- requested_taxonomic_levels[
  requested_taxonomic_levels == "ASV" |
    requested_taxonomic_levels %in% available_taxonomy_ranks
]


missing_requested_levels <- setdiff(
  requested_taxonomic_levels,
  taxonomic_levels
)


if (length(missing_requested_levels) > 0) {
  
  warning(
    paste0(
      "The following requested taxonomic levels are unavailable and will be ",
      "skipped: ",
      paste(missing_requested_levels, collapse = ", ")
    )
  )
  
}


# ==============================================================================
# 7. HELPER FUNCTION: CREATE A UNIQUE TAXON KEY
#
# A hierarchical taxon key avoids incorrectly combining all unclassified taxa
# into one category.
#
# For example, an unclassified genus is distinguished according to its known
# phylum, class, order, and family.
# ==============================================================================

create_taxon_lookup <- function(
    taxonomy_table,
    taxonomic_level
) {
  
  if (taxonomic_level == "ASV") {
    
    lookup <- taxonomy_table %>%
      
      dplyr::transmute(
        
        ASV = ASV,
        
        taxon_id = ASV,
        
        taxon_label = ASV,
        
        taxonomic_level = "ASV"
      )
    
    return(lookup)
    
  }
  
  
  if (!taxonomic_level %in% names(taxonomy_table)) {
    
    stop(
      paste0(
        "Taxonomic level ",
        taxonomic_level,
        " is absent from the taxonomy table."
      )
    )
    
  }
  
  
  rank_position <- match(
    taxonomic_level,
    standard_taxonomy_ranks
  )
  
  
  ranks_for_key <- standard_taxonomy_ranks[
    seq_len(rank_position)
  ]
  
  
  ranks_for_key <- intersect(
    ranks_for_key,
    names(taxonomy_table)
  )
  
  
  lookup <- taxonomy_table
  
  
  lookup$taxon_id <- apply(
    
    lookup[
      ,
      ranks_for_key,
      drop = FALSE
    ],
    
    MARGIN = 1,
    
    FUN = function(x) {
      
      paste(
        x,
        collapse = "|"
      )
      
    }
  )
  
  
  lookup$taxon_label <- lookup[[taxonomic_level]]
  
  lookup$taxonomic_level <- taxonomic_level
  
  
  lookup %>%
    
    dplyr::select(
      ASV,
      taxon_id,
      taxon_label,
      taxonomic_level,
      dplyr::all_of(available_taxonomy_ranks)
    )
  
}


# ==============================================================================
# 8. HELPER FUNCTION: AGGREGATE ASV COUNTS TO A TAXONOMIC LEVEL
# ==============================================================================

aggregate_counts_to_level <- function(
    asv_count_matrix,
    taxonomy_table,
    taxonomic_level
) {
  
  lookup <- create_taxon_lookup(
    taxonomy_table = taxonomy_table,
    taxonomic_level = taxonomic_level
  )
  
  
  taxon_order <- unique(
    lookup$taxon_id
  )
  
  
  # rowsum() aggregates rows, so transpose the matrix temporarily.
  aggregated_taxa_by_sample <- rowsum(
    
    x = t(asv_count_matrix),
    
    group = lookup$taxon_id,
    
    reorder = FALSE
  )
  
  
  aggregated_sample_by_taxon <- t(
    aggregated_taxa_by_sample
  )
  
  
  # Restore a predictable taxon order.
  aggregated_sample_by_taxon <-
    aggregated_sample_by_taxon[
      ,
      taxon_order,
      drop = FALSE
    ]
  
  
  taxon_lookup <- lookup %>%
    
    dplyr::distinct(
      taxon_id,
      taxon_label,
      taxonomic_level
    )
  
  
  if (
    ncol(aggregated_sample_by_taxon) !=
    nrow(taxon_lookup)
  ) {
    
    stop(
      paste0(
        "Taxon aggregation failed at the ",
        taxonomic_level,
        " level."
      )
    )
    
  }
  
  
  list(
    
    count_matrix = aggregated_sample_by_taxon,
    
    taxon_lookup = taxon_lookup,
    
    asv_to_taxon = lookup
  )
  
}


# ==============================================================================
# 9. HELPER FUNCTION: CALCULATE PREVALENCE AND ABUNDANCE
# ==============================================================================

calculate_core_statistics <- function(
    count_matrix,
    taxon_lookup,
    metadata_table,
    taxonomic_level,
    detection_threshold,
    core_threshold
) {
  
  if (!identical(
    rownames(count_matrix),
    metadata_table$sample_id
  )) {
    
    stop(
      paste0(
        "Sample order mismatch while calculating ",
        taxonomic_level,
        " statistics."
      )
    )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Sample-level detection
  # ---------------------------------------------------------------------------
  
  sample_presence_matrix <- count_matrix >= detection_threshold
  
  
  sample_presence_long <- as.data.frame(
    sample_presence_matrix,
    check.names = FALSE
  ) %>%
    
    tibble::rownames_to_column(
      "sample_id"
    ) %>%
    
    tidyr::pivot_longer(
      
      cols = -sample_id,
      
      names_to = "taxon_id",
      
      values_to = "sample_present"
    ) %>%
    
    dplyr::left_join(
      metadata_table,
      by = "sample_id"
    )
  
  
  # ---------------------------------------------------------------------------
  # Collapse sample-level detection to individual-level detection
  # ---------------------------------------------------------------------------
  
  if (individual_presence_rule == "any") {
    
    individual_presence <- sample_presence_long %>%
      
      dplyr::group_by(
        SampleOrigin,
        individual_id,
        taxon_id
      ) %>%
      
      dplyr::summarise(
        
        present = any(
          sample_present,
          na.rm = TRUE
        ),
        
        n_samples_for_individual = dplyr::n(),
        
        .groups = "drop"
      )
    
  } else {
    
    stop(
      paste0(
        "Unsupported individual_presence_rule: ",
        individual_presence_rule
      )
    )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Calculate prevalence within each group
  # ---------------------------------------------------------------------------
  
  prevalence_table <- individual_presence %>%
    
    dplyr::group_by(
      SampleOrigin,
      taxon_id
    ) %>%
    
    dplyr::summarise(
      
      n_individuals = dplyr::n(),
      
      n_individuals_present = sum(
        present,
        na.rm = TRUE
      ),
      
      prevalence =
        n_individuals_present /
        n_individuals,
      
      .groups = "drop"
    ) %>%
    
    dplyr::left_join(
      taxon_lookup,
      by = "taxon_id"
    ) %>%
    
    dplyr::mutate(
      
      taxonomic_level = taxonomic_level,
      
      is_core =
        prevalence >= core_threshold,
      
      core_threshold = core_threshold,
      
      detection_threshold = detection_threshold
    )
  
  
  # ---------------------------------------------------------------------------
  # Calculate relative abundance
  #
  # Relative abundance is calculated using each sample's complete sequencing
  # library as the denominator.
  # ---------------------------------------------------------------------------
  
  full_library_sizes <- rowSums(
    otu_matrix,
    na.rm = TRUE
  )
  
  
  if (any(full_library_sizes <= 0)) {
    
    stop(
      "At least one sample has a sequencing library size of zero."
    )
    
  }
  
  
  relative_abundance_matrix <- sweep(
    
    count_matrix,
    
    MARGIN = 1,
    
    STATS = full_library_sizes,
    
    FUN = "/"
  )
  
  
  relative_abundance_matrix[
    !is.finite(relative_abundance_matrix)
  ] <- 0
  
  
  relative_abundance_long <- as.data.frame(
    relative_abundance_matrix,
    check.names = FALSE
  ) %>%
    
    tibble::rownames_to_column(
      "sample_id"
    ) %>%
    
    tidyr::pivot_longer(
      
      cols = -sample_id,
      
      names_to = "taxon_id",
      
      values_to = "relative_abundance"
    ) %>%
    
    dplyr::left_join(
      metadata_table,
      by = "sample_id"
    )
  
  
  # Average repeated samples within individuals so that each individual
  # contributes equally to the group mean.
  individual_abundance <- relative_abundance_long %>%
    
    dplyr::group_by(
      SampleOrigin,
      individual_id,
      taxon_id
    ) %>%
    
    dplyr::summarise(
      
      individual_mean_relative_abundance = mean(
        relative_abundance,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  
  abundance_table <- individual_abundance %>%
    
    dplyr::group_by(
      SampleOrigin,
      taxon_id
    ) %>%
    
    dplyr::summarise(
      
      n_individuals = dplyr::n(),
      
      mean_relative_abundance = mean(
        individual_mean_relative_abundance,
        na.rm = TRUE
      ),
      
      median_relative_abundance = stats::median(
        individual_mean_relative_abundance,
        na.rm = TRUE
      ),
      
      SD_relative_abundance = stats::sd(
        individual_mean_relative_abundance,
        na.rm = TRUE
      ),
      
      SE_relative_abundance =
        SD_relative_abundance /
        sqrt(n_individuals),
      
      .groups = "drop"
    ) %>%
    
    dplyr::left_join(
      taxon_lookup,
      by = "taxon_id"
    ) %>%
    
    dplyr::mutate(
      
      taxonomic_level = taxonomic_level,
      
      mean_relative_abundance_percent =
        100 * mean_relative_abundance,
      
      median_relative_abundance_percent =
        100 * median_relative_abundance,
      
      SE_relative_abundance_percent =
        100 * SE_relative_abundance
    )
  
  
  combined_table <- prevalence_table %>%
    
    dplyr::left_join(
      
      abundance_table %>%
        
        dplyr::select(
          SampleOrigin,
          taxon_id,
          mean_relative_abundance,
          median_relative_abundance,
          SD_relative_abundance,
          SE_relative_abundance,
          mean_relative_abundance_percent,
          median_relative_abundance_percent,
          SE_relative_abundance_percent
        ),
      
      by = c(
        "SampleOrigin",
        "taxon_id"
      )
    )
  
  
  list(
    
    prevalence = prevalence_table,
    
    abundance = abundance_table,
    
    combined = combined_table,
    
    individual_presence = individual_presence,
    
    individual_abundance = individual_abundance
  )
  
}


# ==============================================================================
# 10. RUN CORE ANALYSIS AT EACH TAXONOMIC LEVEL
# ==============================================================================

core_results <- list()

aggregation_results <- list()


for (current_level in taxonomic_levels) {
  
  message(
    "\nAnalysing taxonomic level: ",
    current_level
  )
  
  
  aggregation_results[[current_level]] <-
    aggregate_counts_to_level(
      
      asv_count_matrix = otu_matrix,
      
      taxonomy_table = taxonomy,
      
      taxonomic_level = current_level
    )
  
  
  core_results[[current_level]] <-
    calculate_core_statistics(
      
      count_matrix =
        aggregation_results[[current_level]]$count_matrix,
      
      taxon_lookup =
        aggregation_results[[current_level]]$taxon_lookup,
      
      metadata_table = metadata_core,
      
      taxonomic_level = current_level,
      
      detection_threshold = detection_threshold,
      
      core_threshold = core_prevalence_threshold
    )
  
  
  # Save level-specific tables.
  
  
  
}


# ==============================================================================
# 11. COMBINE RESULTS ACROSS TAXONOMIC LEVELS
# ==============================================================================

all_core_statistics <- purrr::map_dfr(
  
  core_results,
  
  function(x) {
    x$combined
  }
)




# ==============================================================================
# 12. CREATE WIDE CORE-MEMBERSHIP TABLES
#
# One row represents one taxon.
# One logical column represents whether that taxon is core in each group.
# ==============================================================================

create_core_membership_table <- function(
    core_statistics_table
) {
  
  prevalence_wide <- core_statistics_table %>%
    
    dplyr::select(
      taxon_id,
      taxon_label,
      taxonomic_level,
      SampleOrigin,
      prevalence
    ) %>%
    
    dplyr::mutate(
      
      SampleOrigin = as.character(
        SampleOrigin
      )
    ) %>%
    
    tidyr::pivot_wider(
      
      names_from = SampleOrigin,
      
      values_from = prevalence,
      
      names_prefix = "prevalence_",
      
      values_fill = 0
    )
  
  
  core_wide <- core_statistics_table %>%
    
    dplyr::select(
      taxon_id,
      SampleOrigin,
      is_core
    ) %>%
    
    dplyr::mutate(
      
      SampleOrigin = as.character(
        SampleOrigin
      )
    ) %>%
    
    tidyr::pivot_wider(
      
      names_from = SampleOrigin,
      
      values_from = is_core,
      
      values_fill = FALSE
    )
  
  
  membership <- prevalence_wide %>%
    
    dplyr::left_join(
      core_wide,
      by = "taxon_id"
    )
  
  
  membership %>%
    
    dplyr::mutate(
      
      n_groups_core = rowSums(
        dplyr::across(
          dplyr::all_of(cohort_order)
        )
      ),
      
      core_all_groups =
        n_groups_core == length(cohort_order),
      
      core_unique_to_one_group =
        n_groups_core == 1,
      
      core_pattern = apply(
        
        dplyr::select(
          .,
          dplyr::all_of(cohort_order)
        ),
        
        MARGIN = 1,
        
        FUN = function(x) {
          
          present_groups <- cohort_order[
            as.logical(x)
          ]
          
          if (length(present_groups) == 0) {
            
            return("Not core in any group")
            
          }
          
          paste(
            present_groups,
            collapse = " + "
          )
          
        }
      )
    )
  
}


core_membership_tables <- list()


for (current_level in taxonomic_levels) {
  
  core_membership_tables[[current_level]] <-
    create_core_membership_table(
      core_results[[current_level]]$combined
    )
  
  
  
}


# ==============================================================================
# 13. CORE-SIZE SUMMARY
# ==============================================================================

core_size_by_group <- all_core_statistics %>%
  
  dplyr::filter(
    is_core
  ) %>%
  
  dplyr::count(
    taxonomic_level,
    SampleOrigin,
    name = "n_core_taxa"
  ) %>%
  
  tidyr::complete(
    
    taxonomic_level = taxonomic_levels,
    
    SampleOrigin = factor(
      cohort_order,
      levels = cohort_order
    ),
    
    fill = list(
      n_core_taxa = 0
    )
  ) %>%
  
  dplyr::mutate(
    
    taxonomic_level = factor(
      taxonomic_level,
      levels = taxonomic_levels
    ),
    
    SampleOrigin = factor(
      SampleOrigin,
      levels = cohort_order
    )
  )


core_overlap_summary <- purrr::map_dfr(
  
  taxonomic_levels,
  
  function(current_level) {
    
    membership <-
      core_membership_tables[[current_level]]
    
    
    tibble::tibble(
      
      taxonomic_level = current_level,
      
      n_taxa_core_in_any_group = sum(
        membership$n_groups_core > 0
      ),
      
      n_taxa_core_in_all_groups = sum(
        membership$core_all_groups
      ),
      
      n_taxa_core_in_exactly_one_group = sum(
        membership$n_groups_core == 1
      ),
      
      n_taxa_core_in_exactly_two_groups = sum(
        membership$n_groups_core == 2
      ),
      
      n_taxa_core_in_exactly_three_groups = sum(
        membership$n_groups_core == 3
      )
    )
    
  }
)


unique_core_summary <- purrr::map_dfr(
  
  taxonomic_levels,
  
  function(current_level) {
    
    membership <-
      core_membership_tables[[current_level]]
    
    
    purrr::map_dfr(
      
      cohort_order,
      
      function(group_name) {
        
        other_groups <- setdiff(
          cohort_order,
          group_name
        )
        
        
        unique_to_group <-
          membership[[group_name]] &
          
          rowSums(
            membership[
              ,
              other_groups,
              drop = FALSE
            ]
          ) == 0
        
        
        tibble::tibble(
          
          taxonomic_level = current_level,
          
          SampleOrigin = group_name,
          
          n_unique_core_taxa = sum(
            unique_to_group
          )
        )
        
      }
    )
    
  }
)


core_summary_complete <- core_size_by_group %>%
  
  dplyr::left_join(
    
    unique_core_summary,
    
    by = c(
      "taxonomic_level",
      "SampleOrigin"
    )
  )


readr::write_csv(
  
  core_summary_complete,
  
  file.path(
    table_directory,
    "core_microbiome_size_by_group_and_taxonomic_level.csv"
  )
)


readr::write_csv(
  
  core_overlap_summary,
  
  file.path(
    table_directory,
    "core_microbiome_overlap_by_taxonomic_level.csv"
  )
)


# ==============================================================================
# 14. SAVE SHARED AND GROUP-SPECIFIC CORE TAXA
# ==============================================================================

for (current_level in taxonomic_levels) {
  
  membership <-
    core_membership_tables[[current_level]]
  
  
  shared_all_groups <- membership %>%
    
    dplyr::filter(
      core_all_groups
    )
  
  
  
  
  for (group_name in cohort_order) {
    
    other_groups <- setdiff(
      cohort_order,
      group_name
    )
    
    
    group_specific <- membership %>%
      
      dplyr::filter(
        
        .data[[group_name]],
        
        rowSums(
          dplyr::select(
            .,
            dplyr::all_of(other_groups)
          )
        ) == 0
      )
    
    
    safe_group_name <- stringr::str_replace_all(
      group_name,
      "[^A-Za-z0-9]+",
      "_"
    )
    
    
    
  }
  
}


# ==============================================================================
# 15. FIGURE: NUMBER OF CORE TAXA AT EACH TAXONOMIC LEVEL
# ==============================================================================

core_size_plot <- ggplot2::ggplot(
  
  core_size_by_group,
  
  ggplot2::aes(
    
    x = taxonomic_level,
    
    y = n_core_taxa,
    
    fill = SampleOrigin
  )
) +
  
  ggplot2::geom_col(
    
    position = ggplot2::position_dodge(
      width = 0.8
    ),
    
    width = 0.72
  ) +
  
  ggplot2::scale_fill_manual(
    values = cohort_colours
  ) +
  
  ggplot2::labs(
    
    title = "Core microbiome size across taxonomic levels",
    
    subtitle = paste0(
      "Core taxa were detected at ≥",
      detection_threshold,
      " reads in ≥",
      core_prevalence_threshold * 100,
      "% of individuals within a group"
    ),
    
    x = "Taxonomic level",
    
    y = "Number of core taxa",
    
    fill = "Group"
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    
    panel.grid.minor = ggplot2::element_blank(),
    
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    
    legend.position = "bottom"
  )


print(core_size_plot)


ggplot2::ggsave(
  
  filename = file.path(
    graph_directory,
    "core_microbiome_size_across_taxonomic_levels.pdf"
  ),
  
  plot = core_size_plot,
  
  width = 9,
  height = 6
)


# ==============================================================================
# 16. FIGURE: SHARED VERSUS GROUP-SPECIFIC CORE TAXA
# ==============================================================================

core_shared_plot_data <- core_overlap_summary %>%
  
  dplyr::select(
    
    taxonomic_level,
    
    n_taxa_core_in_all_groups,
    
    n_taxa_core_in_exactly_one_group,
    
    n_taxa_core_in_exactly_two_groups,
    
    n_taxa_core_in_exactly_three_groups
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = -taxonomic_level,
    
    names_to = "overlap_category",
    
    values_to = "n_taxa"
  ) %>%
  
  dplyr::mutate(
    
    taxonomic_level = factor(
      taxonomic_level,
      levels = taxonomic_levels
    ),
    
    overlap_category = dplyr::recode(
      
      overlap_category,
      
      "n_taxa_core_in_all_groups" =
        "Core in all four groups",
      
      "n_taxa_core_in_exactly_three_groups" =
        "Core in three groups",
      
      "n_taxa_core_in_exactly_two_groups" =
        "Core in two groups",
      
      "n_taxa_core_in_exactly_one_group" =
        "Core in one group"
    ),
    
    overlap_category = factor(
      
      overlap_category,
      
      levels = c(
        "Core in all four groups",
        "Core in three groups",
        "Core in two groups",
        "Core in one group"
      )
    )
  )


core_shared_plot <- ggplot2::ggplot(
  
  core_shared_plot_data,
  
  ggplot2::aes(
    
    x = taxonomic_level,
    
    y = n_taxa,
    
    fill = overlap_category
  )
) +
  
  ggplot2::geom_col() +
  
  ggplot2::labs(
    
    title = "Sharing of core taxa among microbiome groups",
    
    x = "Taxonomic level",
    
    y = "Number of core taxa",
    
    fill = "Core membership"
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    
    panel.grid.minor = ggplot2::element_blank(),
    
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    
    legend.position = "bottom"
  )


print(core_shared_plot)


ggplot2::ggsave(
  
  filename = file.path(
    graph_directory,
    "core_microbiome_overlap_across_taxonomic_levels.pdf"
  ),
  
  plot = core_shared_plot,
  
  width = 9,
  height = 6
)


# ==============================================================================
# 17. UPSET PLOTS OF CORE TAXA AT EACH TAXONOMIC LEVEL
# ==============================================================================

create_core_upset_plot <- function(
    membership_table,
    taxonomic_level
) {
  
  upset_input <- membership_table %>%
    
    dplyr::filter(
      n_groups_core > 0
    ) %>%
    
    dplyr::select(
      taxon_id,
      taxon_label,
      dplyr::all_of(cohort_order)
    )
  
  
  if (nrow(upset_input) == 0) {
    
    warning(
      paste0(
        "No core taxa were identified at the ",
        taxonomic_level,
        " level. The UpSet plot was skipped."
      )
    )
    
    return(NULL)
    
  }
  
  
  ComplexUpset::upset(
    
    data = upset_input,
    
    intersect = cohort_order,
    
    name = "Group",
    
    sort_intersections_by = "cardinality",
    
    sort_intersections = "descending",
    
    sort_sets = FALSE,
    
    min_size = 1,
    
    width_ratio = 0.20,
    
    base_annotations = list(
      
      "Intersection size" =
        
        ComplexUpset::intersection_size(
          
          counts = TRUE
        ) +
        
        ggplot2::labs(
          y = "Number of core taxa"
        ) +
        
        ggplot2::theme(
          
          panel.grid.major = ggplot2::element_blank(),
          
          panel.grid.minor = ggplot2::element_blank()
        )
    ),
    
    set_sizes = (
      
      ComplexUpset::upset_set_size(
        
        ggplot2::aes(
          fill = group
        )
      ) +
        
        ggplot2::scale_fill_manual(
          
          values = unname(
            cohort_colours[cohort_order]
          ),
          
          guide = "none"
        ) +
        
        ggplot2::labs(
          x = "Core-set size"
        )
    ),
    
    stripes = ComplexUpset::upset_stripes(
      
      colors = c(
        "white",
        "white"
      )
    )
  ) +
    
    ggplot2::ggtitle(
      
      paste0(
        taxonomic_level,
        "-level core microbiome"
      ),
      
      subtitle = paste0(
        "Core threshold: prevalence ≥ ",
        core_prevalence_threshold * 100,
        "%"
      )
    ) +
    
    ggplot2::theme(
      
      text = ggplot2::element_text(
        size = 11
      ),
      
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )
  
}


core_upset_plots <- list()


for (current_level in taxonomic_levels) {
  
  current_plot <- create_core_upset_plot(
    
    membership_table =
      core_membership_tables[[current_level]],
    
    taxonomic_level =
      current_level
  )
  
  
  core_upset_plots[[current_level]] <- current_plot
  
  
  if (!is.null(current_plot)) {
    
    print(current_plot)
    
    
    ggplot2::ggsave(
      
      filename = file.path(
        graph_directory,
        paste0(
          "Core_UpSet_",
          current_level,
          ".pdf"
        )
      ),
      
      plot = current_plot,
      
      width = 10,
      height = 6
    )
    
  }
  
}


# ==============================================================================
# 18. ASV CORE-PREVALENCE HEATMAP
#
# Includes every ASV that is core in at least one group.
# ==============================================================================

asv_membership <- core_membership_tables[["ASV"]]


core_ASV_ids <- asv_membership %>%
  
  dplyr::filter(
    n_groups_core > 0
  ) %>%
  
  dplyr::pull(
    taxon_id
  )


asv_prevalence_heatmap_data <-
  core_results[["ASV"]]$combined %>%
  
  dplyr::filter(
    taxon_id %in% core_ASV_ids
  ) %>%
  
  dplyr::left_join(
    
    asv_membership %>%
      
      dplyr::select(
        taxon_id,
        n_groups_core,
        core_pattern
      ),
    
    by = "taxon_id"
  ) %>%
  
  dplyr::group_by(
    taxon_id
  ) %>%
  
  dplyr::mutate(
    
    mean_prevalence_across_groups = mean(
      prevalence,
      na.rm = TRUE
    )
  ) %>%
  
  dplyr::ungroup()


asv_row_order <- asv_prevalence_heatmap_data %>%
  
  dplyr::distinct(
    taxon_id,
    n_groups_core,
    core_pattern,
    mean_prevalence_across_groups
  ) %>%
  
  dplyr::arrange(
    
    dplyr::desc(n_groups_core),
    
    core_pattern,
    
    dplyr::desc(mean_prevalence_across_groups),
    
    taxon_id
  ) %>%
  
  dplyr::pull(
    taxon_id
  )


asv_prevalence_heatmap_data <-
  asv_prevalence_heatmap_data %>%
  
  dplyr::mutate(
    
    taxon_id = factor(
      taxon_id,
      levels = rev(asv_row_order)
    ),
    
    SampleOrigin = factor(
      SampleOrigin,
      levels = cohort_order
    )
  )


show_ASV_axis_labels <-
  length(core_ASV_ids) <= maximum_ASV_labels


asv_prevalence_heatmap <- ggplot2::ggplot(
  
  asv_prevalence_heatmap_data,
  
  ggplot2::aes(
    
    x = SampleOrigin,
    
    y = taxon_id,
    
    fill = prevalence
  )
) +
  
  ggplot2::geom_tile(
    
    colour = "white",
    
    linewidth = 0.15
  ) +
  
  ggplot2::scale_fill_viridis_c(
    
    limits = c(0, 1),
    
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      0.90,
      1
    ),
    
    labels = scales::percent_format(
      accuracy = 1
    ),
    
    option = "C"
  ) +
  
  ggplot2::geom_vline(
    
    xintercept = seq(
      1.5,
      length(cohort_order) - 0.5,
      by = 1
    ),
    
    colour = "grey70",
    
    linewidth = 0.25
  ) +
  
  ggplot2::labs(
    
    title = "Prevalence of ASVs belonging to at least one group core",
    
    subtitle = paste0(
      "Core membership requires prevalence ≥",
      core_prevalence_threshold * 100,
      "%"
    ),
    
    x = NULL,
    
    y = "ASV",
    
    fill = "Prevalence"
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    
    panel.grid = ggplot2::element_blank(),
    
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    
    axis.text.y = if (show_ASV_axis_labels) {
      
      ggplot2::element_text(
        size = 5
      )
      
    } else {
      
      ggplot2::element_blank()
      
    },
    
    axis.ticks.y = if (show_ASV_axis_labels) {
      
      ggplot2::element_line()
      
    } else {
      
      ggplot2::element_blank()
      
    },
    
    plot.title = ggplot2::element_text(
      face = "bold"
    )
  )


asv_heatmap_height <- max(
  6,
  min(
    24,
    2.5 + 0.12 * length(core_ASV_ids)
  )
)


print(asv_prevalence_heatmap)


ggplot2::ggsave(
  
  filename = file.path(
    graph_directory,
    "core_ASV_prevalence_heatmap.pdf"
  ),
  
  plot = asv_prevalence_heatmap,
  
  width = 7,
  
  height = asv_heatmap_height,
  
  limitsize = FALSE
)


# ==============================================================================
# 19. ASV MEAN-ABUNDANCE HEATMAP
#
# Includes the union of all ASVs that are core in at least one group.
#
# Values represent group mean relative abundance. The fill is log10-transformed
# for readability, but the underlying values remain percentages.
# ==============================================================================

asv_abundance_heatmap_data <-
  core_results[["ASV"]]$combined %>%
  
  dplyr::filter(
    taxon_id %in% core_ASV_ids
  ) %>%
  
  dplyr::mutate(
    
    abundance_for_plot = log10(
      mean_relative_abundance_percent +
        abundance_pseudocount_percent
    ),
    
    taxon_id = factor(
      taxon_id,
      levels = rev(asv_row_order)
    ),
    
    SampleOrigin = factor(
      SampleOrigin,
      levels = cohort_order
    )
  )


asv_abundance_heatmap <- ggplot2::ggplot(
  
  asv_abundance_heatmap_data,
  
  ggplot2::aes(
    
    x = SampleOrigin,
    
    y = taxon_id,
    
    fill = abundance_for_plot
  )
) +
  
  ggplot2::geom_tile(
    
    colour = "white",
    
    linewidth = 0.15
  ) +
  
  ggplot2::scale_fill_viridis_c(
    
    option = "D",
    
    labels = function(x) {
      
      scales::label_number(
        accuracy = 0.001,
        suffix = "%"
      )(
        10^x
      )
      
    }
  ) +
  
  ggplot2::labs(
    
    title = "Mean relative abundance of core ASVs",
    
    subtitle = paste0(
      "ASVs shown are core in at least one group; fill is displayed on a ",
      "log10 scale"
    ),
    
    x = NULL,
    
    y = "ASV",
    
    fill = "Mean relative\nabundance"
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    
    panel.grid = ggplot2::element_blank(),
    
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    
    axis.text.y = if (show_ASV_axis_labels) {
      
      ggplot2::element_text(
        size = 5
      )
      
    } else {
      
      ggplot2::element_blank()
      
    },
    
    axis.ticks.y = if (show_ASV_axis_labels) {
      
      ggplot2::element_line()
      
    } else {
      
      ggplot2::element_blank()
      
    },
    
    plot.title = ggplot2::element_text(
      face = "bold"
    )
  )


print(asv_abundance_heatmap)


ggplot2::ggsave(
  
  filename = file.path(
    graph_directory,
    "core_ASV_mean_relative_abundance_heatmap.pdf"
  ),
  
  plot = asv_abundance_heatmap,
  
  width = 7,
  
  height = asv_heatmap_height,
  
  limitsize = FALSE
)


# ==============================================================================
# 20. TAXONOMIC COMPOSITION OF ASV-DEFINED CORE MICROBIOMES
#
# This section does NOT independently redefine the core at each broad
# taxonomic level.
#
# Instead:
#
#   1. Identify the ASVs in each group's 90% core.
#   2. Annotate those ASVs by Genus, Family, Order, Class, and Phylum.
#   3. Count how many core ASVs belong to each taxonomic category.
# ==============================================================================

asv_core_long <- core_results[["ASV"]]$combined %>%
  
  dplyr::filter(
    is_core
  ) %>%
  
  dplyr::select(
    
    SampleOrigin,
    
    ASV = taxon_id,
    
    prevalence,
    
    mean_relative_abundance,
    
    mean_relative_abundance_percent
  ) %>%
  
  dplyr::left_join(
    
    taxonomy,
    
    by = "ASV"
  )


composition_levels <- intersect(
  
  c(
    "Genus",
    "Family",
    "Order",
    "Class",
    "Phylum"
  ),
  
  available_taxonomy_ranks
)


ASV_core_composition_tables <- list()

ASV_core_composition_plots <- list()


for (current_level in composition_levels) {
  
  composition_table <- asv_core_long %>%
    
    dplyr::mutate(
      
      taxon_category = .data[[current_level]]
    ) %>%
    
    dplyr::group_by(
      
      SampleOrigin,
      
      taxon_category
    ) %>%
    
    dplyr::summarise(
      
      n_core_ASVs = dplyr::n_distinct(ASV),
      
      summed_mean_relative_abundance =
        sum(
          mean_relative_abundance,
          na.rm = TRUE
        ),
      
      summed_mean_relative_abundance_percent =
        100 * summed_mean_relative_abundance,
      
      .groups = "drop"
    ) %>%
    
    dplyr::group_by(
      SampleOrigin
    ) %>%
    
    dplyr::mutate(
      
      total_core_ASVs = sum(
        n_core_ASVs
      ),
      
      percentage_of_group_core_ASVs =
        100 *
        n_core_ASVs /
        total_core_ASVs
    ) %>%
    
    dplyr::ungroup()
  
  
  ASV_core_composition_tables[[current_level]] <-
    composition_table
  
  
  
  
  # Select the most represented categories across all groups.
  top_categories <- composition_table %>%
    
    dplyr::group_by(
      taxon_category
    ) %>%
    
    dplyr::summarise(
      
      n_core_ASVs_all_groups = sum(
        n_core_ASVs
      ),
      
      .groups = "drop"
    ) %>%
    
    dplyr::slice_max(
      
      order_by = n_core_ASVs_all_groups,
      
      n = maximum_taxa_in_composition_plot,
      
      with_ties = FALSE
    ) %>%
    
    dplyr::pull(
      taxon_category
    )
  
  
  composition_plot_data <- composition_table %>%
    
    dplyr::mutate(
      
      taxon_display = dplyr::if_else(
        
        taxon_category %in% top_categories,
        
        taxon_category,
        
        paste0(
          "Other ",
          stringr::str_to_lower(current_level),
          " taxa"
        )
      )
    ) %>%
    
    dplyr::group_by(
      
      SampleOrigin,
      
      taxon_display
    ) %>%
    
    dplyr::summarise(
      
      n_core_ASVs = sum(
        n_core_ASVs
      ),
      
      .groups = "drop"
    ) %>%
    
    dplyr::group_by(
      SampleOrigin
    ) %>%
    
    dplyr::mutate(
      
      percentage_of_group_core_ASVs =
        100 *
        n_core_ASVs /
        sum(n_core_ASVs)
    ) %>%
    
    dplyr::ungroup()
  
  
  category_order <- composition_plot_data %>%
    
    dplyr::group_by(
      taxon_display
    ) %>%
    
    dplyr::summarise(
      
      total_core_ASVs = sum(
        n_core_ASVs
      ),
      
      .groups = "drop"
    ) %>%
    
    dplyr::arrange(
      total_core_ASVs
    ) %>%
    
    dplyr::pull(
      taxon_display
    )
  
  
  composition_plot_data <- composition_plot_data %>%
    
    dplyr::mutate(
      
      SampleOrigin = factor(
        SampleOrigin,
        levels = cohort_order
      ),
      
      taxon_display = factor(
        taxon_display,
        levels = category_order
      )
    )
  
  
  composition_plot <- ggplot2::ggplot(
    
    composition_plot_data,
    
    ggplot2::aes(
      
      x = SampleOrigin,
      
      y = percentage_of_group_core_ASVs,
      
      fill = taxon_display
    )
  ) +
    
    ggplot2::geom_col(
      
      width = 0.75,
      
      colour = "white",
      
      linewidth = 0.15
    ) +
    
    ggplot2::scale_y_continuous(
      
      labels = function(x) {
        paste0(x, "%")
      },
      
      expand = ggplot2::expansion(
        mult = c(0, 0.02)
      )
    ) +
    
    ggplot2::labs(
      
      title = paste0(
        current_level,
        " composition of ASV-defined core microbiomes"
      ),
      
      subtitle = paste0(
        "Bars show the proportion of each group's core ASVs assigned to each ",
        stringr::str_to_lower(current_level)
      ),
      
      x = NULL,
      
      y = "Percentage of core ASVs",
      
      fill = current_level
    ) +
    
    ggplot2::theme_bw() +
    
    ggplot2::theme(
      
      panel.grid.minor = ggplot2::element_blank(),
      
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        face = "bold"
      ),
      
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      
      legend.position = "right"
    )
  
  
  ASV_core_composition_plots[[current_level]] <-
    composition_plot
  
  
  print(composition_plot)
  
  
  ggplot2::ggsave(
    
    filename = file.path(
      graph_directory,
      paste0(
        "ASV_core_composition_",
        current_level,
        ".pdf"
      )
    ),
    
    plot = composition_plot,
    
    width = 10,
    
    height = 6
  )
  
}


# ==============================================================================
# 21. HEATMAP: NUMBER OF CORE ASVs PER TAXONOMIC CATEGORY
#
# This shows whether a genus or family contributes many core ASVs in one group
# but few or none in another group.
# ==============================================================================

create_core_composition_heatmap <- function(
    composition_table,
    taxonomic_level,
    maximum_categories = 30
) {
  
  top_categories <- composition_table %>%
    
    dplyr::group_by(
      taxon_category
    ) %>%
    
    dplyr::summarise(
      
      total_core_ASVs = sum(
        n_core_ASVs
      ),
      
      .groups = "drop"
    ) %>%
    
    dplyr::slice_max(
      
      order_by = total_core_ASVs,
      
      n = maximum_categories,
      
      with_ties = FALSE
    ) %>%
    
    dplyr::arrange(
      total_core_ASVs
    ) %>%
    
    dplyr::pull(
      taxon_category
    )
  
  
  plot_data <- composition_table %>%
    
    dplyr::filter(
      taxon_category %in% top_categories
    ) %>%
    
    tidyr::complete(
      
      SampleOrigin = factor(
        cohort_order,
        levels = cohort_order
      ),
      
      taxon_category = top_categories,
      
      fill = list(
        n_core_ASVs = 0
      )
    ) %>%
    
    dplyr::mutate(
      
      SampleOrigin = factor(
        SampleOrigin,
        levels = cohort_order
      ),
      
      taxon_category = factor(
        taxon_category,
        levels = top_categories
      )
    )
  
  
  ggplot2::ggplot(
    
    plot_data,
    
    ggplot2::aes(
      
      x = SampleOrigin,
      
      y = taxon_category,
      
      fill = n_core_ASVs
    )
  ) +
    
    ggplot2::geom_tile(
      
      colour = "white",
      
      linewidth = 0.25
    ) +
    
    ggplot2::geom_text(
      
      ggplot2::aes(
        label = n_core_ASVs
      ),
      
      size = 3
    ) +
    
    ggplot2::scale_fill_viridis_c(
      
      option = "B",
      
      direction = -1
    ) +
    
    ggplot2::labs(
      
      title = paste0(
        taxonomic_level,
        " representation within ASV-defined cores"
      ),
      
      x = NULL,
      
      y = taxonomic_level,
      
      fill = "Number of\ncore ASVs"
    ) +
    
    ggplot2::theme_bw() +
    
    ggplot2::theme(
      
      panel.grid = ggplot2::element_blank(),
      
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        face = "bold"
      ),
      
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )
  
}


for (current_level in composition_levels) {
  
  composition_heatmap <- create_core_composition_heatmap(
    
    composition_table =
      ASV_core_composition_tables[[current_level]],
    
    taxonomic_level =
      current_level,
    
    maximum_categories = 30
  )
  
  
  print(composition_heatmap)
  
  
  ggplot2::ggsave(
    
    filename = file.path(
      graph_directory,
      paste0(
        "ASV_core_",
        current_level,
        "_heatmap.pdf"
      )
    ),
    
    plot = composition_heatmap,
    
    width = 8,
    
    height = 9
  )
  
}


# ==============================================================================
# 22. PREVALENCE-THRESHOLD SENSITIVITY ANALYSIS
#
# Reuses the prevalence estimates already calculated above.
#
# It is not necessary to rerun the detection and abundance calculations because
# only the threshold used to define the core changes.
# ==============================================================================

threshold_sensitivity <- purrr::map_dfr(
  
  sensitivity_thresholds,
  
  function(current_threshold) {
    
    all_core_statistics %>%
      
      dplyr::mutate(
        
        core_threshold_tested =
          current_threshold,
        
        is_core_at_tested_threshold =
          prevalence >= current_threshold
      ) %>%
      
      dplyr::group_by(
        
        taxonomic_level,
        
        SampleOrigin,
        
        core_threshold_tested
      ) %>%
      
      dplyr::summarise(
        
        n_core_taxa = sum(
          is_core_at_tested_threshold
        ),
        
        .groups = "drop"
      )
    
  }
) %>%
  
  dplyr::mutate(
    
    taxonomic_level = factor(
      taxonomic_level,
      levels = taxonomic_levels
    ),
    
    SampleOrigin = factor(
      SampleOrigin,
      levels = cohort_order
    )
  )


readr::write_csv(
  
  threshold_sensitivity,
  
  file.path(
    table_directory,
    "core_microbiome_threshold_sensitivity.csv"
  )
)


threshold_sensitivity_plot <- ggplot2::ggplot(
  
  threshold_sensitivity,
  
  ggplot2::aes(
    
    x = core_threshold_tested,
    
    y = n_core_taxa,
    
    colour = SampleOrigin,
    
    group = SampleOrigin
  )
) +
  
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  
  ggplot2::geom_point(
    size = 2
  ) +
  
  ggplot2::facet_wrap(
    
    facets = ggplot2::vars(
      taxonomic_level
    ),
    
    scales = "free_y"
  ) +
  
  ggplot2::scale_colour_manual(
    values = cohort_colours
  ) +
  
  ggplot2::scale_x_continuous(
    
    breaks = sensitivity_thresholds,
    
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  ggplot2::labs(
    
    title = "Sensitivity of core-microbiome size to prevalence threshold",
    
    subtitle = paste0(
      "Detection threshold remained fixed at ≥",
      detection_threshold,
      " reads per sample"
    ),
    
    x = "Minimum prevalence required for core membership",
    
    y = "Number of core taxa",
    
    colour = "Group"
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    
    panel.grid.minor = ggplot2::element_blank(),
    
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    
    legend.position = "bottom"
  )


print(threshold_sensitivity_plot)


ggplot2::ggsave(
  
  filename = file.path(
    graph_directory,
    "core_microbiome_threshold_sensitivity.pdf"
  ),
  
  plot = threshold_sensitivity_plot,
  
  width = 11,
  
  height = 7
)


# ==============================================================================
# 23. THRESHOLD SENSITIVITY OF FOUR-GROUP SHARING
#
# For each threshold and taxonomic level, calculate the number of taxa that are:
#
#   - core in at least one group;
#   - core in all four groups;
#   - unique to exactly one group.
# ==============================================================================

threshold_overlap_sensitivity <- purrr::map_dfr(
  
  taxonomic_levels,
  
  function(current_level) {
    
    current_statistics <-
      core_results[[current_level]]$combined
    
    
    purrr::map_dfr(
      
      sensitivity_thresholds,
      
      function(current_threshold) {
        
        threshold_membership <-
          current_statistics %>%
          
          dplyr::select(
            taxon_id,
            taxon_label,
            SampleOrigin,
            prevalence
          ) %>%
          
          dplyr::mutate(
            
            is_core =
              prevalence >= current_threshold,
            
            SampleOrigin =
              as.character(SampleOrigin)
          ) %>%
          
          dplyr::select(
            taxon_id,
            SampleOrigin,
            is_core
          ) %>%
          
          tidyr::pivot_wider(
            
            names_from = SampleOrigin,
            
            values_from = is_core,
            
            values_fill = FALSE
          )
        
        
        threshold_membership$n_groups_core <-
          rowSums(
            threshold_membership[
              ,
              cohort_order,
              drop = FALSE
            ]
          )
        
        
        tibble::tibble(
          
          taxonomic_level = current_level,
          
          core_threshold = current_threshold,
          
          n_core_in_any_group = sum(
            threshold_membership$n_groups_core > 0
          ),
          
          n_core_in_all_groups = sum(
            threshold_membership$n_groups_core ==
              length(cohort_order)
          ),
          
          n_core_unique_to_one_group = sum(
            threshold_membership$n_groups_core == 1
          )
        )
        
      }
    )
    
  }
)


readr::write_csv(
  
  threshold_overlap_sensitivity,
  
  file.path(
    table_directory,
    "core_microbiome_overlap_threshold_sensitivity.csv"
  )
)


threshold_overlap_plot_data <-
  threshold_overlap_sensitivity %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      n_core_in_any_group,
      n_core_in_all_groups,
      n_core_unique_to_one_group
    ),
    
    names_to = "sharing_category",
    
    values_to = "n_taxa"
  ) %>%
  
  dplyr::mutate(
    
    taxonomic_level = factor(
      taxonomic_level,
      levels = taxonomic_levels
    ),
    
    sharing_category = dplyr::recode(
      
      sharing_category,
      
      "n_core_in_any_group" =
        "Core in ≥1 group",
      
      "n_core_in_all_groups" =
        "Core in all groups",
      
      "n_core_unique_to_one_group" =
        "Unique to one group"
    )
  )


threshold_overlap_plot <- ggplot2::ggplot(
  
  threshold_overlap_plot_data,
  
  ggplot2::aes(
    
    x = core_threshold,
    
    y = n_taxa,
    
    linetype = sharing_category,
    
    shape = sharing_category,
    
    group = sharing_category
  )
) +
  
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  
  ggplot2::geom_point(
    size = 2
  ) +
  
  ggplot2::facet_wrap(
    
    facets = ggplot2::vars(
      taxonomic_level
    ),
    
    scales = "free_y"
  ) +
  
  ggplot2::scale_x_continuous(
    
    breaks = sensitivity_thresholds,
    
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  ggplot2::labs(
    
    title = "Sensitivity of core sharing to the prevalence threshold",
    
    x = "Minimum prevalence required for core membership",
    
    y = "Number of taxa",
    
    linetype = "Core category",
    
    shape = "Core category"
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    
    panel.grid.minor = ggplot2::element_blank(),
    
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    
    legend.position = "bottom"
  )


print(threshold_overlap_plot)


ggplot2::ggsave(
  
  filename = file.path(
    graph_directory,
    "core_microbiome_overlap_threshold_sensitivity.pdf"
  ),
  
  plot = threshold_overlap_plot,
  
  width = 11,
  
  height = 7
)


# ==============================================================================
# 24. HUMAN-READABLE SUMMARY TABLE
# ==============================================================================

human_readable_summary <- core_summary_complete %>%
  
  dplyr::left_join(
    
    core_overlap_summary %>%
      
      dplyr::select(
        taxonomic_level,
        n_taxa_core_in_any_group,
        n_taxa_core_in_all_groups,
        n_taxa_core_in_exactly_one_group
      ),
    
    by = "taxonomic_level"
  ) %>%
  
  dplyr::arrange(
    taxonomic_level,
    SampleOrigin
  )


readr::write_csv(
  
  human_readable_summary,
  
  file.path(
    table_directory,
    "core_microbiome_summary_for_interpretation.csv"
  )
)


print(human_readable_summary)


# ==============================================================================
# 25. SAVE ANALYSIS OBJECTS
#
# Saving the full objects makes it possible to recreate new figures without
# repeating the complete workflow.
# ==============================================================================

saveRDS(
  
  object = list(
    
    settings = list(
      
      detection_threshold =
        detection_threshold,
      
      core_prevalence_threshold =
        core_prevalence_threshold,
      
      individual_column =
        individual_column,
      
      individual_presence_rule =
        individual_presence_rule,
      
      cohort_order =
        cohort_order,
      
      taxonomic_levels =
        taxonomic_levels
    ),
    
    metadata =
      metadata_core,
    
    taxonomy =
      taxonomy,
    
    aggregation_results =
      aggregation_results,
    
    core_results =
      core_results,
    
    core_membership_tables =
      core_membership_tables,
    
    core_size_by_group =
      core_size_by_group,
    
    core_overlap_summary =
      core_overlap_summary,
    
    unique_core_summary =
      unique_core_summary,
    
    threshold_sensitivity =
      threshold_sensitivity,
    
    threshold_overlap_sensitivity =
      threshold_overlap_sensitivity
    
  ),
  
  file = file.path(
    output_directory,
    "Core_microbiome_analysis_objects.rds"
  )
)


# ==============================================================================
# 26. PRINT FINAL CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n",
  "============================================================\n",
  "CORE MICROBIOME ANALYSIS COMPLETE\n",
  "============================================================\n",
  "\n",
  "Detection threshold: ",
  detection_threshold,
  " reads\n",
  sep = ""
)


cat(
  "Core prevalence threshold: ",
  core_prevalence_threshold * 100,
  "%\n",
  sep = ""
)


cat(
  "Prevalence unit: ",
  individual_column,
  "\n",
  sep = ""
)


cat(
  "Taxonomic levels analysed: ",
  paste(
    taxonomic_levels,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat(
  "\nTables saved to:\n  ",
  table_directory,
  "\n",
  sep = ""
)


cat(
  "\nFigures saved to:\n  ",
  graph_directory,
  "\n\n",
  sep = ""
)


print(
  core_overlap_summary
)


cat(
  "\nASVs core in at least one group: ",
  sum(
    asv_membership$n_groups_core > 0
  ),
  "\n",
  sep = ""
)


cat(
  "ASVs core in all four groups: ",
  sum(
    asv_membership$core_all_groups
  ),
  "\n",
  sep = ""
)


for (group_name in cohort_order) {
  
  other_groups <- setdiff(
    cohort_order,
    group_name
  )
  
  
  n_unique <- sum(
    
    asv_membership[[group_name]] &
      
      rowSums(
        asv_membership[
          ,
          other_groups,
          drop = FALSE
        ]
      ) == 0
  )
  
  
  cat(
    "ASVs uniquely core in ",
    group_name,
    ": ",
    n_unique,
    "\n",
    sep = ""
  )
  
}

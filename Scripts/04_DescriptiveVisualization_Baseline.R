# ==============================================================================
# COMPLEXUPSET PLOT WITH FAMILY-LEVEL PIE-CHART INSET
#
# Interpretation:
#
#   UpSet plot
#     - ASV detection threshold: at least 5 reads in a sample.
#     - An ASV is present in a cohort when detected in at least one sample from
#       that cohort.
#     - The highlighted intersection bar represents ASVs detected in all four
#       cohorts.
#
#   Pie-chart inset
#     - Each pie represents a bacterial family containing shared ASVs.
#     - Pie area is proportional to the percentage of shared ASVs assigned to
#       that family.
#     - Pie slices show the distribution of that family's mean relative
#       abundance across the four cohorts.
#
# Required source script:
#   Scripts/00_DataPrep_Baseline.R
#
# Required phyloseq object:
#   ps_all
# ==============================================================================

rm(list = ls())

library(phyloseq)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)

library(ggplot2)
library(scales)
library(scatterpie)
library(patchwork)
library(ComplexUpset)
library(ggtext)
library(ggh4x)


# ==============================================================================
# 3. LOAD DATA
# ==============================================================================

source("Scripts/00_DataPrep_Baseline.R")


stopifnot(
  exists("ps_all"),
  inherits(ps_all, "phyloseq")
)

cohort_colours <- c(
  "Wild:As"      = "#74C476",
  "Wild-like:As" = "#6BAED6",
  "Lab-like:As"  = "#BC80BD",
  "Lab:Mm"       = "#FB8072"
)

# relative abundance visualization ----------------------------------------

## family level

ps_family <- ps_all %>% tax_fix(unknowns = c("Unclassified_Family"))

CFamily <- ps_family %>%
  tax_agg("Family") %>%
  ps_seriate(
    dist = "bray",
    method = "OLO_ward"
  ) %>%
  comp_barplot(
    tax_level = "Family",
    facet_by = "SampleOrigin",
    n_taxa = 15
  ) +
  
  ggh4x::facet_wrap2(
    ~SampleOrigin,
    scales = "free_x",
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = unname(cohort_colours),
        colour = NA
      ),
      text_x = ggh4x::elem_list_text(
        colour = "black",
        face = "bold",
        size = 11
      )
    )
  ) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

print(CFamily)


## genus level
ps_genus <- ps_all %>% tax_fix(unknowns = c("Unclassified_Genus"))

CGenus <- ps_genus %>%
  tax_agg("Genus") %>%
  ps_seriate(
    dist = "bray",
    method = "OLO_ward"
  ) %>%
  comp_barplot(
    tax_level = "Genus",
    facet_by = "SampleOrigin",
    n_taxa = 15
  ) +
  
  ggh4x::facet_wrap2(
    ~SampleOrigin,
    scales = "free_x",
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = unname(cohort_colours),
        colour = NA
      ),
      text_x = ggh4x::elem_list_text(
        colour = "black",
        face = "bold",
        size = 11
      )
    )
  ) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

print(CGenus)

# ==============================================================================
# 5. USER SETTINGS
# ==============================================================================

# Minimum count required for an ASV to be considered detected in a sample.
asv_detection_threshold <- 5


# Maximum number of classified families displayed separately.
#
# Remaining small families are pooled into "Other families".
# Unclassified ASVs remain a separate category.
maximum_families_displayed <- 10


# Size range for the family pies.
#
# Radius is scaled from the square root of the percentage of shared ASVs,
# meaning pie area is approximately proportional to ASV representation.
minimum_pie_radius <- 0.09
maximum_pie_radius <- 0.22


# Position of family pies and labels around the circle.
circle_layout_radius <- 1.00
label_layout_radius <- 1.38


# Final figure dimensions.
final_figure_width <- 12
final_figure_height <- 7


# ==============================================================================
# 6. COHORT ORDER AND COLOURS
# ==============================================================================

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


# Neutral colour for ordinary intersection bars.
default_intersection_colour <- "#737373"


# Colour used only for the four-cohort intersection bar.
#
# This is deliberately distinct from the four cohort colours.
shared_intersection_colour <- "#C44E52"


# ==============================================================================
# 7. EXTRACT AND VALIDATE SAMPLE METADATA
#
# sample_id is created here once and retained as a regular column.
# ==============================================================================

metadata <- data.frame(
  phyloseq::sample_data(ps_all),
  check.names = FALSE
) %>%
  
  tibble::rownames_to_column(
    "sample_id"
  ) %>%
  
  dplyr::mutate(
    
    SampleOrigin = factor(
      as.character(SampleOrigin),
      levels = cohort_order
    )
  )


if (any(is.na(metadata$SampleOrigin))) {
  
  invalid_groups <- unique(
    as.character(
      phyloseq::sample_data(ps_all)$SampleOrigin[
        is.na(metadata$SampleOrigin)
      ]
    )
  )
  
  stop(
    paste0(
      "Some SampleOrigin values did not match cohort_order: ",
      paste(invalid_groups, collapse = ", ")
    )
  )
}


missing_cohorts <- setdiff(
  cohort_order,
  unique(as.character(metadata$SampleOrigin))
)


if (length(missing_cohorts) > 0) {
  
  stop(
    paste0(
      "The following expected cohorts are absent: ",
      paste(missing_cohorts, collapse = ", ")
    )
  )
}


metadata_join <- metadata %>%
  
  dplyr::select(
    sample_id,
    SampleOrigin
  ) %>%
  
  dplyr::distinct()


duplicate_sample_ids <- metadata_join %>%
  
  dplyr::count(
    sample_id,
    name = "n_metadata_rows"
  ) %>%
  
  dplyr::filter(
    n_metadata_rows > 1
  )


if (nrow(duplicate_sample_ids) > 0) {
  
  stop(
    "At least one sample_id occurs more than once in the metadata."
  )
}


# ==============================================================================
# 8. EXTRACT ASV COUNT TABLE
#
# The resulting table must have:
#   rows    = samples
#   columns = ASVs
# ==============================================================================

otu <- as(
  phyloseq::otu_table(ps_all),
  "matrix"
)


if (phyloseq::taxa_are_rows(ps_all)) {
  
  otu <- t(otu)
  
}


otu <- as.data.frame(
  otu,
  check.names = FALSE
)


samples_missing_from_metadata <- setdiff(
  rownames(otu),
  metadata_join$sample_id
)


samples_missing_from_otu <- setdiff(
  metadata_join$sample_id,
  rownames(otu)
)


if (length(samples_missing_from_metadata) > 0) {
  
  stop(
    paste0(
      "Samples in the ASV table but not metadata: ",
      paste(samples_missing_from_metadata, collapse = ", ")
    )
  )
}


if (length(samples_missing_from_otu) > 0) {
  
  stop(
    paste0(
      "Samples in metadata but not the ASV table: ",
      paste(samples_missing_from_otu, collapse = ", ")
    )
  )
}


# Reorder metadata to match the ASV table.
metadata_join <- metadata_join %>%
  
  dplyr::slice(
    match(
      rownames(otu),
      sample_id
    )
  )


stopifnot(
  identical(
    metadata_join$sample_id,
    rownames(otu)
  )
)


# ==============================================================================
# 9. DEFINE ASV PRESENCE
# ==============================================================================

asv_presence <- as.matrix(otu) >=
  asv_detection_threshold


# ==============================================================================
# 10. CREATE ASV SETS FOR EACH COHORT
#
# An ASV is present in a cohort when it reaches the sample-level detection
# threshold in at least one sample from that cohort.
# ==============================================================================

asv_sets <- lapply(
  
  cohort_order,
  
  function(group_name) {
    
    group_samples <- metadata_join %>%
      
      dplyr::filter(
        SampleOrigin == group_name
      ) %>%
      
      dplyr::pull(
        sample_id
      )
    
    
    present_in_group <- colSums(
      
      asv_presence[
        group_samples,
        ,
        drop = FALSE
      ]
      
    ) > 0
    
    
    names(
      present_in_group[
        present_in_group
      ]
    )
  }
)


names(asv_sets) <- cohort_order


# ==============================================================================
# 11. IDENTIFY ASVs DETECTED IN ALL FOUR COHORTS
# ==============================================================================

shared_four_ASVs <- Reduce(
  intersect,
  asv_sets
)


n_shared_four <- length(
  shared_four_ASVs
)


cat(
  "\nASVs detected across all four cohorts: ",
  n_shared_four,
  "\n",
  sep = ""
)


if (n_shared_four != 94) {
  
  warning(
    paste0(
      "The current analysis identified ",
      n_shared_four,
      " shared ASVs rather than 94. ",
      "Confirm that the same samples, filtering and detection threshold were ",
      "used for the original UpSet analysis."
    )
  )
}


# ==============================================================================
# 12. CREATE COMPLEXUPSET INPUT TABLE
#
# Each row represents one ASV.
# Each cohort column indicates whether the ASV is present in that cohort.
# ==============================================================================

all_ASV_ids <- sort(
  unique(
    unlist(
      asv_sets,
      use.names = FALSE
    )
  )
)


upset_data <- tibble::tibble(
  ASV = all_ASV_ids
)


for (group_name in cohort_order) {
  
  upset_data[[group_name]] <-
    upset_data$ASV %in% asv_sets[[group_name]]
  
}


# Confirm the four-way intersection.
four_way_check <- upset_data %>%
  
  dplyr::filter(
    `Wild:As`,
    `Wild-like:As`,
    `Lab-like:As`,
    `Lab:Mm`
  )


stopifnot(
  nrow(four_way_check) ==
    length(shared_four_ASVs)
)


# ==============================================================================
# 13. EXTRACT AND CLEAN TAXONOMY
# ==============================================================================

taxonomy <- as.data.frame(
  phyloseq::tax_table(ps_all),
  stringsAsFactors = FALSE
) %>%
  
  tibble::rownames_to_column(
    "ASV"
  )


if (!"Family" %in% names(taxonomy)) {
  
  stop(
    "The taxonomy table does not contain a Family column."
  )
}


taxonomy_shared <- taxonomy %>%
  
  dplyr::filter(
    ASV %in% shared_four_ASVs
  ) %>%
  
  dplyr::mutate(
    
    Family = as.character(Family),
    
    Family = dplyr::case_when(
      
      is.na(Family) ~
        "Unclassified",
      
      stringr::str_trim(Family) == "" ~
        "Unclassified",
      
      stringr::str_to_lower(
        stringr::str_trim(Family)
      ) %in% c(
        "na",
        "unknown",
        "unclassified",
        "uncultured",
        "uncultured bacterium",
        "unidentified",
        "unassigned"
      ) ~
        "Unclassified",
      
      TRUE ~
        Family
    )
  )


missing_shared_taxonomy <- setdiff(
  shared_four_ASVs,
  taxonomy_shared$ASV
)


if (length(missing_shared_taxonomy) > 0) {
  
  stop(
    paste0(
      length(missing_shared_taxonomy),
      " shared ASVs were absent from the taxonomy table."
    )
  )
}


if (
  dplyr::n_distinct(taxonomy_shared$ASV) !=
  length(shared_four_ASVs)
) {
  
  stop(
    "The taxonomy table does not contain exactly one row per shared ASV."
  )
}


# ==============================================================================
# 14. COUNT SHARED ASVs PER FAMILY
# ==============================================================================

family_ASV_counts <- taxonomy_shared %>%
  
  dplyr::count(
    Family,
    name = "n_shared_ASVs"
  ) %>%
  
  dplyr::mutate(
    
    percentage_of_shared_ASVs =
      100 *
      n_shared_ASVs /
      sum(n_shared_ASVs)
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_shared_ASVs),
    Family
  )


stopifnot(
  sum(family_ASV_counts$n_shared_ASVs) ==
    length(shared_four_ASVs)
)


# ==============================================================================
# 15. POOL SMALL FAMILIES FOR DISPLAY
#
# The most represented classified families are displayed separately.
# Remaining classified families are pooled into "Other families".
# Unclassified ASVs remain separate.
# ==============================================================================

top_family_names <- family_ASV_counts %>%
  
  dplyr::filter(
    Family != "Unclassified"
  ) %>%
  
  dplyr::slice_max(
    order_by = n_shared_ASVs,
    n = maximum_families_displayed,
    with_ties = FALSE
  ) %>%
  
  dplyr::pull(
    Family
  )


taxonomy_shared <- taxonomy_shared %>%
  
  dplyr::mutate(
    
    Family_display = dplyr::case_when(
      
      Family == "Unclassified" ~
        "Unclassified",
      
      Family %in% top_family_names ~
        Family,
      
      TRUE ~
        "Other families"
    )
  )


family_display_counts <- taxonomy_shared %>%
  
  dplyr::count(
    Family_display,
    name = "n_shared_ASVs"
  ) %>%
  
  dplyr::mutate(
    
    percentage_of_shared_ASVs =
      100 *
      n_shared_ASVs /
      sum(n_shared_ASVs),
    
    family_label = paste0(
      Family_display,
      "\n",
      n_shared_ASVs,
      " ASV",
      dplyr::if_else(
        n_shared_ASVs == 1,
        "",
        "s"
      ),
      " (",
      round(
        percentage_of_shared_ASVs,
        1
      ),
      "%)"
    )
  )


# ==============================================================================
# 16. CALCULATE SAMPLE-LEVEL RELATIVE ABUNDANCE
#
# The denominator is each sample's full sequencing library, not only the shared
# ASVs. This preserves the original biological meaning of relative abundance.
# ==============================================================================

sample_library_sizes <- rowSums(
  otu,
  na.rm = TRUE
)


if (any(sample_library_sizes <= 0)) {
  
  stop(
    "At least one sample has a library size of zero."
  )
}


otu_relative_abundance <- sweep(
  
  as.matrix(otu),
  
  MARGIN = 1,
  
  STATS = sample_library_sizes,
  
  FUN = "/"
)


otu_relative_abundance[
  !is.finite(otu_relative_abundance)
] <- 0


# ==============================================================================
# 17. CONVERT SHARED-ASV ABUNDANCE TO LONG FORMAT
# ==============================================================================

shared_ASV_relative_abundance <- as.data.frame(
  
  otu_relative_abundance[
    ,
    shared_four_ASVs,
    drop = FALSE
  ],
  
  check.names = FALSE
) %>%
  
  tibble::rownames_to_column(
    "sample_id"
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = -sample_id,
    
    names_to = "ASV",
    
    values_to = "relative_abundance"
  ) %>%
  
  dplyr::left_join(
    
    taxonomy_shared %>%
      dplyr::select(
        ASV,
        Family,
        Family_display
      ),
    
    by = "ASV"
  ) %>%
  
  dplyr::left_join(
    
    metadata_join,
    
    by = "sample_id"
  )


if (any(is.na(shared_ASV_relative_abundance$Family_display))) {
  
  stop(
    "At least one shared ASV could not be assigned to a display family."
  )
}


if (any(is.na(shared_ASV_relative_abundance$SampleOrigin))) {
  
  stop(
    "At least one sample could not be assigned to a cohort."
  )
}


# ==============================================================================
# 18. SUM SHARED-ASV ABUNDANCE BY FAMILY WITHIN EACH SAMPLE
# ==============================================================================

shared_family_sample_abundance <- shared_ASV_relative_abundance %>%
  
  dplyr::group_by(
    sample_id,
    SampleOrigin,
    Family_display
  ) %>%
  
  dplyr::summarise(
    
    shared_family_relative_abundance = sum(
      relative_abundance,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# ==============================================================================
# 19. CREATE COMPLETE SAMPLE-BY-FAMILY TABLE
#
# Explicit zeros ensure that undetected families contribute zero to cohort means.
# ==============================================================================

shared_family_sample_complete <- tidyr::crossing(
  
  sample_id = metadata_join$sample_id,
  
  Family_display = unique(
    taxonomy_shared$Family_display
  )
  
) %>%
  
  dplyr::left_join(
    
    metadata_join,
    
    by = "sample_id"
  ) %>%
  
  dplyr::left_join(
    
    shared_family_sample_abundance,
    
    by = c(
      "sample_id",
      "SampleOrigin",
      "Family_display"
    )
  ) %>%
  
  dplyr::mutate(
    
    shared_family_relative_abundance =
      tidyr::replace_na(
        shared_family_relative_abundance,
        0
      ),
    
    SampleOrigin = factor(
      SampleOrigin,
      levels = cohort_order
    )
  )


# ==============================================================================
# 20. CALCULATE COHORT MEAN ABUNDANCE FOR EACH FAMILY
# ==============================================================================

family_cohort_abundance <- shared_family_sample_complete %>%
  
  dplyr::group_by(
    Family_display,
    SampleOrigin
  ) %>%
  
  dplyr::summarise(
    
    n_samples = dplyr::n(),
    
    mean_relative_abundance = mean(
      shared_family_relative_abundance,
      na.rm = TRUE
    ),
    
    median_relative_abundance = stats::median(
      shared_family_relative_abundance,
      na.rm = TRUE
    ),
    
    SD_relative_abundance = stats::sd(
      shared_family_relative_abundance,
      na.rm = TRUE
    ),
    
    SE_relative_abundance =
      SD_relative_abundance /
      sqrt(n_samples),
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    
    mean_relative_abundance_percent =
      100 * mean_relative_abundance,
    
    SE_relative_abundance_percent =
      100 * SE_relative_abundance
  )


# ==============================================================================
# 21. CALCULATE PIE FRACTIONS
#
# Within each family, cohort mean relative abundances are normalised to sum to 1.
#
# Thus:
#   - pie area represents the proportion of shared ASVs in the family;
#   - pie slices represent how that family's abundance is distributed among
#     cohorts.
# ==============================================================================

family_pie_values <- family_cohort_abundance %>%
  
  dplyr::group_by(
    Family_display
  ) %>%
  
  dplyr::mutate(
    
    total_mean_abundance_across_cohorts = sum(
      mean_relative_abundance,
      na.rm = TRUE
    ),
    
    pie_fraction = dplyr::if_else(
      
      total_mean_abundance_across_cohorts > 0,
      
      mean_relative_abundance /
        total_mean_abundance_across_cohorts,
      
      0
    )
  ) %>%
  
  dplyr::ungroup()


pie_sum_check <- family_pie_values %>%
  
  dplyr::group_by(
    Family_display
  ) %>%
  
  dplyr::summarise(
    
    pie_sum = sum(
      pie_fraction,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


if (any(abs(pie_sum_check$pie_sum - 1) > 1e-8)) {
  
  stop(
    "At least one family pie does not sum to one."
  )
}


# ==============================================================================
# 22. CONVERT PIE DATA TO WIDE FORMAT
# ==============================================================================

family_pie_wide <- family_pie_values %>%
  
  dplyr::mutate(
    
    cohort_column = dplyr::recode(
      
      as.character(SampleOrigin),
      
      "Wild:As" =
        "pie_Wild",
      
      "Wild-like:As" =
        "pie_Wildlike",
      
      "Lab-like:As" =
        "pie_Lablike",
      
      "Lab:Mm" =
        "pie_Lab"
    )
  ) %>%
  
  dplyr::select(
    Family_display,
    cohort_column,
    pie_fraction
  ) %>%
  
  tidyr::pivot_wider(
    
    names_from = cohort_column,
    
    values_from = pie_fraction,
    
    values_fill = 0
  ) %>%
  
  dplyr::left_join(
    
    family_display_counts,
    
    by = "Family_display"
  )


# ==============================================================================
# 23. SCALE PIE AREA BY SHARED-ASV REPRESENTATION
# ==============================================================================

radius_input <- sqrt(
  family_pie_wide$percentage_of_shared_ASVs
)


if (length(unique(radius_input)) <= 1) {
  
  family_pie_wide$pie_radius <-
    mean(
      c(
        minimum_pie_radius,
        maximum_pie_radius
      )
    )
  
} else {
  
  family_pie_wide$pie_radius <-
    scales::rescale(
      
      radius_input,
      
      to = c(
        minimum_pie_radius,
        maximum_pie_radius
      )
    )
}


# ==============================================================================
# 24. ARRANGE FAMILY PIES AROUND A CIRCLE
# ==============================================================================

family_pie_wide <- family_pie_wide %>%
  
  dplyr::arrange(
    dplyr::desc(n_shared_ASVs),
    Family_display
  ) %>%
  
  dplyr::mutate(
    
    family_order =
      dplyr::row_number(),
    
    # Begin at the top of the circle and continue clockwise.
    angle =
      pi / 2 -
      2 * pi *
      (family_order - 1) /
      dplyr::n(),
    
    x =
      circle_layout_radius *
      cos(angle),
    
    y =
      circle_layout_radius *
      sin(angle),
    
    label_x =
      label_layout_radius *
      cos(angle),
    
    label_y =
      label_layout_radius *
      sin(angle),
    
    label_hjust = dplyr::case_when(
      
      cos(angle) > 0.15 ~
        0,
      
      cos(angle) < -0.15 ~
        1,
      
      TRUE ~
        0.5
    )
  )


circle_guide <- tibble::tibble(
  
  angle = seq(
    0,
    2 * pi,
    length.out = 500
  ),
  
  x =
    circle_layout_radius *
    cos(angle),
  
  y =
    circle_layout_radius *
    sin(angle)
)


# ==============================================================================
# 25. CREATE FAMILY PIE-CHART INSET
#
# The legend is placed in the centre in a 2 × 2 arrangement.
# A white background and border make the panel read as an inset.
# ==============================================================================

family_pie_circle <- ggplot() +
  
  geom_path(
    
    data = circle_guide,
    
    aes(
      x = x,
      y = y
    ),
    
    colour = "grey70",
    
    linewidth = 0.25
  ) +
  
  scatterpie::geom_scatterpie(
    
    data = family_pie_wide,
    
    aes(
      x = x,
      y = y,
      r = pie_radius
    ),
    
    cols = c(
      "pie_Wild",
      "pie_Wildlike",
      "pie_Lablike",
      "pie_Lab"
    ),
    
    colour = "grey25",
    linewidth = 0.25
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 0,
    label = "Shared\nMicrobiome",
    fontface = "bold",
    size = 3.2,
    colour = "grey25",
    lineheight = 0.9
  ) +

  ## Family labels
  geom_text(
    
    data = family_pie_wide,
    
    aes(
      x = label_x,
      y = label_y,
      label = family_label,
      hjust = label_hjust
    ),
    
    size = 2.35,
    lineheight = 0.88
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      "pie_Wild" =
        cohort_colours[["Wild:As"]],
      
      "pie_Wildlike" =
        cohort_colours[["Wild-like:As"]],
      
      "pie_Lablike" =
        cohort_colours[["Lab-like:As"]],
      
      "pie_Lab" =
        cohort_colours[["Lab:Mm"]]
    ),
    
    guide = "none"   # removes legend
  ) +
  
  coord_equal(
    
    xlim = c(-1.72, 1.72),
    ylim = c(-1.64, 1.64),
    
    clip = "off"
  ) +
  
  theme_void() +
  
  theme(
    
    legend.position = "none",
    
    plot.background = element_blank(),
    panel.background = element_blank(),
    
    plot.margin = margin(
      5, 5, 5, 5
    )
  )

print(family_pie_circle)

# ==============================================================================
# COMPLEXUPSET PLOT
#
# Styling:
#   - Set-size bars use the same colours as the cohort pie-chart slices.
#   - Ordinary intersection bars are grey.
#   - The four-cohort shared intersection bar is red.
#   - The four dots and connecting line for that intersection are also red.
#   - Other active dots and lines remain black.
#   - Inactive dots remain grey.
#   - Plotting-area grid lines are removed.
# ==============================================================================

library(ComplexUpset)
library(ggplot2)
library(patchwork)


# ==============================================================================
# COMPLEXUPSET PLOT — LARGER PUBLICATION-READY TEXT
#
# Assumes these objects already exist:
#   upset_data
#   cohort_order
#   cohort_colours
#
# Styling:
#   - cohort-coloured set-size bars;
#   - ordinary intersection bars grey;
#   - four-cohort shared intersection bar, dots and line red;
#   - other active dots black and inactive dots grey;
#   - no panel grid lines;
#   - larger axis titles, axis numbers, cohort labels and bar-count labels.
# ==============================================================================

library(ComplexUpset)
library(ggplot2)


# ==============================================================================
# 1. COHORT ORDER AND COLOURS
# ==============================================================================

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


# ComplexUpset internally indexes the set-size bars numerically.
set_size_colours <- unname(
  cohort_colours[cohort_order]
)


# ==============================================================================
# 2. INTERSECTION COLOURS
# ==============================================================================

default_intersection_colour <- "#737373"

shared_intersection_colour <- "#C44E52"


# ==============================================================================
# 3. TEXT-SIZE SETTINGS
#
# Adjust these values in one place if the exported figure is still too small.
# ==============================================================================

overall_text_size <- 12

axis_title_size <- 13

axis_number_size <- 11

cohort_label_size <- 13

intersection_count_size <- 4.3


# ==============================================================================
# 4. CREATE COMPLEXUPSET PLOT
# ==============================================================================

complex_upset_plot <- ComplexUpset::upset(
  
  data = upset_data,
  
  intersect = cohort_order,
  
  name = "Cohort",
  
  sort_intersections_by = "cardinality",
  
  sort_intersections = "descending",
  
  sort_sets = FALSE,
  
  keep_empty_groups = FALSE,
  
  min_size = 1,
  
  width_ratio = 0.18,
  
  
  # ============================================================================
  # INTERSECTION-SIZE BARS
  # ============================================================================
  
  base_annotations = list(
    
    "Intersection size" =
      
      ComplexUpset::intersection_size(
        
        counts = TRUE,
        
        # Neutral grey for ordinary intersection bars.
        mapping = aes(
          fill = I(default_intersection_colour)
        ),
        
        # Text printed above each intersection bar.
        text = list(
          size = intersection_count_size,
          vjust = -0.25
        )
      ) +
      
      labs(
        y = "Intersection size"
      ) +
      
      theme(
        
        panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        axis.title.y = element_text(
          size = axis_title_size,
          face = "plain"
        ),
        
        axis.text.y = element_text(
          size = axis_number_size
        ),
        
        axis.text.x = element_blank(),
        
        axis.ticks.x = element_blank()
      )
  ),
  
  
  # ============================================================================
  # SET-SIZE BARS
  # ============================================================================
  
  set_sizes = (
    
    ComplexUpset::upset_set_size(
      
      mapping = aes(
        fill = group
      )
    ) +
      
      scale_fill_manual(
        
        values = set_size_colours,
        
        guide = "none"
      ) +
      
      labs(
        x = "Set size"
      ) +
      
      theme(
        
        panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        axis.title.x = element_text(
          size = axis_title_size,
          face = "plain"
        ),
        
        axis.text.x = element_text(
          size = axis_number_size
        ),
        
        axis.text.y = element_blank(),
        
        axis.ticks.y = element_blank()
      )
  ),
  
  
  # ============================================================================
  # HIGHLIGHT THE FOUR-COHORT SHARED INTERSECTION
  #
  # First query:
  #   colours the intersection-size bar.
  #
  # Second query:
  #   colours the corresponding matrix dots and connecting line.
  #
  # `color`, rather than `colour`, is required by upset_query().
  # ============================================================================
  
  queries = list(
    
    ComplexUpset::upset_query(
      
      intersect = cohort_order,
      
      fill = shared_intersection_colour,
      
      color = shared_intersection_colour,
      
      only_components = "Intersection size"
    ),
    
    
    ComplexUpset::upset_query(
      
      intersect = cohort_order,
      
      fill = shared_intersection_colour,
      
      color = shared_intersection_colour,
      
      only_components = "intersections_matrix"
    )
  ),
  
  
  # ============================================================================
  # PANEL-SPECIFIC THEMES
  # ============================================================================
  
  themes = ComplexUpset::upset_modify_themes(
    
    list(
      
      # Set-size panel
      "overall_sizes" = theme(
        
        panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        axis.title.x = element_text(
          size = axis_title_size
        ),
        
        axis.text.x = element_text(
          size = axis_number_size
        )
      ),
      
      
      # Matrix panel containing cohort names and intersection dots
      "intersections_matrix" = theme(
        
        panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        # Wild:As, Wild-like:As, Lab-like:As and Lab:Mm
        axis.text.y = element_text(
          size = cohort_label_size,
          colour = "black"
        ),
        
        axis.title.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.ticks.x = element_blank()
      ),
      
      
      # Intersection-size panel
      "Intersection size" = theme(
        
        panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        axis.title.y = element_text(
          size = axis_title_size
        ),
        
        axis.text.y = element_text(
          size = axis_number_size
        )
      )
    )
  ),
  
  
  # ============================================================================
  # MATRIX BACKGROUND
  # ============================================================================
  
  stripes = ComplexUpset::upset_stripes(
    
    colors = c(
      "white",
      "white"
    )
  )
  
  # No custom matrix layer is supplied.
  #
  # This preserves:
  #   ordinary active dots and lines = black
  #   inactive dots                  = grey
  #
  # The shared four-way column is recoloured by the query above.
) +
  
  theme(
    
    # Default text size inherited by elements without a specific override.
    text = element_text(
      size = overall_text_size,
      colour = "black"
    ),
    
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
  )


# ==============================================================================
# 5. DISPLAY
# ==============================================================================

print(complex_upset_plot)

# ==============================================================================
# 27. ADD PIE GRAPH AS AN INSET
# ==============================================================================

# Wrap the full ComplexUpset composition as one patchwork object
complex_upset_wrapped <- patchwork::wrap_elements(
  full = complex_upset_plot
)

# Add inset relative to the full assembled UpSet figure
complex_upset_with_inset <- complex_upset_wrapped +
  
  patchwork::inset_element(
    family_pie_circle,
    
    left   = 0.48,
    bottom = 0.48,
    right  = 0.99,
    top    = 0.99,
    
    align_to = "full",
    clip = FALSE
  )

print(complex_upset_with_inset)

# ==============================================================================
# COMBINE GENUS COMPOSITION PLOT AND COMPLEXUPSET PLOT
#
# Final layout:
#   A = genus-level relative abundance plot
#   B = ComplexUpset plot with the family pie chart embedded inside it
#
# The pie-chart inset does not receive its own panel letter.
# ==============================================================================

library(cowplot)
library(ggplot2)


# ==============================================================================
# 1. OPTIONAL: CLEAN UP THE GENUS PLOT FOR THE COMBINED FIGURE
# ==============================================================================

CFamily_final <- CFamily +
  
  theme(
    
    # Remove sample labels and ticks along the x-axis
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Keep text readable in the combined figure
    axis.title.y = element_text(
      size = 13
    ),
    
    axis.text.y = element_text(
      size = 11
    ),
    
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ==============================================================================
# 2. CREATE THE FINAL TWO-PANEL FIGURE
#
# CGenus occupies approximately 42% of the width.
# ComplexUpset occupies approximately 58%.
#
# The inset is already part of complex_upset_with_inset, so only the complete
# UpSet panel receives the label "B".
# ==============================================================================

final_combined_plot <- cowplot::ggdraw() +
  
  # ---------------------------------------------------------------------------
# Panel A: genus-level composition
# ---------------------------------------------------------------------------

cowplot::draw_plot(
  
  CFamily_final,
  
  x = 0.00,
  y = 0.00,
  
  width = 0.5,
  height = 1.00
) +
  
  # ---------------------------------------------------------------------------
# Panel B: ComplexUpset plus embedded family-pie inset
# ---------------------------------------------------------------------------

draw_plot(
  complex_upset_with_inset,
  x = 0.485,
  y = 0,
  width = 0.495,
  height = 1
) + 
  
  # ---------------------------------------------------------------------------
# Panel labels
#
# These are added manually so the inset does not receive its own letter.
# ---------------------------------------------------------------------------

cowplot::draw_label(
  
  "A",
  
  x = 0.012,
  y = 0.992,
  
  hjust = 0,
  vjust = 1,
  
  size = 22,
  fontface = "bold"
) +
  
  cowplot::draw_label(
    
    "B",
    
    x = 0.432,
    y = 0.992,
    
    hjust = 0,
    vjust = 1,
    
    size = 22,
    fontface = "bold"
  )


# ==============================================================================
# 3. DISPLAY
# ==============================================================================

print(final_combined_plot)






# ==============================================================================
# 28. SAVE FIGURES
# ==============================================================================

## save final plot
pdf("Output/Graphs/ASV_overlaps.pdf", width = 14, height = 6)
print(final_combined_plot)
dev.off()
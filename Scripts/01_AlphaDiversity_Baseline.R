## ============================================================================
## 01. ALPHA DIVERSITY — BASELINE MICROBIOME
## Repeated-rarefaction primary analysis + unrarefied sensitivity analysis
##
## PRIMARY ANALYSIS
##   - Rarefy every sample repeatedly to 15,857 reads
##   - Calculate Observed richness, Faith's PD, Shannon, and Inverse Simpson
##     diversity after each rarefaction
##   - Average each metric across rarefactions for each sample
##   - Use those sample-level means in all primary models and figures
##
## SENSITIVITY ANALYSIS
##   - Recalculate the same diversity metrics from the unrarefied ASV table
##   - Fit the same mixed models
##   - Compare pairwise estimates from repeated-rarefaction vs unrarefied data
##
## IMPORTANT
##   - Samples are never treated as independent across rarefaction iterations.
##   - The repeated rarefactions are used only to obtain stable sample-level
##     diversity estimates and quantify stochastic rarefaction uncertainty.
## ============================================================================


## ============================================================================
## 0. SETUP
## ============================================================================

rm(list = ls())

# Packages
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(vegan)
library(picante)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)
library(broom.mixed)
library(openxlsx)
library(readr)


## ============================================================================
## 1. ANALYSIS SETTINGS
## ============================================================================

# Number of independent rarefactions.
N_RAREFACTIONS <- 100

# Reproducibility seed.
RAREFACTION_SEED <- 123

# Fixed rarefaction depth used in the manuscript analysis.
RAREFACTION_DEPTH <- 15857

# Plotting order.
sample_origin_plot_levels <- c(
  "Wild:As",
  "Wild-like:As",
  "Lab-like:As",
  "Lab:Mm"
)

# Model reference group: Lab:Mm.
sample_origin_model_levels <- c(
  "Lab:Mm",
  "Lab-like:As",
  "Wild-like:As",
  "Wild:As"
)

# Existing group colours.
sample_origin_colours <- c(
  "Wild:As"      = "#74C476",
  "Wild-like:As" = "#6BAED6",
  "Lab-like:As"  = "#BC80BD",
  "Lab:Mm"       = "#FB8072"
)


## ============================================================================
## 2. INPUT / OUTPUT
## ============================================================================

# Run Scripts/00_DataPrep_Baseline.R first whenever preprocessing changes.
ps_all <- readRDS("Data/Processed/baseline_microbiome_phyloseq.rds")

results_dir <- "Results/01_Alpha_Diversity"
tables_dir <- file.path(results_dir, "Tables")
figures_dir <- file.path(results_dir, "Figures")
objects_dir <- file.path(results_dir, "Objects")


# Output policy:
#   - save ONE multi-sheet workbook containing only manuscript-relevant results;
#   - save only complete reader-facing figures (never individual panels);
#   - save the alpha-panel RDS only because the beta-diversity script needs it
#     to assemble the final combined main-text diversity figure.

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

dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================================
## 3. BASIC DATA CHECKS
## ============================================================================

sample_depths <- phyloseq::sample_sums(ps_all)

# Every retained sample must have at least the selected rarefaction depth.
stopifnot(all(sample_depths >= RAREFACTION_DEPTH))

# Confirm sample names and metadata row names are identical.
metadata_check <- as.data.frame(phyloseq::sample_data(ps_all))

## ============================================================================
## 4. HELPER FUNCTION: CALCULATE ALPHA DIVERSITY
## ============================================================================

# All metrics are calculated from the SAME sample x ASV count matrix.
# This avoids sample-ID inconsistencies between separate diversity functions.

calculate_alpha_metrics <- function(
    ps_object,
    iteration = NA_integer_) {

  sample_ids <- phyloseq::sample_names(ps_object)

  # --------------------------------------------------------------------------
  # Extract ASV count matrix
  # --------------------------------------------------------------------------

  otu <- as(
    phyloseq::otu_table(ps_object),
    "matrix"
  )

  # vegan::diversity() and picante::pd() require:
  # rows = samples
  # columns = taxa
  if (phyloseq::taxa_are_rows(ps_object)) {
    otu <- t(otu)
  }

  # Confirm all sample IDs are retained.
  stopifnot(
    setequal(
      rownames(otu),
      sample_ids
    )
  )

  # Force exact phyloseq sample order.
  otu <- otu[
    sample_ids,
    ,
    drop = FALSE
  ]

  stopifnot(
    identical(
      rownames(otu),
      sample_ids
    )
  )

  # --------------------------------------------------------------------------
  # Observed richness
  # --------------------------------------------------------------------------

  observed <- rowSums(
    otu > 0
  )

  # --------------------------------------------------------------------------
  # Shannon diversity
  # --------------------------------------------------------------------------

  shannon <- vegan::diversity(
    otu,
    index = "shannon"
  )

  # --------------------------------------------------------------------------
  # Inverse Simpson diversity
  # --------------------------------------------------------------------------

  inverse_simpson <- vegan::diversity(
    otu,
    index = "invsimpson"
  )

  # --------------------------------------------------------------------------
  # Faith's phylogenetic diversity
  # --------------------------------------------------------------------------

  faith_pd_results <- picante::pd(
    samp = otu,
    tree = phyloseq::phy_tree(ps_object),
    include.root = TRUE
  )

  stopifnot(
    identical(
      rownames(faith_pd_results),
      sample_ids
    )
  )

  faith_pd <- faith_pd_results$PD

  # --------------------------------------------------------------------------
  # Return one row per sample
  # --------------------------------------------------------------------------

  alpha_metrics <- tibble::tibble(
    SampleID = sample_ids,
    iteration = iteration,
    Observed = as.numeric(observed),
    FaithPD = as.numeric(faith_pd),
    Shannon = as.numeric(shannon),
    InverseSimpson = as.numeric(inverse_simpson)
  )

  # Fail immediately if anything went wrong.
  stopifnot(
    nrow(alpha_metrics) ==
      phyloseq::nsamples(ps_object),

    identical(
      alpha_metrics$SampleID,
      sample_ids
    ),

    !anyNA(alpha_metrics$Observed),
    !anyNA(alpha_metrics$FaithPD),
    !anyNA(alpha_metrics$Shannon),
    !anyNA(alpha_metrics$InverseSimpson)
  )

  alpha_metrics
}


## ============================================================================
## 5. REPEATED RAREFACTION
## ============================================================================

cat(
  "Running ",
  N_RAREFACTIONS,
  " rarefactions at ",
  RAREFACTION_DEPTH,
  " reads/sample.\n",
  sep = ""
)

# Use a different deterministic seed for each rarefaction.
rarefaction_seeds <- RAREFACTION_SEED +
  seq_len(N_RAREFACTIONS) - 1L

alpha_iterations_list <- vector(
  mode = "list",
  length = N_RAREFACTIONS
)

for (i in seq_len(N_RAREFACTIONS)) {

  message(
    "Rarefaction iteration ",
    i,
    " of ",
    N_RAREFACTIONS
  )

  ps_i <- phyloseq::rarefy_even_depth(
    ps_all,
    sample.size = RAREFACTION_DEPTH,
    rngseed = rarefaction_seeds[i],
    replace = FALSE,
    trimOTUs = TRUE,
    verbose = FALSE
  )

  # Rarefaction should preserve all samples and sample names.
  stopifnot(
    phyloseq::nsamples(ps_i) ==
      phyloseq::nsamples(ps_all),

    identical(
      phyloseq::sample_names(ps_i),
      phyloseq::sample_names(ps_all)
    )
  )

  alpha_iterations_list[[i]] <-
    calculate_alpha_metrics(
      ps_object = ps_i,
      iteration = i
    )
}

alpha_iterations <- dplyr::bind_rows(
  alpha_iterations_list
)

# Expected: 239 samples x 100 iterations = 23,900 rows.
stopifnot(
  nrow(alpha_iterations) ==
    phyloseq::nsamples(ps_all) * N_RAREFACTIONS,

  length(unique(alpha_iterations$SampleID)) ==
    phyloseq::nsamples(ps_all),

  setequal(
    unique(alpha_iterations$SampleID),
    phyloseq::sample_names(ps_all)
  )
)

# Confirm no missing values were created.
alpha_missingness_check <- alpha_iterations %>%
  dplyr::summarise(
    Observed_NA = sum(is.na(Observed)),
    FaithPD_NA = sum(is.na(FaithPD)),
    Shannon_NA = sum(is.na(Shannon)),
    InverseSimpson_NA = sum(is.na(InverseSimpson))
  )

print(alpha_missingness_check)

stopifnot(
  all(alpha_missingness_check == 0)
)


## ============================================================================
## 6. SUMMARISE REPEATED RAREFACTIONS PER SAMPLE
## ============================================================================

# The MEAN across rarefactions is used in all primary figures and models.
# SD/min/max quantify stochastic uncertainty due to random subsampling.

alpha_sample_summary <- alpha_iterations %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(

    # IMPORTANT:
    # Calculate SD/min/max BEFORE creating the mean column with the same name.
    # Otherwise dplyr::summarise() can use the newly created scalar mean in
    # subsequent expressions, causing SD to become NA/NaN.

    Observed_SD = stats::sd(
      Observed,
      na.rm = TRUE
    ),
    Observed_min = min(
      Observed,
      na.rm = TRUE
    ),
    Observed_max = max(
      Observed,
      na.rm = TRUE
    ),
    Observed = mean(
      Observed,
      na.rm = TRUE
    ),

    FaithPD_SD = stats::sd(
      FaithPD,
      na.rm = TRUE
    ),
    FaithPD_min = min(
      FaithPD,
      na.rm = TRUE
    ),
    FaithPD_max = max(
      FaithPD,
      na.rm = TRUE
    ),
    FaithPD = mean(
      FaithPD,
      na.rm = TRUE
    ),

    Shannon_SD = stats::sd(
      Shannon,
      na.rm = TRUE
    ),
    Shannon_min = min(
      Shannon,
      na.rm = TRUE
    ),
    Shannon_max = max(
      Shannon,
      na.rm = TRUE
    ),
    Shannon = mean(
      Shannon,
      na.rm = TRUE
    ),

    InverseSimpson_SD = stats::sd(
      InverseSimpson,
      na.rm = TRUE
    ),
    InverseSimpson_min = min(
      InverseSimpson,
      na.rm = TRUE
    ),
    InverseSimpson_max = max(
      InverseSimpson,
      na.rm = TRUE
    ),
    InverseSimpson = mean(
      InverseSimpson,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

stopifnot(
  nrow(alpha_sample_summary) ==
    phyloseq::nsamples(ps_all),

  setequal(
    alpha_sample_summary$SampleID,
    phyloseq::sample_names(ps_all)
  )
)

# Simple numerical summary of rarefaction uncertainty.
alpha_rarefaction_variability <- alpha_sample_summary %>%
  dplyr::summarise(
    mean_Observed_SD =
      mean(Observed_SD, na.rm = TRUE),

    mean_FaithPD_SD =
      mean(FaithPD_SD, na.rm = TRUE),

    mean_Shannon_SD =
      mean(Shannon_SD, na.rm = TRUE),

    mean_InverseSimpson_SD =
      mean(InverseSimpson_SD, na.rm = TRUE)
  )

print(alpha_rarefaction_variability)


## ============================================================================
## 7. ADD METADATA
## ============================================================================

# Extract metadata.
# Convert sample_data explicitly to a plain base data.frame.
# data.frame() strips the phyloseq "sample_data" class more reliably than
# as.data.frame() for downstream dplyr verbs.
ps_meta_rare <- data.frame(
  phyloseq::sample_data(ps_all),
  check.names = FALSE
)

# Use phyloseq sample_names() as the canonical sample ID.
ps_meta_rare$SampleID <- phyloseq::sample_names(ps_all)

# Put ID first.
ps_meta_rare <- ps_meta_rare[
  ,
  c(
    "SampleID",
    setdiff(
      names(ps_meta_rare),
      "SampleID"
    )
  ),
  drop = FALSE
]

# Confirm exact agreement before joining.
stopifnot(
  nrow(ps_meta_rare) ==
    phyloseq::nsamples(ps_all),

  setequal(
    ps_meta_rare$SampleID,
    alpha_sample_summary$SampleID
  )
)

# Join sample-level repeated-rarefaction means to metadata.
ps_meta_rare <- ps_meta_rare %>%
  dplyr::left_join(
    alpha_sample_summary,
    by = "SampleID"
  )

# Final join checks.
stopifnot(
  nrow(ps_meta_rare) ==
    phyloseq::nsamples(ps_all),

  !anyNA(ps_meta_rare$Observed),
  !anyNA(ps_meta_rare$FaithPD),
  !anyNA(ps_meta_rare$Shannon),
  !anyNA(ps_meta_rare$InverseSimpson)
)

cat(
  "Metadata successfully joined to repeated-rarefaction estimates.\n"
)


## ============================================================================
## 8. REPEATED-RAREFACTION OUTPUT POLICY
## ============================================================================

# alpha_iterations and alpha_sample_summary are retained in memory for the
# analysis but are not written as reader-facing files. The manuscript only
# exports results needed for the main text or Supplementary Material.


## ============================================================================
## 9. GROUP SUMMARIES FOR FIGURES
## ============================================================================

ps_meta_plot <- ps_meta_rare %>%
  dplyr::mutate(
    SampleOrigin = factor(
      SampleOrigin,
      levels = sample_origin_plot_levels
    )
  )

summarise_alpha_metric <- function(
    data,
    metric) {

  data %>%
    dplyr::group_by(SampleOrigin) %>%
    dplyr::summarise(
      mean = mean(
        .data[[metric]],
        na.rm = TRUE
      ),

      se = stats::sd(
        .data[[metric]],
        na.rm = TRUE
      ) /
        sqrt(
          sum(
            !is.na(.data[[metric]])
          )
        ),

      sd = stats::sd(
        .data[[metric]],
        na.rm = TRUE
      ),

      min = min(
        .data[[metric]],
        na.rm = TRUE
      ),

      max = max(
        .data[[metric]],
        na.rm = TRUE
      ),

      n = sum(
        !is.na(.data[[metric]])
      ),

      .groups = "drop"
    )
}

alpha_obs_R <- summarise_alpha_metric(
  ps_meta_plot,
  "Observed"
)

alpha_pd_R <- summarise_alpha_metric(
  ps_meta_plot,
  "FaithPD"
)

alpha_shan_R <- summarise_alpha_metric(
  ps_meta_plot,
  "Shannon"
)

alpha_sim_R <- summarise_alpha_metric(
  ps_meta_plot,
  "InverseSimpson"
)

alpha_group_summary <- dplyr::bind_rows(
  dplyr::mutate(
    alpha_obs_R,
    Metric = "Observed richness"
  ),
  dplyr::mutate(
    alpha_pd_R,
    Metric = "Faith's PD"
  ),
  dplyr::mutate(
    alpha_shan_R,
    Metric = "Shannon"
  ),
  dplyr::mutate(
    alpha_sim_R,
    Metric = "Inverse Simpson"
  )
) %>%
  dplyr::select(
    Metric,
    dplyr::everything()
  )



## ============================================================================
## 10. PRIMARY ALPHA-DIVERSITY FIGURES
## ============================================================================

make_alpha_plot <- function(
    data,
    summary_data,
    metric,
    y_label) {

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = SampleOrigin,
      y = .data[[metric]],
      colour = SampleOrigin
    )
  ) +

    ggforce::geom_sina(
      alpha = 1,
      size = 2,
      show.legend = FALSE
    ) +

    ggplot2::geom_point(
      data = summary_data,
      ggplot2::aes(
        x = SampleOrigin,
        y = mean
      ),
      inherit.aes = FALSE,
      colour = "black",
      size = 3
    ) +

    ggplot2::geom_errorbar(
      data = summary_data,
      ggplot2::aes(
        x = SampleOrigin,
        ymin = mean - se,
        ymax = mean + se
      ),
      inherit.aes = FALSE,
      width = 0.3,
      linewidth = 0.5
    ) +

    ggplot2::scale_colour_manual(
      values = sample_origin_colours
    ) +

    ggplot2::scale_x_discrete(
      guide = ggplot2::guide_axis(
        angle = 45
      )
    ) +

    ggplot2::labs(
      x = "",
      y = y_label
    ) +

    ggplot2::theme_bw(
      base_size = 20
    )
}

obs_all_R <- make_alpha_plot(
  ps_meta_plot,
  alpha_obs_R,
  "Observed",
  "Observed richness"
)

pd_all_R <- make_alpha_plot(
  ps_meta_plot,
  alpha_pd_R,
  "FaithPD",
  "Faith's phylogenetic diversity"
)

shan_all_R <- make_alpha_plot(
  ps_meta_plot,
  alpha_shan_R,
  "Shannon",
  "Shannon diversity"
)

sim_all_R <- make_alpha_plot(
  ps_meta_plot,
  alpha_sim_R,
  "InverseSimpson",
  "Inverse Simpson diversity"
)

# --------------------------------------------------------------------------
# Main-text alpha panels for the combined alpha + beta diversity figure
# --------------------------------------------------------------------------

# We do NOT save a standalone main alpha-diversity figure here because the
# manuscript main figure combines alpha- and beta-diversity panels. Instead,
# we save the alpha plotting objects so that 02_BetaDiversity_Baseline can
# assemble the final combined Figure 2.

alpha_main_panels <- list(
  shannon_plot = shan_all_R,
  faith_pd_plot = pd_all_R,
  alpha_plot_data = ps_meta_plot,
  alpha_group_summary = alpha_group_summary,
  sample_origin_colours = sample_origin_colours,
  sample_origin_plot_levels = sample_origin_plot_levels,
  n_rarefactions = N_RAREFACTIONS,
  rarefaction_depth = RAREFACTION_DEPTH
)

saveRDS(
  alpha_main_panels,
  file = file.path(
    objects_dir,
    "alpha_main_panels_for_combined_figure.rds"
  )
)

# --------------------------------------------------------------------------
# Supplementary figure: complementary alpha-diversity metrics only
# --------------------------------------------------------------------------

alpha_supp_metrics <- ggpubr::ggarrange(
  obs_all_R,
  sim_all_R,
  labels = c("A", "B"),
  nrow = 1
)

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Alpha_Diversity_Complementary_Metrics.pdf"
  ),
  plot = alpha_supp_metrics,
  width = 8,
  height = 5
)


## ============================================================================
## 11. RAREFACTION STABILITY QC
## ============================================================================

# Rarefaction variability is retained in memory as quality-control information.
# No separate uncertainty figure/table is exported because it is not necessary
# for interpretation of the manuscript results.

print(alpha_rarefaction_variability)


## ============================================================================
## 12. PREPARE PRIMARY MODEL DATA
## ============================================================================

# SeqPlate is intentionally omitted from the alpha-diversity models because
# sequencing plate is confounded with SampleOrigin in this dataset and
# therefore cannot be estimated cleanly as an independent technical effect.

ps_meta_model <- ps_meta_rare

ps_meta_model$SampleOrigin <- factor(
  ps_meta_model$SampleOrigin,
  levels = sample_origin_model_levels
)

ps_meta_model$Sex <- as.factor(
  ps_meta_model$Sex
)

# Inverse Simpson requires log transformation based on the previous
# diagnostic workflow.
ps_meta_model$Log_InverseSimpson <- log(
  ps_meta_model$InverseSimpson
)


## ============================================================================
## 13. PRIMARY MIXED MODELS
## ============================================================================

# --------------------------------------------------------------------------
# Observed richness
# --------------------------------------------------------------------------

model_observed <- glmmTMB::glmmTMB(
  Observed ~
    SampleOrigin +
    Sex +
    (1 | mouse_ID_dataset),
  data = ps_meta_model
)

summary(model_observed)

emmeans_observed <- emmeans::emmeans(
  model_observed,
  ~ SampleOrigin
)

pairs(
  emmeans_observed,
  adjust = "tukey"
)

car::Anova(
  model_observed,
  type = 3
)

res_observed <- DHARMa::simulateResiduals(
  model_observed
)

plot(res_observed)
DHARMa::testDispersion(res_observed)
DHARMa::testUniformity(res_observed)
DHARMa::testOutliers(res_observed)


# --------------------------------------------------------------------------
# Faith's PD
# --------------------------------------------------------------------------

model_faith <- glmmTMB::glmmTMB(
  FaithPD ~
    SampleOrigin +
    Sex +
    (1 | mouse_ID_dataset),
  data = ps_meta_model
)

summary(model_faith)

emmeans_faith <- emmeans::emmeans(
  model_faith,
  ~ SampleOrigin
)

pairs(
  emmeans_faith,
  adjust = "tukey"
)

car::Anova(
  model_faith,
  type = 3
)

res_faith <- DHARMa::simulateResiduals(
  model_faith
)

plot(res_faith)
DHARMa::testDispersion(res_faith)
DHARMa::testUniformity(res_faith)
DHARMa::testOutliers(res_faith)


# --------------------------------------------------------------------------
# Shannon diversity
# --------------------------------------------------------------------------

model_shannon <- glmmTMB::glmmTMB(
  Shannon ~
    SampleOrigin +
    Sex +
    (1 | mouse_ID_dataset),
  data = ps_meta_model
)

summary(model_shannon)

emmeans_shannon <- emmeans::emmeans(
  model_shannon,
  ~ SampleOrigin
)

pairs(
  emmeans_shannon,
  adjust = "tukey"
)

car::Anova(
  model_shannon,
  type = 3
)

res_shannon <- DHARMa::simulateResiduals(
  model_shannon
)

plot(res_shannon)
DHARMa::testDispersion(res_shannon)
DHARMa::testUniformity(res_shannon)
DHARMa::testOutliers(res_shannon)


# --------------------------------------------------------------------------
# Inverse Simpson diversity
# --------------------------------------------------------------------------

model_inverse <- glmmTMB::glmmTMB(
  Log_InverseSimpson ~
    SampleOrigin +
    Sex +
    (1 | mouse_ID_dataset),
  data = ps_meta_model
)

summary(model_inverse)

emmeans_inverse <- emmeans::emmeans(
  model_inverse,
  ~ SampleOrigin
)

pairs(
  emmeans_inverse,
  adjust = "tukey"
)

car::Anova(
  model_inverse,
  type = 3
)

res_inverse <- DHARMa::simulateResiduals(
  model_inverse
)

plot(res_inverse)
DHARMa::testDispersion(res_inverse)
DHARMa::testUniformity(res_inverse)
DHARMa::testOutliers(res_inverse)


## ============================================================================
## 14. EXPORT PRIMARY MODEL RESULTS
## ============================================================================

alpha_models <- list(
  "Observed richness" = model_observed,
  "Faith PD" = model_faith,
  "Shannon" = model_shannon,
  "log(Inverse Simpson)" = model_inverse
)

# Fixed-effect coefficients.
alpha_GLMM_results <- lapply(
  names(alpha_models),
  function(metric) {

    broom.mixed::tidy(
      alpha_models[[metric]],
      effects = "fixed",
      conf.int = TRUE
    ) %>%
      dplyr::mutate(
        Metric = metric
      )
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::select(
    Metric,
    term,
    estimate,
    std.error,
    statistic,
    p.value,
    conf.low,
    conf.high
  ) %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ ""
    )
  )

# Type III omnibus tests.
tidy_type3_anova <- function(
    model,
    metric) {

  x <- as.data.frame(
    car::Anova(
      model,
      type = 3
    )
  )

  x$Predictor <- rownames(x)
  rownames(x) <- NULL

  x %>%
    dplyr::filter(
      Predictor != "(Intercept)"
    ) %>%
    dplyr::transmute(
      Metric = metric,
      Predictor = Predictor,
      Chisq = Chisq,
      Df = Df,
      p.value = `Pr(>Chisq)`,
      significance = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value < 0.1 ~ ".",
        TRUE ~ ""
      )
    )
}

anova_results <- lapply(
  names(alpha_models),
  function(metric) {

    tidy_type3_anova(
      alpha_models[[metric]],
      metric
    )
  }
) %>%
  dplyr::bind_rows()

# Tukey-adjusted pairwise contrasts.
extract_tukey <- function(
    model,
    metric,
    analysis_label) {

  emmeans::emmeans(
    model,
    ~ SampleOrigin
  ) %>%
    pairs(
      adjust = "tukey"
    ) %>%
    summary(
      infer = TRUE
    ) %>%
    as.data.frame() %>%
    dplyr::mutate(
      Metric = metric,
      Analysis = analysis_label
    )
}

tukey_results <- lapply(
  names(alpha_models),
  function(metric) {

    extract_tukey(
      alpha_models[[metric]],
      metric,
      paste0(
        "Mean of ",
        N_RAREFACTIONS,
        " rarefactions"
      )
    )
  }
) %>%
  dplyr::bind_rows()

# Model metadata.
model_stats <- data.frame(
  Metric = names(alpha_models),
  Transformation = c(
    "None",
    "None",
    "None",
    "Log"
  ),
  Family = rep(
    "Gaussian",
    length(alpha_models)
  ),
  Random_effects = rep(
    "mouse_ID_dataset",
    length(alpha_models)
  ),
  AIC = sapply(
    alpha_models,
    AIC
  ),
  BIC = sapply(
    alpha_models,
    BIC
  ),
  N = sapply(
    alpha_models,
    nobs
  )
)

random_effects <- lapply(
  names(alpha_models),
  function(metric) {

    broom.mixed::tidy(
      alpha_models[[metric]],
      effects = "ran_pars"
    ) %>%
      dplyr::mutate(
        Metric = metric
      )
  }
) %>%
  dplyr::bind_rows()


## ============================================================================
## 15. UNRAREFIED SENSITIVITY ANALYSIS
## ============================================================================

# Use the same model structure as the primary analysis; SeqPlate is omitted here too.

# Calculate all four diversity metrics directly from the unrarefied ASV table.
raw_metrics <- calculate_alpha_metrics(
  ps_object = ps_all,
  iteration = NA_integer_
) %>%
  dplyr::select(
    -iteration
  )

# Metadata with canonical phyloseq sample IDs.
ps_meta_raw <- data.frame(
  phyloseq::sample_data(ps_all),
  check.names = FALSE
)

ps_meta_raw$SampleID <- phyloseq::sample_names(
  ps_all
)

ps_meta_raw <- ps_meta_raw[
  ,
  c(
    "SampleID",
    setdiff(
      names(ps_meta_raw),
      "SampleID"
    )
  ),
  drop = FALSE
]

ps_meta_raw <- dplyr::left_join(
  ps_meta_raw,
  raw_metrics,
  by = "SampleID"
)

stopifnot(
  nrow(ps_meta_raw) ==
    phyloseq::nsamples(ps_all),

  !anyNA(ps_meta_raw$Observed),
  !anyNA(ps_meta_raw$FaithPD),
  !anyNA(ps_meta_raw$Shannon),
  !anyNA(ps_meta_raw$InverseSimpson)
)

ps_meta_raw$SampleOrigin <- factor(
  ps_meta_raw$SampleOrigin,
  levels = sample_origin_model_levels
)

ps_meta_raw$Sex <- as.factor(
  ps_meta_raw$Sex
)

ps_meta_raw$Log_InverseSimpson <- log(
  ps_meta_raw$InverseSimpson
)

# Fit same four models.
raw_models <- list(

  "Observed richness" = glmmTMB::glmmTMB(
    Observed ~
      SampleOrigin +
      Sex +
      (1 | mouse_ID_dataset),
    data = ps_meta_raw
  ),

  "Faith PD" = glmmTMB::glmmTMB(
    FaithPD ~
      SampleOrigin +
      Sex +
      (1 | mouse_ID_dataset),
    data = ps_meta_raw
  ),

  "Shannon" = glmmTMB::glmmTMB(
    Shannon ~
      SampleOrigin +
      Sex +
      (1 | mouse_ID_dataset),
    data = ps_meta_raw
  ),

  "log(Inverse Simpson)" = glmmTMB::glmmTMB(
    Log_InverseSimpson ~
      SampleOrigin +
      Sex +
      (1 | mouse_ID_dataset),
    data = ps_meta_raw
  )
)

raw_tukey_results <- lapply(
  names(raw_models),
  function(metric) {

    extract_tukey(
      raw_models[[metric]],
      metric,
      "Unrarefied"
    )
  }
) %>%
  dplyr::bind_rows()

# Combine primary + sensitivity contrasts.
sensitivity_pairwise <- dplyr::bind_rows(
  tukey_results,
  raw_tukey_results
)



## ============================================================================
## 16. SUPPLEMENTARY FIGURE: PRIMARY VS UNRAREFIED SENSITIVITY
## ============================================================================

sensitivity_plot_data <- sensitivity_pairwise %>%
  dplyr::mutate(

    Metric = factor(
      Metric,
      levels = c(
        "Observed richness",
        "Faith PD",
        "Shannon",
        "log(Inverse Simpson)"
      )
    ),

    Analysis = factor(
      Analysis,
      levels = c(
        paste0(
          "Mean of ",
          N_RAREFACTIONS,
          " rarefactions"
        ),
        "Unrarefied"
      )
    )
  )

sensitivity_plot <- ggplot2::ggplot(
  sensitivity_plot_data,
  ggplot2::aes(
    x = estimate,
    y = contrast,
    shape = Analysis
  )
) +

  ggplot2::geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +

  ggplot2::geom_errorbar(
    ggplot2::aes(
      xmin = asymp.LCL,
      xmax = asymp.UCL
    ),
    position = ggplot2::position_dodge(
      width = 0.5
    ),
    width = 0,
    orientation = "y"
  ) +

  ggplot2::geom_point(
    position = ggplot2::position_dodge(
      width = 0.5
    ),
    size = 2.2
  ) +

  ggplot2::facet_wrap(
    ~ Metric,
    scales = "free_x",
    ncol = 2
  ) +

  ggplot2::labs(
    x = "Estimated pairwise difference (95% CI)",
    y = "",
    shape = "Analysis"
  ) +

  ggplot2::theme_bw(
    base_size = 11
  ) +

  ggplot2::theme(
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Alpha_Diversity_Rarefaction_Sensitivity.pdf"
  ),
  plot = sensitivity_plot,
  width = 10,
  height = 8
)


## ============================================================================
## 17. SAVE ONE STREAMLINED ALPHA-DIVERSITY RESULTS WORKBOOK
## ============================================================================

# This is the ONLY reader-facing table file produced by this script.
# It contains the four result sets needed for the manuscript and supplement.
#
# Sheets:
#   Group Summary
#     Descriptive mean, SE, SD, range, and n for each diversity metric by cohort.
#
#   Omnibus Tests
#     Type III tests for SampleOrigin and Sex from the primary repeated-
#     rarefaction mixed models.
#
#   Pairwise Tests
#     Tukey-adjusted SampleOrigin contrasts from the primary repeated-
#     rarefaction mixed models.
#
#   Rarefaction Sensitivity
#     The same pairwise contrasts from the primary repeated-rarefaction and
#     unrarefied analyses, allowing direct assessment of robustness.

alpha_results_workbook <- list(
  "Group Summary" =
    alpha_group_summary,

  "Omnibus Tests" =
    anova_results,

  "Pairwise Tests" =
    tukey_results,

  "Rarefaction Sensitivity" =
    sensitivity_pairwise
)

openxlsx::write.xlsx(
  alpha_results_workbook,
  file = file.path(
    tables_dir,
    "Alpha_Diversity_results.xlsx"
  ),
  rowNames = FALSE,
  overwrite = TRUE
)

cat(
  "\nSaved alpha-diversity results workbook to: ",
  file.path(
    tables_dir,
    "Alpha_Diversity_results.xlsx"
  ),
  "\n",
  sep = ""
)


## ============================================================================
## 18. OUTPUT SUMMARY
## ============================================================================

cat(
  "\n============================================================\n",
  "ALPHA-DIVERSITY OUTPUTS\n",
  "============================================================\n",
  "Workbook:\n  ",
  file.path(tables_dir, "Alpha_Diversity_results.xlsx"),
  "\n\nReader-facing figures:\n  ",
  file.path(figures_dir, "Alpha_Diversity_Complementary_Metrics.pdf"),
  "\n  ",
  file.path(figures_dir, "Alpha_Diversity_Rarefaction_Sensitivity.pdf"),
  "\n\nDownstream figure object (not a manuscript table/figure):\n  ",
  file.path(objects_dir, "alpha_main_panels_for_combined_figure.rds"),
  "\n============================================================\n",
  sep = ""
)

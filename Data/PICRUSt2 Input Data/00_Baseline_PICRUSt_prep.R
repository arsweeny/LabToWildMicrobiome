# ============================================================
# PREPARE PICRUSt2 INPUT FILES
# ============================================================

library(phyloseq)
library(dplyr)
library(Biostrings)

# Generate ps_all and ASV key
source("Scripts/00_DataPrep_Baseline.R")

# Start from final project phyloseq object
ps_picrust <- ps_all

rm(list = setdiff(ls(), c("asv_key", "ps_picrust")))

# Remove taxa with zero abundance after final sample filtering
ps_picrust <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_picrust) > 0,
  ps_picrust
)

# Remove zero-depth samples as a safety check
ps_picrust <- phyloseq::prune_samples(
  phyloseq::sample_sums(ps_picrust) > 0,
  ps_picrust
)

# ------------------------------------------------------------
# Check ASV mapping
# ------------------------------------------------------------

picrust_asvs <- phyloseq::taxa_names(ps_picrust)

all(picrust_asvs %in% asv_key$ASV)

setdiff(
  picrust_asvs,
  asv_key$ASV
)

# ------------------------------------------------------------
# Make matched sequence key
# ------------------------------------------------------------

sequence_key_picrust <- asv_key %>%
  dplyr::filter(
    ASV %in% picrust_asvs
  ) %>%
  dplyr::mutate(
    ASV = factor(
      ASV,
      levels = picrust_asvs
    )
  ) %>%
  dplyr::arrange(
    ASV
  ) %>%
  dplyr::mutate(
    ASV = as.character(ASV)
  )

# Confirm exact order
stopifnot(
  identical(
    sequence_key_picrust$ASV,
    picrust_asvs
  )
)

# Confirm no missing sequences
stopifnot(
  !any(is.na(sequence_key_picrust$Sequence))
)

# ------------------------------------------------------------
# Build FASTA sequences
# ------------------------------------------------------------

picrust_sequences <- Biostrings::DNAStringSet(
  sequence_key_picrust$Sequence
)

names(picrust_sequences) <- sequence_key_picrust$ASV

# Check
head(picrust_sequences)

# ------------------------------------------------------------
# Export FASTA
# ------------------------------------------------------------

dir.create(
  "PICRUSt2/input",
  recursive = TRUE,
  showWarnings = FALSE
)

Biostrings::writeXStringSet(
  picrust_sequences,
  filepath = "PICRUSt2/input/ASV_sequences.fasta",
  format = "fasta"
)

# ------------------------------------------------------------
# Export raw ASV count table
# ------------------------------------------------------------

asv_counts <- as(
  phyloseq::otu_table(ps_picrust),
  "matrix"
)

if (!phyloseq::taxa_are_rows(ps_picrust)) {
  asv_counts <- t(asv_counts)
}

stopifnot(
  identical(
    rownames(asv_counts),
    picrust_asvs
  )
)

stopifnot(
  all(asv_counts %% 1 == 0)
)

utils::write.table(
  asv_counts,
  file = "PICRUSt2/input/ASV_counts.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# ------------------------------------------------------------
# Final check
# ------------------------------------------------------------

cat(
  "Samples:", ncol(asv_counts), "\n",
  "ASVs:", nrow(asv_counts), "\n",
  "Sequences:", length(picrust_sequences), "\n"
)

cat(
  "ASV IDs match FASTA:",
  identical(
    rownames(asv_counts),
    names(picrust_sequences)
  ),
  "\n"
)

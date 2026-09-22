# data_loader.R
#
# Loads the raw TCGA-COAD + TCGA-READ count matrix, applies the low-count
# filter, and builds the DESeq2 sample/condition table (colData). Same
# structure as the HNSC project's data_loader.R - the only things that
# change per-dataset are the ID lists below and build_coldata()'s call
# site, not the loading/matching logic itself.

library(here)
library(yaml)

load_config <- function(path = here("config", "config.yaml")) {
  yaml::read_yaml(path)
}

load_counts <- function(cfg) {
  raw <- read.delim(here(cfg$paths$raw_counts), header = TRUE, sep = ",")

  # First column is the gene index/id (Ensembl gene ID for this
  # TCGAbiolinks-derived matrix), not a sample column - drop before use.
  count_matrix <- raw[-1]

  if (!all(vapply(count_matrix, is.numeric, logical(1)))) {
    stop("data_loader: non-numeric column found in the count matrix after dropping the index column - check raw_counts formatting.")
  }

  count_matrix
}

filter_low_counts <- function(count_matrix, min_row_mean) {
  count_matrix[rowMeans(count_matrix) > min_row_mean, ]
}

# Patient-level chemotherapy outcome calls, derived from the GDC clinical
# export for TCGA-COAD + TCGA-READ (treatments.treatment_type %in%
# c("Chemotherapy", "Pharmaceutical Therapy, NOS") AND
# treatments.treatment_or_therapy == "yes"), restricted to cases whose
# disease_type is "Adenomas and Adenocarcinomas" (excludes the rarer
# cystic/mucinous/serous neoplasm subtype) and whose chemo treatment
# outcome is a single, unambiguous label across all recorded lines of
# therapy:
#   chemoresistant ("C") = treatments.treatment_outcome == "Progressive Disease"
#   chemosensitive ("S") = treatments.treatment_outcome == "Complete Response"
# Cases with no chemo-outcome record, multiple conflicting outcomes across
# treatment lines, or an intermediate outcome (Partial Response/Stable
# Disease/Treatment Ongoing) are excluded rather than guessed at.
#
# One originally-targeted chemosensitive case (TCGA-AA-3967) had no
# matching Primary Tumor RNA-seq file in GDC and was dropped automatically
# by download_tcga_colorectal_counts.R - already excluded below, so this
# list matches the actual 87-sample count matrix (19 C / 68 S).
chemoresistant_ids <- c(
  "TCGA-5M-AAT6", "TCGA-AA-3680", "TCGA-AA-3844", "TCGA-AA-3930", "TCGA-AA-3972",
  "TCGA-AA-A02K", "TCGA-AD-6964", "TCGA-AG-3584", "TCGA-AG-3999", "TCGA-AG-A016",
  "TCGA-AZ-4682", "TCGA-AZ-4684", "TCGA-AZ-6600", "TCGA-AZ-6606", "TCGA-CI-6620",
  "TCGA-F5-6702", "TCGA-NH-A6GA", "TCGA-NH-A6GB", "TCGA-RU-A8FL"
)
chemosensitive_ids <- c(
  "TCGA-A6-A56B", "TCGA-A6-A5ZU", "TCGA-AA-3517", "TCGA-AA-3542", "TCGA-AA-3548",
  "TCGA-AA-3560", "TCGA-AA-3562", "TCGA-AA-3678", "TCGA-AA-3841", "TCGA-AA-3860",
  "TCGA-AA-3870", "TCGA-AA-3955", "TCGA-AA-3971", "TCGA-AA-3976",
  "TCGA-AA-A00Q", "TCGA-AA-A00U", "TCGA-AA-A010", "TCGA-AA-A01F", "TCGA-AA-A01K",
  "TCGA-AA-A01T", "TCGA-AD-6889", "TCGA-AD-6901", "TCGA-AF-A56L", "TCGA-AF-A56N",
  "TCGA-AG-3591", "TCGA-AG-3593", "TCGA-AG-3600", "TCGA-AG-3609", "TCGA-AG-3611",
  "TCGA-AG-3612", "TCGA-AG-3728", "TCGA-AG-3885", "TCGA-AG-3893", "TCGA-AG-3894",
  "TCGA-AG-3909", "TCGA-AG-4005", "TCGA-AG-4008", "TCGA-AG-4022", "TCGA-AG-A00C",
  "TCGA-AG-A00H", "TCGA-AG-A01L", "TCGA-AG-A01W", "TCGA-AG-A01Y", "TCGA-AG-A02N",
  "TCGA-AG-A036", "TCGA-AH-6643", "TCGA-AY-A8YK", "TCGA-AZ-4308", "TCGA-CA-5254",
  "TCGA-CA-5255", "TCGA-CA-5256", "TCGA-CA-5797", "TCGA-CA-6715", "TCGA-CA-6716",
  "TCGA-D5-6533", "TCGA-DM-A0XF", "TCGA-DY-A0XA", "TCGA-DY-A1DE", "TCGA-F4-6569",
  "TCGA-F4-6805", "TCGA-F4-6807", "TCGA-F5-6864", "TCGA-NH-A50V", "TCGA-QG-A5YV",
  "TCGA-QG-A5YW", "TCGA-QG-A5YX", "TCGA-QG-A5Z1", "TCGA-SS-A7HO"
)

build_coldata <- function(count_matrix, condition_levels) {
  sample_ids <- colnames(count_matrix)

  # read.delim()/data.frame() run column headers through make.names(), which
  # silently turns "-" into "." (e.g. "TCGA-AA-3517" -> "TCGA.AA.3517").
  # Same fix as the HNSC pipeline - normalize both sides before matching
  # rather than relying on the raw file's header formatting.
  normalize_id <- function(x) gsub("[.-]", "-", x)
  sample_ids_norm <- normalize_id(sample_ids)
  chemoresistant_norm <- normalize_id(chemoresistant_ids)
  chemosensitive_norm <- normalize_id(chemosensitive_ids)

  condition_chr <- ifelse(
    sample_ids_norm %in% chemoresistant_norm, "C",
    ifelse(sample_ids_norm %in% chemosensitive_norm, "S", NA_character_)
  )

  unmatched <- sample_ids[is.na(condition_chr)]
  if (length(unmatched) > 0) {
    stop(sprintf(
      "data_loader: %d sample column(s) not found in either chemoresistant_ids or chemosensitive_ids: %s. Update these lists against the GDC clinical export if the cohort has changed.",
      length(unmatched), paste(unmatched, collapse = ", ")
    ))
  }

  condition <- factor(condition_chr, levels = condition_levels)

  if (length(condition) != ncol(count_matrix)) {
    stop(sprintf(
      "data_loader: condition vector length (%d) does not match sample count (%d) in the count matrix - update build_coldata().",
      length(condition), ncol(count_matrix)
    ))
  }

  data.frame(row.names = sample_ids, condition = condition)
}

load_gene_index <- function(cfg) {
  # Row-position -> gene symbol lookup. Matched to deseq2_pipeline.R's
  # results by row index (seq_len), so it must have been generated against
  # the SAME filtered gene set/order as filter_low_counts() produces here -
  # it is not a symbol lookup keyed by a stable gene ID.
  read.csv(here(cfg$paths$gene_index))
}

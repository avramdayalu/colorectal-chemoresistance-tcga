# download_tcga_colorectal_counts.R
#
# Reproducible data-acquisition step for the colorectal chemoresistance
# project: queries, downloads, and assembles the TCGA-COAD + TCGA-READ
# RNA-seq raw count matrix via TCGAbiolinks, restricted to (a) Primary
# Tumor samples, (b) the "Adenomas and Adenocarcinomas" disease-type
# bucket, and (c) exactly the patients with a clean chemoresistant/
# chemosensitive call from the clinical treatment-outcome data (see
# coad_read_id_vectors.R - built from treatments.treatment_type ==
# "Chemotherapy"/"Pharmaceutical Therapy, NOS" + treatments.treatment_or_therapy
# == "yes" + a single unambiguous treatments.treatment_outcome of
# "Progressive Disease" (resistant) or "Complete Response" (sensitive)).
#
# Output is written in the SAME format as the old project's
# data/raw/r_count0007.csv + data/raw/index_assign.csv, so it drops
# straight into the existing data_loader.R with just a config.yaml path
# change and a swap of chemoresistant_ids/chemosensitive_ids.

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("TCGAbiolinks", update = FALSE, ask = FALSE)
}
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  BiocManager::install("SummarizedExperiment", update = FALSE, ask = FALSE)
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(here)

# --- Case lists (from coad_read_id_vectors.R; regenerate that file if the
# clinical classification changes) -------------------------------------
source(here("coad_read_id_vectors.R"))  # defines chemoresistant_ids, chemosensitive_ids
target_cases <- c(chemoresistant_ids, chemosensitive_ids)
message(sprintf(
  "Targeting %d cases (%d chemoresistant, %d chemosensitive).",
  length(target_cases), length(chemoresistant_ids), length(chemosensitive_ids)
))

# --- Query + download ----------------------------------------------------
# One query per project - COAD and READ are separate GDC projects, but
# barcode filtering means we only ever pull the ~88 cases we actually want,
# not the full ~600+ combined cohort.
query_counts <- function(project) {
  GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor",
    barcode = target_cases
  )
}

query_coad <- query_counts("TCGA-COAD")
query_read <- query_counts("TCGA-READ")

GDCdownload(query_coad)
GDCdownload(query_read)

se_coad <- GDCprepare(query_coad)
se_read <- GDCprepare(query_read)

# --- Extract raw (unstranded) counts, drop the STAR summary rows ---------
extract_counts <- function(se) {
  counts <- assay(se, "unstranded")
  counts <- counts[!grepl("^N_", rownames(counts)), , drop = FALSE]
  # Sample barcodes come back at full-aliquot resolution
  # (TCGA-XX-XXXX-01A-...); truncate to the case-level barcode so columns
  # match chemoresistant_ids/chemosensitive_ids exactly.
  colnames(counts) <- substr(colnames(counts), 1, 12)
  counts
}

counts_coad <- extract_counts(se_coad)
counts_read <- extract_counts(se_read)

stopifnot(identical(rownames(counts_coad), rownames(counts_read)))  # same STAR annotation -> same gene order
counts_all <- cbind(counts_coad, counts_read)

# Guard against >1 sample per case (shouldn't happen with sample.type =
# "Primary Tumor" + one aliquot per case, but fail loudly rather than
# silently averaging/duplicating if TCGA ever has a repeat).
dup_cases <- colnames(counts_all)[duplicated(colnames(counts_all))]
if (length(dup_cases) > 0) {
  stop(sprintf(
    "download_tcga_colorectal_counts: multiple samples found for case(s): %s - resolve manually before proceeding.",
    paste(unique(dup_cases), collapse = ", ")
  ))
}

missing_cases <- setdiff(target_cases, colnames(counts_all))
if (length(missing_cases) > 0) {
  warning(sprintf(
    "%d target case(s) had no matching Primary Tumor RNA-seq file in GDC and will be dropped: %s",
    length(missing_cases), paste(missing_cases, collapse = ", ")
  ))
}

# --- Gene symbol lookup (from SummarizedExperiment rowData) --------------
gene_symbols <- rowData(se_coad)[rownames(counts_all), "gene_name"]

# --- Write outputs in the old pipeline's exact raw-data format -----------
dir.create(here("data", "raw"), recursive = TRUE, showWarnings = FALSE)

count_matrix_out <- data.frame(gene = rownames(counts_all), as.data.frame(counts_all), check.names = FALSE)
write.csv(count_matrix_out, here("data", "raw", "r_count_colorectal.csv"), row.names = FALSE)

gene_index_out <- data.frame(Index = seq_len(nrow(counts_all)), `gene symbol` = gene_symbols, check.names = FALSE)
write.csv(gene_index_out, here("data", "raw", "index_assign_colorectal.csv"), row.names = FALSE)

message(sprintf(
  "Wrote data/raw/r_count_colorectal.csv (%d genes x %d samples) and index_assign_colorectal.csv.",
  nrow(counts_all), ncol(counts_all)
))
if (length(missing_cases) > 0) {
  message(sprintf(
    "Final usable cohort: %d cases (%d requested, %d unavailable in GDC).",
    ncol(counts_all), length(target_cases), length(missing_cases)
  ))
}

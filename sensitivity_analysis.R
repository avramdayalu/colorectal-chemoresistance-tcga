# sensitivity_analysis.R
#
# Re-runs DESeq2 + GSEA after excluding the 8 chemoresistant cases that
# received a targeted biologic (bevacizumab/cetuximab/ziv-aflibercept),
# identified from GDC clinical.tsv treatments.therapeutic_agents. Checks
# whether the immune/interferon Hallmark signature (up in C in the full
# 87-case analysis) survives without them, versus the proliferation
# signature (E2F_TARGETS/G2M_CHECKPOINT), which should be more robust to
# this subgroup if it reflects baseline tumor biology rather than treatment
# exposure. See README/discussion for the full rationale.
#
# NOTE: n=11 chemoresistant cases remain after exclusion (vs 19 originally).
# A weaker/noisier result here is also consistent with reduced power, not
# only with "the original signal was an artifact" - read padj changes
# alongside NES magnitude, not padj alone.

library(here)
library(DESeq2)

source(here("src", "data_loader.R"))
source(here("src", "gsea_analysis.R"))  # reuses the existing, working run_gsea()

biologic_treated_C <- c(
  "TCGA-AA-3844", "TCGA-AA-3972", "TCGA-AA-A02K", "TCGA-AG-A016",
  "TCGA-AZ-4682", "TCGA-NH-A6GA", "TCGA-NH-A6GB", "TCGA-RU-A8FL"
)

cfg <- load_config()

# Redirect every results/figures path into a separate subfolder so this run
# cannot overwrite the original 87-case results.
cfg_sens <- cfg
cfg_sens$paths <- lapply(cfg$paths, function(p) {
  sub("^results/", "results/sensitivity_no_biologics/", p)
})
dir.create(here(cfg_sens$paths$figures_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(here(cfg_sens$paths$results_table)), recursive = TRUE, showWarnings = FALSE)

# IMPORTANT: filter_low_counts() runs on the FULL 87-sample counts_raw here,
# exactly as in deseq2_pipeline.R - only AFTER that do we drop columns. This
# keeps the filtered gene set/row order identical to the original analysis,
# so the existing index_assign_colorectal.csv (already realigned earlier)
# stays valid without needing to be rebuilt again.
counts_raw <- load_counts(cfg_sens)
counts_filtered <- filter_low_counts(counts_raw, cfg_sens$filtering$min_row_mean)
coldata_full <- build_coldata(counts_filtered, cfg_sens$design$condition_levels)

normalize_id <- function(x) gsub("[.-]", "-", x)
keep <- !(normalize_id(rownames(coldata_full)) %in% biologic_treated_C)
coldata_sens <- coldata_full[keep, , drop = FALSE]
counts_sens <- counts_filtered[, keep]

message(sprintf(
  "Sensitivity cohort: %d total (%d C, %d S) - excluded %d biologic-treated C case(s)",
  ncol(counts_sens), sum(coldata_sens$condition == "C"),
  sum(coldata_sens$condition == "S"), length(biologic_treated_C)
))

dds_sens <- DESeqDataSetFromMatrix(countData = counts_sens, colData = coldata_sens, design = ~condition)
dds_sens <- DESeq(dds_sens)

# Same contrast convention as the original: S vs C, positive = up in S.
contrast <- unlist(cfg_sens$design$contrast)
res_sens <- results(dds_sens, contrast = contrast)
coef_name <- paste0("condition_", contrast[2], "_vs_", contrast[3])
res_sens_shrink <- lfcShrink(dds_sens, coef = coef_name, type = cfg_sens$shrinkage$method, res = res_sens)

res_df <- as.data.frame(res_sens_shrink)
res_df$stat <- res_sens$stat

gene_index <- load_gene_index(cfg_sens)
if (nrow(res_df) != nrow(gene_index)) {
  stop("Row count mismatch between sensitivity results and gene_index - filtering order was not preserved.")
}
res_df$Index <- seq_len(nrow(res_df))
merged_sens <- merge(res_df, gene_index, by.x = "Index", by.y = "Index", all.x = TRUE)
write.csv(merged_sens, here(cfg_sens$paths$results_table), row.names = FALSE)

pipeline_out_sens <- list(dds = dds_sens, results = merged_sens, cfg = cfg_sens)
gsea_out_sens <- run_gsea(pipeline_out_sens)

# --- Compare NES/padj before vs after, for the pathways that matter ---
gsea_orig <- read.csv(here(cfg$paths$gsea_results))
gsea_sens <- gsea_out_sens$results

compare_pathways <- c(
  "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT"
)

comparison <- merge(
  gsea_orig[gsea_orig$pathway %in% compare_pathways, c("pathway", "NES", "padj")],
  gsea_sens[gsea_sens$pathway %in% compare_pathways, c("pathway", "NES", "padj")],
  by = "pathway", suffixes = c("_original_n87", "_sensitivity_n79"), all = TRUE
)
print(comparison[order(comparison$padj_original_n87), ])
write.csv(comparison, here("results", "tables", "gsea_sensitivity_comparison.csv"), row.names = FALSE)
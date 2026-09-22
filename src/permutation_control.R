# permutation_control.R
#
# Random-subsampling control for sensitivity_analysis.R.
#
# sensitivity_analysis.R dropped 19 -> 11 chemoresistant (C) cases by
# excluding the 8 specifically biologic-treated ones, and several Hallmark
# pathways shifted a lot (E2F_TARGETS/MYC_TARGETS_V1 lost significance,
# the interferon/inflammatory pathways got much MORE significant). That
# result is consistent with two very different explanations:
#
#   (a) those 8 cases specifically carried the proliferation signal and/or
#       were specifically diluting the immune signal, or
#   (b) going from n=19 to n=11 C cases is just a large power/precision
#       change, and ANY random 8-of-19 drop would move these pathways
#       around by a similar amount, with no special role for biologics.
#
# This script can't fully separate these (that would need independent
# replication), but it gives a real empirical answer to "is this specific
# to the biologic-treated cases, or generic to shrinking n?": repeat the
# same drop-8-C-cases-and-rerun-GSEA procedure many times with RANDOMLY
# chosen sets of 8, and see where the real biologic-exclusion result falls
# relative to that random-draw distribution.
#
# Design notes:
#   - GSEA (via run_gsea()) ranks genes by the unshrunk Wald `stat`, not by
#     log2FoldChange, so lfcShrink()/apeglm is skipped here entirely - it
#     only affects LFC display, never the GSEA ranking, and skipping it
#     saves real time across many permutations.
#   - fgsea's own internal permutation testing is seeded via
#     cfg$gsea$seed inside run_gsea() and is IDENTICAL across every
#     iteration here (that's run_gsea()'s existing behavior, untouched) -
#     so the only thing varying between iterations is which 8 C cases were
#     dropped, not fgsea's internal randomness. That's intentional.
#   - Every iteration's gsea_results/gsea_plot files are written to the
#     same scratch path and overwritten - only the numbers we pull into
#     `perm_results` are kept. This is deliberate: keeping N figure sets
#     would just be clutter.
#   - choose(19, 8) = 75,582 possible combinations. n_perm below samples a
#     small fraction of that space at random (without repeating a
#     combination). This is NOT exhaustive - it's a Monte Carlo estimate of
#     "how much do pathway stats move around under a generic 19->11 drop."
#
# RUNTIME: each iteration is a full DESeq2 fit + fgsea run - roughly the
# same as one run of sensitivity_analysis.R. Start with a small n_perm
# (5-10) to sanity-check timing before committing to a larger run.

library(here)
library(DESeq2)

source(here("src", "data_loader.R"))
source(here("src", "gsea_analysis.R"))  # reuses the existing, working run_gsea()

# --- Config -------------------------------------------------------------

n_perm <- 50          # number of random 8-of-19 drops to run. Increase once
# you've confirmed timing on a small run.
perm_seed <- 20260922  # for reproducibility of WHICH combinations get drawn

chemoresistant_ids <- c(
  "TCGA-5M-AAT6", "TCGA-AA-3680", "TCGA-AA-3844", "TCGA-AA-3930", "TCGA-AA-3972",
  "TCGA-AA-A02K", "TCGA-AD-6964", "TCGA-AG-3584", "TCGA-AG-3999", "TCGA-AG-A016",
  "TCGA-AZ-4682", "TCGA-AZ-4684", "TCGA-AZ-6600", "TCGA-AZ-6606", "TCGA-CI-6620",
  "TCGA-F5-6702", "TCGA-NH-A6GA", "TCGA-NH-A6GB", "TCGA-RU-A8FL"
)
biologic_treated_C <- c(
  "TCGA-AA-3844", "TCGA-AA-3972", "TCGA-AA-A02K", "TCGA-AG-A016",
  "TCGA-AZ-4682", "TCGA-NH-A6GA", "TCGA-NH-A6GB", "TCGA-RU-A8FL"
)

compare_pathways <- c(
  "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT"
)

cfg <- load_config()

cfg_perm <- cfg
cfg_perm$paths <- lapply(cfg$paths, function(p) {
  sub("^results/", "results/permutation_control/", p)
})
dir.create(here(cfg_perm$paths$figures_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(here(cfg_perm$paths$gsea_results)), recursive = TRUE, showWarnings = FALSE)

# Same normalize-and-filter approach as sensitivity_analysis.R: filter on the
# full 87-sample matrix BEFORE subsetting columns, so gene row order stays
# identical to index_assign_colorectal.csv across every iteration.
counts_raw <- load_counts(cfg_perm)
counts_filtered <- filter_low_counts(counts_raw, cfg_perm$filtering$min_row_mean)
coldata_full <- build_coldata(counts_filtered, cfg_perm$design$condition_levels)

normalize_id <- function(x) gsub("[.-]", "-", x)
coldata_ids_norm <- normalize_id(rownames(coldata_full))

gene_index <- load_gene_index(cfg_perm)
contrast <- unlist(cfg_perm$design$contrast)

# --- Draw n_perm distinct random 8-of-19 combinations --------------------

set.seed(perm_seed)
used_keys <- character(0)
drawn_combos <- vector("list", n_perm)
for (i in seq_len(n_perm)) {
  repeat {
    candidate <- sort(sample(chemoresistant_ids, length(biologic_treated_C)))
    key <- paste(candidate, collapse = ",")
    if (!(key %in% used_keys)) {
      used_keys <- c(used_keys, key)
      break
    }
  }
  drawn_combos[[i]] <- candidate
}

# --- Run DESeq2 + GSEA once per random drop -------------------------------

perm_results <- list()
t_start <- Sys.time()

for (i in seq_len(n_perm)) {
  dropped_ids <- drawn_combos[[i]]
  n_overlap_with_biologic <- length(intersect(dropped_ids, biologic_treated_C))
  
  keep <- !(coldata_ids_norm %in% dropped_ids)
  n_excluded <- sum(!keep)
  if (n_excluded != length(dropped_ids)) {
    stop(sprintf(
      "permutation_control: iteration %d expected to exclude %d cases but matched %d.",
      i, length(dropped_ids), n_excluded
    ))
  }
  
  coldata_i <- coldata_full[keep, , drop = FALSE]
  counts_i <- counts_filtered[, keep]
  
  dds_i <- DESeqDataSetFromMatrix(countData = counts_i, colData = coldata_i, design = ~condition)
  dds_i <- suppressMessages(DESeq(dds_i, quiet = TRUE))
  
  # GSEA ranks by `stat`, which lfcShrink() never touches - skip shrinkage
  # here entirely (see header note) and pull `stat`/log2FoldChange straight
  # from results().
  res_i <- results(dds_i, contrast = contrast)
  res_df <- as.data.frame(res_i)
  
  if (nrow(res_df) != nrow(gene_index)) {
    stop(sprintf("permutation_control: iteration %d row-count mismatch vs gene_index.", i))
  }
  res_df$Index <- seq_len(nrow(res_df))
  merged_i <- merge(res_df, gene_index, by.x = "Index", by.y = "Index", all.x = TRUE)
  
  pipeline_out_i <- list(dds = dds_i, results = merged_i, cfg = cfg_perm)
  gsea_out_i <- suppressWarnings(suppressMessages(run_gsea(pipeline_out_i)))
  
  gsea_i <- gsea_out_i$results
  hit <- gsea_i[gsea_i$pathway %in% compare_pathways, c("pathway", "NES", "padj")]
  hit$permutation <- i
  hit$n_overlap_with_biologic <- n_overlap_with_biologic
  hit$dropped_ids <- paste(dropped_ids, collapse = "/")
  perm_results[[i]] <- hit
  
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  message(sprintf(
    "[%d/%d] dropped %d case(s) (%d overlap w/ real biologic set) - %.1f min elapsed, ~%.1f min remaining",
    i, n_perm, length(dropped_ids), n_overlap_with_biologic,
    elapsed, elapsed / i * (n_perm - i)
  ))
}

perm_long <- do.call(rbind, perm_results)
rownames(perm_long) <- NULL
dir.create(here("results", "permutation_control", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(perm_long, here("results", "permutation_control", "tables", "permutation_gsea_results.csv"), row.names = FALSE)

# --- Compare the REAL biologic-exclusion result against the random-draw
#     distribution, pathway by pathway ---------------------------------

gsea_orig <- read.csv(here(cfg$paths$gsea_results))
real_sensitivity <- read.csv(here("results", "tables", "gsea_sensitivity_comparison.csv"))

summary_rows <- lapply(compare_pathways, function(pw) {
  orig_row <- gsea_orig[gsea_orig$pathway == pw, c("NES", "padj")]
  real_row <- real_sensitivity[real_sensitivity$pathway == pw, ]
  perm_pw <- perm_long[perm_long$pathway == pw, ]
  
  data.frame(
    pathway = pw,
    NES_original_n87 = orig_row$NES,
    padj_original_n87 = orig_row$padj,
    NES_real_biologic_excl_n79 = real_row$NES_sensitivity_n79,
    padj_real_biologic_excl_n79 = real_row$padj_sensitivity_n79,
    NES_random_median = median(perm_pw$NES, na.rm = TRUE),
    NES_random_range = sprintf("[%.2f, %.2f]", min(perm_pw$NES, na.rm = TRUE), max(perm_pw$NES, na.rm = TRUE)),
    padj_random_median = median(perm_pw$padj, na.rm = TRUE),
    pct_random_sig_padj05 = round(100 * mean(perm_pw$padj <= 0.05, na.rm = TRUE), 1),
    n_random_draws = nrow(perm_pw)
  )
})
summary_df <- do.call(rbind, summary_rows)
print(summary_df)
write.csv(summary_df, here("results", "permutation_control", "tables", "permutation_summary.csv"), row.names = FALSE)

message(sprintf(
  "\nDone: %d random 8-of-19 drops completed in %.1f min. See results/permutation_control/tables/ for full output.",
  n_perm, as.numeric(difftime(Sys.time(), t_start, units = "mins"))
))
message(
  "How to read pct_random_sig_padj05: if it's high (most random 8-drops also keep the ",
  "pathway significant), the REAL biologic-exclusion result losing significance is unremarkable ",
  "- generic to n=11. If it's low (most random drops do NOT lose significance) but the real ",
  "biologic-exclusion result did, that's evidence the biologic-treated cases specifically mattered."
)
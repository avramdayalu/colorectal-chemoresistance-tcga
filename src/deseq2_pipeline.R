# deseq2_pipeline.R
#
# Builds the DESeqDataSet, fits the negative-binomial model, runs QC (PCA),
# extracts results for the S-vs-C contrast, applies LFC shrinkage (apeglm),
# and writes the result tables consumed by make_volcano_plots.R.
#
# Methods: Love MI, Huber W, Anders S (2014). Genome Biology, 15(12), 550.
#   https://doi.org/10.1186/s13059-014-0550-8 (DESeq2)
# Zhu A, Ibrahim JG, Love MI (2019). Bioinformatics, 35(12), 2084-2092.
#   https://doi.org/10.1093/bioinformatics/bty895 (apeglm shrinkage)

library(here)
library(DESeq2)
library(ggplot2)

source(here("src", "data_loader.R"))

run_pipeline <- function(cfg_path = here("config", "config.yaml")) {
  cfg <- load_config(cfg_path)

  counts_raw <- load_counts(cfg)
  counts_filtered <- filter_low_counts(counts_raw, cfg$filtering$min_row_mean)
  coldata <- build_coldata(counts_filtered, cfg$design$condition_levels)

  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = coldata,
    design = ~condition
  )
  dds <- DESeq(dds)

  # --- Quality control: PCA on variance-stabilized counts ---
  vsd <- vst(dds, blind = FALSE)
  pca_plot <- plotPCA(vsd, intgroup = "condition")
  ggsave(here(cfg$paths$pca_plot), pca_plot, width = 6, height = 5)

  # --- Differential expression ---
  contrast <- unlist(cfg$design$contrast)
  res <- results(dds, contrast = contrast)
  # IMPORTANT: do NOT na.omit() here. results(dds) returns exactly one row
  # per input gene, in the same order as counts_filtered - that fixed
  # row order is what index_assign.csv's Index column is built against.
  # Dropping NA rows before assigning Index (as the original scripts did)
  # silently shifts every gene *after* the first NA row by one position,
  # mislabeling it with the previous gene's symbol. Instead, attach Index
  # first (below), and only drop NA rows afterward for the filtered tables.

  # --- LFC shrinkage ---
  # Shrinks noisy, low-confidence log2FoldChange estimates (typically from
  # low-count genes) toward zero, so the fold-change values used for ranking
  # and plotting aren't dominated by low-confidence outliers. Does not alter
  # padj/p-values - those come from the Wald test above, before shrinkage.
  coef_name <- paste0("condition_", contrast[2], "_vs_", contrast[3])
  available_coefs <- resultsNames(dds)
  if (!(coef_name %in% available_coefs)) {
    stop(sprintf(
      "deseq2_pipeline: expected coefficient '%s' not found in resultsNames(dds): %s. Check design$contrast in config.yaml.",
      coef_name, paste(available_coefs, collapse = ", ")
    ))
  }
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = cfg$shrinkage$method, res = res)

  gene_index <- load_gene_index(cfg)
  res_df <- as.data.frame(res_shrunk)
  # Carry the pre-shrinkage Wald statistic through as its own column. It's
  # unused by the volcano/threshold logic (which wants the shrunk log2FC),
  # but it's the standard ranking metric for preranked GSEA (make_gsea.R) -
  # apeglm's shrinkage is meant for effect-size display/ranking-by-magnitude,
  # not for the enrichment test's gene ordering.
  res_df$stat <- res$stat

  if (nrow(res_df) != nrow(gene_index)) {
    stop(sprintf(
      "deseq2_pipeline: results row count (%d) does not match gene_index row count (%d) - the positional Index join is unsafe until these match (e.g. min_row_mean in config.yaml no longer matches how index_assign.csv was built).",
      nrow(res_df), nrow(gene_index)
    ))
  }
  res_df$Index <- seq_len(nrow(res_df))  # safe: attached before any row is dropped
  merged <- merge(res_df, gene_index, by.x = "Index", by.y = "Index", all.x = TRUE)

  write.csv(merged, here(cfg$paths$results_table), row.names = FALSE)

  # Filtered/significant-gene tables drop NA padj rows *after* labeling.
  filtered_05 <- merged[!is.na(merged$padj) & merged$padj <= 0.05, ]
  filtered_01 <- merged[!is.na(merged$padj) & merged$padj <= cfg$thresholds$padj_cutoff, ]
  write.csv(filtered_05, here(cfg$paths$filtered_padj05), row.names = FALSE)
  write.csv(filtered_01, here(cfg$paths$filtered_padj01), row.names = FALSE)

  # Final candidate list: padj AND log2FC criteria together (see config.yaml
  # comment). This is the headline gene list reported in the README and the
  # set labeled on the volcano plot.
  candidate_final <- merged[
    !is.na(merged$padj) &
      merged$padj <= cfg$thresholds$padj_cutoff &
      abs(merged$log2FoldChange) >= cfg$thresholds$log2fc_cutoff,
  ]
  candidate_final <- candidate_final[order(-abs(candidate_final$log2FoldChange)), ]
  write.csv(candidate_final, here(cfg$paths$candidate_genes_final), row.names = FALSE)

  list(dds = dds, results = merged, cfg = cfg)
}

# Allows `Rscript src/deseq2_pipeline.R` to run the pipeline directly, while
# `source()`-ing this file elsewhere (e.g. from make_volcano_plots.R) does not.
if (sys.nframe() == 0) {
  run_pipeline()
}

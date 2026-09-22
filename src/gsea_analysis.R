# gsea_analysis.R
#
# Preranked Gene Set Enrichment Analysis (GSEA) against MSigDB Hallmark gene
# sets. Unlike the padj/log2FC candidate-gene tables, this uses every gene
# in the results table ranked by its Wald statistic - the point of GSEA is
# to detect a pathway shifting coherently even when no single gene in it
# individually clears a significance threshold.
#
# Methods: Subramanian A, Tamayo P, Mootha VK, et al. (2005). PNAS,
#   102(43), 15545-15550. https://doi.org/10.1073/pnas.0506580102 (GSEA)
# Korotkevich G, Sukhov V, Budin N, et al. (2021). bioRxiv 060012.
#   https://doi.org/10.1101/060012 (fgsea implementation - a preprint, not
#   a peer-reviewed journal article, worth citing as such)
# Liberzon A, Birger C, Thorvaldsdottir H, et al. (2015). Cell Systems,
#   1(6), 417-425. https://doi.org/10.1016/j.cels.2015.12.004 (MSigDB Hallmark)

library(here)
library(fgsea)
library(msigdbr)
library(ggplot2)

source(here("src", "deseq2_pipeline.R"))

run_gsea <- function(pipeline_out = run_pipeline()) {
  cfg <- pipeline_out$cfg
  merged <- pipeline_out$results
  gcfg <- cfg$gsea

  # --- Build the ranked gene list ---
  # Rank by the pre-shrinkage Wald statistic (`stat`), not the shrunk
  # log2FoldChange: apeglm shrinkage is meant for effect-size display, and
  # compresses low-confidence genes toward zero in a way that would distort
  # a preranked ordering. `stat` reflects both direction and confidence.
  ranked <- merged[!is.na(merged$stat) & !is.na(merged$gene.symbol) & merged$gene.symbol != "", ]
  if (anyDuplicated(ranked$gene.symbol)) {
    # Keep the row with the largest |stat| for any symbol mapped to more
    # than one row (shouldn't happen with this dataset, but don't silently
    # let duplicate keys corrupt the named vector fgsea expects).
    ranked <- ranked[order(-abs(ranked$stat)), ]
    ranked <- ranked[!duplicated(ranked$gene.symbol), ]
  }
  gene_ranks <- setNames(ranked$stat, ranked$gene.symbol)
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)

  # --- MSigDB Hallmark gene sets (Homo sapiens) ---
  gene_sets <- msigdbr(species = "Homo sapiens", category = gcfg$collection)
  pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)

  # --- Preranked GSEA ---
  set.seed(gcfg$seed)
  gsea_res <- fgsea(
    pathways = pathways,
    stats = gene_ranks,
    minSize = gcfg$min_size,
    maxSize = gcfg$max_size
  )
  gsea_res <- as.data.frame(gsea_res[order(gsea_res$padj), ])
  # `leadingEdge` is a list-column (one character vector per row) - flatten
  # it to a single "/"-delimited string so write.csv doesn't mangle it.
  gsea_res$leadingEdge <- vapply(gsea_res$leadingEdge, paste, character(1), collapse = "/")
  write.csv(gsea_res, here(cfg$paths$gsea_results), row.names = FALSE)

  # --- Plot: significant pathways ranked by NES, colored by direction ---
  top <- gsea_res[!is.na(gsea_res$padj) & gsea_res$padj <= gcfg$padj_cutoff, ]
  if (nrow(top) == 0) {
    warning(sprintf(
      "run_gsea: no Hallmark pathways met padj <= %s; plotting the top %d by padj instead so the script still produces a figure.",
      gcfg$padj_cutoff, gcfg$n_plot
    ))
    top <- gsea_res[seq_len(min(gcfg$n_plot, nrow(gsea_res))), ]
  }
  top <- top[order(top$NES), ]
  top$pathway_label <- gsub("_", " ", sub("^HALLMARK_", "", top$pathway))
  top$pathway_label <- factor(top$pathway_label, levels = top$pathway_label)

  gsea_plot <- ggplot(top, aes(x = NES, y = pathway_label, fill = NES > 0)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(
      values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
      labels = c(`TRUE` = "Up in chemosensitive (S)", `FALSE` = "Up in chemoresistant (C)"),
      name = NULL
    ) +
    labs(
      x = "Normalized Enrichment Score (NES)", y = NULL,
      title = "Hallmark pathway enrichment (preranked GSEA)",
      subtitle = sprintf("Genes ranked by Wald statistic, S vs. C contrast (padj ≤ %s or top %d shown)", gcfg$padj_cutoff, gcfg$n_plot)
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 9, color = "grey30"))

  ggsave(here(cfg$paths$gsea_plot), gsea_plot, width = 9, height = max(4, 0.4 * nrow(top) + 2))

  # `pathways` and `gene_ranks` are returned alongside the results table so
  # leading_edge_analysis.R can reuse them (e.g. for fgsea::plotEnrichment())
  # without re-downloading Hallmark or rebuilding the ranking from scratch.
  invisible(list(results = gsea_res, plot = gsea_plot, pathways = pathways, gene_ranks = gene_ranks, pipeline_out = pipeline_out))
}

if (sys.nframe() == 0) {
  run_gsea()
}

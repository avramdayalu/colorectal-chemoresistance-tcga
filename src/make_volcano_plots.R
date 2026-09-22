# make_volcano_plots.R
#
# Produces the two volcano plots from the pipeline's merged results table:
# one with all genes, one restricted to the final candidate genes (padj <=
# thresholds$padj_cutoff AND |log2FC| >= thresholds$log2fc_cutoff). The
# candidate list is derived directly from the results table rather than
# hardcoded, so it can't drift out of sync with
# results/tables/candidate_genes_final.csv.
#
# Uses: Blighe K, Rana S, Lewis M. EnhancedVolcano (Bioconductor package,
#   no dedicated methods paper - cited as software).
#   https://bioconductor.org/packages/EnhancedVolcano

library(here)
library(EnhancedVolcano)
library(ggplot2)

source(here("src", "deseq2_pipeline.R"))

make_volcano_plots <- function(pipeline_out = run_pipeline()) {
  cfg <- pipeline_out$cfg
  merged <- pipeline_out$results

  y_upper <- max(-log10(merged$padj), na.rm = TRUE) + 1
  # Pad the x-axis beyond the data range so extreme-fold-change genes near
  # the edges have room for their labels instead of being clipped/cramped.
  x_extent <- max(abs(merged$log2FoldChange), na.rm = TRUE)
  x_lim <- c(-x_extent - 1, x_extent + 1)

  plot_all <- EnhancedVolcano(
    merged,
    x = "log2FoldChange", y = "padj", lab = merged$gene.symbol,
    pCutoff = cfg$thresholds$padj_cutoff,
    FCcutoff = cfg$thresholds$log2fc_cutoff,
    labSize = 3,
    xlim = x_lim,
    ylim = c(0, y_upper)
  )
  ggsave(here(cfg$paths$volcano_plot_all), plot_all, width = 10, height = 8)

  candidate_genes <- merged$gene.symbol[
    !is.na(merged$padj) &
      merged$padj <= cfg$thresholds$padj_cutoff &
      abs(merged$log2FoldChange) >= cfg$thresholds$log2fc_cutoff
  ]

  plot_labeled <- EnhancedVolcano(
    merged,
    x = "log2FoldChange", y = "padj", lab = merged$gene.symbol,
    pCutoff = cfg$thresholds$padj_cutoff,
    FCcutoff = cfg$thresholds$log2fc_cutoff,
    labSize = 5, selectLab = candidate_genes, legendLabSize = 12,
    legendPosition = "right",
    xlim = x_lim,
    ylim = c(0, y_upper),
    drawConnectors = TRUE, widthConnectors = 0.5,
    boxedLabels = TRUE, max.overlaps = Inf
  )
  ggsave(here(cfg$paths$volcano_plot_labeled), plot_labeled, width = 11, height = 8)

  invisible(list(all = plot_all, labeled = plot_labeled))
}

if (sys.nframe() == 0) {
  make_volcano_plots()
}

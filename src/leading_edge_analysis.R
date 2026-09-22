# leading_edge_analysis.R
#
# Follow-up diagnostics on the GSEA hits from gsea_analysis.R, for the
# pathways listed in config$leading_edge$focus_pathways:
#
#   1. The classic GSEA running-enrichment-score curve (fgsea::plotEnrichment)
#      - shows the enrichment isn't an artifact of a handful of extreme genes
#      at the tail of the ranking.
#   2. The leading-edge genes themselves (the subset actually driving the
#      score), written out as a table.
#   3. A patient-by-patient heatmap of those leading-edge genes' normalized
#      expression.
#
# IMPORTANT caveat on (3): this is a consistency check, not independent
# validation. Leading-edge genes are selected because they differ between
# chemoresistant and chemosensitive patients, so a heatmap of exactly those
# genes clustering by group is expected almost by construction - it doesn't
# independently confirm the pathway result. What it IS useful for: seeing
# whether the signal holds up broadly across patients, or is actually being
# driven by a handful of outliers - which matters here specifically because
# the whole-transcriptome PCA (see make_volcano_plots.R output / README)
# did NOT show a clean separation.
#
# Methods: enrichment-curve/leading-edge concept and fgsea implementation are
#   the same as gsea_analysis.R (see that file's header for citations).
# Kolde R. pheatmap (CRAN package, no dedicated methods paper - cited as
#   software). https://cran.r-project.org/package=pheatmap

library(here)
library(fgsea)
library(DESeq2)
library(pheatmap)
library(ggplot2)

source(here("src", "gsea_analysis.R"))

run_leading_edge_analysis <- function(gsea_out = run_gsea()) {
  cfg <- gsea_out$pipeline_out$cfg
  dds <- gsea_out$pipeline_out$dds
  gsea_res <- gsea_out$results
  pathways <- gsea_out$pathways
  gene_ranks <- gsea_out$gene_ranks
  focus <- cfg$leading_edge$focus_pathways

  missing_pw <- setdiff(focus, gsea_res$pathway)
  if (length(missing_pw) > 0) {
    stop(sprintf(
      "leading_edge_analysis: pathway(s) not found in GSEA results: %s. Check config$leading_edge$focus_pathways against results/tables/gsea_hallmark_results.csv.",
      paste(missing_pw, collapse = ", ")
    ))
  }

  vsd <- vst(dds, blind = FALSE)
  condition <- colData(dds)$condition
  annotation_col <- data.frame(Group = condition, row.names = colnames(dds))

  # dds/vsd carry no gene-symbol rownames at all - only row *position*
  # (1..N, matching counts_filtered's row order). Gene symbols are attached
  # separately, after the fact, via the Index join in deseq2_pipeline.R. So
  # subsetting assay(vsd) by gene symbol directly silently matches nothing;
  # go through the same Index lookup the rest of the pipeline uses.
  merged <- gsea_out$pipeline_out$results
  index_by_symbol <- setNames(merged$Index, merged$gene.symbol)

  leading_edge_rows <- list()

  for (pw in focus) {
    pw_label <- gsub("_", " ", sub("^HALLMARK_", "", pw))
    row <- gsea_res[gsea_res$pathway == pw, ]

    # --- 1. Enrichment curve ---
    curve <- plotEnrichment(pathways[[pw]], gene_ranks) +
      ggplot2::labs(
        title = sprintf("%s (NES = %.2f, padj = %.2e)", pw_label, row$NES, row$padj),
        subtitle = "Genes ranked by Wald statistic, S vs. C contrast"
      )
    ggsave(here(cfg$paths$figures_dir, sprintf("gsea_enrichment_%s.png", tolower(gsub(" ", "_", pw_label)))),
           curve, width = 8, height = 5)

    # --- 2. Leading-edge genes ---
    # gsea_analysis.R flattens the leadingEdge list-column to a "/"-delimited
    # string before writing gsea_hallmark_results.csv - split it back out.
    le_genes <- strsplit(row$leadingEdge, "/")[[1]]
    leading_edge_rows[[pw]] <- data.frame(pathway = pw, gene_symbol = le_genes)

    # --- 3. Leading-edge heatmap (see caveat in the file header) ---
    gene_idx <- index_by_symbol[le_genes]
    unmatched <- le_genes[is.na(gene_idx)]
    if (length(unmatched) > 0) {
      warning(sprintf(
        "leading_edge_analysis: %d/%d leading-edge gene(s) for %s not found in the results table by symbol (skipped): %s",
        length(unmatched), length(le_genes), pw, paste(unmatched, collapse = ", ")
      ))
    }
    gene_idx <- gene_idx[!is.na(gene_idx)]
    mat <- assay(vsd)[gene_idx, , drop = FALSE]
    rownames(mat) <- names(gene_idx)  # relabel rows with gene symbols for the plot
    pheatmap(
      mat,
      scale = "row",
      annotation_col = annotation_col,
      show_rownames = TRUE,
      main = sprintf("%s - leading-edge genes (per-patient consistency check)", pw_label),
      filename = here(cfg$paths$figures_dir, sprintf("gsea_heatmap_%s.png", tolower(gsub(" ", "_", pw_label)))),
      width = 8, height = max(4, 0.25 * length(le_genes) + 2)
    )
  }

  leading_edge_table <- do.call(rbind, leading_edge_rows)
  rownames(leading_edge_table) <- NULL
  write.csv(leading_edge_table, here(cfg$paths$leading_edge_genes), row.names = FALSE)

  invisible(leading_edge_table)
}

if (sys.nframe() == 0) {
  run_leading_edge_analysis()
}

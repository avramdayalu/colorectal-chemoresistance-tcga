# ora_analysis.R
#
# Complementary Over-Representation Analysis (ORA): a threshold-restricted,
# hypergeometric-test check (clusterProfiler::enrichGO), run separately on
# genes up in chemosensitive (S) patients and genes up in chemoresistant (C)
# patients. If this threshold-restricted method and the threshold-free GSEA
# (gsea_analysis.R) flag the same biological themes, that's a genuine
# cross-check between two different statistical approaches on the same data
# (still not an independent dataset - see README).
#
# Two corrections vs. a naive ORA setup:
#   1. `universe` is set explicitly to every gene actually tested by DESeq2
#      (all rows of the results table), not left at enrichGO()'s default of
#      the whole org.Hs.eg.db annotation. Without this, the hypergeometric
#      test's background is "the whole genome" instead of "the genes we
#      could actually have detected," which inflates significance.
#   2. Genes are split by the *sign* of log2FoldChange under this repo's own
#      contrast (S vs C - see config$design$contrast): positive means up in
#      S, negative means up in C. Getting this backwards is an easy mistake
#      when adapting ORA snippets written against a different contrast.
#
# Methods: Yu G, Wang LG, Han Y, He QY (2012). OMICS, 16(5), 284-287.
#   https://doi.org/10.1089/omi.2011.0118; updated in Wu T, Hu E, Xu S, et al.
#   (2021). The Innovation, 2(3), 100141. (clusterProfiler)
# Ashburner M, Ball CA, Blake JA, et al. (2000). Nature Genetics, 25(1),
#   25-29. https://doi.org/10.1038/75556 (Gene Ontology)

library(here)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

source(here("src", "deseq2_pipeline.R"))

run_ora <- function(pipeline_out = run_pipeline()) {
  cfg <- pipeline_out$cfg
  merged <- pipeline_out$results
  ocfg <- cfg$ora

  universe <- merged$gene.symbol[!is.na(merged$gene.symbol) & merged$gene.symbol != ""]

  sig <- merged[!is.na(merged$padj) & merged$padj <= ocfg$padj_cutoff, ]
  genes_up_in_S <- sig$gene.symbol[sig$log2FoldChange >= ocfg$log2fc_cutoff]
  genes_up_in_C <- sig$gene.symbol[sig$log2FoldChange <= -ocfg$log2fc_cutoff]

  message(sprintf(
    "run_ora: %d genes up in S, %d genes up in C at padj<=%s & |log2FC|>=%s (universe: %d genes)",
    length(genes_up_in_S), length(genes_up_in_C), ocfg$padj_cutoff, ocfg$log2fc_cutoff, length(universe)
  ))

  run_one <- function(genes, direction_label, results_path, plot_path) {
    if (length(genes) < 3) {
      warning(sprintf(
        "run_ora: only %d gene(s) up in %s at this threshold - too few for a meaningful GO test. Loosen ora$padj_cutoff/log2fc_cutoff in config.yaml if you want a result here.",
        length(genes), direction_label
      ))
      return(NULL)
    }
    ora_res <- enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = ocfg$ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
    )
    write.csv(as.data.frame(ora_res), here(results_path), row.names = FALSE)

    if (nrow(as.data.frame(ora_res)) > 0) {
      p <- dotplot(ora_res, showCategory = ocfg$show_categories) +
        ggtitle(sprintf("ORA (GO:%s) - genes up in %s", ocfg$ontology, direction_label))
      ggsave(here(plot_path), p, width = 9, height = max(4, 0.3 * ocfg$show_categories + 2))
    } else {
      warning(sprintf("run_ora: enrichGO returned no significant terms for genes up in %s.", direction_label))
    }
    ora_res
  }

  ora_S <- run_one(genes_up_in_S, "S (chemosensitive)", cfg$paths$ora_up_in_S, cfg$paths$ora_dotplot_up_in_S)
  ora_C <- run_one(genes_up_in_C, "C (chemoresistant)", cfg$paths$ora_up_in_C, cfg$paths$ora_dotplot_up_in_C)

  invisible(list(up_in_S = ora_S, up_in_C = ora_C))
}

if (sys.nframe() == 0) {
  run_ora()
}

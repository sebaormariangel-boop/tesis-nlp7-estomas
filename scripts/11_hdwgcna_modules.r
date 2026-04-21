# 11_hdwgcna_modules.R
# Análisis de módulos hdWGCNA:
#   • ModuleExprScore (UCell) + FeaturePlots
#   • DMEs (FindDMEs) + volcano
#   • ModuleTraitCorrelation
#   • GO enrichment por módulo
#   • Overlap módulo–DEGs con genesectR
#
# Entrada:  qs2/*_hdWGCNA_object.qs2  (salida de 10_hdwgcna_network.R)
# Salida:   figures/ucell_modules.png, plol_dms.png, GO_BP_dotplot.pdf ...
#           tables/trt_comparison_cons.txt, module_scores_*.tsv ...
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

suppressPackageStartupMessages({
  library(hdWGCNA)
  library(WGCNA)
  library(genesectR)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(enrichplot)
})

for (d in c("figures", "tables", "qs2")) dir.create(d, showWarnings = FALSE)

# 1. Cargar objeto qs2 -------------------------------------------------------

coexpr_mode <- "subsetting"   # "full" o "subsetting" — debe coincidir con 10_

in_base <- paste(
  if (grepl("TB$", project)) "ten" else "nat",
  coexpr_mode, "hdWGCNA_object",
  sep = "_"
)
seurat_obj <- qs_read(file.path("qs2", paste0(in_base, ".qs2")))
wgcna_name <- seurat_obj@misc$active_wgcna
info(seurat_obj)

TRT <- trt_map(tech_mode)
lev1 <- TRT$trt
lev2 <- TRT$ctrl

# 2. Module scores (UCell, por célula) ----------------------------------------

seurat_obj <- ModuleExprScore(seurat_obj, n_genes = 25, method = "UCell")

scores_mat <- GetModuleScores(seurat_obj)
score_cols <- paste0("score_", colnames(scores_mat))
colnames(scores_mat) <- score_cols
seurat_obj <- AddMetaData(seurat_obj, metadata = scores_mat)

# Feature plots (hMEs)
reduction_plot <- if (coexpr_mode == "subsetting") "umap.ct" else "umap.harmony"

plot_list <- ModuleFeaturePlot(
  seurat_obj, reduction = reduction_plot,
  features = "hMEs", order = TRUE
)
p_feats <- patchwork::wrap_plots(plot_list, ncol = 6)
save_plot_px(p_feats, file.path("figures", "ucell_modules.png"), w = 3500)

plot_list2 <- ModuleFeaturePlot(
  seurat_obj, reduction = reduction_plot,
  features = "scores", order = "shuffle", ucell = TRUE
)
p_feats2 <- patchwork::wrap_plots(plot_list2, ncol = 6)
save_plot_px(p_feats2, file.path("figures", "ucell_hubgenes_modules.png"), w = 3500)

# Guardar scores
save_tsv(
  GetModuleScores(seurat_obj),
  file.path("tables", paste0("module_scores_", coexpr_mode, "_", project, ".tsv"))
)

# 3. DotPlot de módulos (hMEs) ------------------------------------------------

MEs      <- GetMEs(seurat_obj, harmonized = TRUE)
modules  <- GetModules(seurat_obj)
mods     <- levels(modules$module)
mods     <- setdiff(mods, "grey")

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p_dot <- DotPlot(seurat_obj, features = mods, group.by = "annotation") +
  RotatedAxis() +
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")

save_plot_px(p_dot, file.path("figures", "dotplot_modules.png"), w = 3500)

# RadarPlot por tipo celular
p_radar <- ModuleRadarPlot(
  seurat_obj,
  group.by          = "annotation",
  barcodes          = colnames(seurat_obj),
  axis.label.size   = 2,
  grid.label.size   = 4,
  axis.label.offset = 0.8,
  plot.extent.x.sf  = 1.15
) + coord_cartesian(clip = "off") +
    theme(plot.margin = margin(10, 25, 10, 25))

save_plot_px(p_radar, file.path("figures", "radarplot_modules.png"))

# 4. Comparación por tratamiento — enfoque conservador (por muestra) ----------

cell_meta <- seurat_obj[[]]
MEs_cells <- GetMEs(seurat_obj, harmonized = TRUE)

df_me <- cbind(
  cell_meta[, c("Sample", "Treatment"), drop = FALSE],
  as.data.frame(MEs_cells[rownames(cell_meta), , drop = FALSE])
)

df_sample <- df_me |>
  dplyr::group_by(Sample, Treatment) |>
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean, na.rm = TRUE),
                   .groups = "drop")

me_cols <- setdiff(colnames(df_sample), c("Sample", "Treatment"))

res_sample <- tibble::tibble(
  module    = me_cols,
  mean_lev1 = sapply(me_cols, \(m) mean(df_sample[[m]][df_sample$Treatment == lev1], na.rm = TRUE)),
  mean_lev2 = sapply(me_cols, \(m) mean(df_sample[[m]][df_sample$Treatment == lev2], na.rm = TRUE))
) |>
  dplyr::mutate(diff = mean_lev1 - mean_lev2) |>
  dplyr::filter(!grepl("grey", module, ignore.case = TRUE)) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    p_t = tryCatch(
      t.test(df_sample[[module]] ~ df_sample$Treatment)$p.value,
      error = \(e) NA_real_
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(p_adj = p.adjust(p_t, method = "BH")) |>
  dplyr::arrange(p_adj, dplyr::desc(abs(diff)))

print(head(res_sample, 15))
save_tsv(res_sample, file.path("tables", "trt_comparison_cons.tsv"))

# 5. DME analysis (estilo hdWGCNA) -------------------------------------------

barcodes1 <- rownames(cell_meta)[cell_meta$Treatment == lev1]
barcodes2 <- rownames(cell_meta)[cell_meta$Treatment == lev2]

DMEs <- FindDMEs(seurat_obj, barcodes1 = barcodes1, barcodes2 = barcodes2,
                 test.use = "wilcox")

p_lol <- PlotDMEsLollipop(seurat_obj, DMEs,
                           wgcna_name = wgcna_name, pvalue = "p_val_adj")
save_plot_px(p_lol, file.path("figures", "dmes_lollipop.png"), h = 1500)

p_vol <- PlotDMEsVolcano(seurat_obj, DMEs, wgcna_name = wgcna_name)
save_plot_px(p_vol, file.path("figures", "dmes_volcano.png"))

# 6. Module-trait correlation -------------------------------------------------

seurat_obj$Treatment_num <- as.integer(seurat_obj$Treatment == lev1)

seurat_obj$celltype <- droplevels(factor(seurat_obj$celltype))

seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits   = c("Treatment_num", "nFeature_RNA", "nCount_RNA"),
  features = "hMEs",
  cor_meth = "spearman",
  group.by = seurat_obj$celltype
)

p_mtc <- PlotModuleTraitCorrelation(
  seurat_obj,
  label        = "fdr",
  label_symbol = "stars",
  text_size    = 2,
  text_digits  = 2,
  text_color   = "white",
  high_color   = "yellow",
  mid_color    = "black",
  low_color    = "purple",
  plot_max     = 0.2,
  combine      = TRUE
)
save_plot_px(p_mtc, file.path("figures", "module_trait_correlation.png"),
             w = 2000, h = 2500)

# 7. GO enrichment por módulo (BP) -------------------------------------------

modules_tbl <- GetModules(seurat_obj) |> dplyr::filter(module != "grey")
gene_sets   <- split(modules_tbl$gene_name, modules_tbl$module)
universe    <- seurat_obj@misc[[wgcna_name]]$wgcna_genes

ego_bp <- clusterProfiler::compareCluster(
  geneCluster   = gene_sets,
  fun           = "enrichGO",
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  universe      = universe,
  pAdjustMethod = "BH"
)

p_go <- enrichplot::dotplot(ego_bp, showCategory = 10, font.size = 8)
ggsave(file.path("figures", "GO_BP_dotplot.pdf"), p_go, width = 14, height = 30)

save_tsv(modules_tbl,
  file.path("tables", paste0(coexpr_mode, "_modules_nogrey.tsv")))
save_tsv(GetModules(seurat_obj),
  file.path("tables", paste0(coexpr_mode, "_modules_full.tsv")))

# 8. Overlap módulo–DEGs con genesectR ----------------------------------------
#
#  Para poder ejecutar este bloque se requieren los archivos de DEGs:
#    deg_leaf_path  → tabla TSV de DEGs en hoja   (make_deg_df)
#    deg_stoma_path → tabla TSV de DEGs en estomas (make_deg_df)

deg_leaf_path  <- file.path(bulk_leaf_dir,
  "DE_all_Genotipo_en__Well_watered_nlp7_vs_WT.tsv")
deg_stoma_path <- file.path(bulk_stoma_dir,
  "DE_all_ESTOMA__Genotipo_en_Estoma_SinN_nlp7_vs_WT.tsv")

make_deg_df <- function(path, contrast_name, baseMean_min = 10) {
  readr::read_tsv(path, show_col_types = FALSE) |>
    dplyr::transmute(
      gene      = gene,
      avg_log2FC = log2FoldChange,
      p_val_adj  = padj,
      baseMean   = baseMean,
      contrast   = contrast_name
    ) |>
    dplyr::filter(!is.na(p_val_adj), baseMean >= baseMean_min)
}

if (file.exists(deg_leaf_path) && file.exists(deg_stoma_path)) {
  deg_leaf  <- make_deg_df(deg_leaf_path,  "leaf_WW_nlp7_vs_WT")
  deg_stoma <- make_deg_df(deg_stoma_path, "stoma_SinN_nlp7_vs_WT")
  deg_df    <- dplyr::bind_rows(deg_leaf, deg_stoma)

  fc <- 0.58
  deg_sets <- list(
    leaf_UP    = with(subset(deg_df, contrast == "leaf_WW_nlp7_vs_WT"  & p_val_adj <= 0.05 & avg_log2FC >=  fc), unique(gene)),
    leaf_DOWN  = with(subset(deg_df, contrast == "leaf_WW_nlp7_vs_WT"  & p_val_adj <= 0.05 & avg_log2FC <= -fc), unique(gene)),
    stoma_UP   = with(subset(deg_df, contrast == "stoma_SinN_nlp7_vs_WT" & p_val_adj <= 0.05 & avg_log2FC >=  fc), unique(gene)),
    stoma_DOWN = with(subset(deg_df, contrast == "stoma_SinN_nlp7_vs_WT" & p_val_adj <= 0.05 & avg_log2FC <= -fc), unique(gene))
  )

  wgcna_universe <- unique(seurat_obj@misc[[wgcna_name]]$wgcna_genes)
  bulk_universe  <- unique(deg_df$gene)
  master_set     <- intersect(bulk_universe, wgcna_universe)

  sets_all <- c(gene_sets, deg_sets)
  sets_all <- lapply(sets_all, \(x) intersect(unique(as.character(x)), master_set))

  gs <- genesectR::gs_import(sets_all, master_set)
  gs <- genesectR::gs_compute_matricies(gs)

  png(file.path("figures", "overlap_genesectR_modules_vs_degs.png"),
      width = 6000, height = 5000, res = 500)
  genesectR::gs_multi_plot(gs, names(gene_sets))
  dev.off()

  png(file.path("figures", "overlap_genesectR_degsets.png"),
      width = 6000, height = 5000, res = 500)
  genesectR::gs_multi_plot(gs, names(deg_sets))
  dev.off()
} else {
  message("Archivos de DEGs no encontrados — se omite el análisis genesectR.")
  message("  ", deg_leaf_path)
  message("  ", deg_stoma_path)
}

# 9. Guardar objeto actualizado -----------------------------------------------

save_obj(seurat_obj, file.path("qs2", in_base), ext = "qs2")
message("Objeto guardado: qs2/", in_base, ".qs2")

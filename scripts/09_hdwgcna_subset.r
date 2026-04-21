# 09_hdwgcna_subset.R
# Preparación del objeto para hdWGCNA:
#   • Elegir dataset (Tenorio / Natanella) y modo de coexpresión
#   • Cargar objeto anotado desde 08_annotation.R
#   • Subset del tipo celular de interés
#   • Normalización, PCA, Harmony y clustering del subset
#   • Proporciones por tratamiento
#
# Entrada:  *_final_annotated.qs2 / .rds  (salida de 08_annotation.R)
# Salida:   *_<mode>_subclustered.qs2
#           figures/<ct>_clustree_<project>.png
#           figures/<ct>_umap_<project>.png
#           tables/gc_cluster_props_*.csv
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

suppressPackageStartupMessages({
  library(hdWGCNA)
  library(WGCNA)
})
enableWGCNAThreads(nThreads = 56)

for (d in c("figures", "tables", "qs2")) dir.create(d, showWarnings = FALSE)

# 1. Selección interactiva ----------------------------------------------------

{
  tech_opt <- c("Single-cell | Tenorio Berrío et al.",
                "Single-nucleus | Illouz-Eliaz et al.")
  tech_ch  <- menu(tech_opt, title = bold("Seleccionar tecnología:"))
  if (tech_ch == 0) stop("Cancelado.")
  tech_mode <- c("sc", "sn")[tech_ch]

  mode_opt <- c("Dataset completo", "Subset de tipo celular")
  mode_ch  <- menu(mode_opt, title = bold("Modo de coexpresión:"))
  if (mode_ch == 0) stop("Cancelado.")
  coexpr_mode <- c("full", "subsetting")[mode_ch]

  message(bold(red("==> ")), tech_opt[tech_ch],
          " | modo: ", mode_opt[mode_ch])
}

# 2. Cargar objeto anotado ----------------------------------------------------
#
#  Ajusta el path al archivo producido por 08_annotation.R.
#  Si tienes find_files(), úsalo; si no, asigna directamente:

obj_path <- if (tech_mode == "sc") {
  file.path(qs2_dir, "tenorio_final_annotated.qs2")
} else {
  file.path(qs2_dir, "natanella_final_annotated.qs2")
}

seurat_obj <- if (grepl("\\.qs2$", obj_path)) {
  qs_read(obj_path)
} else {
  readRDS(obj_path)
}

DefaultAssay(seurat_obj) <- "RNA"
info(seurat_obj)

# 3. Elegir tipo celular de interés -------------------------------------------

if (interactive_mode) {
  message("Tipos celulares disponibles:")
  print(unique(seurat_obj$celltype))
  ct <- trimws(readline(prompt = "Tipo celular: "))
} else {
  ct <- "Guard Cells"   # ← CAMBIAR si no es interactivo
}

stopifnot(ct %in% seurat_obj$celltype)

ct_subset <- subset(seurat_obj, subset = celltype == ct)
DefaultAssay(ct_subset) <- "RNA"
info(ct_subset)

# 4. Preprocesar objeto completo (opcional, solo para modo "full") -------------

if (coexpr_mode == "full") {
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- NormalizeData(seurat_obj) |>
    FindVariableFeatures(nfeatures = 2000, selection.method = "vst") |>
    ScaleData()
  info(seurat_obj)
}

# 5. Normalización e integración del subset -----------------------------------

ct_subset <- NormalizeData(ct_subset) |>
  FindVariableFeatures(nfeatures = 2000, selection.method = "vst") |>
  ScaleData()

ct_subset <- RunPCA(ct_subset,
  features = VariableFeatures(ct_subset), npcs = 50
)

ct_subset <- RunHarmony(
  ct_subset,
  group.by.vars = "Sample",
  reduction.use = "pca",
  dims.use      = DIMS_USE
)

# 6. Exploración de resoluciones (clustree) -----------------------------------

DefaultAssay(ct_subset) <- "RNA"

ct_subset <- FindNeighbors(
  ct_subset,
  reduction  = "harmony",
  dims       = DIMS_USE,
  graph.name = c("ct_nn", "ct_snn")
)

res_grid <- c(1e-6, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

p_clustree <- clustree_plot(
  ct_subset, res_grid,
  prefix     = "ct_clusters.",
  graph.name = "ct_snn"
)
save_plot_px(
  p_clustree,
  file.path("figures", paste0(ct, "_clustree_", project, ".png"))
)

# 7. Clustering y UMAP final del subset ---------------------------------------

CLUST_RES_CT <- switch(tech_mode, sc = 0.7, sn = 0.6)

ct_subset <- FindClusters(
  ct_subset,
  resolution   = CLUST_RES_CT,
  cluster.name = "ct_clusters",
  algorithm    = 4,
  random.seed  = SEED,
  graph.name   = "ct_snn"
) |>
  RunUMAP(
    reduction      = "harmony",
    dims           = DIMS_USE,
    seed.use       = SEED,
    reduction.name = "umap.ct",
    umap.method    = "umap-learn",
    metric         = "correlation"
  )

p_trt <- DimPlot(ct_subset, reduction = "umap.ct",
                 group.by = "Treatment", shuffle = TRUE) + umap_theme()
p_clu <- DimPlot(ct_subset, reduction = "umap.ct",
                 group.by = "ct_clusters", shuffle = TRUE) +
  umap_theme() + ggtitle(paste0(ct, " UMAP"))

print(p_trt + p_clu)
save_plot_px(p_trt + p_clu,
  file.path("figures", paste0(ct, "_umap_", project, ".png")),
  w = 2500, h = 1250
)

# 8. Proporciones por tratamiento + muestra ------------------------------------

TRT <- trt_map(tech_mode)

sample_trt <- ct_subset@meta.data |>
  dplyr::distinct(Sample, Treatment)

df_prop <- ct_subset@meta.data |>
  dplyr::mutate(cluster = Idents(ct_subset)) |>
  dplyr::count(Sample, cluster, name = "n") |>
  tidyr::complete(Sample, cluster, fill = list(n = 0)) |>
  dplyr::left_join(sample_trt, by = "Sample") |>
  dplyr::group_by(Sample) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::ungroup()

wide <- df_prop |>
  dplyr::group_by(Treatment, cluster) |>
  dplyr::summarise(mean_prop = mean(prop), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Treatment, values_from = mean_prop,
                     values_fill = 0)

eps <- 1e-6
wide2 <- wide |>
  dplyr::mutate(
    delta_pp = 100 * (.data[[TRT$trt]] - .data[[TRT$ctrl]]),
    log2FC   = log2((.data[[TRT$trt]] + eps) / (.data[[TRT$ctrl]] + eps)),
    FC       = 2^log2FC
  ) |>
  dplyr::arrange(dplyr::desc(delta_pp))

save_csv(df_prop, file.path("tables", "ct_cluster_props_by_sample.csv"))
save_csv(wide,   file.path("tables", "ct_cluster_props_wide.csv"))
save_csv(wide2,  file.path("tables", "ct_cluster_props_contrasts.csv"))

# 9. Añadir anotación al objeto principal -------------------------------------

ct_subset$annotation <- paste(ct, Idents(ct_subset), sep = "_")

if (coexpr_mode == "full") {
  seurat_obj$annotation <- NA_character_
  seurat_obj$annotation[colnames(ct_subset)] <- ct_subset$annotation
} else {
  seurat_obj <- ct_subset
}

info(seurat_obj)

# 10. Guardar -----------------------------------------------------------------

out_base <- paste(
  if (grepl("TB$", project)) "ten" else "nat",
  coexpr_mode, "subclustered",
  sep = "_"
)
save_obj(seurat_obj, file.path("qs2", out_base), ext = "qs2")
message("Objeto guardado: qs2/", out_base, ".qs2")

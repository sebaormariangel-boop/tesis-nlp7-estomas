# 05_integration.R
# Normalización del ensayo "bulkOK", selección de HVGs, PCA,
# integración por Harmony, exploración de dimensiones y UMAP final.
#
# Entrada:  act_bulkOK.qs2
# Salida:   act_integrated.qs2
#           harmony_dims_10_50_step2.pdf   (diagnóstico de dims)
#           clustree_<project>.png
#           umap_integrated_<project>.png
#           dotplot_<project>.png
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

# Parámetros (ajustar tras inspección visual) ---------------------------------

DIMS_FINAL <- 1:32    # dimensiones para UMAP y clustering definitivo
RES_FINAL  <- 0.8     # resolución de clustering final (0.8 nat; 0.7 ten)

# 1. Cargar objeto ------------------------------------------------------------

act <- qs_read("act_bulkOK.qs2")
DefaultAssay(act) <- "bulkOK"
act[["bulkOK"]]   <- split(act[["bulkOK"]], f = act$Sample)
info(act)

# 2. Normalización, HVGs y escalado -------------------------------------------

act <- NormalizeData(act)
act <- FindVariableFeatures(act, nfeatures = 2000, selection.method = "vst")

p_var <- VariableFeaturePlot(act)
save_plot_px(p_var,
  file.path(fig_dir, paste0("hvgs_", project, ".png")),
  w = 3500, h = 3000
)

act <- ScaleData(act)
info(act)

# 3. PCA y exploración --------------------------------------------------------

act <- RunPCA(act, npcs = 50, verbose = FALSE)

png(file.path(fig_dir, paste0("dimheatmap_pcs28_50_", project, ".png")),
    width = 2000, height = 3000, res = 150)
DimHeatmap(act, dims = 28:50, cells = 500, balanced = TRUE)
dev.off()

p_elbow <- ElbowPlot(act, ndims = 50)
save_plot_px(p_elbow,
  file.path(fig_dir, paste0("elbowplot_", project, ".png")),
  w = 1800, h = 1200
)

# 4. UMAP sin integrar (referencia) -------------------------------------------

act <- FindNeighbors(act, dims = 1:50, reduction = "pca")
act <- FindClusters(act, resolution = 1, cluster.name = "unintegrated_clusters")
act <- RunUMAP(
  act,
  reduction = "pca", dims = 1:50,
  reduction.name = "umap.unintegrated",
  umap.method = "umap-learn", seed.use = SEED, metric = "correlation"
)

p_unint <- DimPlot(
  act, reduction = "umap.unintegrated",
  group.by = c("Treatment", "Sample", "unintegrated_clusters"),
  ncol = 3, label = TRUE, alpha = 0.3, pt.size = 0.3, shuffle = TRUE
)
save_plot_px(p_unint,
  file.path(fig_dir, paste0("umap_unintegrated_", project, "_50pcs.png")),
  w = 6500, h = 2000
)

# 5. Integración con Harmony --------------------------------------------------

DefaultAssay(act) <- "bulkOK"

capas_sep <- (
  length(grep("data\\.",   Layers(act), value = TRUE)) > 1 &&
  length(grep("counts\\.", Layers(act), value = TRUE)) > 1
)

if (capas_sep) {
  message("Capas separadas — usando IntegrateLayers()")
  act <- IntegrateLayers(
    act,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    theta          = 4,
    lambda         = 0.5,
    sigma          = 0.1,
    verbose        = FALSE
  )
} else {
  message("Capas unidas — usando RunHarmony()")
  act <- RunHarmony(
    act, group.by.vars = "Sample",
    reduction.use = "pca", theta = 2,
    dims.use = DIMS_FINAL,
    .options = harmony_options(
      max.iter.cluster = 100,
      epsilon.cluster  = 1e-4,
      epsilon.harmony  = 1e-4
    )
  )
}

# 6. Exploración de dimensiones (dims_grid → PNG por cada d → PDF) -----------
#    Genera un UMAP + DotPlot para cada número de dims en dims_grid.
#    Útil para elegir DIMS_FINAL antes de correr el bloque definitivo.

dims_grid <- seq(10, 50, by = 2)
out_dim_dir <- file.path(fig_dir, paste0("png_harmony_dims_", project))
dir.create(out_dim_dir, showWarnings = FALSE, recursive = TRUE)

for (d in dims_grid) {
  dims_use <- 1:d

  act <- FindNeighbors(
    act, reduction = "harmony", dims = dims_use,
    graph.name = c("harmony_nn", "harmony_snn")
  )
  act <- FindClusters(
    act, graph.name = "harmony_snn",
    resolution = 1, cluster.name = "harmony_clusters",
    algorithm = 4, random.seed = SEED
  )
  act <- RunUMAP(
    act, reduction = "harmony", dims = dims_use,
    seed.use = SEED, reduction.name = "umap.harmony",
    umap.method = "umap-learn", metric = "correlation",
    reduction.key = "UMAPH_"
  )

  DefaultAssay(act) <- "RNA"
  if (!"data" %in% Layers(act[["RNA"]])) act <- NormalizeData(act, assay = "RNA")

  p_u <- DimPlot(
    act, reduction = "umap.harmony",
    group.by = c("Treatment", "Sample", "harmony_clusters"),
    ncol = 3, label = TRUE, alpha = 0.3, pt.size = 0.3, shuffle = TRUE
  ) +
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3))) +
    ggtitle(paste0("Harmony dims 1:", d))

  p_d <- dots_markers(act)

  ggsave(
    filename = file.path(out_dim_dir, sprintf("harmony_dims_%02d.png", d)),
    plot     = p_u / p_d,
    width    = 30, height = 12, units = "in",
    dpi = 300, device = ragg::agg_png, bg = "white"
  )

  DefaultAssay(act) <- "bulkOK"
}

# Compilar PNGs en un PDF (requiere magick)
png_files <- sort(list.files(out_dim_dir, pattern = "\\.png$", full.names = TRUE))
imgs <- magick::image_read(png_files)
magick::image_write(
  imgs,
  path   = file.path(fig_dir, paste0("harmony_dims_10_50_step2_", project, ".pdf")),
  format = "pdf"
)
message("PDF de dims guardado en: ", fig_dir)

# 7. (Opcional) Sensibilidad a la semilla de UMAP ----------------------------
#
# Descomenta para verificar estabilidad del UMAP con distintas semillas.
# Usa el objeto `act` con DIMS_FINAL ya definido.

# act <- FindNeighbors(act, reduction = "harmony", dims = DIMS_FINAL,
#                      graph.name = c("harmony_nn", "harmony_snn"))
# act <- FindClusters(act, graph.name = "harmony_snn", resolution = RES_FINAL,
#                     cluster.name = "harmony_clusters", algorithm = 4, random.seed = SEED)
#
# seed_plots <- lapply(1:5, function(s) {
#   act <- RunUMAP(act, reduction = "harmony", dims = DIMS_FINAL,
#                  seed.use = s, reduction.name = paste0("umap32_s", s),
#                  umap.method = "umap-learn", metric = "correlation")
#   DimPlot(act, reduction = paste0("umap32_s", s), label = TRUE,
#           alpha = 0.3, pt.size = 0.3) + ggtitle(paste0("seed = ", s))
# })
# save_plot_px(wrap_plots(seed_plots, ncol = 3),
#   file.path(fig_dir, paste0("umap_seeds_", project, ".png")),
#   w = 9000, h = 6000)

# 8. (Opcional) Diagnóstico de artefactos de UMAP ----------------------------
#
# Inspeccionnar vecinos de clústeres sospechosos para detectar artefactos.

# nn_name <- names(act@neighbors)[1]
# nn_idx  <- act@neighbors[[nn_name]]@nn.idx
#
# cluster_sospechoso <- "19"   # ← cambiar al clúster de interés
# c_cells <- WhichCells(act, idents = cluster_sospechoso)
# i_cells <- match(c_cells, colnames(act))
# i_cells <- i_cells[!is.na(i_cells)]
#
# nb_idx   <- nn_idx[i_cells, , drop = FALSE]
# nb_cells <- colnames(act)[as.vector(nb_idx)]
# message("Distribución de vecinos del clúster ", cluster_sospechoso, ":")
# print(prop.table(sort(table(Idents(act)[nb_cells]), decreasing = TRUE)))

# 9. Clustree con DIMS_FINAL --------------------------------------------------

DefaultAssay(act) <- "bulkOK"

act <- FindNeighbors(
  act, reduction = "harmony", dims = DIMS_FINAL,
  graph.name = c("harmony_nn", "harmony_snn")
)

res_grid <- c(1e-6, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

p_clustree <- clustree_plot(
  act,
  res        = res_grid,
  prefix     = "harmony_clusters.",
  graph.name = "harmony_snn"
)
save_plot_px(p_clustree,
  file.path(fig_dir, paste0("clustree_", project, ".png")),
  w = 4000, h = 5000
)

# 10. UMAP y clustering definitivo --------------------------------------------

message("Usando dims ", paste(range(DIMS_FINAL), collapse = ":"),
        " | resolución ", RES_FINAL)

act <- FindClusters(
  act, graph.name = "harmony_snn",
  resolution   = RES_FINAL,
  cluster.name = "harmony_clusters",
  algorithm    = 4,
  random.seed  = SEED
)
act <- RunUMAP(
  act, reduction = "harmony", dims = DIMS_FINAL,
  seed.use = SEED, reduction.name = "umap.harmony",
  umap.method = "umap-learn", metric = "correlation"
)

# Diagnóstico de integración
message("\nProporciones de muestra por clúster:")
print(prop.table(table(act$Sample, act$harmony_clusters, act$Treatment), 2))

message("\nCélulas por clúster:")
print(
  act@meta.data |>
    mutate(cluster = Idents(act)) |>
    group_by(cluster) |>
    summarise(n_cells = n(), .groups = "drop")
)

# 11. UMAP + DotPlot ----------------------------------------------------------

DefaultAssay(act) <- "RNA"
if (!"data" %in% Layers(act[["RNA"]])) act <- NormalizeData(act, assay = "RNA")

p_umap <- DimPlot(
  act, reduction = "umap.harmony",
  group.by = c("Treatment", "Sample", "harmony_clusters"),
  ncol = 3, label = TRUE, alpha = 0.3, pt.size = 0.3, shuffle = TRUE
) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3))) +
  theme(legend.position = "right")

p_markers <- dots_markers(act)

save_plot_px(p_umap,
  file.path(fig_dir, paste0("umap_integrated_", project,
                             "_", max(DIMS_FINAL), "pcs.png")),
  w = 6500, h = 2000
)
save_plot_px(p_markers,
  file.path(fig_dir, paste0("dotplot_", project,
                             "_", max(DIMS_FINAL), "pcs.png")),
  w = 4000, h = 3000
)
ggsave(
  file.path(fig_dir, paste0("umap_3cols_", project, ".png")),
  p_umap / p_markers,
  width = 30, height = 12, units = "in",
  dpi = 300, device = ragg::agg_png, bg = "white"
)

# 12. Guardar -----------------------------------------------------------------

qs_save(act, file.path(qs2_dir, paste0(project, "_fixed_umap.qs2")))
qs_save(act, "act_integrated.qs2")   # alias estándar para el pipeline
message("Objeto guardado en: act_integrated.qs2")
info(act)

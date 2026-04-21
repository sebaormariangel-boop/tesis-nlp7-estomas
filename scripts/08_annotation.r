# 08_annotation.R
# Anotación manual de clústeres, enriquecimiento GO y UMAP final publicable.
#
# Entrada:  act_integrated.qs2  (o act_digest_scored.qs2)
# Salida:   *_final_annotated.qs2
#           UMAP_celltypes_*.png / .pdf
#           GO_BP_dotplot.pdf
#           best_cluster_per_signature_<project>.csv
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

# 1. Cargar objeto -------------------------------------------------------------

in_file <- if (file.exists("act_digest_scored.qs2")) {
  "act_digest_scored.qs2"
} else {
  "act_integrated.qs2"
}
act <- qs_read(in_file)

DefaultAssay(act) <- "RNA"
if (!"data" %in% Layers(act[["RNA"]])) act <- NormalizeData(act, assay = "RNA")
info(act)

# 2. Marcadores diferenciales por clúster ------------------------------------
#
#  Necesarios para asignar identidades y para GO enrichment.
#  presto::wilcoxauc es mucho más rápido que FindAllMarkers con test.use="wilcox".

all_markers <- presto::wilcoxauc(act, group_by = "harmony_clusters", assay = "data")
all_markers <- all_markers |>
  dplyr::filter(logFC > 0.25, pct_in > 10, padj < 0.05) |>
  dplyr::arrange(group, desc(auc))

# Top 20 por clúster (diagnóstico)
top20 <- all_markers |>
  group_by(group) |>
  slice_head(n = 20) |>
  ungroup()

message("Marcadores diferenciales calculados: ", nrow(all_markers), " filas")

# Exportar tabla completa
write.csv(all_markers,
  file      = paste0("findallmarkers_", project, ".csv"),
  row.names = FALSE
)

# 3. Enriquecimiento GO (Biological Process) por clúster ---------------------

gene_list_by_cluster <- split(all_markers$feature, all_markers$group)

# Convertir AGI a Entrez
gene_list_entrez <- lapply(gene_list_by_cluster, function(g) {
  tryCatch(
    clusterProfiler::bitr(
      g, fromType = "TAIR", toType = "ENTREZID",
      OrgDb = org.At.tair.db
    )$ENTREZID,
    error = function(e) character(0)
  )
})

# Eliminar clústeres sin conversión
gene_list_entrez <- Filter(function(x) length(x) > 0, gene_list_entrez)

go_BP <- clusterProfiler::compareCluster(
  geneCluster = gene_list_entrez,
  fun         = "enrichGO",
  OrgDb       = org.At.tair.db,
  keyType     = "ENTREZID",
  ont         = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

p_go <- enrichplot::dotplot(go_BP, showCategory = 5) +
  ggtitle("GO - Biological Process")

ggsave(
  file.path(fig_dir, paste0("GO_BP_dotplot_", project, ".pdf")),
  p_go, width = 14, height = 30
)

message("GO enrichment guardado.")

# 4. Mapa de anotación --------------------------------------------------------
#
#  ⚠  Ajusta los números de clúster a los que aparezcan en tu UMAP.
#     Inspecciona: top20, best_cluster_per_signature_*.csv y dotplots.

if (tech_mode == "sc") {

  ct_map <- c(
    `1`  = "Mesophyll",
    `2`  = "Mesophyll",
    `3`  = "Mesophyll",
    `4`  = "Mesophyll",
    `5`  = "S phase",
    `6`  = "Mesophyll",
    `7`  = "Bundle Sheath",
    `8`  = "Mesophyll",
    `9`  = "Abaxial Epidermis",
    `10` = "Cambium",
    `11` = "Mitotic/Cytokinesis (G2/M)",
    `12` = "Cycling (S-phase licensing)",
    `13` = "Adaxial Epidermis",
    `14` = "Bundle Sheath",
    `15` = "Guard Cells",
    `16` = "Companion Cells",
    `17` = "Hydathode",
    `18` = "Myrosin Idioblasts",
    `19` = "Mesophyll"
  )

  ct_order <- c(
    "Mesophyll", "Bundle Sheath", "Cambium",
    "Companion Cells", "Guard Cells",
    "Adaxial Epidermis", "Abaxial Epidermis",
    "Hydathode", "Myrosin Idioblasts",
    "S phase", "Cycling (S-phase licensing)", "Mitotic/Cytokinesis (G2/M)"
  )

  ct_cols <- c(
    "Mesophyll"                      = "#00A84F",
    "Bundle Sheath"                  = "#F0E442",
    "Cambium"                        = "#8B4513",
    "Companion Cells"                = "#FF8C00",
    "Guard Cells"                    = "#6F2DBD",
    "Adaxial Epidermis"              = "#37474F",
    "Abaxial Epidermis"              = "#00BFC4",
    "Hydathode"                      = "#1F77B4",
    "Myrosin Idioblasts"             = "#FF4DA6",
    "S phase"                        = "#FFA500",
    "Cycling (S-phase licensing)"    = "#FFDD00",
    "Mitotic/Cytokinesis (G2/M)"     = "#E31A1C"
  )

} else if (tech_mode == "sn") {

  ct_map <- c(
    `1`  = "Mesophyll",
    `2`  = "Mesophyll",
    `3`  = "Epidermis",
    `4`  = "Mesophyll",
    `5`  = "Mesophyll",
    `6`  = "Mesophyll",
    `7`  = "Mesophyll",
    `8`  = "Mesophyll",
    `9`  = "Mesophyll",
    `10` = "Guard Cells",
    `11` = "Vasculature",
    `12` = "Epidermis",
    `13` = "Cambium",
    `14` = "Cycling",
    `15` = "Phloem",
    `16` = "Unknown 1",
    `17` = "Trichomes",
    `18` = "Guard Cells",
    `19` = "Unknown 2",
    `20` = "Myrosin Idioblasts",
    `21` = "Unknown 3"
  )

  ct_order <- c(
    "Mesophyll", "Epidermis", "Trichomes",
    "Guard Cells", "Vasculature", "Phloem",
    "Cambium", "Myrosin Idioblasts",
    "Cycling", "Unknown 1", "Unknown 2", "Unknown 3"
  )

  ct_cols <- c(
    "Mesophyll"      = "#2FBF71",
    "Epidermis"      = "#4E79A7",
    "Trichomes"      = "#17BEBB",
    "Guard Cells"    = "#7B2CBF",
    "Vasculature"    = "#1D4ED8",
    "Phloem"         = "#FF8C00",
    "Cambium"        = "#8B4513",
    "Myrosin Idioblasts" = "#D81B60",
    "Cycling"        = "#F2C300",
    "Unknown 1"      = "#374151",
    "Unknown 2"      = "#9CA3AF",
    "Unknown 3"      = "#6B7280"
  )

} else {
  stop("tech_mode no reconocido: ", tech_mode)
}

# Filtrar a los tipos que realmente existen en el objeto
ct_map_use <- ct_map[as.character(ct_map) %in% ct_order]

# 5. Aplicar anotación --------------------------------------------------------

act$celltype <- factor(
  unname(ct_map[as.character(act$harmony_clusters)]),
  levels = ct_order
)
Idents(act) <- act$celltype

message("\nDistribución por tipo celular:")
print(table(act$celltype))

# 6. UMAP final ---------------------------------------------------------------

## 6.1 Estándar (con leyenda de colores)
p_umap_leg <- DimPlot(
  act, reduction = "umap.harmony",
  group.by = "celltype", cols = ct_cols,
  pt.size = 0.01
) + umap_theme() + ggtitle("")

## 6.2 Con etiquetas repelidas
p_umap_lab <- LabelClusters(
  plot = DimPlot(
    act, reduction = "umap.harmony",
    group.by = "celltype", cols = ct_cols,
    pt.size = 0.12, alpha = 0.85, raster = FALSE
  ),
  id = "celltype", repel = TRUE, size = 3.5
) + umap_theme() + theme(legend.position = "none")

## 6.3 Fancy UMAP numerado (requiere helpers.R:fancy_umap)
ct_cols_palette <- RColorBrewer::brewer.pal(
  min(12, length(ct_order)), "Paired"
)
# Si hay más tipos que colores en la paleta, usa rainbow
if (length(ct_order) > 12) {
  ct_cols_palette <- rainbow(length(ct_order))
}
ct_cols_named <- setNames(
  ct_cols_palette[seq_along(ct_order)],
  ct_order
)

res_fancy <- fancy_umap(
  act,
  reduction           = "umap.harmony",
  ct_cols             = ct_cols,          # usa colores manuales definidos arriba
  circle.size         = 6.0,
  num.size            = 2.5,
  legend_title        = "Cell Type",
  legend_title_size   = 4.0,
  legend_row_spacing  = 0.05,
  legend_title_pad    = 0.04,
  legend_dot_size     = 6.0,
  legend_num_circle_size = 6.0,
  legend_num_text_size   = 2.5,
  legend_label_size   = 3.5,
  legend_x_dot        = 0.00,
  legend_x_num        = 0.15,
  legend_x_label      = 0.25,
  legend_xlim_right   = 2.0,
  legend_width        = 0.4
)

## 6.4 DotPlot con marcadores canónicos
p_dot <- dots_markers(act, group.by = "celltype")

# 7. Guardar figuras ----------------------------------------------------------

ggsave(
  file.path(fig_dir, "UMAP_celltypes_legend.pdf"),
  p_umap_leg, width = 7, height = 5, useDingbats = FALSE
)
save_plot_px(p_umap_leg,
  file.path(fig_dir, "UMAP_celltypes_legend.png"),
  w = 3500, h = 2500
)
save_plot_px(p_umap_lab,
  file.path(fig_dir, "UMAP_celltypes_labeled.png"),
  w = 3500, h = 2500
)
save_plot_px(res_fancy$plot,
  file.path(fig_dir, paste0("umap_fancy_", project, ".png")),
  w = 4000
)
save_plot_px(p_dot,
  file.path(fig_dir, paste0("dotplot_final_", project, ".png")),
  w = 4000, h = 3000
)
ggsave(
  file.path(fig_dir, paste0("umap_final_panel_", project, ".png")),
  p_umap_leg / p_dot,
  width = 20, height = 12, units = "in",
  dpi = 300, device = ragg::agg_png, bg = "white"
)

# 8. Guardar objeto final -----------------------------------------------------

out_name <- paste0(
  if (tech_mode == "sc") "tenorio" else "natanella",
  "_final_annotated"
)

qs_save(act, file.path(qs2_dir, paste0(out_name, ".qs2")))
saveRDS(act, paste0(out_name, ".rds"))

message("Objeto final guardado: ", out_name)
info(act)

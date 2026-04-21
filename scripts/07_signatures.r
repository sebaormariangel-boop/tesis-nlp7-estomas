# 07_signatures.R
# Puntajes UCell con lista combinada de marcadores (cell_cycle, epidermis,
# gc, vasculatura, mesófilo) y genes de ciclo celular.
#
# Entrada:  act_digest_scored.qs2  (o act_integrated.qs2)
# Salida:   best_cluster_per_signature_<project>.csv
#           UCell_featureplots_<project>.pdf
#           dotplot_<project>_ucell.png
#           cellcycle_*.png
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

N_CORES <- min(56L, parallel::detectCores(logical = FALSE))

# 1. Cargar objeto -------------------------------------------------------------

in_file <- if (file.exists("act_digest_scored.qs2")) {
  "act_digest_scored.qs2"
} else {
  "act_integrated.qs2"
}
act <- qs_read(in_file)
act <- JoinLayers(act)

DefaultAssay(act) <- "RNA"
if (!"data" %in% Layers(act[["RNA"]])) act <- NormalizeData(act, assay = "RNA")
info(act)

# 2. Lista maestra de marcadores (master list externa, opcional) --------------

if (file.exists(file.path(raw_path, "marker_master_list.txt"))) {
  marker_tbl <- readr::read_delim(
    file.path(raw_path, "marker_master_list.txt"),
    show_col_types = FALSE
  ) |> dplyr::select(-Name)

  pat_excl <- paste(
    c("sepal","petal","flower","stamen","carpel","ovule",
      "pollen","seed coat","stem vasculature","sam ","meristem"),
    collapse = "|"
  )

  marker_leaf <- marker_tbl |>
    dplyr::filter(is.na(Organ) | Organ == "Leaf") |>
    dplyr::select(-Organ) |>
    dplyr::filter(!grepl(pat_excl, `Cell type`, ignore.case = TRUE))

  marker_list_ext <- marker_leaf |>
    group_by(`Cell type`) |>
    summarise(m = list(unique(`Gene ID`)), .groups = "drop") |>
    tibble::deframe()
  names(marker_list_ext) <- paste0(names(marker_list_ext), "_B")
} else {
  marker_list_ext <- list()
  message("marker_master_list.txt no encontrado — se usa solo la lista interna.")
}

# 3. Firmas de células de guarda (archivos externos, opcional) ----------------

gc_files <- list.files(
  file.path(raw_path, "gc_signatures"),
  pattern    = "*(\\D*)gcgenes\\.txt$",
  full.names = FALSE
)
gc_signatures <- setNames(
  lapply(gc_files, function(f)
    read.table(file.path(raw_path, "gc_signatures", f))$V1),
  stringr::str_match(gc_files, "^(.*?)_gcgenes")[, 2]
)

# 4. Lista combinada principal de marcadores ----------------------------------
#
#  Combina todas las listas definidas en 00_setup.R.
#  `markers_ucell` = objeto enviado a ScoreSignatures_UCell.

markers_ucell <- c(
  cell_cycle,
  epidermis_markers,
  gc_markers,
  vasc_markers,
  mes_markers,
  marker_list_ext,
  gc_signatures
)

message("Total de firmas UCell: ", length(markers_ucell))

# 5. Marcadores extendidos de mesófilo (visualización exploratoria) -----------

mesophyll_markers <- sort(unique(c(
  # Rubisco / Calvin
  "AT1G67090", "AT5G38430", "AT5G38420", "AT5G38410",
  "AT3G01500",
  # Traducción cloroplástica
  "AT3G63140", "AT3G52150",
  # Antenas LHC (LHCA)
  "AT3G54890", "AT3G61470", "AT1G61520", "AT3G47470", "AT1G45474", "AT1G19150",
  # Antenas LHC (LHCB)
  "AT1G29920", "AT1G29910", "AT1G29930", "AT2G34430", "AT2G34420",
  "AT2G05100", "AT2G05070", "AT3G27690",
  "AT5G54270", "AT5G01530", "AT3G08940", "AT2G40100", "AT4G10340", "AT1G15820",
  # Otros LHC / NPQ
  "AT1G44575", "AT3G22840", "AT4G14690", "AT4G34190", "AT4G17600", "AT1G76570", "AT5G28450",
  # Mesófilo (señal / capas)
  "AT4G12970", "AT2G45190", "AT3G17185", "AT4G23060", "AT4G23600",
  # Otros marcadores de atlas scRNA-seq
  "AT1G12090", "AT5G20630"
)))

# DotPlot exploratorio de marcadores de mesófilo
p_mes <- dots_markers(act, features = mesophyll_markers)
save_plot_px(p_mes,
  file.path(fig_dir, paste0("dotplot_mesophyll_", project, ".png")),
  w = 8000, h = 4000, dpi = 300
)

# 6. Calcular rankings UCell --------------------------------------------------

m      <- median(act$nFeature_RNA, na.rm = TRUE)
q75    <- as.numeric(quantile(act$nFeature_RNA, 0.75, na.rm = TRUE))
maxRank <- round(min(max(m, 1200), q75, 3000), -2)
maxRank <- max(maxRank, 1200L)
message("maxRank UCell: ", maxRank)

ranks2 <- UCell::StoreRankings_UCell(
  LayerData(act, assay = "RNA", layer = "data"),
  maxRank = maxRank,
  ncores  = N_CORES
)

# 7. Puntajes UCell -----------------------------------------------------------

s2 <- UCell::ScoreSignatures_UCell(
  precalc.ranks = ranks2,
  features      = markers_ucell,
  ncores        = N_CORES
)

message("Firmas puntuadas: ", ncol(s2))

# Añadir a metadata
cells_common <- intersect(colnames(act), rownames(s2))
seu_gc       <- subset(act, cells = cells_common)
s2_aligned   <- s2[cells_common, , drop = FALSE]
stopifnot(identical(rownames(s2_aligned), colnames(seu_gc)))
seu_gc <- Seurat::AddMetaData(seu_gc, metadata = as.data.frame(s2_aligned))

ucell_cols <- grep("_UCell$", colnames(seu_gc@meta.data), value = TRUE)
message("Columnas UCell añadidas: ", length(ucell_cols))

# 8. Feature plots en PDF -----------------------------------------------------

pdf(
  file.path(fig_dir, paste0("UCell_featureplots_", project, ".pdf")),
  width = 14, height = 10
)
chunk_size <- 12
for (i in seq(1, length(ucell_cols), by = chunk_size)) {
  feats <- ucell_cols[i:min(i + chunk_size - 1, length(ucell_cols))]
  print(
    FeaturePlot(seu_gc, reduction = "umap.harmony",
                features = feats, ncol = 3,
                min.cutoff = "q1", max.cutoff = "q99")
  )
}
dev.off()

# DotPlot general de todas las firmas
p_dp <- dots_markers(act, features = markers_ucell)
save_plot_px(p_dp,
  file.path(fig_dir, paste0("dotplot_", project, "_ucell.png")),
  w = 12000, h = 4000, dpi = 300
)

# 9. Promedio de puntaje UCell por clúster ------------------------------------

cluster_col <- "harmony_clusters"

ord_cluster <- function(x) {
  xn <- suppressWarnings(as.numeric(as.character(x)))
  ifelse(is.na(xn), Inf, xn)
}

ucluster_mean <- seu_gc@meta.data |>
  dplyr::select(all_of(cluster_col), all_of(ucell_cols)) |>
  group_by(.data[[cluster_col]]) |>
  summarise(across(all_of(ucell_cols), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

# Mejor firma por clúster
best_sig <- ucluster_mean |>
  pivot_longer(-all_of(cluster_col), names_to = "signature", values_to = "mean_score") |>
  group_by(.data[[cluster_col]]) |>
  slice_max(mean_score, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(ord_cluster(.data[[cluster_col]]))

message("\nMejor firma por clúster:")
print(n = Inf, best_sig)

# Top 5 firmas por clúster
top5_sig <- ucluster_mean |>
  pivot_longer(-all_of(cluster_col), names_to = "signature", values_to = "mean_score") |>
  group_by(.data[[cluster_col]]) |>
  slice_max(mean_score, n = 5, with_ties = FALSE) |>
  ungroup() |>
  arrange(ord_cluster(.data[[cluster_col]]), desc(mean_score))

message("\nTop 5 firmas por clúster:")
print(n = Inf, top5_sig)

write.csv(
  best_sig,
  file      = paste0("best_cluster_per_signature_", project, ".csv"),
  row.names = FALSE
)

# 10. Ciclo celular -----------------------------------------------------------

DefaultAssay(act) <- "RNA"

cc_agi <- names(unlist(cell_cycle))   # todos los genes de las sublistas
cc_agi <- unique(cc_agi[cc_agi %in% rownames(act)])

message("Genes de ciclo celular presentes: ", length(cc_agi), " / ",
        length(unique(names(unlist(cell_cycle)))))

p_cc_feat <- FeaturePlot(act, reduction = "umap.harmony",
                          features = cc_agi, pt.size = 0.1)
save_plot_px(p_cc_feat,
  file.path(fig_dir, paste0("cellcycle_featureplot_", project, ".png")),
  w = 3200, h = 3200
)

p_cc_dot <- DotPlot(act, features = cc_agi) + RotatedAxis()
save_plot_px(p_cc_dot,
  file.path(fig_dir, paste0("cellcycle_dotplot_", project, ".png")),
  w = 2000, h = 1500
)

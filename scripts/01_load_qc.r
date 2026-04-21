# 01_load_qc.R
# Lectura de matrices CellRanger, filtrado QC, visualización y merge.
# Salida: nat_prefilt.qs2  (o ten_prefilt.qs2 según tech_mode)
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

# 1. Leer matrices y aplicar QC -----------------------------------------------

options(cli.progress_show_after = 0)

all_seu   <- vector("list", length(samples))
names(all_seu) <- samples

top_dirs <- list.dirs(sample_path, full.names = TRUE, recursive = FALSE)

cli::cli_progress_bar(
  "Leyendo matrices y aplicando QC Seurat",
  total = length(samples),
  clear = TRUE
)
cli_progress_update(set = 0, force = TRUE)

for (k in seq_along(samples)) {
  s       <- samples[k]
  mtx_dir <- find_10x_dir(s, top_dirs)
  mtx     <- Read10X(mtx_dir, gene.column = 1)

  seu <- CreateSeuratObject(
    counts       = mtx,
    project      = s,
    min.cells    = 0,
    min.features = 0
  )
  seu$Sample <- s

  cli::cli_progress_output("Células iniciales en {s}: {bold(red(ncol(seu)))}")

  # Metadatos de réplica y tratamiento
  if (s %in% names(map_repeat)) seu$Repeat    <- factor(unname(map_repeat[s]), levels = unique(map_repeat))
  if (s %in% names(map_treat))  seu$Treatment <- factor(unname(map_treat[s]),  levels = unique(map_treat))

  # Porcentaje mitocondrial y cloroplástico
  seu[["percent.mt"]]  <- PercentageFeatureSet(seu, pattern = "^ATMG")
  seu[["percent.chl"]] <- PercentageFeatureSet(seu, pattern = "^ATCG")

  # Filtro QC
  seu <- subset(
    seu,
    subset =
      nCount_RNA   >= qc_cfg$min_UMI      &
      nFeature_RNA >= qc_cfg$min_features  &
      percent.mt   <= qc_cfg$max_mt        &
      percent.chl  <= qc_cfg$max_chl
  )

  cli::cli_progress_output("Células tras filtrado en {s}: {bold(red(ncol(seu)))}")

  all_seu[[s]] <- seu
  cli_progress_update()
}

cli_progress_done(result = "Muestras cargadas y filtradas.")
invisible(lapply(all_seu, info))

# 2. Filtrar genes con expresión en < 5 células --------------------------------

all_seu <- lapply(all_seu, function(x) {
  DefaultAssay(x) <- "RNA"
  mat        <- GetAssayData(x, assay = "RNA", layer = "counts")
  keep_genes <- Matrix::rowSums(mat > 0) >= 5
  x          <- subset(x, features = rownames(mat)[keep_genes])

  message(
    "Se eliminaron ", sum(!keep_genes),
    " genes expresados en < 5 células en ", unique(x$Sample)
  )
  x
})

invisible(lapply(all_seu, info))

# 3. Visualización QC por muestra ----------------------------------------------

for (k in seq_along(all_seu)) {
  s <- samples[[k]]

  # Violin plots de métricas QC
  p_vln <- VlnPlot(
    object   = all_seu[[k]],
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"),
    ncol     = 4,
    pt.size  = 0
  )
  save_plot_px(p_vln,
    file   = paste0("vlnplot_", s, "_", project, ".png"),
    width  = 4500, height = 2000
  )

  # Scatter de correlaciones
  p_scatter <- (
    FeatureScatter(all_seu[[k]], "nCount_RNA",   "nFeature_RNA") +
    FeatureScatter(all_seu[[k]], "nCount_RNA",   "percent.chl")  +
    FeatureScatter(all_seu[[k]], "nFeature_RNA", "percent.chl")  +
    FeatureScatter(all_seu[[k]], "nCount_RNA",   "percent.mt")   +
    FeatureScatter(all_seu[[k]], "nFeature_RNA", "percent.mt")
  )
  save_plot_px(p_scatter,
    file   = paste0(s, "_", project, "_feature_scatter.png"),
    width  = 6500, height = 2000
  )
}

# 4. Merge y guardado ----------------------------------------------------------

act <- merge(
  all_seu[[1]],
  y            = all_seu[-1],
  add.cell.ids = samples,
  project      = project
)

out_file <- paste0(tolower(sub("Leaf_", "", project)), "_prefilt.qs2")
qs_save(act, out_file)
message("Objeto guardado en: ", out_file)

info(act)

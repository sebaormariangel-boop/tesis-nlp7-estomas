# 01_preprocessing.R --------------------------------------------------------
# Contenidos:
#  - Lectura de muestras
#  - Filtro por número de genes, de moléculas de RNA, %mitocondrial, %plastidial
#  - Remoción de dobletes

# 1 | Selección de dataset --------------------------------------------------

dataset_options <- c(
  "Single-cell | Tenorio Berrío et al.",
  "Single-nucleus | Illouz-Eliaz et al."
)

choice <- menu(dataset_options, title = "Seleccionar tecnología")

if (choice == 0) {
  stop("Selección cancelada por el usuario.", call. = FALSE)
}

tech_mode <- c("sc", "sn")[choice]

message(crayon::bold(crayon::red("==== ")),
        crayon::bold(crayon::red(dataset_options[choice])),
        crayon::bold(crayon::red(" ====")))

get_dataset_config <- function(tech_mode) {
  switch(
    tech_mode,
    sc = list(
      samples = c("MD_F_rep1", "MD_F_rep2", "WW_F_rep1", "WW_F_rep2"),
      qc = list(
        min_UMI = 1250L,
        min_features = 1000L,
        max_mt = 5,
        max_chl = 40
      ),
      sample_path = file.path(RAW_DIR, "tenorio"),
      results_path = file.path(RESULTS_DIR, "tenorio"),
      project = "Leaf_MD_WW_TB",
      repeat_map = c(
        MD_F_rep1 = "R1",
        MD_F_rep2 = "R2",
        WW_F_rep1 = "R1",
        WW_F_rep2 = "R2"
      ),
      treat_map = c(
        MD_F_rep1 = "MD",
        MD_F_rep2 = "MD",
        WW_F_rep1 = "WW",
        WW_F_rep2 = "WW"
      )
    ),
    sn = list(
      samples = c("DR_T0_rep1", "DR_T0_rep2", "WW_T0_rep1", "WW_T0_rep2"),
      qc = list(
        min_UMI = 300L,
        min_features = 300L,
        max_mt = 1,
        max_chl = 40
      ),
      sample_path = file.path(RAW_DIR, "natanella"),
      results_path = file.path(RESULTS_DIR, "natanella"),
      project = "Leaf_DR_WW_IE",
      repeat_map = c(
        DR_T0_rep1 = "R1",
        DR_T0_rep2 = "R2",
        WW_T0_rep1 = "R1",
        WW_T0_rep2 = "R2"
      ),
      treat_map = c(
        DR_T0_rep1 = "DR",
        DR_T0_rep2 = "DR",
        WW_T0_rep1 = "WW",
        WW_T0_rep2 = "WW"
      )
    ),
    stop("tech_mode debe ser 'sc' o 'sn'", call. = FALSE)
  )
}

cfg <- get_dataset_config(tech_mode)

samples <- cfg$samples
qc_cfg <- cfg$qc
sample_path <- cfg$sample_path
project <- cfg$project
map_repeat <- cfg$repeat_map
map_treat <- cfg$treat_map

setwd(cfg$results_path)

message(crayon::green("Muestras: "), paste(samples, collapse = ", "))
message(
  crayon::green("Umbrales QC: "),
  sprintf(
    "min_UMI = %d, min_features = %d, max_mt = %.1f%%, max_chl = %.1f%%",
    qc_cfg$min_UMI, qc_cfg$min_features, qc_cfg$max_mt, qc_cfg$max_chl
  )
)
message(crayon::green("Directorio: "), cfg$results_path)
message(crayon::green("Proyecto: "), project)

# 2 | Lectura de matrices y QC inicial --------------------------------------

read_and_filter_sample <- function(sample_name, top_dirs, qc_cfg, map_repeat, map_treat) {
  mtx_dir <- find_10x_dir(sample_name, top_dirs)
  mtx <- Seurat::Read10X(mtx_dir, gene.column = 1)

  seu <- Seurat::CreateSeuratObject(
    counts = mtx,
    project = sample_name,
    min.cells = 0,
    min.features = 0
  )

  seu$Sample <- sample_name
  seu$Repeat <- factor(unname(map_repeat[sample_name]), levels = unique(map_repeat))
  seu$Treatment <- factor(unname(map_treat[sample_name]), levels = unique(map_treat))

  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^ATMG")
  seu[["percent.chl"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^ATCG")

  seu <- subset(
    seu,
    subset =
      nCount_RNA >= qc_cfg$min_UMI &
      nFeature_RNA >= qc_cfg$min_features &
      percent.mt <= qc_cfg$max_mt &
      percent.chl <= qc_cfg$max_chl
  )

  seu
}

options(cli.progress_show_after = 0)

all_seu <- vector("list", length(samples))
names(all_seu) <- samples

top_dirs <- list.dirs(sample_path, full.names = TRUE, recursive = FALSE)

cli::cli_progress_bar(
  "Leyendo matrices y aplicando QC Seurat",
  total = length(samples),
  clear = TRUE
)

for (i in seq_along(samples)) {
  s <- samples[i]

  mtx_dir <- find_10x_dir(s, top_dirs)
  mtx <- Seurat::Read10X(mtx_dir, gene.column = 1)

  seu_raw <- Seurat::CreateSeuratObject(
    counts = mtx,
    project = s,
    min.cells = 0,
    min.features = 0
  )

  cli::cli_progress_output("Células iniciales en {s}: {bold(red(ncol(seu_raw)))}")

  seu <- read_and_filter_sample(
    sample_name = s,
    top_dirs = top_dirs,
    qc_cfg = qc_cfg,
    map_repeat = map_repeat,
    map_treat = map_treat
  )

  cli::cli_progress_output("Células tras el filtrado en {s}: {bold(red(ncol(seu)))}")

  all_seu[[s]] <- seu
  cli::cli_progress_update()
}

cli::cli_progress_done(result = "Muestras cargadas y filtradas.")
invisible(lapply(all_seu, info))

# 3 | Filtrar genes expresados en menos de 5 células ------------------------

filter_genes_by_detection <- function(seu, min_cells = 5) {
  Seurat::DefaultAssay(seu) <- "RNA"
  mat <- SeuratObject::GetAssayData(seu, assay = "RNA", layer = "counts")

  keep_genes <- Matrix::rowSums(mat > 0) >= min_cells
  removed_genes <- rownames(mat)[!keep_genes]

  seu <- subset(seu, features = rownames(mat)[keep_genes])

  message(
    "Se eliminaron ", length(removed_genes),
    " genes expresados en < ", min_cells,
    " células en ", unique(seu$Sample)
  )

  seu
}

all_seu <- lapply(all_seu, filter_genes_by_detection, min_cells = 5)
invisible(lapply(all_seu, info))

# 4 | Plots QC por muestra --------------------------------------------------

save_qc_plots <- function(seu, sample_name, project, suffix = "") {
  vln <- Seurat::VlnPlot(
    object = seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chl"),
    ncol = 4,
    pt.size = 0
  )

  ggplot2::ggsave(
    filename = paste0("vlnplot_", sample_name, "_", project, suffix, ".png"),
    plot = vln,
    width = 4500,
    height = 2000,
    units = "px",
    dpi = 300
  )

  plot1 <- Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.chl")
  plot2 <- Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- Seurat::FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.chl")
  plot4 <- Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot5 <- Seurat::FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.mt")

  qc_scatter <- plot2 + plot1 + plot3 + plot4 + plot5

  ggplot2::ggsave(
    filename = paste0(sample_name, "_", project, "_feature_scatter", suffix, ".png"),
    plot = qc_scatter,
    width = 6500,
    height = 2000,
    units = "px",
    dpi = 300
  )
}

plot_suffix <- if (exists("dobletes_filtrados") && isTRUE(dobletes_filtrados)) "_nodbls" else ""

for (i in seq_along(all_seu)) {
  save_qc_plots(
    seu = all_seu[[i]],
    sample_name = samples[i],
    project = project,
    suffix = plot_suffix
  )
}

# 5 | Merge global ----------------------------------------------------------

act <- merge(
  x = all_seu[[1]],
  y = all_seu[-1],
  add.cell.ids = samples,
  project = project
)

qs2::qs_save(act, paste0(project, "_prefilt.qs2"))
info(act)

# 6 | Detección y remoción de dobletes --------------------------------------

Seurat::DefaultAssay(act) <- "RNA"

sce <- SeuratObject::JoinLayers(act) |>
  SingleCellExperiment::as.SingleCellExperiment()

sce <- scds::cxds_bcds_hybrid(sce, estNdbl = TRUE, verb = TRUE)

act$hybrid_score <- SummarizedExperiment::colData(sce)$hybrid_score
act$hybrid_call <- SummarizedExperiment::colData(sce)$hybrid_call

doublets <- colnames(subset(act, subset = hybrid_call))
act <- subset(act, subset = !hybrid_call)

message("Dobletes eliminados: ", length(doublets))

saveRDS(doublets, "scds_dbls.rds")
dobletes_filtrados <- TRUE

print(table(SummarizedExperiment::colData(sce)$Sample,
            SummarizedExperiment::colData(sce)$hybrid_call))

info(act)

# 7 | Reaplicar filtros opcionales desde archivos ---------------------------

apply_optional_filters <- function(seu) {
  if (file.exists("scds_dbls.rds")) {
    doublets <- readRDS("scds_dbls.rds")
    seu <- subset(seu, cells = setdiff(colnames(seu), doublets))
    message("Dobletes removidos: ", length(doublets))
    dobletes_filtrados <<- TRUE
    info(seu)
  }

  if (file.exists("responsive_cells_ucell.rds")) {
    responsive_cells <- readRDS("responsive_cells_ucell.rds")
    seu <- subset(seu, cells = setdiff(colnames(seu), responsive_cells))
    message("Células responsivas a digestión removidas: ", length(responsive_cells))
    no_digestion <<- TRUE
    info(seu)
  }

  if (file.exists("final_cells.rds") && file.exists("final_genes.rds")) {
    final_cells <- readRDS("final_cells.rds")
    final_genes <- readRDS("final_genes.rds")

    seu <- subset(
      seu,
      cells = intersect(colnames(seu), final_cells),
      features = intersect(rownames(seu), final_genes)
    )

    no_digestion <<- TRUE
    info(seu)
  }

  seu
}

act <- apply_optional_filters(act)

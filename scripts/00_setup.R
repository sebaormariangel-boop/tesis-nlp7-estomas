# 00_setup.R ---------------------------------------------------------------
## Contenido:
##  - Variables globales
##  - Helpers
##  - Librerías
##  - Rutas

# 0 | Configuración global -------------------------------------------------

SEED <- 123
DIMS_USE <- 1:32
N_THREADS <- 56

ASSAY_USE <- "RNA"
LAYER_USE <- "data"   # Seurat v5
SAMPLE_COL <- "Sample"
TREAT_COL <- "Treatment"
REDUCTION_USE <- "harmony"

INTERACTIVE_MODE <- TRUE

# Barrido de fraction para hdWGCNA
FRAC_GRID <- c(0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15)
MIN_GENES <- 4000
MAX_GENES <- 15000

# Activar paquetes pesados de regulación/motifs solo cuando los necesites
LOAD_OPTIONAL_REGULON_PKGS <- FALSE

set.seed(SEED)

# 1 | Paquetes -------------------------------------------------------------

BASE_PKGS <- c(
  "Seurat", "SeuratObject", "hdWGCNA", "WGCNA",
  "Matrix", "dplyr", "tidyr", "tibble", "purrr",
  "ggplot2", "cowplot", "patchwork", "clustree",
  "stringr", "readr", "readxl", "qs2",
  "UCell", "scds", "SingleCellExperiment",
  "harmony", "edgeR",
  "clusterProfiler", "org.At.tair.db", "enrichplot",
  "genesectR", "presto", "tidytext",
  "digest", "reticulate",
  "cli", "crayon", "mclust", "ggrepel",
  "RColorBrewer", "AnnotationDbi"
)

OPTIONAL_REGULON_PKGS <- c(
  "TFBSTools", "motifmatchr", "GenomicRanges", "GenomeInfoDb",
  "ensembldb", "rtracklayer", "RSQLite", "JASPAR2024",
  "BSgenome.Athaliana.TAIR.TAIR9", "xgboost"
)

load_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Faltan paquetes por instalar:\n- ",
      paste(missing_pkgs, collapse = "\n- "),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages(
    invisible(lapply(
      pkgs,
      function(pkg) library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    ))
  )
}

load_packages(BASE_PKGS)

if (isTRUE(LOAD_OPTIONAL_REGULON_PKGS)) {
  load_packages(OPTIONAL_REGULON_PKGS)
}

ggplot2::theme_set(cowplot::theme_cowplot())

# 2 | Reticulate -----------------------------------------------------------

RETICULATE_PYTHON_PATH <- file.path(
  path.expand("~"),
  "miniconda3", "envs", "seurat", "bin", "python"
)

Sys.setenv(RETICULATE_PYTHON = RETICULATE_PYTHON_PATH)

PY_CONFIG <- reticulate::py_config()
UMAP_PY_AVAILABLE <- reticulate::py_module_available("umap")

message("RETICULATE_PYTHON = ", Sys.getenv("RETICULATE_PYTHON"))
message("Python activo: ", PY_CONFIG$python)
message("umap disponible en Python: ", UMAP_PY_AVAILABLE)

# hilos para WGCNA
WGCNA::enableWGCNAThreads(nThreads = N_THREADS)

# 3 | Rutas y archivos -----------------------------------------------------

HOME_DIR <- path.expand("~")
DOCS_DIR <- file.path(HOME_DIR, "Documentos")

# Ajusta este bloque si cambia tu estructura
PROJ_DIR <- file.path(DOCS_DIR, "projects", "sortiz", "2do_semestre", "singlecell")
RAW_DIR <- file.path(PROJ_DIR, "final_samples")
RESULTS_DIR <- file.path(PROJ_DIR, "results")

TABLES_DROUGHT_DIR <- file.path(DOCS_DIR, "tables_DROUGHT")
TABLES_STOMATA_DIR <- file.path(DOCS_DIR, "tables_STOMATA")

AT_SYMBOL_PATH <- file.path(HOME_DIR, "atsymbol_nate.txt")
PP_GENES_PATH <- file.path(RAW_DIR, "protoplasting_leaves.txt")
FIND_FILES_PATH <- file.path(HOME_DIR, "projects", "2do_semestre", "scripts", "find_files.R")

# Cargar archivo de símbolos
map_df <- if (file.exists(AT_SYMBOL_PATH)) {
  readr::read_tsv(AT_SYMBOL_PATH, col_types = readr::cols()) %>%
    dplyr::filter(!is.na(gene), !is.na(symbol)) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
} else {
  warning("No encontré el archivo de símbolos: ", AT_SYMBOL_PATH)
  tibble::tibble(gene = character(), symbol = character())
}

gene2sym <- setNames(map_df$symbol, map_df$gene)

# Genes de digestión / protoplasting
pp.genes <- if (file.exists(PP_GENES_PATH)) {
  scan(PP_GENES_PATH, what = "character", quiet = TRUE)
} else {
  warning("No encontré el archivo de genes de digestión: ", PP_GENES_PATH)
  character(0)
}

# Cargar find_files() si existe
if (file.exists(FIND_FILES_PATH)) {
  source(FIND_FILES_PATH)
}

# 4 | Marcadores y abreviaturas --------------------------------------------

markers <- list(
  Epidermis          = c("AT2G42840", "AT3G16370"),
  Vasculature        = c("AT1G64700", "AT5G03610", "AT1G80520", "AT3G14990"),
  Mesophyll          = c("AT4G26530", "AT2G34430", "AT4G12970"),
  Bundle_Sheath      = c("AT1G25230", "AT4G13770"),
  Cambium            = c("AT2G39700"),
  Companion_Cells    = c("AT4G19840", "AT5G18600", "AT1G64370"),
  Sieve_Element      = c("AT1G05760", "AT5G04890"),
  Phloem_Parenchyma  = c("AT3G11930", "AT5G24800"),
  Xylem              = c("AT3G10080", "AT5G60490"),
  Hydathode          = c("AT1G22900"),
  Myrosin_Idioblasts = c("AT5G26000", "AT5G25980", "AT3G16400"),
  Guard_Cells        = c("AT1G04800", "AT1G08810", "AT1G12480", "AT4G33950"),
  Pavement_Cells     = c("AT5G63180")
)

abbr <- c(
  Epidermis          = "E",
  Vasculature        = "V",
  Mesophyll          = "M",
  Bundle_Sheath      = "BS",
  Cambium            = "Cam",
  Companion_Cells    = "CC",
  Sieve_Element      = "SE",
  Phloem_Parenchyma  = "PP",
  Xylem              = "Xyl",
  Hydathode          = "Hyd",
  Myrosin_Idioblasts = "MI",
  Guard_Cells        = "GC",
  Pavement_Cells     = "Pav"
)

marker_names_old <- names(markers)
names(markers) <- ifelse(
  marker_names_old %in% names(abbr),
  abbr[marker_names_old],
  marker_names_old
)

cap_text <- stringr::str_wrap(
  paste(names(abbr), abbr, sep = " = ", collapse = "; "),
  width = 90
)

# 5 | Helpers ---------------------------------------------------------------

quiet <- function(expr) {
  invisible(
    suppressWarnings(
      suppressMessages(
        capture.output(eval.parent(substitute(expr)))
      )
    )
  )
}

style_text <- function(x, fg = NULL, bold = FALSE) {
  x <- as.character(x)

  if (!requireNamespace("crayon", quietly = TRUE)) {
    return(x)
  }

  if (!is.null(fg)) {
    f <- switch(
      fg,
      red = crayon::red,
      green = crayon::green,
      blue = crayon::blue,
      x
    )
    x <- f(x)
  }

  if (bold) {
    x <- crayon::bold(x)
  }

  x
}

info <- function(obj,
                 sample_col = SAMPLE_COL,
                 feature_col = "nFeature_RNA") {
  md <- obj[[]]

  medians <- "—"
  if (all(c(sample_col, feature_col) %in% names(md))) {
    med_tbl <- md |>
      dplyr::group_by(.data[[sample_col]]) |>
      dplyr::summarise(
        median_genes = median(.data[[feature_col]], na.rm = TRUE),
        .groups = "drop"
      )

    medians <- paste0(
      med_tbl[[sample_col]], ": ", med_tbl$median_genes
    ) |>
      paste(collapse = ", ")
  }

  vf <- length(Seurat::VariableFeatures(obj))
  genes_line <- paste0(
    nrow(obj),
    if (vf > 0) paste0(" (", vf, ")") else ""
  )

  lines <- c(
    style_text(Project(obj), fg = "red", bold = TRUE),
    paste0(style_text("Genes:", fg = "green"), " ", genes_line),
    paste0(style_text("Células:", fg = "green"), " ", ncol(obj)),
    paste0(style_text("Capas:", fg = "green"), " ", paste(SeuratObject::Layers(obj), collapse = ", ")),
    paste0(style_text("Ensayos:", fg = "green"), " ", paste(Seurat::Assays(obj), collapse = ", ")),
    paste0(style_text("Ensayo activo:", fg = "green"), " ", Seurat::DefaultAssay(obj)),
    paste0(style_text("Reducciones:", fg = "green"), " ", paste(Seurat::Reductions(obj), collapse = ", ")),
    paste0(style_text("Mediana de genes (por muestra):", fg = "green"), " ", medians)
  )

  message(paste(lines, collapse = "\n"))
  invisible(NULL)
}

dots_markers <- function(obj, group.by = "harmony_clusters") {
  Seurat::DotPlot(
    obj,
    features = markers,
    assay = "RNA",
    group.by = group.by
  ) +
    ggplot2::labs(caption = cap_text) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 10),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(size = 10),
      panel.spacing = grid::unit(1, "lines"),
      axis.title = ggplot2::element_blank(),
      plot.caption.position = "plot",
      plot.caption = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(t = 5, r = 12, b = 30, l = 12)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_color_gradient2(
      low = "navy",
      mid = "white",
      high = "firebrick3",
      midpoint = 0,
      limits = c(-2.5, 2.5)
    )
}

umap_theme <- function(base_size = 11) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
}

find_10x_dir <- function(sample_name, top_dirs) {
  hits <- top_dirs[grepl(sample_name, basename(top_dirs), fixed = TRUE)]

  if (length(hits) == 0) {
    stop("No hay carpeta para: ", sample_name, call. = FALSE)
  }
  if (length(hits) > 1) {
    stop(
      "Más de una carpeta para ", sample_name, ":\n",
      paste(basename(hits), collapse = "\n"),
      call. = FALSE
    )
  }

  base <- hits[1]
  candidates <- c(
    base,
    file.path(base, "filtered_feature_bc_matrix"),
    file.path(base, "raw_feature_bc_matrix")
  )

  ok <- candidates[
    file.exists(file.path(candidates, "matrix.mtx.gz")) |
      file.exists(file.path(candidates, "matrix.mtx"))
  ]

  if (length(ok) == 0) {
    stop(
      "Existe la carpeta base, pero no matrix.mtx(.gz) para: ",
      sample_name, "\nBase: ", base,
      call. = FALSE
    )
  }

  ok[1]
}

clustree_plot <- function(obj,
                          res,
                          prefix = "clusters.",
                          graph.name = "harmony_snn") {
  n <- length(res)
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  on.exit(close(pb), add = TRUE)

  message("Clustering: probando ", n, " resoluciones.")

  for (k in seq_along(res)) {
    r <- res[[k]]
    obj <- Seurat::FindClusters(
      obj,
      graph.name = graph.name,
      resolution = r,
      cluster.name = paste0(prefix, r),
      algorithm = 4,
      random.seed = SEED,
      verbose = FALSE
    )
    utils::setTxtProgressBar(pb, k)
  }

  message("Clustering terminado.")
  clustree::clustree(obj, prefix = prefix)
}

clust <- function(obj, res) {
  options(
    cli.dynamic = TRUE,
    cli.num_ansi_colors = 256,
    cli.progress_show_after = 0,
    cli.progress_clear = TRUE
  )

  fmt_eta <- function(start_time, current, total) {
    done <- max(current, 1L)
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- elapsed / done
    remaining <- max(0, rate * (total - done))
    mm <- floor(remaining / 60)
    ss <- round(remaining %% 60)
    sprintf("%02d:%02d", mm, ss)
  }

  old_opt <- options(cli.progress_show_after = 0)
  on.exit(options(old_opt), add = TRUE)

  n <- length(res)
  t0 <- Sys.time()

  id <- cli::cli_progress_bar(
    "Buscando clusters...",
    total = n,
    clear = FALSE,
    format = "{current}/{total} {cli_progress_bar} {percent}% | ETA {eta}",
    format_done = "{symbol$tick} Clustering terminado | {elapsed}"
  )

  for (k in seq_along(res)) {
    r <- res[k]

    quiet({
      obj <- Seurat::FindClusters(
        obj,
        resolution = r,
        cluster.name = paste0("harmony_clusters.", r),
        algorithm = 4,
        random.seed = SEED,
        verbose = FALSE
      )
    })

    cli::cli_progress_message(sprintf(
      "Paso %d/%d | res=%.3g | ETA %s",
      k, n, r, fmt_eta(t0, k, n)
    ))

    cli::cli_progress_update(id, inc = 1)
  }

  elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  mm <- floor(elapsed_sec / 60)
  ss <- round(elapsed_sec %% 60)

  final <- sprintf("%d/%d (100%%) | tiempo: %02d:%02d", n, n, mm, ss)
  cli::cli_alert_success(final)
  cli::cli_progress_done(id)

  clustree::clustree(obj, prefix = "harmony_clusters.")
}

hash_vec <- function(x, algo = "xxhash64") {
  x <- as.character(x)
  digest::digest(x, algo = algo, serialize = TRUE)
}

hash_list <- function(x, algo = "xxhash64") {
  digest::digest(x, algo = algo, serialize = TRUE)
}

fraction_sweep <- function(
  gc_obj,
  frac_grid,
  group_by,
  tmp_prefix = "tmp_frac_",
  assay = ASSAY_USE,
  strict = TRUE
) {
  if (length(frac_grid) < 2) {
    stop("fraction_sweep: frac_grid debe tener >= 2 valores.", call. = FALSE)
  }

  frac_grid <- sort(as.numeric(frac_grid))
  if (anyNA(frac_grid)) {
    stop("fraction_sweep: frac_grid tiene valores no numéricos.", call. = FALSE)
  }
  if (anyDuplicated(frac_grid)) {
    stop("fraction_sweep: frac_grid tiene fracciones duplicadas.", call. = FALSE)
  }

  out <- vector("list", length(frac_grid))
  pb <- utils::txtProgressBar(min = 0, max = length(frac_grid), style = 3)
  on.exit(close(pb), add = TRUE)

  prev_n <- NULL
  prev_hash <- NULL

  for (i in seq_along(frac_grid)) {
    f <- frac_grid[[i]]
    tmp_name <- paste0(tmp_prefix, gsub("\\.", "_", as.character(f)))

    tmp <- gc_obj

    tmp <- hdWGCNA::SetupForWGCNA(
      seurat_obj = tmp,
      assay = assay,
      gene_select = "fraction",
      fraction = f,
      group.by = group_by,
      wgcna_name = tmp_name
    )

    genes <- tmp@misc[[tmp_name]]$wgcna_genes
    if (is.null(genes) || length(genes) == 0) {
      stop("fraction_sweep: wgcna_genes es NULL/0 para fraction = ", f, call. = FALSE)
    }

    genes <- unique(as.character(genes))
    ghash <- hash_vec(sort(genes))
    n_genes <- length(genes)

    out[[i]] <- tibble::tibble(
      fraction = f,
      n_genes = n_genes,
      genes_hash = ghash
    )

    if (strict && !is.null(prev_n)) {
      if (n_genes > prev_n) {
        stop(
          "fraction_sweep: n_genes aumentó al subir fraction (",
          frac_grid[[i - 1]], " -> ", f, "): ", prev_n, " -> ", n_genes,
          call. = FALSE
        )
      }
      if (n_genes == prev_n) {
        stop(
          "fraction_sweep: plateau de n_genes (",
          frac_grid[[i - 1]], " y ", f, " dan ", n_genes, " genes).",
          call. = FALSE
        )
      }
      if (identical(ghash, prev_hash)) {
        stop(
          "fraction_sweep: set de genes idéntico entre fracciones consecutivas (",
          frac_grid[[i - 1]], " y ", f, ").",
          call. = FALSE
        )
      }
    }

    prev_n <- n_genes
    prev_hash <- ghash
    utils::setTxtProgressBar(pb, i)
  }

  df <- dplyr::bind_rows(out)

  if (strict) {
    if (anyDuplicated(df$n_genes)) {
      stop("fraction_sweep: hay n_genes repetidos.", call. = FALSE)
    }
    if (anyDuplicated(df$genes_hash)) {
      stop("fraction_sweep: hay genes_hash repetidos.", call. = FALSE)
    }
  }

  df
}

trt_map <- function(tech_mode) {
  switch(
    tech_mode,
    sc = list(ctrl = "WW", trt = "MD"),
    sn = list(ctrl = "WW", trt = "DR"),
    stop("tech_mode debe ser 'sc' o 'sn'")
  )
}

choose_one <- function(x, title = "Selecciona una opción:") {
  x <- sort(unique(x))
  ch <- utils::menu(x, title = title)
  if (ch == 0) {
    stop("Selección cancelada.", call. = FALSE)
  }
  x[ch]
}

save_txt <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, con = path)
}

save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = path)
}

save_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, file = path)
}

save_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(df, file = path)
}

save_plot_px <- function(p,
                         filename,
                         w = 2500,
                         h = 2500,
                         dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = filename,
    plot = p,
    width = w,
    height = h,
    units = "px",
    dpi = dpi
  )
}

save_obj <- function(obj, path, overwrite_ext = TRUE) {
  save_options <- c("qs2", "rds")

  save_choice <- utils::menu(
    save_options,
    title = "Seleccionar formato de guardado:"
  )

  if (save_choice == 0) {
    stop("Selección cancelada por el usuario.", call. = FALSE)
  }

  ext <- save_options[[save_choice]]

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (overwrite_ext) {
    path_final <- sub("\\.[^.]*$", paste0(".", ext), path)
    if (!grepl("\\.[^.]*$", path)) {
      path_final <- paste0(path, ".", ext)
    }
  } else {
    path_final <- if (grepl("\\.[^.]*$", path)) path else paste0(path, ".", ext)
  }

  if (ext == "qs2") {
    qs2::qs_save(obj, file = path_final)
  } else if (ext == "rds") {
    saveRDS(obj, file = path_final)
  } else {
    stop("Formato inesperado.", call. = FALSE)
  }

  invisible(path_final)
}

to_sym <- function(x) {
  y <- unname(gene2sym[x])
  ifelse(is.na(y), x, y)
}

to_sym_agi <- function(x) {
  sym <- unname(gene2sym[x])
  ifelse(is.na(sym) | sym == "", x, paste0(sym, " (", x, ")"))
}

# 6 | Resumen al cargar -----------------------------------------------------

message("PROJ_DIR = ", PROJ_DIR)
message("RAW_DIR = ", RAW_DIR)
message("RESULTS_DIR = ", RESULTS_DIR)
message("Genes de digestión cargados: ", length(pp.genes))
message("Mapa gene -> symbol cargado: ", nrow(map_df), " entradas")

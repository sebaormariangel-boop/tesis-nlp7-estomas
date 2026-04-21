# helpers.R
# Funciones auxiliares reutilizables para el pipeline de scRNA/snRNA-seq.
# Requiere que 00_setup.R haya sido sourced previamente.
# ─────────────────────────────────────────────────────────────────────────────

# 1. Diagnóstico de objetos Seurat --------------------------------------------

#' Imprime un resumen formateado del objeto Seurat
#'
#' @param obj         Objeto Seurat
#' @param sample_col  Columna de metadatos con el nombre de muestra
#' @param feature_col Columna de metadatos con nFeature
info <- function(obj,
                 sample_col  = "Sample",
                 feature_col = "nFeature_RNA") {
  md <- obj[[]]

  style <- function(x, fg = NULL, bold = FALSE) {
    x <- as.character(x)
    if (!requireNamespace("crayon", quietly = TRUE)) return(x)
    if (!is.null(fg)) {
      f <- switch(fg, red = crayon::red, green = crayon::green, identity)
      x <- f(x)
    }
    if (bold) x <- crayon::bold(x)
    x
  }

  medians <- "—"
  if (all(c(sample_col, feature_col) %in% names(md))) {
    med_tbl <- md |>
      dplyr::group_by(.data[[sample_col]]) |>
      dplyr::summarise(
        median_genes = median(.data[[feature_col]], na.rm = TRUE),
        .groups = "drop"
      )
    medians <- paste0(med_tbl[[sample_col]], ": ", med_tbl$median_genes) |>
      paste(collapse = ", ")
  }

  vf <- length(VariableFeatures(obj))
  genes_line <- paste0(nrow(obj), if (vf > 0) paste0(" (", vf, ")") else "")

  lines <- c(
    style(Project(obj), fg = "red", bold = TRUE),
    paste0(style("Genes:",         fg = "green"), " ", genes_line),
    paste0(style("Células:",       fg = "green"), " ", ncol(obj)),
    paste0(style("Capas:",         fg = "green"), " ", paste(Layers(obj),         collapse = ", ")),
    paste0(style("Ensayos:",       fg = "green"), " ", paste(Seurat::Assays(obj), collapse = ", ")),
    paste0(style("Ensayo activo:", fg = "green"), " ", DefaultAssay(obj)),
    paste0(style("Reducciones:",   fg = "green"), " ", paste(Reductions(obj),     collapse = ", ")),
    paste0(style("Mediana genes (por muestra):", fg = "green"), " ", medians)
  )

  message(paste(lines, collapse = "\n"))
  invisible(NULL)
}

# 2. Clustering exploratorio (clustree) ----------------------------------------

#' Corre FindClusters sobre un vector de resoluciones y devuelve el árbol
#'
#' @param obj        Objeto Seurat (con FindNeighbors ya ejecutado)
#' @param res        Vector numérico de resoluciones
#' @param prefix     Prefijo para los nombres de columna en metadata
#' @param graph.name Nombre del grafo SNN a usar (ej. "harmony_snn")
#' @return Plot de clustree
clustree_plot <- function(obj,
                          res,
                          prefix     = "clusters.",
                          graph.name = "harmony_snn") {
  n  <- length(res)
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  on.exit(close(pb), add = TRUE)

  message("Clustering: probando ", n, " resoluciones...")

  for (k in seq_along(res)) {
    r   <- res[[k]]
    obj <- Seurat::FindClusters(
      obj,
      graph.name   = graph.name,
      resolution   = r,
      cluster.name = paste0(prefix, r),
      algorithm    = 4,    # Leiden
      random.seed  = SEED,
      verbose      = FALSE
    )
    utils::setTxtProgressBar(pb, k)
  }

  message("\nClustering terminado.")
  clustree::clustree(obj, prefix = prefix)
}

# 3. Visualización -------------------------------------------------------------

#' DotPlot estandarizado con etiquetas de símbolo en el eje X
#'
#' @param obj      Objeto Seurat
#' @param features Lista de marcadores (por defecto usa `markers` del entorno global)
#' @param group.by Columna de agrupación
#' @return ggplot
dots_markers <- function(obj,
                          features = NULL,
                          group.by = "harmony_clusters") {
  if (is.null(features)) {
    features <- get("markers", envir = parent.env(environment()))
    if (is.null(features)) features <- get("markers", envir = globalenv())
  }

  dp  <- DotPlot(obj, features = features, assay = "RNA", group.by = group.by)
  ord <- levels(dp$data$features.plot)

  dp +
    scale_x_discrete(
      breaks = ord,
      labels = function(x) to_sym(sub("\n.*$", "", x))
    ) +
    theme(
      axis.text.y           = element_text(size = 8),
      axis.text.x           = element_text(angle = 90, hjust = 1, size = 10),
      strip.placement       = "outside",
      strip.text            = element_text(size = 10),
      panel.spacing         = unit(1, "lines"),
      axis.title            = element_blank(),
      plot.caption.position = "plot",
      plot.caption          = element_text(hjust = 0.5),
      plot.margin           = margin(t = 5, r = 12, b = 30, l = 12)
    ) +
    coord_cartesian(clip = "off") +
    scale_color_gradient2(
      low      = "navy",
      mid      = "white",
      high     = "firebrick3",
      midpoint = 0,
      limits   = c(-2.5, 2.5)
    )
}

#' Tema minimalista para UMAPs
umap_theme <- function() {
  theme_void(base_size = 12) +
    theme(
      legend.title      = element_blank(),
      legend.text       = element_text(size = 8),
      legend.key.height = unit(0.5, "cm"),
      plot.title        = element_text(face = "plain", size = 14, hjust = 0.5)
    )
}

# 4. I/O -----------------------------------------------------------------------

#' Guarda un ggplot en píxeles con ragg (más rápido que el backend base)
#'
#' @param plot Objeto ggplot
#' @param file Ruta del archivo de salida
#' @param w    Ancho en píxeles  (default 3000)
#' @param h    Alto  en píxeles  (default 2000)
#' @param dpi  Resolución        (default 300)
save_plot_px <- function(plot, file, w = 3000, h = 2000, dpi = 300) {
  ggplot2::ggsave(
    filename = file,
    plot     = plot,
    width    = w,
    height   = h,
    units    = "px",
    dpi      = dpi,
    device   = ragg::agg_png,
    bg       = "white"
  )
}

# 5. Búsqueda de directorios de matrices 10x ----------------------------------

#' Localiza el directorio de matrices CellRanger para una muestra dada
find_10x_dir <- function(s, top_dirs) {
  hits <- top_dirs[grepl(s, basename(top_dirs), fixed = TRUE)]

  if (length(hits) == 0) stop("No se encontró carpeta para: ", s)
  if (length(hits) > 1)  stop(
    "Más de una carpeta coincide con '", s, "':\n",
    paste(basename(hits), collapse = "\n")
  )

  candidates <- c(
    hits[1],
    file.path(hits[1], "filtered_feature_bc_matrix"),
    file.path(hits[1], "raw_feature_bc_matrix")
  )
  ok <- candidates[
    file.exists(file.path(candidates, "matrix.mtx.gz")) |
    file.exists(file.path(candidates, "matrix.mtx"))
  ]
  if (length(ok) == 0) stop("Sin matrix.mtx(.gz) para '", s, "': ", hits[1])
  ok[1]
}

# 6. Supresión de output -------------------------------------------------------

#' Silencia mensajes, advertencias y stdout de una expresión
quiet <- function(expr) {
  invisible(suppressWarnings(suppressMessages(
    capture.output(eval.parent(substitute(expr)))
  )))
}

# ============================================================
# FUNCIONES ADICIONALES — hdWGCNA (scripts 09–14)
# ============================================================

# 8. Conversión gen ↔ símbolo (variante con AGI entre paréntesis) ------------

#' Devuelve "Símbolo (AGI)" si hay símbolo, o sólo "AGI" si no lo hay
to_sym_agi <- function(x) {
  sym <- unname(gene2sym[x])
  ifelse(is.na(sym) | sym == "", x, paste0(sym, " (", x, ")"))
}

# 9. Hash de vectores / listas (para control de cambios) ----------------------

#' Hash rápido de un vector de cadenas
hash_vec <- function(x, algo = "xxhash64") {
  digest::digest(as.character(x), algo = algo, serialize = TRUE)
}

#' Hash rápido de cualquier objeto R
hash_list <- function(x, algo = "xxhash64") {
  digest::digest(x, algo = algo, serialize = TRUE)
}

# 10. Mapa de tratamientos por tech_mode --------------------------------------

#' Devuelve lista(ctrl, trt) para la tecnología dada
#'
#' @param tech_mode "sc" (MD/WW) o "sn" (DR/WW)
trt_map <- function(tech_mode) {
  switch(
    tech_mode,
    sc = list(ctrl = "WW", trt = "MD"),
    sn = list(ctrl = "WW", trt = "DR"),
    stop("tech_mode debe ser 'sc' o 'sn'")
  )
}

# 11. Funciones de guardado seguro --------------------------------------------

#' Escribe un vector de texto a un archivo (crea carpeta si no existe)
save_txt <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, con = path)
}

#' Guarda un objeto R como .rds (crea carpeta si no existe)
save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = path)
}

#' Escribe un data.frame como CSV con readr (crea carpeta si no existe)
save_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, file = path)
}

#' Escribe un data.frame como TSV con readr (crea carpeta si no existe)
save_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(df, file = path)
}

#' Guarda un objeto Seurat como .qs2 o .rds (menú interactivo)
#'
#' En modo no interactivo, usa ext = "qs2" o ext = "rds" directamente.
#'
#' @param obj          Objeto a guardar
#' @param path         Ruta base (sin extensión, o con extensión que se sobreescribirá)
#' @param overwrite_ext Reemplazar extensión existente (default TRUE)
#' @param ext          Extensión fija ("qs2" / "rds"); si NULL se pregunta interactivamente
save_obj <- function(obj, path, overwrite_ext = TRUE, ext = NULL) {
  save_options <- c("qs2", "rds")

  if (is.null(ext)) {
    save_choice <- menu(save_options,
                        title = "Seleccionar formato de guardado:")
    if (save_choice == 0) stop("Selección cancelada.", call. = FALSE)
    ext <- save_options[[save_choice]]
  }

  stopifnot(ext %in% save_options)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  path_final <- if (overwrite_ext) {
    sub("\\.[^.]*$", paste0(".", ext), path)
  } else {
    if (!grepl("\\.[^.]*$", path)) paste0(path, ".", ext) else path
  }

  if (ext == "qs2") {
    qs2::qs_save(obj, file = path_final)
  } else {
    saveRDS(obj, file = path_final)
  }

  message("Guardado en: ", path_final)
  invisible(path_final)
}

# 12. Fraction sweep para SetupForWGCNA ---------------------------------------

#' Explora distintas fracciones de detección y devuelve n_genes por fracción
#'
#' Útil para elegir la fracción óptima antes de llamar a SetupForWGCNA.
#'
#' @param gc_obj    Objeto Seurat
#' @param frac_grid Vector numérico de fracciones a probar (decreasing)
#' @param group_by  Variable de agrupación (ej. "celltype" o "Sample")
#' @param tmp_prefix Prefijo para nombres de wgcna temporales
#' @param assay     Ensayo a usar
#' @param strict    Si TRUE, lanza error ante plateau o reordenamiento
#' @return tibble con columnas fraction, n_genes, genes_hash
fraction_sweep <- function(gc_obj,
                            frac_grid,
                            group_by,
                            tmp_prefix = "tmp_frac_",
                            assay      = "RNA",
                            strict     = TRUE) {
  if (length(frac_grid) < 2)
    stop("fraction_sweep: frac_grid debe tener >= 2 valores.", call. = FALSE)

  frac_grid <- sort(as.numeric(frac_grid))
  if (anyNA(frac_grid))
    stop("fraction_sweep: frac_grid tiene valores no numéricos.", call. = FALSE)
  if (anyDuplicated(frac_grid))
    stop("fraction_sweep: frac_grid tiene fracciones duplicadas.", call. = FALSE)

  out <- vector("list", length(frac_grid))
  pb  <- utils::txtProgressBar(min = 0, max = length(frac_grid), style = 3)
  on.exit(close(pb), add = TRUE)

  prev_n <- NULL; prev_hash <- NULL

  for (i in seq_along(frac_grid)) {
    f        <- frac_grid[[i]]
    tmp_name <- paste0(tmp_prefix, gsub("\\.", "_", as.character(f)))
    tmp      <- gc_obj

    tmp <- hdWGCNA::SetupForWGCNA(
      seurat_obj  = tmp,
      assay       = assay,
      gene_select = "fraction",
      fraction    = f,
      group.by    = group_by,
      wgcna_name  = tmp_name
    )

    genes  <- unique(as.character(tmp@misc[[tmp_name]]$wgcna_genes))
    if (length(genes) == 0)
      stop("fraction_sweep: wgcna_genes es 0 para fraction=", f, call. = FALSE)

    ghash  <- hash_vec(sort(genes))
    n_genes <- length(genes)

    if (strict && !is.null(prev_n)) {
      if (n_genes > prev_n)
        stop("fraction_sweep: n_genes aumentó al subir fraction (",
             frac_grid[[i-1]], " -> ", f, "): ", prev_n, " -> ", n_genes, call. = FALSE)
      if (n_genes == prev_n)
        stop("fraction_sweep: plateau en n_genes (", f, ")", call. = FALSE)
      if (identical(ghash, prev_hash))
        stop("fraction_sweep: set idéntico entre fracciones consecutivas", call. = FALSE)
    }

    out[[i]] <- tibble::tibble(fraction = f, n_genes = n_genes, genes_hash = ghash)
    prev_n   <- n_genes; prev_hash <- ghash
    utils::setTxtProgressBar(pb, i)
  }

  df <- dplyr::bind_rows(out)
  if (strict && (anyDuplicated(df$n_genes) || anyDuplicated(df$genes_hash)))
    stop("fraction_sweep: hay n_genes o hashes repetidos.", call. = FALSE)
  df
}

# 13. Utilidades hdWGCNA ------------------------------------------------------

#' Extrae los MEs (o los calcula si no existen)
#'
#' @param seurat_obj Objeto Seurat con hdWGCNA
#' @param harmonized TRUE = hMEs, FALSE = MEs crudos
#' @param wgcna_name Nombre del análisis hdWGCNA
#' @param group.by.vars Variable de batch para armonización (ej. "Sample")
#' @return Lista: $obj (Seurat actualizado), $MEs (matriz celdas × módulos)
get_mes_safe <- function(seurat_obj,
                          harmonized    = TRUE,
                          wgcna_name    = NULL,
                          group.by.vars = NULL) {
  out <- tryCatch(
    hdWGCNA::GetMEs(seurat_obj, harmonized = harmonized, wgcna_name = wgcna_name),
    error = function(e) NULL
  )
  if (is.null(out) || ncol(out) == 0) {
    message("No encontré ", if (harmonized) "hMEs" else "MEs",
            ". Corriendo ModuleEigengenes()...")
    seurat_obj <- hdWGCNA::ModuleEigengenes(
      seurat_obj,
      group.by.vars = group.by.vars,
      wgcna_name    = wgcna_name
    )
    out <- hdWGCNA::GetMEs(seurat_obj, harmonized = harmonized,
                           wgcna_name = wgcna_name)
  }
  list(obj = seurat_obj, MEs = out)
}

#' Extrae los sets de genes por módulo (sin grey)
#'
#' @param seurat_obj Objeto Seurat con hdWGCNA
#' @param wgcna_name Nombre del análisis hdWGCNA
#' @param prefix     Prefijo opcional para los nombres de módulo
#' @return Lista nombrada gene → módulo
get_module_sets <- function(seurat_obj,
                             wgcna_name,
                             prefix = "") {
  mods <- hdWGCNA::GetModules(seurat_obj, wgcna_name = wgcna_name) |>
    dplyr::filter(!module %in% c("grey", "gray"))

  gene_col <- dplyr::case_when(
    "gene"      %in% colnames(mods) ~ "gene",
    "gene_name" %in% colnames(mods) ~ "gene_name",
    TRUE ~ NA_character_
  )
  if (is.na(gene_col))
    stop("GetModules() sin columna gene / gene_name.")

  sets <- split(as.character(mods[[gene_col]]), as.character(mods$module))

  # orden numérico
  mod_num <- suppressWarnings(
    as.integer(gsub("^.*?(\\d+)$", "\\1", names(sets)))
  )
  sets <- sets[order(mod_num, names(sets), na.last = TRUE)]

  if (nzchar(prefix)) names(sets) <- paste0(prefix, "::", names(sets))
  sets
}

#' Estandariza columnas de una tabla DE para uso en fgsea / ORA
#'
#' Compatible con DESeq2 (stat_raw, log2FoldChange_shr) y edgeR/otros.
#'
#' @param df   Data.frame con resultados DE
#' @param which Etiqueta del contraste (ej. "leaf::WW_nlp7_vs_WT")
#' @return tibble con columnas: gene, padj, stat, logFC, which
standardize_deg <- function(df, which) {
  nms <- colnames(df)

  if (!("gene" %in% nms) && !is.null(rownames(df))) {
    df  <- tibble::rownames_to_column(df, var = "gene"); nms <- colnames(df)
  }

  gene_col <- dplyr::case_when(
    "gene"    %in% nms ~ "gene",
    "Gene"    %in% nms ~ "Gene",
    "gene_id" %in% nms ~ "gene_id",
    TRUE ~ NA_character_
  )
  padj_col <- dplyr::case_when(
    "padj_raw" %in% nms ~ "padj_raw",
    "padj"     %in% nms ~ "padj",
    "p_val_adj" %in% nms ~ "p_val_adj",
    "FDR"      %in% nms ~ "FDR",
    TRUE ~ NA_character_
  )
  stat_col <- dplyr::case_when(
    "stat_raw"  %in% nms ~ "stat_raw",
    "stat"      %in% nms ~ "stat",
    "statistic" %in% nms ~ "statistic",
    TRUE ~ NA_character_
  )
  lfc_col <- dplyr::case_when(
    "log2FoldChange_shr" %in% nms ~ "log2FoldChange_shr",
    "log2FoldChange_raw" %in% nms ~ "log2FoldChange_raw",
    "log2FoldChange"     %in% nms ~ "log2FoldChange",
    "logFC"              %in% nms ~ "logFC",
    "avg_log2FC"         %in% nms ~ "avg_log2FC",
    TRUE ~ NA_character_
  )

  if (is.na(gene_col) || is.na(padj_col) || (is.na(stat_col) && is.na(lfc_col)))
    stop("No pude detectar columnas clave en:\n", paste(nms, collapse = ", "))

  df |>
    dplyr::transmute(
      gene  = trimws(sub("\\.\\d+$", "", as.character(.data[[gene_col]]))),
      padj  = as.numeric(.data[[padj_col]]),
      stat  = if (!is.na(stat_col)) as.numeric(.data[[stat_col]]) else NA_real_,
      logFC = if (!is.na(lfc_col))  as.numeric(.data[[lfc_col]])  else NA_real_,
      which = which
    ) |>
    dplyr::filter(!is.na(gene), gene != "", !is.na(padj),
                  is.finite(padj), padj >= 0)
}

#' Construye un vector de ranking para fgsea a partir de una tabla DE estandarizada
#'
#' @param df          Salida de standardize_deg()
#' @param cap         Máximo de -log10(padj) cuando se usa logFC × -log10(padj)
#' @param jitter_ties Añadir ruido pequeño para romper empates
make_rank_vec <- function(df, cap = 50, jitter_ties = TRUE) {
  padj_clip <- pmax(df$padj, 1e-300)
  score <- if (!all(is.na(df$stat))) {
    df$stat
  } else {
    if (all(is.na(df$logFC))) stop("No hay 'stat' ni 'logFC' para rankear.")
    df$logFC * pmin(-log10(padj_clip), cap)
  }

  rk <- tibble::tibble(gene = df$gene, score = score) |>
    dplyr::filter(!is.na(gene), !is.na(score), is.finite(score)) |>
    dplyr::group_by(gene) |>
    dplyr::summarise(score = score[which.max(abs(score))], .groups = "drop")

  if (jitter_ties)
    rk <- dplyr::mutate(rk, score = score + rnorm(dplyr::n(), 0, 1e-9))

  rk |> dplyr::arrange(dplyr::desc(score)) |>
    (\(x) setNames(x$score, x$gene))()
}


#' UMAP con círculos numerados en los centroides y leyenda compacta lateral
#'
#' @param seurat_obj    Objeto Seurat
#' @param celltype_col  Columna de metadatos con el tipo celular
#' @param reduction     Nombre de la reducción (ej. "umap.harmony")
#' @param ct_cols       Vector nombrado de colores por tipo celular
#' @param ...           Parámetros estéticos opcionales (ver cuerpo de la función)
#' @return Lista con $plot (completo), $umap y $legend
fancy_umap <- function(
    seurat_obj,
    celltype_col        = "celltype",
    reduction           = "umap.harmony",
    ct_cols,
    pt.size             = 0.15,
    pt.alpha            = 0.80,
    circle.size         = 4.8,
    circle.fill         = "white",
    circle.stroke       = 0.9,
    circle.color        = "black",
    num.size            = 2.0,
    num.color           = "black",
    legend_title        = "Cell Type",
    legend_title_size   = 5.2,
    legend_row_spacing  = 0.62,
    legend_title_pad    = 0.55,
    legend_dot_size     = 2.6,
    legend_num_circle_size = 3.0,
    legend_num_text_size   = 3.0,
    legend_label_size   = 3.4,
    legend_num_circle_lwd  = 0.9,
    legend_x_dot        = 0.00,
    legend_x_num        = 0.38,
    legend_x_label      = 0.72,
    legend_xlim_right   = 2.35,
    legend_width        = 0.38
) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  if (!celltype_col %in% colnames(seurat_obj@meta.data))
    stop("No existe '", celltype_col, "' en metadata.")
  if (!reduction %in% names(seurat_obj@reductions))
    stop("No existe la reducción '", reduction, "'.")
  if (missing(ct_cols) || is.null(ct_cols))
    stop("Debes pasar 'ct_cols'.")

  emb <- Seurat::Embeddings(seurat_obj, reduction = reduction)
  df  <- as.data.frame(emb[, 1:2, drop = FALSE])
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  df[[celltype_col]] <- seurat_obj@meta.data[[celltype_col]]

  # Orden de tipos celulares
  if (!is.null(names(ct_cols)) && all(nzchar(names(ct_cols)))) {
    ct_levels <- intersect(names(ct_cols), unique(as.character(df[[celltype_col]])))
    if (length(ct_levels) == 0)
      stop("names(ct_cols) no coincide con valores de '", celltype_col, "'.")
    ct_cols_use <- ct_cols[ct_levels]
  } else {
    ct_levels   <- sort(unique(as.character(df[[celltype_col]])))
    ct_cols_use <- ct_cols
    names(ct_cols_use) <- ct_levels
  }

  df[[celltype_col]] <- factor(as.character(df[[celltype_col]]), levels = ct_levels)
  ct_num <- setNames(seq_along(ct_levels), ct_levels)
  df$num <- unname(ct_num[as.character(df[[celltype_col]])])

  centers <- df |>
    dplyr::group_by(.data[[celltype_col]]) |>
    dplyr::summarise(
      UMAP_1 = median(.data$UMAP_1),
      UMAP_2 = median(.data$UMAP_2),
      ct_chr = dplyr::first(as.character(.data[[celltype_col]])),
      .groups = "drop"
    ) |>
    dplyr::mutate(num = unname(ct_num[ct_chr])) |>
    dplyr::select(-ct_chr)

  # UMAP principal
  p_umap <- ggplot(df, aes(x = .data$UMAP_1, y = .data$UMAP_2)) +
    geom_point(aes(color = .data[[celltype_col]]),
               size = pt.size, alpha = pt.alpha, show.legend = FALSE) +
    geom_point(data = centers, inherit.aes = FALSE,
               aes(x = .data$UMAP_1, y = .data$UMAP_2),
               shape = 21, size = circle.size,
               fill = circle.fill, color = circle.color, stroke = circle.stroke) +
    geom_text(data = centers, inherit.aes = FALSE,
              aes(x = .data$UMAP_1, y = .data$UMAP_2, label = .data$num),
              size = num.size, color = num.color, fontface = "bold") +
    scale_color_manual(values = ct_cols_use) +
    theme_void()

  # Leyenda manual compacta
  leg <- data.frame(
    celltype = factor(ct_levels, levels = ct_levels),
    num      = unname(ct_num[ct_levels]),
    y        = rev(seq_along(ct_levels)) * legend_row_spacing
  )
  y_top <- max(leg$y)
  y_min <- min(leg$y)

  p_leg <- ggplot(leg, aes(y = .data$y)) +
    annotate("text",
             x = legend_x_dot, y = y_top + legend_title_pad,
             label = legend_title, hjust = 0,
             fontface = "bold", size = legend_title_size) +
    geom_point(aes(x = legend_x_dot, color = .data$celltype),
               shape = 16, size = legend_dot_size) +
    geom_point(aes(x = legend_x_num),
               shape = 21, size = legend_num_circle_size,
               fill = "white", color = "black", stroke = legend_num_circle_lwd) +
    geom_text(aes(x = legend_x_num, label = .data$num),
              fontface = "bold", size = legend_num_text_size) +
    geom_text(aes(x = legend_x_label, label = .data$celltype),
              hjust = 0, size = legend_label_size) +
    scale_color_manual(values = ct_cols_use, guide = "none") +
    coord_cartesian(
      xlim = c(legend_x_dot - 0.05, legend_xlim_right),
      ylim = c(y_min - 0.35 * legend_row_spacing,
               y_top + legend_title_pad + 0.35),
      clip = "off"
    ) +
    theme_void()

  p <- p_umap + p_leg + patchwork::plot_layout(widths = c(1, legend_width))

  list(plot = p, umap = p_umap, legend = p_leg)
}

# 04_bulk_corr.R
# Correlación de Spearman entre pseudobulk y bulk RNA-seq.
#
# Mejora respecto a v1:
#   • Correlación separada por condición (WW vs WW, DR/MD vs DR/MD)
#   • Pairwise entre todas las réplicas de bulk × pseudobulk
#   • Outliers consenso: genes en top-1% de diferencia de rangos en >50% de pares
#   • Genes outliers = unión entre condiciones → más conservador
#
# Entrada:  act_filtered.qs2
# Salida:   act_bulkOK.qs2
#           <project>_spearman_outliers.txt
#           spearman_bulk_vs_pseudobulk_<project>_<cond>.png
#           genes_outliers_spearman_<project>.csv
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

# 1. Cargar objeto -------------------------------------------------------------

act <- qs_read("act_filtered.qs2")
DefaultAssay(act) <- "RNA"
info(act)

# 2. Pseudobulk por muestra ---------------------------------------------------

pb <- AggregateExpression(
  act,
  group.by     = "Sample",
  assays       = "RNA",
  slot         = "counts",
  return.seurat = FALSE
)
counts_pb <- as.matrix(pb$RNA)

# Metadata de pseudobulk
cond_levels <- if (grepl(".*TB$", project)) c("WW", "MD") else c("WW", "DR")

meta_pb <- data.frame(
  sample    = colnames(counts_pb),
  condition = sub("-.*$", "", colnames(counts_pb))   # "MD", "WW", "DR"
) |>
  mutate(condition = factor(condition, levels = cond_levels))
rownames(meta_pb) <- meta_pb$sample

dge_pb <- DGEList(
  counts  = counts_pb,
  group   = meta_pb$condition,
  samples = meta_pb
)
design_pb <- model.matrix(~ 0 + condition, data = dge_pb$samples)
dge_pb    <- dge_pb[filterByExpr(dge_pb, design_pb), , keep.lib.sizes = FALSE]
dge_pb    <- calcNormFactors(dge_pb, method = "TMM")
logCPM_pb <- cpm(dge_pb, log = TRUE, prior.count = 1)

message("Pseudobulk: ", nrow(logCPM_pb), " genes × ", ncol(logCPM_pb), " muestras")

# 3. Bulk RNA-seq --------------------------------------------------------------

counts_bulk <- switch(
  tech_mode,
  sc = {
    read_excel(
      file.path(raw_path, "tenorio/bulk/geo/tenorio_bulk.xlsx"),
      n_max = 32833
    ) |>
      dplyr::select(1, 12:15) |>
      column_to_rownames(var = "Gene") |>
      as.matrix()
  },
  sn = {
    read.delim(file.path(raw_path, "natanella/bulk/Run1.Counts.txt"))[,
      c("D30.1", "D30.2", "D30.3", "ww.0.1", "ww.0.2", "ww.0.3")]
  }
)

metadata_bulk <- switch(
  tech_mode,
  sc = read.table(file.path(raw_path, "tenorio/bulk/geo/metadata.txt"),
                  header = TRUE, row.names = 1),
  sn = read.table(file.path(raw_path, "natanella/bulk/metadata.txt"),
                  header = TRUE, row.names = 1)
)
metadata_bulk <- metadata_bulk[colnames(counts_bulk), , drop = FALSE]

metadata_bulk <- switch(
  tech_mode,
  sc = metadata_bulk |>
    mutate(treatment = factor(treatment, levels = c("well_watered", "mild_drought")),
           replicate = factor(replicate, levels = c("1", "2"))),
  sn = metadata_bulk |>
    mutate(treatment = factor(treatment, levels = c("well_watered", "drought")),
           replicate = factor(replicate, levels = c("1", "2", "3")))
)

stopifnot(all(colnames(counts_bulk) == rownames(metadata_bulk)))

y      <- DGEList(counts = counts_bulk, group = metadata_bulk$treatment,
                  samples = metadata_bulk)
design <- model.matrix(~ 0 + treatment, data = y$samples)
y      <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
y      <- calcNormFactors(y, method = "TMM")
logCPM_bulk <- cpm(y, log = TRUE, prior.count = 1)

message("Bulk: ", nrow(logCPM_bulk), " genes × ", ncol(logCPM_bulk), " muestras")

# 4. Función auxiliar: outliers consenso por par de condición -----------------

#' Para una condición dada, calcula outliers de Spearman en todos los pares
#' bulk × pseudobulk.
#'
#' @param bulk_mat  Matriz logCPM bulk (genes × réplicas de la condición)
#' @param pb_mat    Matriz logCPM pseudobulk (genes × réplicas de la condición)
#' @param p_outlier Percentil de corte para considerar outlier (default 0.99)
#' @param min_prop  Proporción mínima de pares en que debe aparecer (default 0.50)
#' @return Vector de genes outlier consenso
spearman_outliers_cond <- function(bulk_mat, pb_mat,
                                   p_outlier = 0.99, min_prop = 0.50) {
  common <- intersect(rownames(bulk_mat), rownames(pb_mat))
  message("  Genes comunes: ", length(common))

  bm <- bulk_mat[common, , drop = FALSE]
  pm <- pb_mat[common, , drop = FALSE]

  # Correlación pairwise (diagnóstico)
  cor_mat <- outer(
    colnames(bm), colnames(pm),
    Vectorize(\(a, b) cor(bm[, a], pm[, b],
                          method = "spearman", use = "pairwise.complete.obs"))
  )
  dimnames(cor_mat) <- list(colnames(bm), colnames(pm))
  message("  Correlaciones Spearman:")
  print(round(cor_mat, 3))

  # Rangos
  bulk_ranks <- apply(bm, 2, rank, ties.method = "average", na.last = "keep")
  pb_ranks   <- apply(pm, 2, rank, ties.method = "average", na.last = "keep")

  pairs <- tidyr::crossing(bulk = colnames(bulk_ranks), pseudo = colnames(pb_ranks))

  outlier_tbl <- purrr::pmap_dfr(pairs, function(bulk, pseudo) {
    diff <- abs(bulk_ranks[, bulk] - pb_ranks[, pseudo])
    thr  <- quantile(diff, p_outlier, na.rm = TRUE)
    tibble(
      bulk           = bulk,
      pseudo         = pseudo,
      genes_outliers = list(common[diff >= thr]),
      n              = sum(diff >= thr, na.rm = TRUE)
    )
  })

  # Consenso: genes en > min_prop de los pares
  cons <- outlier_tbl |>
    dplyr::filter(lengths(genes_outliers) > 0) |>
    tidyr::unnest(genes_outliers) |>
    dplyr::count(genes_outliers, sort = TRUE) |>
    dplyr::rename(gene = genes_outliers, n_pairs = n) |>
    dplyr::mutate(prop_pairs = n_pairs / nrow(outlier_tbl))

  cons |> dplyr::filter(prop_pairs > min_prop) |> dplyr::pull(gene)
}

# 5. Separar bulk y pseudobulk por condición ----------------------------------

if (tech_mode == "sc") {
  ww_bulk_cols <- grep("_ww_", colnames(logCPM_bulk), value = TRUE)
  md_bulk_cols <- grep("_md_", colnames(logCPM_bulk), value = TRUE)
  bulk_WW <- logCPM_bulk[, ww_bulk_cols, drop = FALSE]
  bulk_MD <- logCPM_bulk[, md_bulk_cols, drop = FALSE]

  pb_WW <- logCPM_pb[, grep("^WW", colnames(logCPM_pb), value = TRUE), drop = FALSE]
  pb_MD <- logCPM_pb[, grep("^MD", colnames(logCPM_pb), value = TRUE), drop = FALSE]

  message("\n--- Condición WW ---")
  outliers_WW <- spearman_outliers_cond(bulk_WW, pb_WW)
  message("  Outliers WW: ", length(outliers_WW))

  message("\n--- Condición MD ---")
  outliers_MD <- spearman_outliers_cond(bulk_MD, pb_MD)
  message("  Outliers MD: ", length(outliers_MD))

} else {   # sn / natanella
  bulk_WW <- logCPM_bulk[, grep("^ww", colnames(logCPM_bulk), value = TRUE), drop = FALSE]
  bulk_DR <- logCPM_bulk[, grep("^D30", colnames(logCPM_bulk), value = TRUE), drop = FALSE]

  pb_WW <- logCPM_pb[, grep("^WW", colnames(logCPM_pb), value = TRUE), drop = FALSE]
  pb_DR <- logCPM_pb[, grep("^DR", colnames(logCPM_pb), value = TRUE), drop = FALSE]

  message("\n--- Condición WW ---")
  outliers_WW <- spearman_outliers_cond(bulk_WW, pb_WW)
  message("  Outliers WW: ", length(outliers_WW))

  message("\n--- Condición DR ---")
  outliers_DR <- spearman_outliers_cond(bulk_DR, pb_DR)
  message("  Outliers DR: ", length(outliers_DR))

  outliers_MD <- outliers_DR   # alias para el bloque de unión
}

# 6. Unión de outliers entre condiciones --------------------------------------

final_outliers <- union(outliers_WW, outliers_MD)

message("\nResumen outliers Spearman:")
message("  WW:    ", length(outliers_WW))
message("  MD/DR: ", length(outliers_MD))
message("  Intersección: ", length(intersect(outliers_WW, outliers_MD)))
message("  Unión (final): ", length(final_outliers))

# Exportar lista de outliers
write.table(
  final_outliers,
  file      = paste0(project, "_spearman_outliers.txt"),
  quote     = FALSE, col.names = FALSE, row.names = FALSE
)

# 7. Crear ensayo "bulkOK" -----------------------------------------------------

act      <- JoinLayers(act)
genes_ok <- setdiff(rownames(act), final_outliers)

act[["bulkOK"]] <- CreateAssayObject(
  counts = GetAssayData(act, assay = "RNA", layer = "counts")[genes_ok, ]
)

message("Genes retenidos en bulkOK: ", length(genes_ok), " / ", nrow(act))

DefaultAssay(act) <- "RNA"
act <- NormalizeData(act)

qs_save(act, "act_bulkOK.qs2")
message("Objeto guardado en: act_bulkOK.qs2")
info(act)

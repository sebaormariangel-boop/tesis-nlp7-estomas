# 13_hdwgcna_enrichment.R
# Enriquecimiento funcional centrado en los módulos hdWGCNA:
#   • fgsea (módulos hdWGCNA nat_gc + ten_all) vs bulk hoja y estomas
#   • ORA + GSEA (firmas de tipo celular sc) vs bulk
#   • Análisis de hub genes (tabla enriquecida con bulk + overlap)
#   • NLP7 (AT4G24020): pseudobulk en Guard Cells vs hME del módulo
#
# Entradas:
#   qs2/*_hdWGCNA_object.qs2      (módulos hdWGCNA)
#   tables_DROUGHT / tables_STOMATA (archivos DE_combined_*.tsv)
#   tables/*_markers_.tsv          (FindAllMarkers de 08_annotation.R)
#   tables/overlap_genes_by_cell.tsv (generado externamente)
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

suppressPackageStartupMessages({
  library(hdWGCNA)
  library(WGCNA)
  library(fgsea)
  library(genesectR)
})

for (d in c("figures", "tables")) dir.create(d, showWarnings = FALSE)

# Parámetros globales ---------------------------------------------------------

ALPHA_BULK    <- 0.05
ALPHA_MARKERS <- 0.05
TOP_N_MARKERS <- 200
MIN_SIZE      <- 15
MAX_SIZE      <- 1000
RANK_CAP      <- 50
JITTER_TIES   <- TRUE
FGSEA_EPS     <- 0

# 1. Cargar objetos hdWGCNA ---------------------------------------------------

nat_gc  <- qs_read(file.path(qs2_dir, "nat_subsetting_hdWGCNA_object.qs2"))
ten_all <- qs_read(file.path(qs2_dir, "ten_full_hdWGCNA_object.qs2"))

# 2. Leer archivos bulk (DE_combined_*.tsv) -----------------------------------

pick_de_files <- function(dir_path, run_all = TRUE, preferred = NULL, label = "") {
  if (!dir.exists(dir_path)) stop("No existe: ", dir_path, " (", label, ")")
  files <- list.files(dir_path, pattern = "^DE_combined_.*\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop("Sin DE_combined_*.tsv en: ", dir_path)
  if (run_all) return(files)
  hit <- if (!is.null(preferred)) files[grepl(preferred, basename(files))] else files[1]
  if (length(hit) == 0) files[1] else hit[1]
}

leaf_files  <- pick_de_files(bulk_leaf_dir,  label = "LEAF")
stoma_files <- pick_de_files(bulk_stoma_dir, label = "STOMA")

message("Leaf:  ", length(leaf_files),  " archivo(s)")
message("Stoma: ", length(stoma_files), " archivo(s)")

read_bulk_file <- function(fp, bulk_label) {
  df  <- readr::read_tsv(fp, show_col_types = FALSE)
  tag <- sub("^DE_combined_", "", sub("\\.tsv$", "", basename(fp)))
  list(
    std  = standardize_deg(df, which = paste0(bulk_label, "::", tag)),
    tag  = tag,
    bulk = bulk_label,
    file = fp
  )
}

bulk_list <- c(
  lapply(leaf_files,  read_bulk_file, bulk_label = "leaf"),
  lapply(stoma_files, read_bulk_file, bulk_label = "stoma")
)

# 3. fgsea: módulos hdWGCNA vs bulk ------------------------------------------

module_info <- list(
  nat_gc  = list(
    pathways   = get_module_sets(nat_gc,  nat_gc@misc$active_wgcna),
    wgcna_name = nat_gc@misc$active_wgcna
  ),
  ten_all = list(
    pathways   = get_module_sets(ten_all, ten_all@misc$active_wgcna),
    wgcna_name = ten_all@misc$active_wgcna
  )
)

run_fgsea_bulk <- function(pathways, ranks, source, wgcna_name, bulk, tag) {
  universe <- names(ranks)
  pw_use   <- lapply(pathways, intersect, universe)
  pw_use   <- pw_use[lengths(pw_use) >= MIN_SIZE & lengths(pw_use) <= MAX_SIZE]
  if (length(pw_use) == 0) return(tibble::tibble())

  res <- fgsea::fgseaMultilevel(
    pathways = pw_use, stats = ranks,
    minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = FGSEA_EPS
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(module_source = source, wgcna_name = wgcna_name,
                  bulk = bulk, contrast_tag = tag) |>
    dplyr::arrange(padj, dplyr::desc(abs(NES))) |>
    dplyr::mutate(
      leadingEdge_n   = lengths(leadingEdge),
      leadingEdge_str = vapply(leadingEdge, \(x) paste(head(x, 50), collapse = ";"), ""),
      .keep = "all"
    ) |>
    dplyr::select(-leadingEdge)

  res
}

all_fgsea <- purrr::map_dfr(names(module_info), function(ms) {
  pw <- module_info[[ms]]$pathways
  w  <- module_info[[ms]]$wgcna_name

  purrr::map_dfr(bulk_list, function(bk) {
    rk <- make_rank_vec(bk$std, cap = RANK_CAP, jitter_ties = JITTER_TIES)
    run_fgsea_bulk(pw, rk, source = ms, wgcna_name = w,
                   bulk = bk$bulk, tag = bk$tag)
  })
})

save_tsv(all_fgsea, file.path("tables", "fgsea_hdwgcna_modules_vs_bulk.tsv"))

top20_fgsea <- all_fgsea |>
  dplyr::group_by(module_source, bulk, contrast_tag) |>
  dplyr::arrange(padj, dplyr::desc(abs(NES))) |>
  dplyr::slice_head(n = 20) |>
  dplyr::ungroup()

save_tsv(top20_fgsea, file.path("tables", "fgsea_top20_hdwgcna_modules_vs_bulk.tsv"))

p_fgsea <- all_fgsea |>
  dplyr::mutate(mlog10 = -log10(pmax(padj, 1e-300)), sig = padj < 0.05) |>
  ggplot(aes(x = NES, y = mlog10)) +
  geom_point(aes(alpha = sig), size = 2) +
  facet_grid(module_source ~ bulk, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "NES", y = "-log10(FDR)",
       title = "fgsea: módulos hdWGCNA vs bulk") +
  theme_classic() + theme(legend.position = "none")

ggsave(file.path("figures", "fgsea_hdwgcna_modules_vs_bulk.png"),
       p_fgsea, width = 10, height = 6, dpi = 300)

# 4. ORA + GSEA: firmas de tipo celular sc vs bulk ---------------------------

markers_path <- file.path("tables", paste0(project, "_markers_.tsv"))
if (!file.exists(markers_path)) {
  message("Tabla de marcadores no encontrada: ", markers_path)
  message("Ejecuta 08_annotation.R (FindAllMarkers) y copia el .tsv aquí.")
} else {
  mk <- readr::read_tsv(markers_path, show_col_types = FALSE) |>
    dplyr::mutate(gene = trimws(sub("\\.\\d+$", "", gene)))

  marker_sets <- mk |>
    dplyr::filter(p_val_adj < ALPHA_MARKERS) |>
    dplyr::arrange(cluster, dplyr::desc(avg_log2FC)) |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n = TOP_N_MARKERS) |>
    dplyr::summarise(genes = list(unique(gene)), .groups = "drop") |>
    dplyr::filter(lengths(genes) >= MIN_SIZE) |>
    (\(x) setNames(x$genes, x$cluster))()

  message("Firmas de tipo celular: ", length(marker_sets), " clusters")

  # ORA (Fisher) por bulk-file
  fisher_ora <- function(target, markers, universe) {
    target  <- intersect(unique(target), universe)
    markers <- intersect(unique(markers), universe)
    x <- length(intersect(target, markers))
    m <- length(target); k <- length(markers); N <- length(universe)
    mat <- matrix(c(x, m-x, k-x, N-m-k+x), 2)
    ft  <- fisher.test(mat, alternative = "greater")
    tibble::tibble(overlap = x, size_target = m, size_markers = k,
                   OR = unname(ft$estimate), pval = ft$p.value)
  }

  ora_all <- purrr::map_dfr(bulk_list, function(bk) {
    std      <- bk$std
    universe <- unique(std$gene)
    dirs     <- std |>
      dplyr::mutate(dir = dplyr::case_when(
        padj < ALPHA_BULK & logFC > 0 ~ "UP",
        padj < ALPHA_BULK & logFC < 0 ~ "DOWN",
        TRUE ~ "NS"
      )) |>
      dplyr::filter(dir != "NS") |>
      dplyr::group_by(dir) |>
      dplyr::summarise(genes = list(unique(gene)), .groups = "drop")

    if (nrow(dirs) == 0) return(tibble::tibble())

    tidyr::expand_grid(dir = dirs$dir, cluster = names(marker_sets)) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        tgt  = list(dirs$genes[[match(dir, dirs$dir)]]),
        mkrs = list(marker_sets[[cluster]]),
        stats = list(fisher_ora(tgt, mkrs, universe))
      ) |>
      tidyr::unnest(stats) |>
      dplyr::ungroup() |>
      dplyr::mutate(padj = p.adjust(pval, method = "BH"),
                    bulk = bk$bulk, contrast_tag = bk$tag,
                    set  = paste0(bk$which, "_", dir))
  })

  save_tsv(ora_all, file.path("tables", "ORA_bulk_vs_celltype_signatures.tsv"))

  # GSEA de tipo celular
  gsea_ct <- purrr::map_dfr(bulk_list, function(bk) {
    rk      <- make_rank_vec(bk$std, cap = RANK_CAP, jitter_ties = JITTER_TIES)
    universe <- names(rk)
    pw_use  <- lapply(marker_sets, intersect, universe)
    pw_use  <- pw_use[lengths(pw_use) >= MIN_SIZE & lengths(pw_use) <= MAX_SIZE]
    if (length(pw_use) == 0) return(tibble::tibble())

    fgsea::fgseaMultilevel(pw_use, rk, minSize = MIN_SIZE, maxSize = MAX_SIZE,
                           eps = FGSEA_EPS) |>
      tibble::as_tibble() |>
      dplyr::mutate(bulk = bk$bulk, contrast_tag = bk$tag,
                    leadingEdge_n   = lengths(leadingEdge),
                    leadingEdge_str = vapply(leadingEdge,
                                            \(x) paste(head(x,50), collapse = ";"), "")) |>
      dplyr::select(-leadingEdge)
  })

  save_tsv(gsea_ct, file.path("tables", "GSEA_bulk_vs_celltype_signatures.tsv"))
  message("ORA + GSEA de tipo celular guardados.")
}

# 5. NLP7 (AT4G24020): pseudobulk en Guard Cells vs hME ----------------------

obj_for_nlp7 <- nat_gc   # cambiar a ten_all si es necesario
wgcna_nlp7   <- obj_for_nlp7@misc$active_wgcna
gc_label      <- "Guard Cells"

tmp_mes <- get_mes_safe(obj_for_nlp7, harmonized = TRUE,
                         wgcna_name = wgcna_nlp7, group.by.vars = "Sample")
obj_for_nlp7 <- tmp_mes$obj
hMEs_nlp7    <- tmp_mes$MEs

# Módulo que contiene NLP7 (o elegir el mayor ΔGC)
mods_tbl_nlp7 <- hdWGCNA::GetModules(obj_for_nlp7, wgcna_name = wgcna_nlp7)
nlp7_mod <- mods_tbl_nlp7 |>
  dplyr::filter(gene_name == nlp7_gene) |>
  dplyr::pull(module) |>
  unique()

mod_for_nlp7 <- if (length(nlp7_mod) == 1 && nlp7_mod %in% colnames(hMEs_nlp7)) {
  message("NLP7 en módulo: ", nlp7_mod)
  nlp7_mod
} else {
  all_mods <- setdiff(colnames(hMEs_nlp7), "grey")
  message("NLP7 no asignado. Eligiendo módulo con mayor varianza...")
  all_mods[which.max(apply(hMEs_nlp7[, all_mods, drop = FALSE], 2, var))]
}

cells_gc <- colnames(obj_for_nlp7)[obj_for_nlp7$celltype == gc_label]
stopifnot(length(cells_gc) > 0)

# Expresión NLP7 + hME del módulo por célula
cell_df <- tibble::tibble(
  cell        = cells_gc,
  Sample      = obj_for_nlp7$Sample[cells_gc],
  nlp7_data   = as.numeric(
    GetAssayData(obj_for_nlp7, assay = "RNA", layer = "data")[nlp7_gene, cells_gc]
  ),
  nlp7_counts = as.numeric(
    GetAssayData(obj_for_nlp7, assay = "RNA", layer = "counts")[nlp7_gene, cells_gc]
  ),
  module_hME  = as.numeric(hMEs_nlp7[cells_gc, mod_for_nlp7])
)

# Pseudobulk por muestra
pb_nlp7 <- cell_df |>
  dplyr::group_by(Sample) |>
  dplyr::summarise(
    n_cells           = dplyr::n(),
    mean_nlp7         = mean(nlp7_data, na.rm = TRUE),
    mean_module       = mean(module_hME, na.rm = TRUE),
    prop_detect_nlp7  = mean(nlp7_counts > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(nlp7_detected = prop_detect_nlp7 > 0)

print(pb_nlp7)
save_tsv(pb_nlp7,
  file.path("tables", paste0("pseudobulk_GC_NLP7_vs_hME_", mod_for_nlp7, ".tsv")))

cor_nlp7 <- cor.test(pb_nlp7$mean_nlp7, pb_nlp7$mean_module,
                      method = "spearman", exact = FALSE)

p_scatter <- ggplot(pb_nlp7, aes(mean_nlp7, mean_module)) +
  geom_point(aes(size = n_cells, alpha = prop_detect_nlp7)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title    = paste0("Guard Cells pseudobulk: NLP7 vs hME(", mod_for_nlp7, ")"),
    subtitle = paste0("Spearman rho = ", signif(unname(cor_nlp7$estimate), 3),
                      ", p = ", signif(cor_nlp7$p.value, 3)),
    x = paste0("Mean ", nlp7_gene, " (data)"),
    y = paste0("Mean hME(", mod_for_nlp7, ")")
  ) + theme_classic()

ggsave(file.path("figures", paste0("scatter_GC_NLP7_hME_", mod_for_nlp7, ".png")),
       p_scatter, width = 7, height = 5, dpi = 300)

message("Script 13 completado.")
writeLines(capture.output(sessionInfo()),
           file.path("tables", "sessionInfo_enrichment.txt"))

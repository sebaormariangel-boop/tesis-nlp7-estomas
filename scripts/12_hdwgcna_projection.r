# 12_hdwgcna_projection.R
# Proyección de módulos entre datasets y análisis de preservación.
#   • ProjectModules (REF → QUERY)
#   • ModulePreservation (Z-summary)
#   • genesectR: solapamiento entre módulos de ambos datasets
#   • Heatmap z-score de scores por celltype × tratamiento
#   • Violin del módulo más específico para Guard Cells (ΔGC)
#
# Entrada:  qs2/*_hdWGCNA_object.qs2  (uno por dataset)
#           qs2/*_final_annotated.qs2 (objeto completo anotado)
# Salida:   figures/projection_*/..., tables/projection_*/...
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

suppressPackageStartupMessages({
  library(hdWGCNA)
  library(WGCNA)
  library(genesectR)
})
enableWGCNAThreads(nThreads = 56)

for (d in c("figures", "tables", "qs2")) dir.create(d, showWarnings = FALSE)

# 1. Cargar objetos -----------------------------------------------------------
#  Ajusta los paths a los archivos generados en pasos anteriores.

nat_gc  <- qs_read(file.path(qs2_dir, "nat_subsetting_hdWGCNA_object.qs2"))
ten_all <- qs_read(file.path(qs2_dir, "ten_full_hdWGCNA_object.qs2"))
nat_all <- qs_read(file.path(qs2_dir, "natanella_final_annotated.qs2"))

# 2. Selección interactiva de REF y QUERY ------------------------------------

ref_opts <- c("nat_gc (Guard Cells natanella)", "ten_all (completo tenorio)")
ref_ch   <- menu(ref_opts, title = "Objeto REF (contiene módulos hdWGCNA):")
if (ref_ch == 0) stop("Cancelado.")
ref_obj  <- list(nat_gc, ten_all)[[ref_ch]]
ref_name <- c("nat_gc", "ten_all")[[ref_ch]]

qry_opts <- c("nat_all (completo natanella)", "ten_all (completo tenorio)")
qry_ch   <- menu(qry_opts, title = "Objeto QUERY (destino):")
if (qry_ch == 0) stop("Cancelado.")
query_obj  <- list(nat_all, ten_all)[[qry_ch]]
query_name <- c("nat_all", "ten_all")[[qry_ch]]

message("REF: ", ref_name, " | QUERY: ", query_name)

wgcna_ref  <- ref_obj@misc$active_wgcna
wgcna_proj <- paste0("proj_from_", ref_name)

meta_celltype <- "celltype"
gc_label      <- "Guard Cells"
reduction_use <- "umap.harmony"

OUT_FIG <- file.path("figures", paste0("projection_", ref_name, "_to_", query_name))
OUT_TAB <- file.path("tables",  paste0("projection_", ref_name, "_to_", query_name))
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# 3. Proyectar módulos -------------------------------------------------------

DefaultAssay(ref_obj)   <- "RNA"
DefaultAssay(query_obj) <- "RNA"

query_obj <- ProjectModules(
  seurat_obj      = query_obj,
  seurat_ref      = ref_obj,
  wgcna_name      = wgcna_ref,
  wgcna_name_proj = wgcna_proj,
  group.by.vars   = NULL
)

# kME en query (puede dar "TOM not found" — es normal en proyectado)
query_obj <- tryCatch(
  ModuleConnectivity(query_obj, group.by = meta_celltype,
                     group_name = gc_label, wgcna_name = wgcna_proj),
  error = function(e) {
    message("kME en query falló (TOM no encontrado): ", conditionMessage(e))
    query_obj
  }
)

# UCell scores
query_obj <- ModuleExprScore(query_obj, method = "UCell",
                              wgcna_name = wgcna_proj)

# Feature plots de hMEs
plot_list <- ModuleFeaturePlot(
  query_obj, wgcna_name = wgcna_proj,
  reduction = reduction_use, features = "hMEs",
  order_points = TRUE, restrict_range = FALSE
)
p_proj <- patchwork::wrap_plots(plot_list, ncol = 4)
save_plot_px(
  p_proj,
  file.path(OUT_FIG, paste0("umap_projected_hMEs_", ref_name, "_to_", query_name, ".png")),
  w = 3500, h = 2500
)

# 4. ΔGC: ranking de módulos por especificidad en Guard Cells ----------------

MEs <- GetMEs(query_obj, wgcna_name = wgcna_proj)
df  <- dplyr::bind_cols(
  dplyr::select(query_obj@meta.data, all_of(meta_celltype)),
  as.data.frame(MEs)
) |> dplyr::rename(celltype = all_of(meta_celltype))

med <- df |>
  dplyr::group_by(celltype) |>
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), median, na.rm = TRUE),
                   .groups = "drop")

spec <- med |>
  tidyr::pivot_longer(-celltype, names_to = "module", values_to = "med_hME") |>
  dplyr::group_by(module) |>
  dplyr::summarise(
    gc        = med_hME[celltype == gc_label],
    other_max = max(med_hME[celltype != gc_label], na.rm = TRUE),
    delta_GC  = gc - other_max,
    .groups   = "drop"
  ) |>
  dplyr::filter(module != "grey") |>
  dplyr::arrange(dplyr::desc(delta_GC))

save_tsv(spec, file.path(OUT_TAB, "delta_GC_projected_modules.tsv"))
message("Módulo con mayor ΔGC: ", spec$module[1])

# 5. Heatmap z-score (celltype × tratamiento) --------------------------------

score_prefix <- "modscore_"
scores_raw   <- GetModuleScores(query_obj, wgcna_name = wgcna_proj)
scores_raw   <- scores_raw[, setdiff(colnames(scores_raw), "grey"), drop = FALSE]
colnames(scores_raw) <- paste0(score_prefix, colnames(scores_raw))
query_obj <- AddMetaData(query_obj, metadata = scores_raw)

plot_df <- query_obj@meta.data |>
  tibble::rownames_to_column("cell") |>
  dplyr::select(cell, Sample, Treatment, all_of(meta_celltype),
                dplyr::starts_with(score_prefix)) |>
  tidyr::pivot_longer(dplyr::starts_with(score_prefix),
                      names_to = "module", values_to = "score") |>
  dplyr::mutate(module = sub(paste0("^", score_prefix), "", module)) |>
  dplyr::group_by(Sample, Treatment, .data[[meta_celltype]], module) |>
  dplyr::summarise(score_median = median(score, na.rm = TRUE),
                   .groups = "drop")

heat_df <- plot_df |>
  dplyr::group_by(Treatment, .data[[meta_celltype]], module) |>
  dplyr::summarise(mean_score = mean(score_median, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(group = paste(.data[[meta_celltype]], Treatment, sep = " | ")) |>
  dplyr::select(module, group, mean_score) |>
  tidyr::pivot_wider(names_from = group, values_from = mean_score)

mat <- heat_df |>
  tibble::column_to_rownames("module") |>
  as.matrix()

mat_z <- t(apply(mat, 1, function(x) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
}))

plot_z <- as.data.frame(mat_z) |>
  tibble::rownames_to_column("module") |>
  tidyr::pivot_longer(-module, names_to = "group", values_to = "zscore") |>
  tidyr::separate(group, into = c("cell_type", "treatment"), sep = " \\| ")

p_heat <- ggplot(plot_z, aes(x = treatment, y = module, fill = zscore)) +
  geom_tile() +
  facet_wrap(~cell_type, scales = "free_x") +
  scale_fill_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick3",
                       midpoint = 0, name = "Row z-score") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_plot_px(p_heat,
  file.path(OUT_FIG, paste0("heatmap_zscore_", ref_name, "_to_", query_name, ".png")),
  w = 3500, h = 2000
)

# Heatmap solo Guard Cells
plot_z_gc <- dplyr::filter(plot_z, trimws(cell_type) == gc_label)

p_heat_gc <- ggplot(plot_z_gc, aes(x = trimws(treatment), y = module, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick3",
                       midpoint = 0, name = "z-score\npor fila") +
  labs(x = "Tratamiento", y = "Módulo") +
  theme_bw(base_size = 10) + theme(panel.grid = element_blank())

save_plot_px(p_heat_gc,
  file.path(OUT_FIG, paste0("heatmap_zscore_GC_", ref_name, "_to_", query_name, ".png")),
  w = 1200, h = 2000
)

# 6. Violin del módulo con mayor ΔGC -----------------------------------------

mod_star <- spec$module[1]

if (!is.null(mod_star) && mod_star %in% colnames(MEs)) {
  query_obj[[mod_star]] <- MEs[, mod_star]

  dfv <- Seurat::FetchData(query_obj, vars = c(meta_celltype, mod_star)) |>
    dplyr::rename(celltype = 1, score = 2)

  ct_ord <- dfv |>
    dplyr::group_by(celltype) |>
    dplyr::summarise(med = median(score), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(med)) |>
    dplyr::pull(celltype)

  query_obj[[meta_celltype]] <- factor(query_obj[[meta_celltype]][, 1], levels = ct_ord)

  p_vln <- VlnPlot(query_obj, features = mod_star,
                   group.by = meta_celltype, pt.size = 0) +
    labs(title = paste0("Módulo proyectado: ", mod_star,
                        " (", ref_name, " → ", query_name, ")"),
         x = NULL, y = mod_star) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  save_plot_px(p_vln,
    file.path(OUT_FIG, paste0("violin_", mod_star, "_", ref_name, "_to_", query_name, ".png")),
    w = 3000, h = 1800
  )
}

# 7. Module preservation (Z-summary) -----------------------------------------

pres_name <- paste0(ref_name, "_to_", query_name)

query_obj <- ModulePreservation(
  query_obj,
  seurat_ref      = ref_obj,
  name            = pres_name,
  n_permutations  = 250,
  verbose         = 3,
  wgcna_name      = wgcna_proj,
  wgcna_name_ref  = wgcna_ref
)

mod_pres_Z   <- GetModulePreservation(query_obj, pres_name)$Z
mod_pres_obs <- GetModulePreservation(query_obj, pres_name)$obs

save_tsv(as.data.frame(mod_pres_Z),
  file.path(OUT_TAB, "module_preservation_Z.tsv"))
save_tsv(as.data.frame(mod_pres_obs),
  file.path(OUT_TAB, "module_preservation_obs.tsv"))

for (stat in c("summary", "rank", "all")) {
  plist <- PlotModulePreservation(query_obj, name = pres_name, statistics = stat)
  p_pres <- patchwork::wrap_plots(plist, ncol = min(length(plist), 6))
  ggsave(
    file.path(OUT_FIG, paste0("module_preservation_", stat, ".png")),
    p_pres,
    width  = if (stat == "all") 25 else 10,
    height = if (stat == "rank") 10 else 6,
    dpi = 300
  )
}

# 8. genesectR: solapamiento entre módulos REF y QUERY -----------------------

sets_ref <- get_module_sets(ref_obj,   wgcna_ref,  prefix = toupper(ref_name))
sets_qry <- get_module_sets(query_obj, wgcna_proj, prefix = toupper(query_name))

u_ref    <- GetWGCNAGenes(ref_obj,   wgcna_name = wgcna_ref)
u_qry    <- GetWGCNAGenes(query_obj, wgcna_name = wgcna_proj)
master   <- intersect(u_ref, u_qry)

sets_ref <- lapply(sets_ref, intersect, master)
sets_qry <- lapply(sets_qry, intersect, master)

gs <- genesectR::gs_import(c(sets_ref, sets_qry), master)
gs <- genesectR::gs_compute_matricies(gs, mc = FALSE)

png(file.path(OUT_FIG, "genesect_modules_ref_vs_query.png"),
    width = 5000, height = 4500, res = 300)
par(mar = c(2, 5, 10, 10))
genesectR::gs_multi_plot(gs, expand = names(sets_ref))
dev.off()

message("Script 12 completado.")

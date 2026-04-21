# 10_hdwgcna_network.R
# Construcción de la red de coexpresión con hdWGCNA:
#   fraction sweep → SetupForWGCNA → metacélulas (QC) →
#   SetDatExpr → soft-power → ConstructNetwork → eigengenes → kME → hub genes
#
# Entrada:  qs2/*_subclustered.qs2  (salida de 09_hdwgcna_subset.R)
# Salida:   qs2/*_with_metacells.qs2
#           qs2/*_hdWGCNA_object.qs2
#           hdWGCNA_TOM/<name>_TOM.rda
#           tables/power_table.txt, hub_df_*.txt, metacells_*.txt
#           figures/wgcna_fraction_sweep.png, soft_power_sweep_*.png, dendrogram_*.png
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

suppressPackageStartupMessages({
  library(hdWGCNA)
  library(WGCNA)
})
enableWGCNAThreads(nThreads = 56)

for (d in c("figures", "tables", "qs2", "hdWGCNA_TOM")) {
  dir.create(d, showWarnings = FALSE)
}

# 1. Cargar objeto ------------------------------------------------------------

#  Ajusta el nombre al archivo producido por 09_hdwgcna_subset.R.
#  coexpr_mode y ct deben coincidir con lo elegido en el paso anterior.

coexpr_mode <- "subsetting"   # "full" o "subsetting"
ct          <- "Guard Cells"  # tipo celular elegido en 09_

in_base <- paste(
  if (grepl("TB$", project)) "ten" else "nat",
  coexpr_mode, "subclustered",
  sep = "_"
)
seurat_obj <- qs_read(file.path("qs2", paste0(in_base, ".qs2")))
info(seurat_obj)

# 2. Nombre del análisis WGCNA -----------------------------------------------

ct_words  <- tolower(strsplit(ct, "[ _]+")[[1]])
wgcna_name <- paste(c(ct_words, "data", "v1"), collapse = "_")
message("wgcna_name: ", wgcna_name)

# 3. Fraction sweep ----------------------------------------------------------

group_var <- if ("Sample" %in% colnames(seurat_obj@meta.data)) "Sample" else "orig.ident"

df_frac <- fraction_sweep(
  seurat_obj,
  FRAC_GRID,
  group_by   = if (coexpr_mode == "full") "celltype" else group_var,
  tmp_prefix = "tmp_frac_"
)

p_frac <- ggplot(df_frac, aes(fraction, n_genes)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = c(MIN_GENES, MAX_GENES), linetype = 2) +
  scale_x_continuous(breaks = FRAC_GRID) +
  labs(title = "SetupForWGCNA: fraction vs n_genes")

save_plot_px(p_frac, file.path("figures", "wgcna_fraction_sweep.png"),
             w = 2200, h = 1400)

# 4. SetupForWGCNA -----------------------------------------------------------
#
#  ⚠ Elige la fracción que dé entre MIN_GENES y MAX_GENES genes.

frac_use <- switch(tech_mode, sc = 0.20, sn = 0.10)
message("Fracción elegida: ", frac_use,
        " → n_genes ≈ ", df_frac$n_genes[which.min(abs(df_frac$fraction - frac_use))])

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction    = frac_use,
  group.by    = if (coexpr_mode == "full") "celltype" else group_var,
  wgcna_name  = wgcna_name,
  assay       = assay_use
)
seurat_obj <- SetActiveWGCNA(seurat_obj, wgcna_name)

wgcna_genes <- GetWGCNAGenes(seurat_obj)
message("Genes WGCNA: ", length(wgcna_genes))

# Fracción de detección (diagnóstico)
mat_wg  <- GetAssayData(seurat_obj, assay = assay_use, layer = layer_use)[wgcna_genes, ]
det_frac <- Matrix::rowMeans(mat_wg > 0)
message("Detección — resumen:")
print(summary(det_frac))

# 5. Metacélulas (por Sample) -------------------------------------------------

message("\nCélulas por muestra:")
print(table(seurat_obj$Sample))

met_params <- if (coexpr_mode == "full") {
  switch(tech_mode,
    sc = list(k_met = 20, max_shared = 4,  min_cells = 100),
    sn = list(k_met = 40, max_shared = 12, min_cells = 200)
  )
} else {
  switch(tech_mode,
    sc = list(k_met = 12, max_shared = 2, min_cells = 25),
    sn = list(k_met = 15, max_shared = 2, min_cells = 25)
  )
}
list2env(met_params, envir = environment())

removed_groups <- character()

seurat_obj <- withCallingHandlers(
  MetacellsByGroups(
    seurat_obj     = seurat_obj,
    group.by       = c(if (coexpr_mode == "full") "celltype" else "annotation", "Sample"),
    ident.group    = "Sample",
    reduction      = reduction_use,
    dims           = DIMS_USE,
    target_metacells = 400,
    k              = k_met,
    max_shared     = max_shared,
    min_cells      = min_cells,
    layer          = layer_use,
    mode           = "average",
    verbose        = FALSE
  ),
  warning = function(w) {
    msg <- gsub("\n", " ", conditionMessage(w))
    pat <- "Removing the following groups that did not meet min_cells:\\s*(.*)$"
    if (grepl(pat, msg)) {
      grps <- trimws(strsplit(sub(pat, "\\1", msg), ",\\s*")[[1]])
      removed_groups <<- unique(c(removed_groups, grps))
      invokeRestart("muffleWarning")
    }
  }
)
message("Grupos eliminados: ", paste(removed_groups, collapse = ", "))

seurat_obj    <- NormalizeMetacells(seurat_obj)
metacell_obj  <- GetMetacellObject(seurat_obj)

message("\nMetacélulas por muestra:")
print(table(metacell_obj$Sample))

# 5a. QC metacélulas ---------------------------------------------------------

md_meta         <- metacell_obj[[]]
md_meta$n_cells <- lengths(strsplit(md_meta$cells_merged, ","))

sizes <- md_meta |>
  dplyr::group_by(Sample) |>
  dplyr::summarise(
    n_metacells   = dplyr::n(),
    min_cells     = min(n_cells),
    median_cells  = median(n_cells),
    max_cells     = max(n_cells),
    .groups = "drop"
  )
print(sizes)
save_csv(sizes, file.path("tables", "metacells_sizes.csv"))

# Cobertura de células únicas
cells_df <- md_meta |>
  tibble::rownames_to_column("metacell") |>
  dplyr::transmute(metacell, Sample,
                   cell = stringr::str_split(cells_merged, ",\\s*")) |>
  tidyr::unnest(cell) |>
  dplyr::mutate(cell = trimws(cell)) |>
  dplyr::filter(cell != "")

stopifnot(all(cells_df$cell %in% colnames(seurat_obj)))

tot_cells <- table(seurat_obj[["Sample"]][, 1])
cov <- cells_df |>
  dplyr::group_by(Sample) |>
  dplyr::reframe(
    unique_cells = dplyr::n_distinct(cell),
    total_cells  = as.integer(tot_cells[unique(Sample)]),
    coverage     = unique_cells / total_cells
  )
print(cov)
save_csv(cov, file.path("tables", "metacells_coverage.csv"))

# 5b. Visualizar metacélulas (opcional) -------------------------------------

seurat_obj <- seurat_obj |>
  ScaleMetacells(features = VariableFeatures(seurat_obj)) |>
  RunPCAMetacells(features = VariableFeatures(seurat_obj)) |>
  RunHarmonyMetacells(group.by.vars = "Sample") |>
  RunUMAPMetacells(
    reduction = "harmony", dims = DIMS_USE,
    seed.use = SEED, umap.method = "umap-learn", metric = "correlation"
  )

p1 <- DimPlotMetacells(seurat_obj, group.by = "Sample") +
  ggtitle("Metacells: Sample") + umap_theme()

save_plot_px(p1,
  file.path("figures", paste0("metacell_umap_", project, ".png")),
  w = 2500, h = 2000
)

# 6. SetDatExpr --------------------------------------------------------------

metacell_obj$WGCNA_group <- "all"
seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj)

kept_groups  <- setdiff(
  unique(paste0(seurat_obj$celltype, "#", seurat_obj$Sample)),
  removed_groups
)
valid_celltypes <- unique(sub("#.*$", "", kept_groups))

seurat_obj <- SetDatExpr(
  seurat_obj,
  group.by      = "WGCNA_group",
  group_name    = "all",
  use_metacells = TRUE,
  assay         = assay_use,
  layer         = layer_use
)

datExpr <- GetDatExpr(seurat_obj)
message("dim(datExpr) [muestras × genes]: ", paste(dim(datExpr), collapse = " × "))

mad_g <- apply(as.matrix(datExpr), 2, mad)
message("Genes con MAD=0: ", sum(mad_g == 0),
        " (", round(mean(mad_g == 0) * 100, 3), "%)")

# 7. Soft-power sweep --------------------------------------------------------

seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")
plot_list  <- PlotSoftPowers(seurat_obj)
print(patchwork::wrap_plots(plot_list, ncol = 2))

save_plot_px(
  patchwork::wrap_plots(plot_list, ncol = 2),
  file.path("figures", paste0("soft_power_sweep_", project, ".png"))
)

power_tbl <- GetPowerTable(seurat_obj)
save_csv(power_tbl, file.path("tables", paste0(project, "_power_table.csv")))
print(head(power_tbl, 10))

# Elegir el menor power con SFT.R.sq >= 0.8
cand       <- dplyr::filter(power_tbl, SFT.R.sq >= 0.8) |> dplyr::arrange(Power)
soft_power <- if (nrow(cand) > 0) cand$Power[1] else {
  message("⚠ SFT.R.sq nunca alcanza 0.8 — usando power=10 por defecto.")
  10
}
message("Soft power elegido: ", soft_power)

# 8. Construir red (TOM en disco) --------------------------------------------

tom_name <- paste0(
  paste(ct_words, collapse = ""),
  "_TOM_p", soft_power
)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name      = tom_name,
  tom_outdir    = "hdWGCNA_TOM",
  overwrite_tom = TRUE
)

# Dendrograma
png(file.path("figures",
              paste0("dendrogram_p", soft_power, "_", project, ".png")),
    width = 3200, height = 1800, res = 300)
PlotDendrogram(seurat_obj,
               main = paste0(ct, " hdWGCNA dendrogram"))
dev.off()

modules <- GetModules(seurat_obj)
message("\nMódulos:")
print(table(modules$module))

# 9. Module eigengenes + conectividad (kME) ----------------------------------

seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "Sample")

hMEs <- GetMEs(seurat_obj)
MEs  <- GetMEs(seurat_obj, harmonized = FALSE)
message("dim(MEs): ", paste(dim(MEs), collapse = " × "))

save_rds(hMEs, file.path("qs2", paste0("hmes_", project, ".rds")))
save_rds(MEs,  file.path("qs2", paste0("mes_",  project, ".rds")))

seurat_obj$WGCNA_group <- "all"
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by   = "WGCNA_group",
  group_name = "all"
)

# Renombrar módulos
new_mod_prefix <- paste(c(ct_words, "M"), collapse = "-")
seurat_obj <- ResetModuleNames(seurat_obj, new_name = new_mod_prefix)

# 10. Hub genes + kME plot ---------------------------------------------------

p_kme <- PlotKMEs(seurat_obj, n_hubs = 25, ncol = 5)
save_plot_px(p_kme,
  file.path("figures", paste0("kme_", project, ".png")),
  w = 5000
)

hub_df <- GetHubGenes(seurat_obj, n_hubs = 25)
hub_df <- dplyr::mutate(hub_df, symbol = to_sym(gene_name))
message("\nTop hub genes:")
print(head(hub_df, 20))
save_tsv(hub_df, file.path("tables", paste0("hub_df_", project, ".tsv")))

# 11. Guardar ----------------------------------------------------------------

out_base <- paste(
  if (grepl("TB$", project)) "ten" else "nat",
  coexpr_mode, "hdWGCNA_object",
  sep = "_"
)
save_obj(seurat_obj, file.path("qs2", out_base), ext = "qs2")
message("Objeto guardado: qs2/", out_base, ".qs2")
info(seurat_obj)

# 03_filter.R
# Filtrado opcional (modular): aplica los .rds que existan en el directorio.
#   - scds_dbls.rds          → elimina dobletes
#   - responsive_cells_ucell.rds → elimina células responsivas a digestión
#   - final_cells.rds + final_genes.rds → subset a células/genes finales
#
# Entrada:  act_postdbls.qs2  (o cualquier qs2/rds reciente del proyecto)
# Salida:   act_filtered.qs2
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

# 1. Cargar objeto -------------------------------------------------------------

# Prioriza el archivo post-dobletes; si no existe, cae sobre el pre-filtro
candidate_files <- c("act_postdbls.qs2",
                     paste0(tolower(sub("Leaf_", "", project)), "_prefilt.qs2"))

in_file <- candidate_files[file.exists(candidate_files)][1]
if (is.na(in_file)) stop("No se encontró ningún archivo de entrada.")

act <- qs_read(in_file)
message("Cargado desde: ", in_file)
info(act)

# 2. Aplicar filtros disponibles -----------------------------------------------

## 2.1 Dobletes ----------------------------------------------------------------
if (file.exists("scds_dbls.rds")) {
  doublets <- readRDS("scds_dbls.rds")
  act      <- subset(act, cells = setdiff(colnames(act), doublets))
  message("Dobletes removidos: ", length(doublets))
  info(act)
} else {
  message("scds_dbls.rds no encontrado — se omite filtro de dobletes.")
}

## 2.2 Células responsivas a digestión ----------------------------------------
if (file.exists("responsive_cells_ucell.rds")) {
  responsive_cells <- readRDS("responsive_cells_ucell.rds")
  act <- subset(act, cells = setdiff(colnames(act), responsive_cells))
  message("Células responsivas a digestión removidas: ", length(responsive_cells))
  info(act)
} else {
  message("responsive_cells_ucell.rds no encontrado — se omite filtro de digestión.")
}

## 2.3 Subset a células/genes finales -----------------------------------------
if (file.exists("final_cells.rds") && file.exists("final_genes.rds")) {
  final_cells <- readRDS("final_cells.rds")
  final_genes <- readRDS("final_genes.rds")

  act <- subset(
    act,
    cells    = intersect(colnames(act), final_cells),
    features = intersect(rownames(act), final_genes)
  )
  message("Subset a células y genes finales aplicado.")
  info(act)
} else {
  message("final_cells.rds / final_genes.rds no encontrados — se omite subset final.")
}

# 3. Guardar ------------------------------------------------------------------

qs_save(act, "act_filtered.qs2")
message("Objeto filtrado guardado en: act_filtered.qs2")
info(act)

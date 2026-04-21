# 02_doublets.R
# Detección de dobletes con scds (cxds + bcds + hybrid score).
# Entrada:  *_prefilt.qs2
# Salida:   scds_dbls.rds   (vector de barcodes de dobletes)
#           act_postdbls.qs2
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

# 1. Cargar objeto pre-filtrado ------------------------------------------------

in_file <- paste0(tolower(sub("Leaf_", "", project)), "_prefilt.qs2")
stopifnot(file.exists(in_file))

act <- qs_read(in_file)
info(act)

# 2. Convertir a SingleCellExperiment y estimar dobletes ----------------------

DefaultAssay(act) <- "RNA"
sce <- JoinLayers(act) |> as.SingleCellExperiment()

sce <- cxds_bcds_hybrid(sce, estNdbl = TRUE, verb = TRUE)

# Tabla de dobletes por muestra (útil para diagnóstico)
message("\nDobletes por muestra:")
print(table(colData(sce)$Sample, colData(sce)$hybrid_call))

# 3. Transferir scores al objeto Seurat y filtrar ------------------------------

act$hybrid_score <- colData(sce)$hybrid_score
act$hybrid_call  <- colData(sce)$hybrid_call   # TRUE = doblete

doublets <- colnames(act)[act$hybrid_call]
message("\nDobletes detectados: ", length(doublets))
message("Dobletes / total: ",
        round(100 * length(doublets) / ncol(act), 1), "%")

act <- subset(act, subset = !hybrid_call)
message("Células retenidas: ", ncol(act))

# 4. Guardar ------------------------------------------------------------------

saveRDS(doublets, "scds_dbls.rds")
qs_save(act, "act_postdbls.qs2")

info(act)
message("Dobletes guardados en: scds_dbls.rds")
message("Objeto guardado en:    act_postdbls.qs2")

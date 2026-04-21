# 14_hdwgcna_tfnetwork.R
# Red de regulación transcripcional con hdWGCNA:
#   • JASPAR 2024 (plantas) → motif matching en promotores
#   • ConstructTFNetwork (xgboost) + parche para xgboost v2/3
#   • AssignTFRegulons + RegulonScores
#   • Visualización de regulones y network plot
#
# Nota: requiere Bioconductor (TFBSTools, motifmatchr, BSgenome, ensembldb).
#
# Entrada:  qs2/*_hdWGCNA_object.qs2  (salida de 10_hdwgcna_network.R)
# Salida:   qs2/*_with_motifs.qs2, qs2/*_tfnet.qs2
#           figures/regulons.png, tables/hub_centric_*.tsv
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

suppressPackageStartupMessages({
  library(hdWGCNA)
  library(WGCNA)
  # Bioconductor
  library(TFBSTools)
  library(motifmatchr)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ensembldb)
  library(rtracklayer)
  library(RSQLite)
  library(JASPAR2024)
  library(BSgenome.Athaliana.TAIR.TAIR9)
  library(AnnotationDbi)
  library(org.At.tair.db)
  library(xgboost)
})
enableWGCNAThreads(nThreads = 56)

for (d in c("figures", "tables", "qs2", "RDS")) dir.create(d, showWarnings = FALSE)

# Parámetros ------------------------------------------------------------------

GTF_PATH       <- "/home/administrador/Documentos/projects/sortiz/databases/TAIR10_62/Arabidopsis_thaliana.TAIR10.62.gtf"
ENSDB_OUT      <- "EnsDb.Athaliana.TAIR10.62.sqlite"
N_THREADS      <- 56
REG_THRESH_PCT <- 0.90   # percentil de Gain para definir umbral de regulón
N_TFS_PER_GENE <- 6      # estrategia A: top TFs por gen

# 1. Cargar objeto ------------------------------------------------------------

coexpr_mode <- "subsetting"
ct          <- "Guard Cells"

in_base <- paste(
  if (grepl("TB$", project)) "ten" else "nat",
  coexpr_mode, "hdWGCNA_object",
  sep = "_"
)
gc_obj <- qs_read(file.path("qs2", paste0(in_base, ".qs2")))
wgcna_name <- gc_obj@misc$active_wgcna
info(gc_obj)

# 2. JASPAR 2024 PFMs (plantas, Arabidopsis thaliana taxon = 3702) ------------

jaspar_db <- JASPAR2024()
sq        <- DBI::dbConnect(RSQLite::SQLite(), TFBSTools::db(jaspar_db))
pfm_at    <- TFBSTools::getMatrixSet(
  x    = sq,
  opts = list(collection = "CORE", tax_group = "plants",
              species = 3702, all_versions = FALSE)
)
DBI::dbDisconnect(sq)
message("Motivos JASPAR cargados: ", length(pfm_at))

# 3. EnsDb desde GTF (offline) ------------------------------------------------

genome_at <- BSgenome.Athaliana.TAIR.TAIR9

if (!file.exists(ENSDB_OUT)) {
  message("Construyendo EnsDb desde GTF...")
  gtf_gr <- rtracklayer::import(GTF_PATH)

  seqmap <- c("1"="Chr1","2"="Chr2","3"="Chr3","4"="Chr4","5"="Chr5",
               "Mt"="ChrM","Pt"="ChrC")
  common <- intersect(names(seqmap), GenomeInfoDb::seqlevels(gtf_gr))
  if (length(common) > 0)
    gtf_gr <- GenomeInfoDb::renameSeqlevels(gtf_gr, seqmap[common])

  gtf_gr <- GenomeInfoDb::keepSeqlevels(
    gtf_gr,
    intersect(GenomeInfoDb::seqlevels(gtf_gr), GenomeInfoDb::seqlevels(genome_at)),
    pruning.mode = "coarse"
  )
  GenomeInfoDb::seqinfo(gtf_gr) <- GenomeInfoDb::seqinfo(genome_at)[
    GenomeInfoDb::seqlevels(gtf_gr)]

  ensDbFromGRanges(
    gtf_gr,
    organism      = "Arabidopsis thaliana",
    genomeVersion = "TAIR10",
    version       = 62,
    outfile       = ENSDB_OUT
  )
}

edb_at <- EnsDb(ENSDB_OUT)

# 4. Promotores (2 kb upstream) -----------------------------------------------

genes_gr <- ensembldb::genes(edb_at)
genes_gr <- GenomeInfoDb::keepSeqlevels(
  genes_gr,
  intersect(GenomeInfoDb::seqlevels(genes_gr), GenomeInfoDb::seqlevels(genome_at)),
  pruning.mode = "coarse"
)
GenomeInfoDb::seqinfo(genes_gr) <- GenomeInfoDb::seqinfo(genome_at)[
  GenomeInfoDb::seqlevels(genes_gr)]

genes_gr <- genes_gr[mcols(genes_gr)$gene_id %in% rownames(gc_obj)]
gene_ids <- mcols(genes_gr)$gene_id

prom_gr        <- GenomicRanges::promoters(genes_gr, upstream = 2000, downstream = 0)
prom_gr        <- GenomicRanges::trim(prom_gr)
names(prom_gr) <- gene_ids

# Remapear seqlevels si vienen como "1","2"…
lvl <- GenomeInfoDb::seqlevels(prom_gr)
if (any(lvl %in% c("1","2","3","4","5","Mt","Pt"))) {
  seqmap2 <- c("1"="Chr1","2"="Chr2","3"="Chr3","4"="Chr4","5"="Chr5",
                "Mt"="ChrM","Pt"="ChrC")
  com2 <- intersect(names(seqmap2), lvl)
  prom_gr <- GenomeInfoDb::renameSeqlevels(prom_gr, seqmap2[com2])
}
prom_gr <- GenomeInfoDb::keepSeqlevels(
  prom_gr,
  intersect(GenomeInfoDb::seqlevels(prom_gr), GenomeInfoDb::seqlevels(genome_at)),
  pruning.mode = "coarse"
)
GenomeInfoDb::seqinfo(prom_gr) <- GenomeInfoDb::seqinfo(genome_at)[
  GenomeInfoDb::seqlevels(prom_gr)]

stopifnot(length(prom_gr) > 0)
message("Promotores: ", length(prom_gr))

# 5. Motif matching ----------------------------------------------------------

mm       <- motifmatchr::matchMotifs(pfm_at, prom_gr, genome_at, out = "matches")
m_mat    <- Matrix::Matrix(motifmatchr::motifMatches(mm), sparse = TRUE)
rownames(m_mat) <- names(prom_gr)

motif_ids <- vapply(pfm_at, TFBSTools::ID,   character(1))
tf_symbol <- vapply(pfm_at, TFBSTools::name, character(1))
colnames(m_mat) <- motif_ids

# Mapear motivo → TF (TAIR ID)
tair_id <- AnnotationDbi::mapIds(
  org.At.tair.db,
  keys     = tf_symbol,
  keytype  = "SYMBOL",
  column   = "TAIR",
  multiVals = "first"
)

motif_df <- tibble::tibble(
  motif_ID  = motif_ids,
  gene_name = unname(tair_id)
) |>
  dplyr::filter(!is.na(gene_name), gene_name %in% rownames(gc_obj)) |>
  dplyr::distinct(motif_ID, .keep_all = TRUE)

m_mat2 <- m_mat[, motif_df$motif_ID, drop = FALSE]
message("Motivos válidos (con TF en el objeto): ", ncol(m_mat2))

motif_targets <- lapply(seq_len(ncol(m_mat2)),
  function(j) rownames(m_mat2)[m_mat2[, j] > 0])
names(motif_targets) <- colnames(m_mat2)

# 6. Guardar motivos en gc_obj ------------------------------------------------

gc_obj <- SetMotifMatrix(gc_obj, m_mat2)
gc_obj <- SetMotifs(gc_obj, motif_df)
gc_obj <- SetMotifTargets(gc_obj, motif_targets)
gc_obj <- SetPFMList(gc_obj, pfm_at)
stopifnot(nrow(GetMotifMatrix(gc_obj)) > 0)

save_obj(gc_obj,
  file.path("qs2", paste0(in_base, "_with_motifs")), ext = "qs2")

# 7. ConstructTFNetwork (xgboost) — con parche para xgboost v2/3 --------------

model_params <- list(
  objective = "reg:squarederror",
  max_depth = 1,
  eta       = 0.1,
  nthread   = N_THREADS,
  alpha     = 0.5
)

# Parche: convierte modelos guardados por callback a xgb.Booster
as_booster <- function(mod) {
  if (inherits(mod, "xgb.Booster")) return(mod)
  if (is.raw(mod))  return(xgboost::xgb.load.raw(mod))
  if (is.list(mod) && !is.null(mod$raw) && is.raw(mod$raw))
    return(xgboost::xgb.load.raw(mod$raw))
  if (is.list(mod) && !is.null(mod$model)) return(as_booster(mod$model))
  stop("No se puede convertir a xgb.Booster (tipo: ",
       paste(class(mod), collapse=", "), ")")
}

construct_tf_network_safe <- function(seurat_obj, model_params,
                                       wgcna_name = NULL, nfold = 5) {
  if (is.null(wgcna_name)) wgcna_name <- seurat_obj@misc$active_wgcna
  hdWGCNA::CheckWGCNAName(seurat_obj, wgcna_name)

  motif_matrix <- GetMotifMatrix(seurat_obj)
  motif_df_loc <- GetMotifs(seurat_obj)

  cb      <- list(xgboost::xgb.cb.cv.predict(save_models = TRUE))
  datExpr <- as.matrix(GetDatExpr(seurat_obj, wgcna_name = wgcna_name))
  genes_use <- intersect(colnames(datExpr), rownames(motif_matrix))

  importance_df <- data.frame()
  eval_df       <- data.frame()
  pb <- utils::txtProgressBar(min = 0, max = length(genes_use), style = 3)

  for (i in seq_along(genes_use)) {
    cur_gene  <- genes_use[i]
    cur_motifs <- colnames(motif_matrix)[which(motif_matrix[cur_gene, ] != 0)]
    cur_tfs   <- unique(motif_df_loc$gene_name[motif_df_loc$motif_ID %in% cur_motifs])
    cur_tfs   <- setdiff(cur_tfs[cur_tfs %in% genes_use], cur_gene)
    if (length(cur_tfs) < 2) { utils::setTxtProgressBar(pb, i); next }

    x_vars <- datExpr[, cur_tfs, drop = FALSE]
    y_var  <- as.numeric(datExpr[, cur_gene])
    if (all(y_var == 0)) { utils::setTxtProgressBar(pb, i); next }

    tf_cor <- as.numeric(cor(x = as.matrix(x_vars), y = y_var))
    names(tf_cor) <- cur_tfs

    xgb_cv <- xgboost::xgb.cv(
      params    = model_params,
      data      = xgboost::xgb.DMatrix(x_vars, label = y_var),
      nrounds   = 100,
      nfold     = nfold,
      showsd    = FALSE,
      callbacks = cb,
      verbose   = 0
    )

    models <- xgb_cv$cv_predict$models %||% xgb_cv$models
    if (is.null(models)) { utils::setTxtProgressBar(pb, i); next }

    imp <- Reduce("+", lapply(seq_len(nfold), function(f) {
      cur <- xgboost::xgb.importance(colnames(x_vars), as_booster(models[[f]]))
      ix  <- match(colnames(x_vars), as.character(cur$Feature))
      mat <- as.matrix(cur[ix, -1, drop = FALSE])
      mat[is.na(mat)] <- 0
      mat
    })) / nfold

    imp_df <- as.data.frame(imp)
    imp_df$tf   <- colnames(x_vars)
    imp_df$gene <- cur_gene
    imp_df$Cor  <- as.numeric(tf_cor)

    importance_df <- rbind(importance_df,
                            dplyr::select(imp_df, tf, gene, Gain, Cover, Frequency, Cor))
    eval_df       <- rbind(eval_df,
                            cbind(as.data.frame(xgb_cv$evaluation_log), variable = cur_gene))
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  seurat_obj <- SetTFNetwork(seurat_obj, importance_df, wgcna_name = wgcna_name)
  seurat_obj <- SetTFEval(seurat_obj,    eval_df,        wgcna_name = wgcna_name)
  seurat_obj
}

# Operador %||% (null-coalescing)
`%||%` <- function(a, b) if (!is.null(a)) a else b

gc_obj <- construct_tf_network_safe(gc_obj, model_params,
                                     wgcna_name = wgcna_name)

results    <- GetTFNetwork(gc_obj)
reg_thresh <- as.numeric(quantile(results$Gain, REG_THRESH_PCT, na.rm = TRUE))
message("Umbral de Gain (P", round(REG_THRESH_PCT*100), "): ", round(reg_thresh, 5))

save_obj(gc_obj, file.path("qs2", paste0(in_base, "_tfnet")), ext = "qs2")

# 8. Definir regulones --------------------------------------------------------

gc_obj <- AssignTFRegulons(gc_obj,
  strategy   = "A",
  reg_thresh = reg_thresh,
  n_tfs      = N_TFS_PER_GENE
)

reg <- GetTFRegulons(gc_obj)

# Targets de NLP7
targets_nlp7 <- reg |>
  dplyr::filter(tf == nlp7_gene) |>
  dplyr::pull(gene)

message("Targets de NLP7 (", to_sym(nlp7_gene), "): ", length(targets_nlp7))
message("Símbolos: ", paste(to_sym(targets_nlp7), collapse = ", "))

# Plot de regulón NLP7
p_reg <- RegulonBarPlot(gc_obj, selected_tf = nlp7_gene)
p_reg$data$gene_lab <- to_sym_agi(p_reg$data$gene)

# Reemplazar mapping de label (capa geom_text)
i_txt <- which(vapply(p_reg$layers, \(l) inherits(l$geom, "GeomText"), logical(1)))[1]
if (!is.na(i_txt)) {
  p_reg$layers[[i_txt]]$mapping$label <- rlang::quo(gene_lab)
}
p_reg <- p_reg + ggtitle(paste0(to_sym_agi(nlp7_gene), " — predicted targets"))

save_plot_px(p_reg, file.path("figures", "regulon_NLP7.png"),
             w = 5000, h = 5000)

# 9. Regulon scores -----------------------------------------------------------

gc_obj <- RegulonScores(gc_obj, target_type = "positive", ncores = 8)
gc_obj <- RegulonScores(gc_obj, target_type = "negative",
                         cor_thresh = -0.05, ncores = 8)

pos_scores <- GetRegulonScores(gc_obj, target_type = "positive")
neg_scores <- GetRegulonScores(gc_obj, target_type = "negative")

message("Regulon scores positivos: ", ncol(pos_scores), " TFs")
message("Regulon scores negativos: ", ncol(neg_scores), " TFs")

# Network plot (si TFNetworkPlot está disponible)
tf_interest <- to_agi("NLP7")
if (!is.null(tf_interest) && tf_interest %in% reg$tf) {
  p_net <- TFNetworkPlot(gc_obj,
                          selected_tfs = tf_interest,
                          target_type  = "positive",
                          depth        = 1)
  save_plot_px(p_net, file.path("figures", "tfnetwork_NLP7.png"))
}

# 10. Guardar objeto final ----------------------------------------------------

save_obj(gc_obj, file.path("qs2", paste0(in_base, "_tfnet_final")), ext = "qs2")
message("Script 14 completado.")
writeLines(capture.output(sessionInfo()),
           file.path("tables", "sessionInfo_tfnetwork.txt"))

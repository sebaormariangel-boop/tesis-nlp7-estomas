# 00_setup.R
# Configuración global: librerías, rutas, marcadores y parámetros de QC.
# Este archivo debe ser sourced al inicio de cada script del pipeline.
# ─────────────────────────────────────────────────────────────────────────────

# 0. Entorno -------------------------------------------------------------------

rm(list = ls())
gc()
set.seed(123)
SEED <- 123L   # constante usada en FindClusters / RunUMAP / permutaciones

# 1. Python / UMAP -------------------------------------------------------------

Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/envs/seurat/bin/python")
library(reticulate)
py_config()
stopifnot(py_module_available("umap"))

# 2. Librerías -----------------------------------------------------------------

suppressPackageStartupMessages({
  # Datos y manipulación
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(readxl)
  library(readr)
  library(tidytext)        # reorder_within

  # scRNA-seq
  library(Seurat)
  library(SeuratObject)
  library(sctransform)
  library(harmony)
  library(clustree)
  library(UCell)
  library(scds)
  library(SingleCellExperiment)
  library(hdWGCNA)
  library(WGCNA)
  library(presto)          # FindAllMarkers rápido

  # Expresión diferencial / normalización
  library(edgeR)

  # Enriquecimiento funcional
  library(clusterProfiler)
  library(org.At.tair.db)
  library(enrichplot)

  # Visualización
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(patchwork)
  library(ragg)            # backend rápido para ggsave
  library(RColorBrewer)
  library(gridExtra)
  library(magick)          # combinar imágenes / exportar PDF

  # Utilidades
  library(Matrix)
  library(digest)
  library(qs2)
  library(cli)
  library(crayon)
  library(mclust)
  library(genesectR)
})

theme_set(theme_cowplot())

# 2b. Librerías adicionales — hdWGCNA + análisis de red ----------------------
#     Se cargan solo cuando se van a correr los scripts 09–14.
#     Descomenta si los scripts anteriores ya han sido ejecutados.

# suppressPackageStartupMessages({
#   library(hdWGCNA)
#   library(WGCNA)
#   library(fgsea)
#   # Motif / TF network (requiere Bioconductor)
#   library(TFBSTools)
#   library(motifmatchr)
#   library(GenomicRanges)
#   library(GenomeInfoDb)
#   library(ensembldb)
#   library(rtracklayer)
#   library(RSQLite)
#   library(JASPAR2024)
#   library(BSgenome.Athaliana.TAIR.TAIR9)
#   library(AnnotationDbi)
#   library(xgboost)
# })
#
# # Hilos WGCNA (ajustar según servidor)
# WGCNA::enableWGCNAThreads(nThreads = 56)

# 3. Rutas del proyecto --------------------------------------------------------

proj_dir <- "/home/administrador/Documentos/projects/sortiz/2do_semestre/singlecell"
raw_path <- file.path(proj_dir, "final_samples")

# Genes de protoplasteo (para score de digestión)
pp.genes <- scan(
  file.path(raw_path, "protoplasting_leaves.txt"),
  what  = "character",
  quiet = TRUE
)

# 4. Mapeo gen ↔ símbolo (Arabidopsis) ----------------------------------------

map_df <- readr::read_tsv(
  "~/atsymbol_nate.txt",
  col_types = cols(),
  show_col_types = FALSE
) |>
  dplyr::filter(!is.na(gene), !is.na(symbol)) |>
  dplyr::distinct(gene, .keep_all = TRUE)

gene2sym <- setNames(map_df$symbol, map_df$gene)
sym2gene <- setNames(map_df$gene,   map_df$symbol)

#' Convierte AGI a símbolo (devuelve AGI si no hay símbolo)
to_sym <- function(x) {
  y <- unname(gene2sym[x])
  ifelse(is.na(y), x, y)
}

#' Convierte símbolo a AGI (devuelve el input si no se encuentra)
to_agi <- function(x) {
  y <- unname(sym2gene[x])
  ifelse(is.na(y), x, y)
}

# 5. Selección de tecnología / dataset ----------------------------------------
#
#  Asigna directamente:
#    tech_mode <- "sc"   # Tenorio Berrío et al. (scRNA-seq, BD Rhapsody)
#    tech_mode <- "sn"   # Illouz-Eliaz et al.   (snRNA-seq, 10x)
#
#  Para selección interactiva, comenta la línea siguiente y descomenta
#  el bloque de menu().

tech_mode <- "sc"

# ## Alternativa interactiva:
# opciones  <- c("Single-cell | Tenorio Berrío et al.",
#                "Single-nucleus | Illouz-Eliaz et al.")
# choice    <- menu(opciones, title = "Seleccionar tecnología")
# if (choice == 0) stop("Selección cancelada.")
# tech_mode <- c("sc", "sn")[choice]

stopifnot(tech_mode %in% c("sc", "sn"))

# 6. Parámetros dependientes de tecnología ------------------------------------

samples <- switch(
  tech_mode,
  sc = c("MD_F_rep1", "MD_F_rep2", "WW_F_rep1", "WW_F_rep2"),
  sn = c("DR_T0_rep1", "DR_T0_rep2", "WW_T0_rep1", "WW_T0_rep2")
)

qc_cfg <- switch(
  tech_mode,
  sc = list(min_UMI = 1250L, min_features = 1000L, max_mt = 5,  max_chl = 40),
  sn = list(min_UMI =  300L, min_features =  300L, max_mt = 1,  max_chl = 40)
)

sample_path <- switch(
  tech_mode,
  sc = file.path(raw_path, "tenorio"),
  sn = file.path(raw_path, "natanella")
)

project <- switch(tech_mode, sc = "Leaf_MD_WW_TB", sn = "Leaf_DR_WW_IE")

map_repeat <- switch(
  tech_mode,
  sc = c(MD_F_rep1 = "R1", MD_F_rep2 = "R2", WW_F_rep1 = "R1", WW_F_rep2 = "R2"),
  sn = c(DR_T0_rep1 = "R1", DR_T0_rep2 = "R2", WW_T0_rep1 = "R1", WW_T0_rep2 = "R2")
)

map_treat <- switch(
  tech_mode,
  sc = c(MD_F_rep1 = "MD", MD_F_rep2 = "MD", WW_F_rep1 = "WW", WW_F_rep2 = "WW"),
  sn = c(DR_T0_rep1 = "DR", DR_T0_rep2 = "DR", WW_T0_rep1 = "WW", WW_T0_rep2 = "WW")
)

# 7. Marcadores de referencia (DotPlot canónico) ------------------------------

markers_ref <- list(
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

abbr_ref <- c(
  Epidermis="E", Vasculature="V", Mesophyll="M", Bundle_Sheath="BS",
  Cambium="Cam", Companion_Cells="CC", Sieve_Element="SE",
  Phloem_Parenchyma="PP", Xylem="Xyl", Hydathode="Hyd",
  Myrosin_Idioblasts="MI", Guard_Cells="GC", Pavement_Cells="Pav"
)

# `markers` = alias abreviado de markers_ref (usado en dots_markers por defecto)
markers <- markers_ref
names(markers) <- ifelse(
  names(markers) %in% names(abbr_ref),
  abbr_ref[names(markers)],
  names(markers)
)

cap_text <- stringr::str_wrap(
  paste(names(abbr_ref), abbr_ref, sep = " = ", collapse = "; "),
  width = 90
)

# 8. Listas detalladas de marcadores (UCell / DotPlots exploratorios) ---------

## Ciclo celular
cell_cycle <- list(
  S_Lic                    = c("AT2G29680", "AT2G07690"),          # CDC6A, MCM5
  S_Rep                    = c("AT1G07370", "AT2G29570"),          # PCNA1, PCNA2
  G2M                      = c("AT4G37490", "AT3G54180",
                                "AT4G33270", "AT4G33260"),         # CYCB1;1, CDKB1;1, CDC20
  Cytokinesis_Phragmoplast = c("AT1G08560", "AT1G18370",
                                "AT3G17360", "AT3G19050"),         # KNOLLE, HINKEL, POK1/2
  Chromatin_H2A            = c("AT5G54640", "AT4G27230",
                                "AT5G59870", "AT1G51060",
                                "AT3G20670"),
  DNA_Dmg_Resp             = c("AT1G08880"),                       # H2AXA (γ-H2AX)
  Chromatin_HMG            = c("AT5G23420", "AT4G23800"),
  Epi_Maint                = c("AT5G49160")                        # MET1
)

## Epidermis
epidermis_markers <- list(
  Protoderm_L1 = c("AT2G42840", "AT4G21750"),
  Cuticle      = c("AT2G26250", "AT1G27950", "AT1G51500", "AT1G68530"),
  PC           = c("AT1G04040", "AT5G28630", "AT2G32990", "AT5G63180"),
  Adaxial_PC   = c("AT5G62210"),
  Patt_Trich   = c("AT2G30424"),
  Trich        = c("AT5G40330"),
  Flavonoid_Ad = c("AT5G13930"),
  Leaf_Abaxial = c("AT2G26580", "AT4G00180", "AT5G44430", "AT5G16560"),
  Exp          = c("AT2G37640", "AT3G29030"),
  Aux_Ad       = c("AT1G04250"),
  Phototrop    = c("AT1G18810")
)

## Células de guarda
gc_markers <- list(
  Mature_GC        = c("AT1G04800", "AT1G08810", "AT1G12480", "AT4G33950",
                        "AT1G22690", "AT4G14480", "AT1G79840", "AT3G16370"),
  Stomatal_Lineage = c("AT5G53210", "AT3G06120", "AT3G24140",
                        "AT1G80080", "AT1G34245")
)

## Vasculatura
vasc_markers <- list(
  BS       = c("AT1G25230", "AT4G13770", "AT1G77990", "AT2G43150", "AT1G12010"),
  Hyd      = c("AT1G22900", "AT1G62510", "AT3G16660", "AT3G16670"),
  MI       = c("AT5G26000", "AT5G25980", "AT3G16400", "AT2G47260"),
  Vasc     = c("AT1G64700", "AT5G03610", "AT1G80520", "AT3G14990"),
  Phl      = c("AT1G22150", "AT3G19710", "AT5G09220", "AT5G23010",
               "AT5G02380", "AT1G05760", "AT5G04890"),
  CC       = c("AT4G19840", "AT5G18600", "AT1G64370", "AT1G22710",
               "AT1G79430", "AT1G52190", "AT4G24450", "AT5G15950",
               "AT5G57350", "AT5G63850"),
  PP       = c("AT3G11930", "AT5G24800", "AT3G48740"),
  PP_CC    = c("AT5G23660", "AT5G54660"),
  SE       = c("AT3G01670", "AT3G01680"),
  Xyl      = c("AT3G10080", "AT5G60490", "AT1G67560", "AT4G35350",
               "AT4G32880", "AT5G19530", "AT1G20850", "AT5G62380"),
  Procam   = c("AT3G25710", "AT1G46480", "AT2G27230", "AT1G19850"),
  Cam      = c("AT2G39700", "AT5G61480"),
  Protophl = c("AT3G09070", "AT1G31880", "AT2G37590",
               "AT5G02460", "AT3G52490", "AT4G29920", "AT5G57130")
)

## Mesófilo
mes_markers <- list(
  Mes      = c("AT4G26530", "AT2G34430", "AT4G12970",
               "AT1G29910", "AT1G55670", "AT2G06520",
               "AT5G17170", "AT2G05100"),
  Palisade = c("AT5G05690", "AT1G70410", "AT4G23060"),
  Spongy   = c("AT4G23600")
)

# 9. Directorio de resultados --------------------------------------------------

results_dir <- switch(
  tech_mode,
  sc = file.path(proj_dir, "results/tenorio"),
  sn = file.path(proj_dir, "results/natanella")
)

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
setwd(results_dir)

# 10. Directorios de figuras y objetos ----------------------------------------

fig_dir <- file.path("Documentos", "figures")
qs2_dir <- file.path("Documentos", "qs2")
for (d in c(fig_dir, qs2_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# 11. Reporte ------------------------------------------------------------------

# 10b. Constantes para hdWGCNA (scripts 09–14) --------------------------------

## Reducción y capas
assay_use     <- "RNA"
layer_use     <- "data"       # log-normalizado (Seurat v5)
reduction_use <- "harmony"
DIMS_USE      <- 1:32         # dimensiones Harmony para clustering del subset

## Fraction sweep
FRAC_GRID <- c(0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15)
MIN_GENES <- 4000
MAX_GENES <- 15000

## Modo interactivo (TRUE = usa readline/menu; FALSE = asigna ct directamente)
interactive_mode <- TRUE

## Rutas de datos externos (ajustar según servidor)
bulk_leaf_dir   <- "/home/sortiz/Documentos/tables_DROUGHT"
bulk_stoma_dir  <- "/home/sortiz/Documentos/tables_STOMATA"
nlp7_gene       <- "AT4G24020"  # NLP7 (AT4G24020)
myb60_gene      <- "AT1G08810"  # MYB60

message(green("tech_mode:     "), bold(tech_mode))
message(green("Muestras:      "), paste(samples, collapse = ", "))
message(
  green("QC:            "),
  sprintf("min_UMI=%d  min_features=%d  max_mt=%.0f%%  max_chl=%.0f%%",
          qc_cfg$min_UMI, qc_cfg$min_features, qc_cfg$max_mt, qc_cfg$max_chl)
)
message(green("Resultados en: "), results_dir)
message(green("Proyecto:      "), project)

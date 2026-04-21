# 06_digestion.R
# Puntaje de digestión (protoplasteo) con UCell y modelo nulo por permutaciones.
# Identifica células significativamente afectadas por el protocolo de digestión.
#
# Entrada:  act_integrated.qs2
# Salida:   act_digest_scored.qs2
#           responsive_cells_ucell.rds
#           final_cells.rds / final_genes.rds
#           digest_summary_by_sample.csv
#           vln_digresp_ucell.png
#           dig_resp_umap_*.png
#           sig_digest_barplot.png
#           sig_digest_heatmap.png
#           digestion_supl.png
# ─────────────────────────────────────────────────────────────────────────────

source("00_setup.R")
source("helpers.R")

N_PERM   <- 500   # permutaciones para modelo nulo
N_CORES  <- 56    # núcleos para UCell

# 1. Cargar objeto -------------------------------------------------------------

act <- qs_read("act_integrated.qs2")
act <- JoinLayers(act)

DefaultAssay(act) <- "RNA"
if (!"data" %in% Layers(act[["RNA"]])) act <- NormalizeData(act, assay = "RNA")
info(act)

# 2. Calcular score UCell real -------------------------------------------------

m      <- median(act$nFeature_RNA)
q75    <- quantile(act$nFeature_RNA, 0.75)
maxRank <- round(min(max(m, 1200), q75, 3000), -2)
message("maxRank UCell: ", maxRank)

ranks2 <- UCell::StoreRankings_UCell(
  LayerData(act, assay = "RNA", layer = "data"),
  maxRank = maxRank,
  ncores  = N_CORES
)

s2 <- UCell::ScoreSignatures_UCell(
  precalc.ranks = ranks2,
  features      = list(Digestion = pp.genes),
  ncores        = N_CORES
)

score_real <- as.numeric(s2[, 1])

# 3. Modelo nulo por permutaciones (bin de detección × expresión) --------------

mat  <- LayerData(act, assay = "RNA", layer = "data")
univ <- intersect(rownames(ranks2), rownames(mat))
firma <- intersect(pp.genes, univ)

det   <- Matrix::rowMeans(mat[univ, , drop = FALSE] > 0)
mu    <- Matrix::rowMeans(mat[univ, , drop = FALSE])

d_bin <- cut(det, breaks = unique(quantile(det, seq(0, 1, 0.1), na.rm = TRUE)),
             include.lowest = TRUE, labels = FALSE)
m_bin <- cut(mu,  breaks = unique(quantile(mu,  seq(0, 1, 0.1), na.rm = TRUE)),
             include.lowest = TRUE, labels = FALSE)

bin2d <- paste(d_bin, m_bin, sep = "_")
names(bin2d) <- univ

firma_bins <- bin2d[firma]

set.seed(123)
perm_sets <- replicate(N_PERM, {
  unlist(lapply(firma_bins, function(b) {
    pool <- names(bin2d)[bin2d == b]
    sample(pool, 1)
  }), use.names = FALSE)
}, simplify = FALSE)
names(perm_sets) <- sprintf("null_%04d", seq_len(N_PERM))

s_null <- UCell::ScoreSignatures_UCell(
  precalc.ranks = ranks2,
  features      = perm_sets,
  ncores        = N_CORES
)

null_mat <- as.matrix(s_null)
ge_counts <- rowSums(null_mat >= score_real)
p_emp     <- (ge_counts + 1) / (N_PERM + 1)
fdr       <- p.adjust(p_emp, method = "BH")

null_mu <- rowMeans(null_mat)
null_sd <- apply(null_mat, 1, sd)

act$Digest_score        <- score_real
act$Digest_FDR          <- fdr
act$Digest_resp         <- fdr < 0.05
act$Digest_resp_strict  <- (fdr < 0.05) & (score_real > null_mu + 2 * null_sd)
act$Digest_mlog10fdr    <- pmin(-log10(act$Digest_FDR), 50)

# 4. Resumen por muestra -------------------------------------------------------

meta_sum <- act@meta.data |>
  mutate(sample = Sample, cluster = Idents(act)) |>
  group_by(sample) |>
  summarise(
    n_cells       = n(),
    n_resp        = sum(Digest_resp),
    n_resp_strict = sum(Digest_resp_strict),
    pct_resp      = 100 * n_resp / n_cells,
    .groups       = "drop"
  )

print(meta_sum)
write.csv(meta_sum, "digest_summary_by_sample.csv", row.names = FALSE)

# 5. Violin plot: score por muestra x estado ----------------------------------

levs <- act@meta.data |>
  dplyr::distinct(Sample, .data[["Digest_resp"]]) |>
  mutate(
    sample_resp = paste0(Sample, "__", ifelse(Digest_resp, "resp", "unresp")),
    estado      = ifelse(Digest_resp, 1, 0)
  ) |>
  arrange(Sample, estado) |>
  pull(sample_resp)

act$sample_resp <- factor(
  paste0(act$Sample, "__", ifelse(act$Digest_resp, "resp", "unresp")),
  levels = levs
)

cols_resp <- setNames(
  ifelse(grepl("__resp$", levs), "#D55E00", "#7F7F7F"),
  levs
)

p_vln <- VlnPlot(
  act, features = "Digest_score",
  group.by = "sample_resp", pt.size = 0, cols = cols_resp
) +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               fill = "white", color = "black", linewidth = 0.25) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.title.x  = element_blank(),
    axis.title.y  = element_blank(),
    legend.position = "none"
  )

save_plot_px(p_vln, "vln_digresp_ucell.png", width = 2400, height = 1800)

# Violín facetado (más limpio)
plot_df <- FetchData(act, vars = c("Digest_score", "Sample", "Digest_resp")) |>
  mutate(estado = ifelse(Digest_resp, "resp", "unresp"))

p_vln2 <- ggplot(plot_df, aes(x = estado, y = Digest_score, fill = estado)) +
  geom_violin(scale = "width", trim = FALSE, color = "grey20", linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  facet_wrap(~Sample, nrow = 1) +
  scale_fill_manual(values = c(unresp = "#B0B0B0", resp = "#2A6F97")) +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        axis.title.x = element_blank(), legend.position = "none") +
  labs(y = "Puntaje de digestión")

# 6. UMAPs de digestión -------------------------------------------------------

# UMAP discreto (resp / no resp)
p_umap_resp <- DimPlot(
  act, reduction = "umap.harmony", group.by = "Digest_resp"
) + ggtitle("") + umap_theme()

save_plot_px(p_umap_resp, "dig_resp_umap.png", width = 3000, height = 1800)

# UMAP continuo por -log10(FDR)
emb <- as.data.frame(Embeddings(act, "umap.harmony"))
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
emb$mlog10q <- act$Digest_mlog10fdr
emb$sig     <- act$Digest_resp

p_sig_umap <- ggplot(emb, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = mlog10q, alpha = sig), size = 0.5) +
  scale_alpha_manual(values = c(`FALSE` = 0.15, `TRUE` = 1)) +
  scale_color_viridis_c(option = "viridis") +
  coord_equal() +
  theme_void() +
  labs(color = "-log10(FDR)", alpha = "FDR ≤ 0.05")

save_plot_px(p_sig_umap, "dig_resp_umap_fdr.png", width = 3000, height = 1800)

# 7. Barras: % células responsivas por clúster + test Fisher ------------------

cl   <- Idents(act)
flag <- as.logical(act$Digest_resp)
p0   <- mean(flag, na.rm = TRUE)

bar_df <- tibble(cluster = cl, flag = flag) |>
  group_by(cluster) |>
  summarise(n = n(), k = sum(flag, na.rm = TRUE),
            pct = 100 * k / n, .groups = "drop") |>
  mutate(cluster = factor(cluster, levels = levels(cl)))

efecto_min_pp <- 0.05
OR_min        <- 1.25

res_fisher <- do.call(rbind, lapply(levels(cl), function(k) {
  idx <- which(cl == k)
  a <- sum(flag[idx], na.rm = TRUE);  b <- length(idx) - a
  j <- setdiff(seq_along(flag), idx)
  cc <- sum(flag[j], na.rm = TRUE);   d <- length(j) - cc
  ft <- fisher.test(matrix(c(a, b, cc, d), 2, 2), alternative = "greater")
  data.frame(
    cluster   = k,
    prop      = a / (a + b),
    prop_rest = cc / (cc + d),
    OR        = (a / b) / (cc / d),
    p_fisher  = ft$p.value
  )
}))
res_fisher$q_fisher <- p.adjust(res_fisher$p_fisher, "BH")
res_fisher$enriched <- (res_fisher$q_fisher < 0.05) &
  (res_fisher$OR >= OR_min) &
  ((res_fisher$prop - res_fisher$prop_rest) >= efecto_min_pp)

bar_plot_df <- bar_df |>
  left_join(res_fisher[, c("cluster", "q_fisher", "OR", "enriched")],
            by = "cluster") |>
  mutate(cluster = fct_relevel(cluster, levels(cl)))

p_bar <- ggplot(bar_plot_df, aes(x = cluster, y = pct, fill = enriched)) +
  geom_col(width = 0.75, color = "grey25", linewidth = 0.2) +
  geom_text(aes(label = ifelse(enriched, "*", "")), vjust = -0.4, size = 5) +
  scale_fill_manual(values = c(`TRUE` = "#2A6F97", `FALSE` = "#B0B0B0"),
                    guide = "none") +
  scale_y_continuous(name = "% células positivas",
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Clúster") +
  theme_classic(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)))

save_plot_px(p_bar, "sig_digest_barplot.png", width = 2700, height = 1800)

# 8. Heatmap: % células responsivas por muestra × clúster -------------------

tab     <- with(act@meta.data, table(Sample, Idents(act), Digest_resp))
tot_mat <- apply(tab, c(1, 2), sum)
pct_mat <- 100 * (tab[, , "TRUE"] / tot_mat)
pct_mat[is.nan(pct_mat)] <- NA_real_
pct_mat <- as.matrix(pct_mat)[, levels(Idents(act)), drop = FALSE]

heat_df <- as.data.frame.matrix(pct_mat) |>
  rownames_to_column("Sample") |>
  pivot_longer(-Sample, names_to = "Cluster", values_to = "Percent") |>
  mutate(
    Sample  = factor(Sample,  levels = rownames(pct_mat)),
    Cluster = factor(Cluster, levels = colnames(pct_mat))
  )

p_heatmap <- ggplot(heat_df, aes(x = Cluster, y = Sample, fill = Percent)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(Percent), "", round(Percent, 0))), size = 3) +
  scale_fill_gradientn(
    colours  = c("white", "#9ecae1", "#3182bd", "#08519c"),
    limits   = c(0, 100),
    na.value = "grey95",
    name     = "% células\nresponsivas"
  ) +
  labs(title = "% células responsivas", x = "Clúster", y = "Muestra") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid    = element_blank(),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    axis.title    = element_text(face = "bold"),
    plot.title    = element_text(face = "bold", hjust = 0.5)
  )

save_plot_px(p_heatmap, "sig_digest_heatmap.png", width = 2400, height = 1800)

# 9. Figura compuesta (suplementaria) -----------------------------------------

p_umap_clusters <- DimPlot(
  act,
  reduction = "umap.harmony",
  raster    = TRUE,
  shuffle   = TRUE,
  label     = TRUE
) + theme(legend.position = "none") + umap_theme()

top_row <- (
  (p_umap_clusters + ggtitle("A") + theme(plot.margin = unit(c(0, 15, 0, 0), "pt"))) |
  (p_sig_umap      + ggtitle("B") + theme(plot.margin = unit(c(0,  0, 0, 12), "pt")))
) + plot_layout(widths = c(1.8, 2))

bottom_row <- (
  (p_vln2    + ggtitle("C")) |
  (p_bar     + ggtitle("D")) |
  (p_heatmap + ggtitle("E"))
) + plot_layout(widths = c(1.35, 0.9, 1.8))

dig_supl <- (top_row / bottom_row) +
  plot_layout(heights = c(1, 1.15)) &
  theme(
    plot.title          = element_text(hjust = 0, face = "bold"),
    plot.title.position = "plot"
  )

ggsave("digestion_supl.png", dig_supl,
       width = 4500, height = 2500, units = "px", dpi = 300)

# 10. Guardar objeto y células responsivas ------------------------------------

qs_save(act, "act_digest_scored.qs2")

responsive_cells <- colnames(act)[act$Digest_resp]
message("Células responsivas: ", length(responsive_cells))
saveRDS(responsive_cells, "responsive_cells_ucell.rds")

final_cells <- setdiff(colnames(act), responsive_cells)
saveRDS(final_cells, "final_cells.rds")
saveRDS(rownames(act), "final_genes.rds")

message("Guardado: act_digest_scored.qs2")
message("Guardado: responsive_cells_ucell.rds  (", length(responsive_cells), " células)")
message("Guardado: final_cells.rds             (", length(final_cells), " células)")

##A&B
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(readxl)
library(irGSEA)

sce.all <- readRDS("data/zenodo_final.rds")

gene_list <- readxl::read_excel("Ribosis.xlsx")
genes <- as.character(gene_list$gene)
custom_geneset <- list(Ribosome_Genes = genes)

sce.all <- irGSEA.score(
  object = sce.all, 
  assay = "RNA", 
  slot = "data", 
  seeds = 123, 
  ncores = 4, 
  min.cells = 3, 
  custom = TRUE, 
  geneset = custom_geneset, 
  msigdb = FALSE, 
  species = "Homo sapiens", 
  method = c("singscore","ssgsea", "viper"), 
  kcdf = 'Gaussian'
)

singscore_data <- GetAssayData(sce.all, assay = "singscore", slot = "data")
sce.all$meta.data$singscore <- as.numeric(singscore_data["Ribosome-Genes", ])

ssgsea_data <- GetAssayData(sce.all, assay = "ssgsea", slot = "data")
sce.all$meta.data$ssgsea <- as.numeric(ssgsea_data["Ribosome-Genes", ])

viper_data <- GetAssayData(sce.all, assay = "viper", slot = "data")
sce.all$meta.data$viper <- as.numeric(viper_data["Ribosome-Genes", ])

data <- data.frame(
  cell_type = sce.all$meta.data$celltype,
  singscore = sce.all$meta.data$singscore,
  ssgsea = sce.all$meta.data$ssgsea,
  viper = sce.all$meta.data$viper
)

normalize_with_negatives <- function(x) {
  (x - mean(x)) / sd(x)
}

data[, c("singscore", "ssgsea", "viper")] <- lapply(
  data[, c("singscore", "ssgsea", "viper")], 
  normalize_with_negatives
)

melted_data <- tidyr::pivot_longer(
  data,
  cols = c(singscore, ssgsea, viper),
  names_to = "method",
  values_to = "score"
)

summary_data <- melted_data %>%
  group_by(cell_type, method) %>%
  summarise(
    average_expression = mean(score),
    percent_expressed = sum(score > 0) / n() * 100,
    .groups = "drop"
  )

summary_data$method <- factor(
  summary_data$method,
  levels = c("singscore", "ssgsea", "viper")
)

p_a <- ggplot(summary_data, aes(x = method, y = cell_type)) +
  geom_point(
    aes(size = percent_expressed, color = average_expression),
    alpha = 0.8,
    shape = 16
  ) +
  scale_size_continuous(
    name = "Percent Expressed (%)",
    range = c(2, 12),
    breaks = seq(20, 100, 20),
    limits = c(0, 100)
  ) +
  scale_color_gradient2(
    name = "Average Expression",
    low = "#2166ac",
    mid = "#f7f7f7",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-2, 2)
  ) +
  labs(x = "Scoring Method", y = "Cell Type") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    coord_fixed(ratio = 0.8)
  )

df_scores <- sce.all$meta.data[, c("singscore", "ssgsea", "viper")]
scale_scores <- scale(df_scores) %>% as.data.frame()
pca_res <- prcomp(scale_scores, center = TRUE, scale. = TRUE)
scale_scores$combined <- pca_res$x[, 1]
sce.all$meta.data$combined_score <- scale_scores$combined
sce.all$meta.data$score_group <- ifelse(
  sce.all$meta.data$combined_score > median(sce.all$meta.data$combined_score),
  "High_score", "Low_score"
)

p_b1 <- DimPlot(sce.all, group.by = "celltype", reduction = "umap", pt.size = 0.5) +
  NoLegend() +
  ggtitle("Cell Type")

p_b2 <- DimPlot(sce.all, group.by = "score_group", reduction = "umap", pt.size = 0.5) +
  NoLegend() +
  ggtitle("Score Group")

p_b <- plot_grid(p_b1, p_b2, ncol = 2, labels = c("B1", "B2"))

ggsave("plots/fig3/fig3A_bubble.pdf", p_a, width = 8, height = 6)
ggsave("plots/fig3/fig3B_umap.pdf", p_b, width = 12, height = 6)

print(p_a)
print(p_b)
##C
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(SCP)
library(grid)

rm(list = ls())

sce.all <- readRDS("AML_zenodo.rds")

Idents(sce.all) <- 'celltype'

dat <- sce.all@meta.data
data <- dat[, c("celltype", "score_group")]

R_oe <- calTissueDist(data,
                      byPatient = FALSE,
                      colname.cluster = "celltype",
                      colname.tissue = "score_group",
                      method = "chisq",
                      min.rowSum = 0)

Roe <- R_oe

col_breaks <- c(0, 0.5, 1, 1.5, 2)
col_palette <- c("#4575B4", "#91BFDB", "#FFFFBF", "#FEE090", "#FC8D59")
col_fun <- colorRamp2(col_breaks, col_palette)

celltype_colors <- c(
  "HSC" = "#66C5CC",
  "T/NK" = "#F6CF71",
  "B" = "#F89C74",
  "Erythroid" = "#DCB0F2",
  "GMP" = "#87C55F",
  "Mono/Marco" = "#9EB9F3",
  "MEP" = "#FE88B1",
  "Plasma" = "#C9DB74",
  "pDC" = "#8BE0A4",
  "Neutrophils" = "#B1A5EB"
)

row_ha <- rowAnnotation(
  CellType = rownames(Roe),
  col = list(CellType = celltype_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

hm <- Heatmap(
  as.matrix(Roe),
  name = "Ro/e",
  col = col_fun,
  rect_gp = gpar(col = "grey90", lwd = 0.5),
  border = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_side = "left",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "Ro/e Index",
    at = col_breaks,
    labels = col_breaks
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- Roe[i, j]
    if (!is.na(val)) {
      label <- case_when(
        val > 1.5 ~ "+++",
        val > 1 ~ "++",
        val > 0.5 ~ "+",
        val > 0 ~ "+/-",
        val == 0 ~ "-",
        TRUE ~ ""
      )
      rgb_col <- col2rgb(fill)
      brightness <- (0.299 * rgb_col[1] + 0.587 * rgb_col[2] + 0.114 * rgb_col[3]) / 255
      text_color <- ifelse(brightness > 0.6, "black", "white")
      grid.text(label, x, y, gp = gpar(fontsize = 10, col = text_color))
    }
  },
  left_annotation = row_ha
)

lgd_list <- list(
  Legend(
    labels = c("+++ (>1.5)", "++ (>1)", "+ (>0.5)", "+/- (>0)", "- (0)"),
    title = "Significance Levels",
    type = "points",
    pch = c(rep("+",4), "-"),
    background = "transparent"
  ),
  Legend(
    labels = names(celltype_colors),
    title = "Cell Types",
    type = "points",
    pch = 15,
    legend_gp = gpar(col = celltype_colors)
  )
)

pdf("ROe.pdf", width = 7, height = 9, bg = "white")
draw(hm, annotation_legend_list = lgd_list)
dev.off()
##D
library(Seurat)
library(infercnv)
library(AnnoProbe)
library(future)
library(data.table)

expr_matrix <- as.matrix(sce.all@assays$RNA@counts)
cell_annot <- data.frame(
  cell_type = sce.all$celltype,
  row.names = colnames(sce.all)
)

reference_group <- "T/NK"
observation_groups <- c("HSC", "GMP", "MEP")

cells_use <- colnames(sce.all)[sce.all$celltype %in% c(reference_group, observation_groups)]
dat <- as.data.frame(expr_matrix[, cells_use])
groupinfo <- data.frame(v1 = colnames(dat), v2 = sce.all$celltype[cells_use])

geneInfor <- annoGene(rownames(dat), "SYMBOL", "human")
geneInfor <- geneInfor[with(geneInfor, order(chr, start)), c(1, 4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]
dat <- dat[rownames(dat) %in% geneInfor[, 1], ]
dat <- dat[match(geneInfor[, 1], rownames(dat)), ]

write.table(groupinfo, file = "groupFiles.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(dat, file = "expFile.txt", sep = "\t", quote = FALSE)
write.table(geneInfor, file = "geneFile.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

infercnv_obj <- CreateInfercnvObject(
  delim = "\t",
  raw_counts_matrix = "expFile.txt",
  annotations_file = "groupFiles.txt",
  gene_order_file = "geneFile.txt",
  ref_group_names = reference_group
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "inferCNV/",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  plot_steps = TRUE,
  analysis_mode = "subclusters",
  HMM = FALSE,
  num_threads = 6
)

plot_cnv(
  infercnv_obj,
  obs_title = "Observations (Cells)",
  ref_title = "References (Cells)",
  cluster_by_groups = TRUE,
  plot_chr_scale = TRUE,
  x.center = 1,
  x.range = "auto",
  hclust_method = "ward.D",
  custom_color_pal = color.palette(c("#0071B2", "white", "#C3250A"), c(2, 2)),
  color_safe_pal = FALSE,
  output_filename = "infercnv_D",
  output_format = "pdf",
  png_res = 600,
  dynamic_resize = 0
)
##E
library(Seurat)
library(GSVA)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(RColorBrewer)
library(circlize)
library(BiocParallel)
library(SingleCellExperiment)

sce.all <- readRDS("data/zenodo_final.rds")

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

meta_data <- sce.all@meta.data %>%
  tibble::rownames_to_column("cell_id")

stratified_sample <- meta_data %>%
  group_by(celltype, score_group) %>%
  sample_frac(0.1)

sce_subset <- subset(sce.all, cells = stratified_sample$cell_id)

expr_matrix <- LayerData(sce_subset, assay = "RNA", layer = "data")

params <- gsvaParam(
  exprData = as.matrix(expr_matrix),
  geneSets = hallmark_list,
  minSize = 10,
  maxSize = 500,
  kcdf = "Gaussian",
  tau = 0.25,
  maxDiff = TRUE,
  absRanking = FALSE
)

register(MulticoreParam(workers = 4))

GSVA_hall <- gsva(
  params,
  verbose = TRUE,
  BPPARAM = MulticoreParam(workers = 4, progressbar = TRUE)
)

gsva_scores <- t(GSVA_hall)
colnames(gsva_scores) <- paste0("HALLMARK_", colnames(gsva_scores))

rownames(gsva_scores) <- colnames(sce_subset)
sce_subset@meta.data <- cbind(sce_subset@meta.data, gsva_scores)

df <- sce_subset@meta.data
hallmark_cols <- grep("^HALLMARK_", colnames(df), value = TRUE)

avg_df <- aggregate(
  df[, colnames(df) %in% hallmark_cols],
  by = list(score_group = df$score_group, celltype = df$celltype),
  FUN = mean
)

custom_celltype_order <- c("HSC", "Mono/Marco", "T/NK", "GMP", "B", "MEP","pDC", "Erythroid", "Plasma", "Neutrophils")
score_group_order <- c("High_score", "Low_score")

avg_df$celltype <- factor(avg_df$celltype, levels = custom_celltype_order)
avg_df$score_group <- factor(avg_df$score_group, levels = score_group_order)
avg_df <- avg_df[order(avg_df$score_group, avg_df$celltype), ]

group_combo <- paste(avg_df$score_group, avg_df$celltype, sep = "_")
exp <- scale(as.matrix(avg_df[, -(1:2)]))
rownames(exp) <- group_combo

column_ha <- HeatmapAnnotation(
  Group = avg_df$score_group,
  CellType = avg_df$celltype,
  col = list(
    Group = c("High_score" = "#FF7F00", "Low_score" = "#1F78B4"),
    CellType = setNames(brewer.pal(10, "Paired"), custom_celltype_order)
  ),
  annotation_name_side = "left"
)

h_gsva <- Heatmap(
  t(exp),
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c('#1A5592', 'white', "#B83D3D")),
  top_annotation = column_ha,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_split = avg_df$score_group,
  cluster_column_slices = FALSE,
  column_title = c("High Score Group", "Low Score Group"),
  column_names_rot = 45,
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = "gray80", lwd = 0.5),
  heatmap_legend_param = list(title_position = "lefttop-rot")
)

pdf("plots/fig3/fig3E_gsva_heatmap.pdf", width = 14, height = 10)
draw(h_gsva, heatmap_legend_side = "right")
dev.off()
##F
library(scMetabolism)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)
library(circlize)
library(Seurat)

sce.all <- readRDS("data/zenodo_final.rds")

set.seed(123)
sampled_cells <- sample(colnames(sce.all), size = floor(0.20 * length(colnames(sce.all))))
sce_sampled <- sce.all[, sampled_cells]

sce.metab <- sc.metabolism.Seurat(
  obj = sce_sampled,
  method = "AUCell",
  imputation = FALSE,
  ncores = 1,
  metabolism.type = "KEGG"
)

score <- sce.metab@assays$METABOLISM$score
colnames(score) <- gsub("\\.", "-", colnames(score))
sce.metab@meta.data <- cbind(sce.metab@meta.data, t(score))

ribo_related <- c(
  "Amino sugar and nucleotide sugar metabolism",
  "Purine metabolism",
  "Pyrimidine metabolism",
  "Alanine, aspartate and glutamate metabolism",
  "Glycine, serine and threonine metabolism",
  "Cysteine and methionine metabolism",
  "Glycolysis / Gluconeogenesis",
  "Citrate cycle (TCA cycle)",
  "Oxidative phosphorylation",
  "Pyruvate metabolism",
  "Glutathione metabolism",
  "Nicotinate and nicotinamide metabolism",
  "Pantothenate and CoA biosynthesis"
)

selected_pathways <- c(
  ribo_related,
  "Fatty acid biosynthesis",
  "Fatty acid degradation",
  "Steroid biosynthesis",
  "Glycerophospholipid metabolism",
  "Sphingolipid metabolism",
  "Arachidonic acid metabolism",
  "Lysine degradation",
  "Arginine and proline metabolism",
  "Tryptophan metabolism",
  "N-Glycan biosynthesis",
  "Vitamin B6 metabolism",
  "Folate biosynthesis",
  "One carbon pool by folate",
  "Porphyrin and chlorophyll metabolism",
  "Drug metabolism - cytochrome P450",
  "Pentose phosphate pathway",
  "Propanoate metabolism",
  "Butanoate metabolism",
  "Inositol phosphate metabolism",
  "Sulfur metabolism",
  "Synthesis and degradation of ketone bodies",
  "Glycerolipid metabolism",
  "Ether lipid metabolism",
  "Linoleic acid metabolism",
  "Valine, leucine and isoleucine degradation",
  "Histidine metabolism",
  "Tyrosine metabolism",
  "beta-Alanine metabolism",
  "Taurine and hypotaurine metabolism",
  "Riboflavin metabolism",
  "Biotin metabolism",
  "Lipoic acid metabolism",
  "Retinol metabolism",
  "Ubiquinone and other terpenoid-quinone biosynthesis"
)

selected_pathways <- unique(intersect(selected_pathways, rownames(sce.metab@assays$METABOLISM$score)))
selected_pathways <- selected_pathways[1:50]

df <- sce.metab@meta.data
avg_df <- aggregate(
  df[, colnames(df) %in% selected_pathways],
  by = list(score_group = df$score_group, celltype = df$celltype),
  FUN = mean
)

custom_order <- c("HSC", "Mono/Marco", "T/NK", "GMP", "B", "MEP", "pDC", "Erythroid", "Plasma", "Neutrophils")
avg_df$celltype <- factor(avg_df$celltype, levels = custom_order)
avg_df$score_group <- factor(avg_df$score_group, levels = c("High_score", "Low_score"))
avg_df <- avg_df[order(avg_df$score_group, avg_df$celltype), ]

group_combo <- paste(avg_df$score_group, avg_df$celltype, sep = "_")
exp <- scale(as.matrix(avg_df[, -(1:2)]))
rownames(exp) <- group_combo

column_ha <- HeatmapAnnotation(
  Group = avg_df$score_group,
  CellType = avg_df$celltype,
  col = list(
    Group = c("High_score" = "#FF7F00", "Low_score" = "#1F78B4"),
    CellType = setNames(brewer.pal(10, "Paired"), custom_order)
  ),
  annotation_name_side = "left"
)

row_anno <- rowAnnotation(
  Ribosome_related = ifelse(rownames(t(exp)) %in% ribo_related, "Yes", "No"),
  col = list(Ribosome_related = c("Yes" = "gold", "No" = "transparent")),
  show_legend = TRUE,
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 10)
)

h_combined <- Heatmap(
  t(exp),
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c('#1A5592', 'white', "#B83D3D")),
  top_annotation = column_ha,
  right_annotation = row_anno,
  cluster_columns = FALSE,
  column_split = avg_df$score_group,
  cluster_column_slices = FALSE,
  column_title = c("High Score Group", "Low Score Group"),
  column_names_rot = 45,
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  row_title = "Metabolic Pathways",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  rect_gp = gpar(col = "gray80", lwd = 0.5),
  border = TRUE,
  heatmap_legend_param = list(
    title_position = "lefttop-rot",
    legend_height = unit(4, "cm")
  )
)

pdf("plots/figX/figXE_metabolic_heatmap.pdf", width = 14, height = 10)
draw(h_combined, heatmap_legend_side = "right")
dev.off()

print(h_combined)

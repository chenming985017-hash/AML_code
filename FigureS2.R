##A&B
rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(future)

# Data loading and Seurat object creation
mydata <- Read10X_h5("AML_GSE116256_expression.h5")
phe <- data.table::fread("AML_GSE116256_CellMetainfo_table.tsv")
phe <- as.data.frame(phe)
rownames(phe) <- phe$Cell
pbmc <- CreateSeuratObject(counts = mydata, project = "GSE116256", min.cells = 5, min.features = 300)
pbmc <- AddMetaData(pbmc, phe, col.name = colnames(phe))
pbmc <- pbmc[, pbmc$Source %in% "AML"]

# Normalization, variable features, PCA, Harmony, UMAP
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunHarmony(pbmc, "Patient")
pbmc <- RunUMAP(pbmc, dims = 1:15, reduction = "harmony")

# Assign celltype and plot UMAP (Figure A)
pbmc$celltype <- pbmc$Celltype..minor.lineage.
mycolors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3',
              '#E95C59', '#E59CC4', '#AB3282',  '#BD956A', 
              '#9FA3A8', '#E0D4CA',  '#C5DEBA', '#F7F398',
              '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
              '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
              '#968175')
pdf("Supplementary_Figure_S2A_UMAP.pdf", height = 8, width = 10)
DimPlot(pbmc, reduction = "umap", group.by = "celltype", cols = mycolors, label = TRUE, repel = TRUE, raster = FALSE)
dev.off()

# Find marker genes and plot heatmap (Figure B)
pbmc <- SetIdent(pbmc, value = "celltype")
plan("multicore", workers = 8)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("Supplementary_Figure_S2B_Heatmap.pdf", height = 10, width = 12)
DoHeatmap(pbmc, features = top_markers$gene, group.by = "celltype", size = 3) + NoLegend()
dev.off()
save(pbmc,file = "AML.RData")
##C&D
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(readxl)
library(irGSEA)

# 
load("AML.RData")  
sce.all <- pbmc  

# 
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

ggsave("bubble.pdf", p_a, width = 8, height = 6)
ggsave("umap.pdf", p_b, width = 12, height = 6)

print(p_a)
print(p_b)
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

load("AML.RData")  
sce.all <- pbmc 

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

load("AML.RData")  
sce.all <- pbmc

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

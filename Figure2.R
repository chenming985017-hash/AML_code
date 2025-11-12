library(Seurat)
library(SCP)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)
library(tidydr)
library(cols4all)
library(mascarade)
library(ggtext)
library(RColorBrewer)
library(ggh4x)
library(ggalluvial)
library(ggthemes)
library(ClusterGVis)
library(org.Hs.eg.db)
library(jjAnno)
sce.all <- readRDS("data/zenodo_final.rds")
##A
CellDimPlot(sce.all,group.by="celltype",reduction="UMAP")
##B
elltype <- Idents(sce.all)
df <- FetchData(object=sce.all, vars=c("umap_1","umap_2","RNA_snn_res.0.2","celltype","seurat_clusters"))
head(df)
str(df)
dim(df)
maskTable <- generateMask(dims=df[,1:2], cluster=df$celltype, minDensity = 1.5, smoothSigma = 0.05)
class(maskTable)
dim(maskTable)
head(maskTable)
UMAP <- as.data.frame(sce.all@reductions$umap@cell.embeddings)
UMAP <- cbind(UMAP, celltype)
p <- ggplot(df, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=celltype), size = 1) + 
  geom_path(data=maskTable, aes(group=group), linewidth=0.3, linetype = 2) +
  coord_fixed() + 
  theme_classic()
p
p4 <- p +
  theme_dr(xlength = 0.2, 
           ylength = 0.2, 
           arrow = grid::arrow(length = unit(0.1, "inches"), 
                               ends = 'last', type = "closed")) + 
  theme(panel.grid = element_blank())
p4
label <- UMAP %>%
  group_by(celltype) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))
head(label)
p5 <- p4 +
  geom_text(data = label,
            aes(x = umap_1, y = umap_2, label = celltype),
            fontface = "bold",
            color = 'black', size = 3)
p5
c4a_gui()
mycol <- c4a('light', 10)
color2 <- mycol
CairoPDF("umap_celltype_label.pdf", width=12, height=9)
p7 <- p5 +
  scale_color_manual(values = color2) +
  scale_fill_manual(values = color2)
dev.off()
##C
genes_to_check <- list(
  HSC = c("CD34","HOXA9"),
  "T/NK" = c("CD3D","CD3G"),
  B = c("CD19","CD79A"),
  Erythroid = c("HBB","HBA1"),
  GMP = c("AZU1","CTSG","MPO"),
  "Mono/Marco"  = c("LYZ","CD68","FTL"),
  MEP = c("GATA1","NFE2"),
  Plasma = c("SDC1", "MZB1"),
  pDC = c("CD1C", "CLEC10A","IRF8" ),
  Neutrophils = c("CD177","FCGR3B","CMTM2")
)

all_genes <- unlist(genes_to_check)
valid_genes <- all_genes[all_genes %in% rownames(sce.all)]
if(length(valid_genes) < length(all_genes)) {
  missing_genes <- all_genes[!all_genes %in% rownames(sce.all)]
  warning(paste("以下基因在数据集中不存在:", paste(missing_genes, collapse=", ")))
}

counts_data <- GetAssayData(sce.all, slot = "counts")

if (!"celltype" %in% colnames(sce.all@meta.data)) {
  stop("sce.all 对象的元数据中不存在 'celltype' 列，请检查并添加。")
}

result_data <- data.frame(
  Function = character(),
  Gene = character(),
  Group = character(),
  AvgExpr = numeric(),
  PctExpr = numeric()
)

for(func_name in names(genes_to_check)) {
  func_genes <- genes_to_check[[func_name]]
  func_genes <- func_genes[func_genes %in% rownames(counts_data)]
  
  for(group in levels(sce.all$celltype)) {
    cells_in_group <- WhichCells(sce.all, expression = celltype == group)
    
    for(gene in func_genes) {
      gene_expr <- counts_data[gene, cells_in_group]
      avg_expr <- mean(gene_expr)
      pct_expr <- sum(gene_expr > 0) / length(gene_expr) * 100
      
      result_data <- rbind(result_data, data.frame(
        Function = func_name,
        Gene = gene,
        Group = group,
        AvgExpr = avg_expr,
        PctExpr = pct_expr
      ))
    }
  }
}

for(gene in unique(result_data$Gene)) {
  gene_idx <- result_data$Gene == gene
  max_expr <- max(result_data$AvgExpr[gene_idx])
  if(max_expr > 0) {
    result_data$AvgExpr[gene_idx] <- result_data$AvgExpr[gene_idx] / max_expr
  }
}

result_data$Function <- factor(result_data$Function, levels = names(genes_to_check))
result_data$Gene <- factor(result_data$Gene, levels = rev(all_genes))
result_data$Group <- factor(result_data$Group, 
                            levels = c("HSC", "T/NK", "B", "Erythroid", 
                                       "GMP","Mono/Marco", "MEP", "Plasma", "pDC", 
                                       "Neutrophils"))

cell_cluster_colors <- c(
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

result_data$Function <- factor(
  result_data$Function,
  levels = names(cell_cluster_colors)
)
function_colors <- cell_cluster_colors[levels(result_data$Function)]

p <- ggplot(result_data, aes(y = Group, x = Gene, color = AvgExpr, size = PctExpr)) +
  geom_point() +
  scale_color_gradientn(
    colours = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")),
    name = "Relative\nExpression"
  ) +
  scale_size_continuous(name = "Percentage\nExpressed", range = c(0, 6)) +
  ggh4x::facet_grid2(
    . ~ Function,
    scales = "free_x",
    space = "free_x",
    strip = strip_themed(background_x = elem_list_rect(fill = function_colors))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(color = "black")
  )

CairoPDF("plots/fig3/fig3_bubble_plot.pdf", width=12, height=8)
print(p)
dev.off()
##D&E
meta <- sce.all@meta.data
meta <- meta[, c("orig.ident", "celltype")]
colnames(meta) <- c("Sample", "CellType")

sankey_data <- meta %>%
  group_by(Sample, CellType) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count))

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

p_d <- ggplot(sankey_data, aes(
  x = Sample,
  stratum = CellType,
  alluvium = CellType,
  y = Proportion,
  fill = CellType
)) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = celltype_colors, name = "Cell Type") +
  geom_flow() +
  geom_stratum(width = 0.6, colour = "grey40") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title = element_blank()
  )

meta$Group <- ifelse(grepl("AML", meta$Sample), "AML", "HD")

bar_data <- meta %>%
  group_by(Group, CellType) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Proportion = Count / sum(Count))

p_e <- ggplot(bar_data, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = celltype_colors, name = "Cell Type") +
  theme_few() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(x = "", y = "Relative proportion(%)")

ggsave("plots/figX/figXD_sankey.pdf", p_d, width = 12, height = 8)
ggsave("plots/figX/figXE_stacked_bar.pdf", p_e, width = 6, height = 6)

print(p_d)
print(p_e)
##F
pbmc.markers.all <- FindAllMarkers(sce.all,
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25)

filtered_markers <- pbmc.markers.all[!grepl("\\.", rownames(pbmc.markers.all)) & 
                                       !grepl("^MT-", rownames(pbmc.markers.all)), ]

pbmc.markers <- filtered_markers %>%
  group_by(cluster) %>%
  top_n(n = 40, wt = avg_log2FC)

st.data <- prepareDataFromscRNA(object = sce.all,
                                diffData = pbmc.markers,
                                showAverage = TRUE)

enrich.go <- enrichCluster(object = st.data,
                           OrgDb = org.Hs.eg.db,
                           type = "BP",
                           organism = "hsa",
                           pvalueCutoff = 0.5,
                           topn = 5,
                           seed = 5201314)

enrich.kegg <- enrichCluster(object = st.data,
                             OrgDb = org.Hs.eg.db,
                             type = "KEGG",
                             organism = "hsa",
                             pvalueCutoff = 0.5,
                             topn = 5,
                             seed = 5201314)

markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),60,
                                             replace = F)]

go.col = rep(useMyCol("calm",n = 10),each = 5)
kegg.col = rep(useMyCol("stallion",n = 10),each = 5)

pdf('plots/figF/figF_gokegg.pdf',height = 10,width = 20,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich.go,
           go.col = go.col,
           add.bar = T,
           by.go = "anno_block",
           annoKegg.data = enrich.kegg,
           kegg.col = kegg.col,
           add.kegg.bar = T,
           by.kegg = "anno_block",
           word_wrap = F,
           add_new_line = F,
           line.side = "left",
           cluster.order = c(1:10)
)
dev.off()



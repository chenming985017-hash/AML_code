##A
rm(list = ls())
library(limma)
library(ggplot2)
library(ggrepel)

expFile <- "geneexpr.txt"
clusterFile <- "Cluster.txt"

rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
exp <- rt
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
data <- data[rowMeans(data) > 0, ]
data <- t(data)

data.pca <- prcomp(data, scale. = TRUE)
pcaPredict <- predict(data.pca)
write.table(pcaPredict, file = "newTab.xls", quote = FALSE, sep = "\t")

cluster <- read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cluster <- as.factor(cluster[,1])

bioCol <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
CluCol <- bioCol[1:nlevels(cluster)]

PCA <- data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2], Cluster = cluster)
PCA.mean <- aggregate(PCA[,1:2], list(Cluster = PCA$Cluster), mean)

p <- ggplot(PCA, aes(PC1, PC2)) +
  stat_ellipse(aes(color = Cluster), type = "t", level = 0.95, size = 0.8, linetype = "dashed", alpha = 0.5) +
  geom_point(aes(color = Cluster, fill = Cluster), size = 3.5, shape = 21, alpha = 0.8) +
  geom_label_repel(data = PCA.mean, aes(label = Cluster, fill = Cluster), color = "white", size = 5, fontface = "bold", box.padding = 0.8, show.legend = FALSE) +
  scale_color_manual(values = CluCol) +
  scale_fill_manual(values = CluCol) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey30", fill = NA, size = 0.5)
  ) +
  labs(title = "PCA Analysis with Cluster Visualization",
       x = paste0("PC1 (", round(summary(data.pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(data.pca)$importance[2,2]*100, 1), "%)")) +
  coord_fixed(ratio = 1) +
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0, color = "grey30", arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("segment", y = -Inf, yend = Inf, x = 0, xend = 0, color = "grey30", arrow = arrow(length = unit(0.3, "cm")))

ggsave("fig7A_PCA.pdf", plot = p, width = 10, height = 8, dpi = 300)
ggsave("fig7A_PCA.png", plot = p, width = 10, height = 8, dpi = 300)
##B
rm(list = ls())
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
library(tidyverse)

expFile <- "merged_GSE_expression.csv"
gmtFile <- "immune.gmt"
clusterFile <- "Cluster.txt"

outTab <- read.csv(expFile, row.names = 1, check.names = FALSE)
data_exp <- as.matrix(outTab)

geneSets <- getGmt(gmtFile, geneIdType = SymbolIdentifier())
ssgseaScore <- gsva(data_exp, geneSets, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

normalize <- function(x){(x - min(x))/(max(x) - min(x))}
ssgseaScore <- normalize(ssgseaScore)

ssgseaOut <- rbind(id = colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut, file = "ssGSEA.result.txt", sep = "\t", quote = FALSE, col.names = FALSE)

cluster <- read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

ssgseaScore <- t(ssgseaScore)
sameSample <- intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore <- ssgseaScore[sameSample, , drop = FALSE]
cluster <- cluster[sameSample, , drop = FALSE]
scoreCluster <- cbind(ssgseaScore, cluster)

data_melt <- melt(scoreCluster, id.vars = c("cluster"))
colnames(data_melt) <- c("cluster", "Immune", "Fraction")

bioCol <- c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol <- bioCol[1:length(levels(factor(data_melt[,"cluster"])))]

data_melt$Immune <- str_replace_all(data_melt$Immune, "na", "")

p <- ggplot(data_melt, aes(x = Immune, y = Fraction, fill = cluster)) +
  geom_boxplot(width = 0.7, outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.5, alpha = 0.7) +
  geom_jitter(aes(color = cluster), position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.4) +
  scale_fill_manual(values = bioCol) +
  scale_color_manual(values = bioCol) +
  stat_compare_means(
    aes(group = cluster),
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    label = "p.signif",
    size = 4
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  labs(
    title = "Immune Infiltration by Cluster",
    x = "",
    y = "Immune Infiltration Score",
    fill = "Cluster",
    color = "Cluster"
  )

ggsave(
  "fig7B_immune_infiltration.pdf",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)
##C
rm(list = ls())
library(IOBR)
library(tidyverse)
library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
library(pheatmap)
library(ggpubr)
library(reshape2)

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

# Read data
expFile <- "merged_GSE_expression.csv"
clusterFile <- "Cluster.txt"
outTab <- read.csv(expFile, row.names = 1, check.names = FALSE)
eset_stad <- as.data.frame(outTab)
rs_merge <- read.table(clusterFile, header = TRUE, sep = "\t", row.names = 1)

# TME deconvolution
im_mcpcounter <- deconvo_tme(eset = eset_stad, method = "mcpcounter")
im_epic <- deconvo_tme(eset = eset_stad, method = "epic", arrays = FALSE)
im_xcell <- deconvo_tme(eset = eset_stad, method = "xcell", arrays = FALSE)
im_cibersort <- deconvo_tme(eset = eset_stad, method = "cibersort", arrays = FALSE, perm = 1000)
im_ips <- deconvo_tme(eset = eset_stad, method = "ips", plot = FALSE)
im_quantiseq <- deconvo_tme(eset = eset_stad, method = "quantiseq", scale_mrna = TRUE)
im_estimate <- deconvo_tme(eset = eset_stad, method = "estimate")
im_estimate$ID <- str_replace_all(im_estimate$ID, "-", ".")
im_timer <- deconvo_tme(eset = eset_stad, method = "timer", group_list = rep("gbm", dim(eset_stad)[2]))

# Combine TME scores
tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by = "ID") %>% 
  inner_join(im_xcell, by = "ID") %>% 
  inner_join(im_cibersort, by = "ID") %>% 
  inner_join(im_ips, by = "ID") %>% 
  inner_join(im_quantiseq, by = "ID") %>% 
  inner_join(im_estimate, by = "ID") %>% 
  inner_join(im_timer, by = "ID")

write.csv(tme_combine, file = "tme.score.csv", row.names = FALSE)
save(tme_combine, file = "tme.score.RDATA")

# Prepare data for heatmap
load("tme.score.RDATA")
tme_combine <- as.data.frame(tme_combine)
tme_combine <- tme_combine[, colSums(is.na(tme_combine)) == 0]
rownames(tme_combine) <- tme_combine$ID
tme_combine <- tme_combine[, -1]

aaa <- data.frame(ID = colnames(tme_combine), Type = NA)
aaa$Type <- sapply(strsplit(as.character(aaa$ID), "_"), function(x) tail(x, n = 1))
dataa <- as.data.frame(t(tme_combine))
genename1 <- intersect(aaa$ID, rownames(dataa))
dataa <- dataa[genename1, ]
aaa <- dplyr::filter(aaa, ID %in% genename1)

rownames(rs_merge) <- rs_merge$ID
sameid <- intersect(colnames(dataa), rownames(rs_merge))
dataa <- dataa[, sameid]
rs_merge <- rs_merge[sameid, ]
hmdat <- as.data.frame(t(dataa))

# Annotation
annCol <- data.frame(
  cluster = rs_merge$cluster,
  row.names = rownames(rs_merge),
  stringsAsFactors = FALSE
)

# Statistical test for gene names
sigVec <- c()
sampleType <- annCol$cluster
for(i in colnames(hmdat)){
  test <- wilcox.test(hmdat[, i] ~ sampleType)
  pvalue <- test$p.value
  Sig <- ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")))
  sigVec <- c(sigVec, paste0(i, Sig))
}
colnames(hmdat) <- sigVec

# Standardization function
standarize.fun <- function(indata = NULL, halfwidth = NULL, centerFlag = TRUE, scaleFlag = TRUE) {  
  outdata <- t(scale(t(indata), center = centerFlag, scale = scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata > halfwidth] <- halfwidth
    outdata[outdata < (-halfwidth)] <- -halfwidth
  }
  return(outdata)
}

# Annotation colors
Type.col <- brewer.pal(n = length(unique(aaa$Type)), name = "Paired")
annRow <- data.frame(Type = factor(aaa$Type, levels = unique(aaa$Type)),
                     row.names = colnames(hmdat),
                     stringsAsFactors = FALSE)

annColors <- list(
  Type = setNames(Type.col, unique(aaa$Type)),
  cluster = c("A" = "red", "B" = "blue")
)

# Sample order
samorder <- rownames(rs_merge[order(rs_merge$cluster), ])
indata <- as.data.frame(t(hmdat))
indata[indata == 0] <- 0.00000000000000000000000000001
plotdata <- standarize.fun(indata, halfwidth = 2)

# Plot heatmap
pheatmap::pheatmap(
  mat = as.matrix(plotdata[, samorder]),
  border_color = NA,
  color = colorRampPalette(c(rep("purple", 5), "white", rep("orange", 5)))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annCol[samorder, , drop = FALSE],
  annotation_row = annRow,
  annotation_colors = annColors,
  gaps_row = cumsum(table(annRow$Type)),
  cellwidth = 1.2,
  cellheight = 10,
  filename = "fig7C_immune_heatmap.pdf"
)

# Prepare data for boxplots
hmdat$ID <- rownames(hmdat)
annCol$ID <- rownames(annCol)
annCol <- dplyr::select(annCol, ID, everything())
newdata <- inner_join(annCol, hmdat, by = "ID")

colnames(newdata) <- gsub("\\*", "", colnames(newdata))
genename2 <- str_replace_all(genename1, "-", "")
genename2 <- str_replace_all(genename2, "\\+", "")
genename2 <- str_replace_all(genename2, "\\(", "")
genename2 <- str_replace_all(genename2, "\\)", "")
colnames(newdata) <- str_replace_all(colnames(newdata), "-", "")
colnames(newdata) <- str_replace_all(colnames(newdata), "\\+", "")
colnames(newdata) <- str_replace_all(colnames(newdata), "\\(", "")
colnames(newdata) <- str_replace_all(colnames(newdata), "\\)", "")

# Boxplots
bioCol <- c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol <- bioCol[1:length(levels(factor(newdata[,"cluster"])))]

for (i in 1:length(genename1)) {
  immunecell <- genename2[i]
  if (immunecell %in% colnames(newdata)) {
    p <- ggboxplot(newdata, x = "cluster", y = immunecell, color = "cluster",
                   ylab = immunecell,
                   xlab = "",
                   legend.title = "cluster",
                   palette = bioCol)
    
    p <- p + stat_compare_means(method = "kruskal.test",
                                label = "p.signif",
                                hide.ns = FALSE,
                                label.y = max(newdata[, immunecell], na.rm = TRUE) + 0.0001,
                                label.x = "B")
    
    pdf(file = paste0("fig7C_", immunecell, "_boxplot.pdf"), width = 4, height = 4) 
    print(p)
    dev.off()
    print(i)
  }
}
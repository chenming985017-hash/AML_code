rm(list = ls())
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(tidyverse)

expFile <- "merged_GSE_expression.csv"
clusterFile <- "Cluster.txt"
gmtFile <- "h.all.v2023.1.Hs.symbols.gmt"

outTab <- read.csv(expFile, row.names = 1, check.names = FALSE)
exp <- outTab
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

geneSets <- getGmt(gmtFile, geneIdType = SymbolIdentifier())
gsvaResult <- gsva(data, 
                   geneSets, 
                   min.sz = 10, 
                   max.sz = 500, 
                   verbose = TRUE,
                   parallel.sz = 1)
gsvaOut <- rbind(id = colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file = "gsvaOut.txt", sep = "\t", quote = FALSE, col.names = FALSE)

cluster <- read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

gsvaResult <- t(gsvaResult)
sameSample <- intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult <- gsvaResult[sameSample, , drop = FALSE]
cluster <- cluster[sameSample, , drop = FALSE]
gsvaCluster <- cbind(gsvaResult, cluster)
Project <- "AML"
gsvaCluster <- cbind(gsvaCluster, Project)
colnames(gsvaCluster)[which(colnames(gsvaCluster) == "cluster")] <- "cluster"

adj.P.Val.Filter <- 0.05
allType <- as.vector(gsvaCluster$cluster)
comp <- combn(levels(factor(allType)), 2)
fig_labels <- c("A", "B", "C", "D", "E")
fig_idx <- 1

for(i in 1:ncol(comp)){
  treat <- gsvaCluster[gsvaCluster$cluster == comp[2,i], ]
  con <- gsvaCluster[gsvaCluster$cluster == comp[1,i], ]
  data_bind <- rbind(con, treat)
  Type <- as.vector(data_bind$cluster)
  ann <- data_bind[, c(ncol(data_bind), (ncol(data_bind)-1))]
  data_mat <- t(data_bind[, -c((ncol(data_bind)-1), ncol(data_bind))])
  design <- model.matrix(~0 + factor(Type))
  colnames(design) <- levels(factor(Type))
  fit <- lmFit(data_mat, design)
  contrast <- paste0(comp[2,i], "-", comp[1,i])
  cont.matrix <- makeContrasts(contrast, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  allDiff <- topTable(fit2, adjust = 'fdr', number = 200000)
  allDiffOut <- rbind(id = colnames(allDiff), allDiff)
  write.table(allDiffOut, file = paste0(contrast, ".all.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  diffSig <- allDiff[with(allDiff, (abs(logFC) > 0.1 & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut <- rbind(id = colnames(diffSig), diffSig)
  write.table(diffSigOut, file = paste0(contrast, ".diff.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  bioCol <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
  ann_colors <- list()
  CluCol <- bioCol[1:length(levels(factor(allType)))]
  names(CluCol) <- levels(factor(allType))
  ann_colors[["cluster"]] <- CluCol[c(comp[1,i], comp[2,i])]
  
  termNum <- 20
  diffTermName <- as.vector(rownames(diffSig))
  diffLength <- length(diffTermName)
  if(diffLength < termNum){termNum <- diffLength}
  if(termNum > 0){
    hmGene <- diffTermName[1:termNum]
    hmExp <- data_mat[hmGene, ]
    if(fig_idx <= length(fig_labels)){
      pdf(file = paste0("fig6", fig_labels[fig_idx], "_heatmap.pdf"), height = 7, width = 10)
      pheatmap(hmExp, 
               annotation = ann,
               annotation_colors = ann_colors,
               color = colorRampPalette(c(rep("purple",2), "white", rep("orange",2)))(50),
               cluster_cols = FALSE,
               show_colnames = FALSE,
               gaps_col = as.vector(cumsum(table(Type))),
               scale = "row",
               fontsize = 10,
               fontsize_row = 7,
               fontsize_col = 10)
      dev.off()
      fig_idx <- fig_idx + 1
    }
  }
}
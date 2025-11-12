library(Seurat)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(gridExtra)

load("AML.RData")  
sce.all <- pbmc  

scRNA.list <- SplitObject(sce.all, split.by = "score_group")

sc.high <- scRNA.list[["High_score"]]
data.input.high <- GetAssayData(sc.high, layer = "data")
meta.high <- sc.high@meta.data[, "celltype", drop = FALSE]

CellChatDB <- CellChatDB.human
cellchat.high <- createCellChat(object = data.input.high, meta = meta.high, group.by = "celltype")
cellchat.high@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat.high <- subsetData(cellchat.high)
cellchat.high <- identifyOverExpressedGenes(cellchat.high)
cellchat.high <- identifyOverExpressedInteractions(cellchat.high)
cellchat.high <- computeCommunProb(cellchat.high, raw.use = TRUE, population.size = TRUE)
cellchat.high <- filterCommunication(cellchat.high, min.cells = 5)
cellchat.high <- computeCommunProbPathway(cellchat.high)
cellchat.high <- aggregateNet(cellchat.high)
cellchat.high <- netAnalysis_computeCentrality(cellchat.high)

sc.low <- scRNA.list[["Low_score"]]
data.input.low <- GetAssayData(sc.low, layer = "data")
meta.low <- sc.low@meta.data[, "celltype", drop = FALSE]

cellchat.low <- createCellChat(object = data.input.low, meta = meta.low, group.by = "celltype")
cellchat.low@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat.low <- subsetData(cellchat.low)
cellchat.low <- identifyOverExpressedGenes(cellchat.low)
cellchat.low <- identifyOverExpressedInteractions(cellchat.low)
cellchat.low <- computeCommunProb(cellchat.low, raw.use = TRUE, population.size = TRUE)
cellchat.low <- filterCommunication(cellchat.low, min.cells = 5)
cellchat.low <- computeCommunProbPathway(cellchat.low)
cellchat.low <- aggregateNet(cellchat.low)
cellchat.low <- netAnalysis_computeCentrality(cellchat.low)

cellchat.list <- list(Low = cellchat.low, High = cellchat.high)
cellchat <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))

p1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
p2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
combined_plot <- p1 + p2
ggsave("fig4A_B.pdf", plot = combined_plot, width = 12, height = 6)

par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

p5 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
p6 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
combined_rank <- p5 + p6
ggsave("fig4C_D.pdf", plot = combined_rank, width = 12, height = 6)

i <- 1
pathway.union <- union(cellchat.list[[i]]@netP$pathways, cellchat.list[[i+1]]@netP$pathways)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(cellchat.list)[i], width = 5, height = 6)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(cellchat.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ggsave("fig4_outgoing_heatmap.pdf", plot = grid.grab(), width = 10, height = 6)

ht3 <- netAnalysis_signalingRole_heatmap(cellchat.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(cellchat.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht4 <- netAnalysis_signalingRole_heatmap(cellchat.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(cellchat.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
ggsave("fig4_incoming_heatmap.pdf", plot = grid.grab(), width = 10, height = 6)

ht5 <- netAnalysis_signalingRole_heatmap(cellchat.list[[i]], pattern = "all", signaling = pathway.union, title = names(cellchat.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht6 <- netAnalysis_signalingRole_heatmap(cellchat.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(cellchat.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))
ggsave("fig4_overall_heatmap_E_F.pdf", plot = grid.grab(), width = 10, height = 6)

p7 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:10), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in High", angle.x = 45, remove.isolate = T)
p8 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:10), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in High", angle.x = 45, remove.isolate = T)
combind_bubble <- p7 + p8
ggsave("fig4_GMP_bubble.pdf", plot = combind_bubble, width = 11, height = 7)
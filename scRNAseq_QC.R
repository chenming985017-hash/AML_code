setwd("zenodo")
rm(list=ls())
options(stringsAsFactors = F)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  BiocManager::install("DoubletFinder", version = "3.16")
}
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(ggplot2)
source('scRNA_scripts/lib.R')

sce.raw <- readRDS("scRNA_merged_raw.rds")

sce.raw <- subset(
  sce.raw,
  nFeature_RNA > 50 &
    nCount_RNA > 100
)

sce.doublet <- sce.raw %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(features = VariableFeatures(.), verbose = FALSE, npcs = 30)

sweep.res.list <- paramSweep(sce.doublet, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

nExp <- round(ncol(sce.doublet) * 0.07)

sce.doublet <- doubletFinder_v3(
  seurat_object = sce.doublet,
  PCs = 1:10,
  pK = optimal.pK,
  nExp = nExp,
  reuse.pANN = FALSE,
  sct = FALSE
)

doublet.col <- colnames(sce.doublet@meta.data)[grep("DF.*doublets", colnames(sce.doublet@meta.data))]
sce.singlet <- subset(sce.doublet, subset = !!sym(doublet.col) == "Singlet")

saveRDS(sce.singlet, file = "scRNA_singlet_merged.rds")

sce.all <- readRDS("scRNA_singlet_merged.rds")

names(sce.all@assays$RNA@layers)
sce.all[["RNA"]]$counts 

LayerData(sce.all, assay = "RNA", layer = "counts")
sce.all <- JoinLayers(sce.all)
dim(sce.all[["RNA"]]$counts )

as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident)

sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all[["percent.rp"]] <- PercentageFeatureSet(sce.all, pattern = "^RP[SL]")
sce.all[["percent.hb"]] <- PercentageFeatureSet(sce.all, pattern = "^HB[^(P)]")

VlnPlot(sce.all, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0.01, group.by = "orig.ident")

sce.all <- subset(
  sce.all,
  percent.mt < 20 &
    nCount_RNA < 50000 & nCount_RNA > 1000 &
    nFeature_RNA < 6000 & nFeature_RNA > 200
)
table(sce.all@meta.data$orig.ident)

sc <- sce.all

f <- "zenodo.Rdata"
library(harmony)
if(!file.exists(f)){
  sc <- sc %>% 
    NormalizeData() %>%  
    FindVariableFeatures() %>%  
    ScaleData(features = VariableFeatures(.)) %>%  
    RunPCA(pc.genes = VariableFeatures(.))  %>%
    RunHarmony("orig.ident") %>%
    FindNeighbors(dims = 1:20, reduction = "harmony") %>% 
    FindClusters(resolution = 0.2) %>% 
    RunUMAP(dims = 1:20, reduction = "harmony") %>% 
    RunTSNE(dims = 1:20, reduction = "harmony")
  save(sce.all, file = f)
}
load("zenodo.Rdata")

ElbowPlot(sc)
p1 <- DimPlot(sc, reduction = "umap", label = T, pt.size = 0.5) + NoLegend(); p1
p2 <- DimPlot(sc, reduction = "tsne", label = T, pt.size = 0.5) + NoLegend(); p2

genes_to_check <- list(
  HSC = c("CD34","HOXA9","SMIM24"),
  "T/NK" = c("CD3D","CD3G"),
  B = c("CD19","CD79A"),
  Erythroid = c("HBB","HBA1"),
  GMP = c("AZU1","CTSG","MPO","SOX4"),
  "Mono/Marco" = c("LYZ","CD68","FTL"),
  MEP = c("GATA1","NFE2"),
  Plasma = c("SDC1", "MZB1"),
  pDC = c("CD1C", "CLEC10A","IRF8" ),
  Neutrophils = c("CD177", "FCGR3B")
) 
p_all_markers <- DotPlot(sc, features = genes_to_check, assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1), axis.title = element_blank()) 
p_all_markers

print(sc@active.ident)
sc <- RenameIdents(
  sc,
  "0" = "HSC", "1" = "Mono/Marco", "2" = "T/NK", "3" = "GMP",
  "4" = "GMP", "5" = "GMP", "6" = "T/NK",
  "7" = "B", "8" = "MEP", "9" = "pDC", "10" = "Erythroid",
  "11" = "Plasma", "12" = "Neutrophils", "13" = "Mono/Marco",
  "14" = "HSC"
)

saveRDS(sc, file = "zenodo.rds")
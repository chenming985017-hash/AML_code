##A
rm(list = ls())
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata <- t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth] <- halfwidth
    outdata[outdata<(-halfwidth)] <- -halfwidth
  }
  return(outdata)
}

load("tme.score.RDATA")
load("RS_Risk_clin.RDATA")

clin <- OSdata
tme_combine <- as.data.frame(tme_combine)
tme_combine <- tme_combine[, colSums(is.na(tme_combine)) == 0]
rownames(tme_combine) <- tme_combine$ID
tme_combine <- tme_combine[, -1]
tme_combine <- tme_combine[clin$ID, ]

aaa <- data.frame(ID = colnames(tme_combine), Type = NA)
aaa$Type <- sapply(strsplit(as.character(aaa$ID), "_"), function(x) tail(x, n = 1))
dataa <- as.data.frame(t(tme_combine))
genename1 <- intersect(aaa$ID, rownames(dataa))
dataa <- dataa[genename1, ]
aaa <- dplyr::filter(aaa, ID %in% genename1)

risk <- clin
rownames(risk) <- risk$ID
risk <- risk[, -1]
sameid <- intersect(colnames(dataa), rownames(risk))
dataa <- dataa[, sameid]
risk <- risk[sameid, ]
hmdat <- as.data.frame(t(dataa))

annRow <- data.frame(Type = factor(aaa$Type, levels = unique(aaa$Type)),
                     row.names = aaa$ID,
                     stringsAsFactors = FALSE)

annCol <- data.frame(RiskType = risk$Risk_group,
                     OS = risk$OS,
                     RS = scale(risk$RS),
                     row.names = rownames(risk),
                     stringsAsFactors = FALSE)

sigVec <- c()
sampleType <- annCol$RiskType
for (i in colnames(hmdat)) {
  test <- wilcox.test(hmdat[, i] ~ sampleType)
  pvalue <- test$p.value
  Sig <- ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")))
  sigVec <- c(sigVec, paste0(i, Sig))
}
colnames(hmdat) <- sigVec

indata <- as.data.frame(t(hmdat))
indata[indata == 0] <- 1e-20
plotdata <- standarize.fun(indata, halfwidth = 2)

Type.col <- brewer.pal(n = length(unique(annRow$Type)), name = "Paired")
annColors <- list(
  Type = setNames(Type.col, unique(annRow$Type)),
  RiskType = setNames(c("blue", "red"), c("Low_risk", "High_risk")),
  RS = greenred(64)
)

samorder <- rownames(annCol[order(annCol$RS), ])

pheatmap(
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
  gaps_col = sum(annCol$RiskType == "Low_risk"),
  gaps_row = cumsum(table(annRow$Type)),
  cellwidth = 1.5,
  cellheight = 12,
  filename = "Figure11A_immune_heatmap.pdf"
)
##B
rm(list = ls())
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata <- t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth] <- halfwidth
    outdata[outdata<(-halfwidth)] <- -halfwidth
  }
  return(outdata)
}

outTab <- read.csv("merge.csv", row.names = 1, check.names = FALSE)
load("RS_Risk_clin.RDATA")
dataa <- as.data.frame(outTab)

dataa <- dataa[, OSdata$ID]
aaa <- read.csv("Cytokine.csv", header = T, sep = ",")
genename1 <- intersect(aaa$ID, rownames(dataa))
dataa <- dataa[genename1, ]
aaa <- dplyr::filter(aaa, ID %in% genename1)

risk <- OSdata
rownames(risk) <- risk$ID
sameid <- intersect(colnames(dataa), rownames(risk))
dataa <- dataa[, sameid]
risk <- risk[sameid, ]
hmdat <- as.data.frame(t(dataa))

type <- data.frame(Type = aaa$Type, row.names = aaa$ID)
Typeid <- type$Type

annCol <- data.frame(score = scale(risk$RS),
                     RiskType = risk$Risk_group,
                     row.names = rownames(risk),
                     stringsAsFactors = F)

sampleType <- annCol$RiskType
sigVec <- c()
for (i in colnames(hmdat)) {
  test <- wilcox.test(hmdat[, i] ~ sampleType)
  pvalue <- test$p.value
  Sig <- ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")))
  sigVec <- c(sigVec, paste0(i, Sig))
}
colnames(hmdat) <- sigVec

annRow <- data.frame(Type = factor(Typeid, levels = unique(Typeid)),
                     row.names = colnames(hmdat),
                     stringsAsFactors = F)

Type.col <- brewer.pal(n = length(unique(Typeid)), name = "Paired")
annColors <- list(Type = setNames(Type.col, unique(Typeid)),
                  "RS" = greenred(64), 
                  "RiskType" = setNames(c("red", "blue"), c("High_risk", "Low_risk")))

samorder <- rownames(risk[order(risk$RS), ])

indata <- t(hmdat)
indata <- indata[, colSums(indata) > 0]
plotdata <- standarize.fun(indata, halfwidth = 2)

pheatmap(
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
  gaps_col = table(annCol$RiskType)[2],
  gaps_row = cumsum(table(annRow$Type)),
  cellwidth = 1.2,
  cellheight = 10,
  filename = "FigureB_immune_heatmap.pdf"
)
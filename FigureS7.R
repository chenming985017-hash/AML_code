rm(list = ls())
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(tidyverse)

exp <- readRDS(file = "GDSC2_Expr (RMA Normalized and Log Transformed).rds")
dim(exp)

drug <- readRDS("GDSC2_Res.rds")
drug <- exp(drug)
dim(drug)

identical(rownames(drug), colnames(exp))

outTab <- read.csv("merge.csv", row.names = 1, check.names = FALSE)
load("RS_Risk_clin.RDATA")
clin <- OSdata

outTab <- outTab[, clin$ID]
outTab <- as.matrix(outTab)
calcPhenotype(trainingExprData = exp,
              trainingPtype = drug,
              testExprData = outTab,
              batchCorrect = 'eb',
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              printOutput = TRUE,
              removeLowVaringGenesFrom = 'rawData')

drugdata <- read.csv("calcPhenotype_Output/DrugPredictions.csv", header = T, sep = ",")
colnames(drugdata)[1] <- "ID"
rownames(clin) <- clin$ID
drugdata$Risk_group <- clin[drugdata$ID, "Risk_group"]
drugdata <- dplyr::select(drugdata, ID, Risk_group, everything())

Index <- colnames(drugdata)[3:200]

drugdata$Risk_group <- factor(drugdata$Risk_group, levels = c('High_risk', 'Low_risk'))
my_comparisons <- list(c("High_risk", "Low_risk"))

str(drugdata)
for (i in 1:198) {
  p <- ggboxplot(drugdata, x = "Risk_group", y = Index[i],
                 fill = "Risk_group", legend = F, palette = c("#E7B800", "#00AFBB"), bxp.errorbar = T) +
    theme(legend.position = 'none') +
    ylab(label = paste0(Index[i], " IC50")) +
    xlab(label = "Risk_group") +
    stat_compare_means(comparisons = my_comparisons, method = 't.test', aes(label = ..p.signif..))
  ggsave(p, filename = paste(Index[i], "IC50_difference.pdf", sep = "_"), width = 5.6, height = 4.22)
  print(i)
}
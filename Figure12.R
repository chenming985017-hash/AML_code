library(pheatmap)
library(tidyverse)
library(plyr)
library(ggpubr)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

outTab <- read.csv("merge.csv", row.names = 1, check.names = FALSE)
load("RS_Risk_clin.RDATA")
clin <- OSdata
outTab <- as.data.frame(outTab)

TIDE <- outTab[, clin$ID]
TIDE <- sweep(TIDE, 2, apply(TIDE, 2, median, na.rm=T))
TIDE <- sweep(TIDE, 1, apply(TIDE, 1, median, na.rm=T))
write.table(TIDE, "TIDE_input.self_subtract", sep = "\t", row.names = T, col.names = NA, quote = F)

rm(list = ls())
outTab <- read.csv("merge.csv", row.names = 1, check.names = FALSE)
load("RS_Risk_clin.RDATA")
clin <- OSdata

TIDE <- outTab[, clin$ID]
TIDE <- as.data.frame(TIDE)
ann <- clin
rownames(ann) <- ann$ID
TIDE.res <- read.csv("TIDE_output.csv", header = T, row.names = 1, check.names = F, stringsAsFactors = F)
ann$Response <- TIDE.res[rownames(ann), "Responder"]
ann$TIDE.score <- TIDE.res[rownames(ann), "TIDE"]
print(table(ann$Response, ann$Risk_group))
print(fisher.test(table(ann$Response, ann$Risk_group)))

trait <- "Response"
rt <- ann

bioCol <- c("#0066FF", "#FF0000", "#FF9900", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
bioCol <- bioCol[1:length(unique(rt[, trait]))]

rt1 <- rt[, c(trait, "Risk_group")]
rt1 <- na.omit(rt1)
colnames(rt1) <- c("trait", "Risk_group")
table(rt1$trait)
df <- as.data.frame(table(rt1))
df <- ddply(df, .(Risk_group), transform, percent = Freq/sum(Freq) * 100)
df <- ddply(df, .(Risk_group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label <- paste0(sprintf("%.0f", df$percent), "%")
df$Risk_group <- factor(df$Risk_group, levels = c("Low_risk", "High_risk"))

p <- ggplot(df, aes(x = factor(Risk_group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values = bioCol) +
  xlab("score") + ylab("Percent weight") + guides(fill = guide_legend(title = trait)) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_bw()

pdf(file = paste0(trait, "_TIDE_barplot.pdf"), width = 4, height = 5)
print(p)
dev.off()

colnames(TIDE.res)
index <- "Dysfunction"
rt[, index] <- TIDE.res[rownames(rt), index]
rt2 <- rt[, c("Risk_group", index)]
rt2 <- na.omit(rt2)
colnames(rt2) <- c("Risk_group", "score")
type <- levels(factor(rt2[, "Risk_group"]))
rt2$Risk_group <- factor(rt2$Risk_group, levels = type)
comp <- combn(type, 2)
my_comparisons <- list()
for (i in 1:ncol(comp)) {my_comparisons[[i]] <- comp[, i]}

boxplot <- ggboxplot(rt2, x = "Risk_group", y = "score", fill = "Risk_group",
                     xlab = "Risk_group",
                     ylab = index,
                     legend.title = "Risk_group",
                     palette = bioCol) + 
  stat_compare_means(comparisons = my_comparisons)

pdf(file = paste0(index, "_TIDE_boxplot.pdf"), width = 5, height = 4.5)
print(boxplot)
dev.off()

colnames(TIDE.res)
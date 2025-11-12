##A
rm(list = ls())
library(limma)
library(survival)
library(survminer)
library(igraph)
library(psych)
library(reshape2)
library(RColorBrewer)
library(readxl)
library(tidyverse)

expData <- read.csv("prognostic_gene_expression_matrix.csv", row.names = 1, check.names = FALSE)
ribo_genes <- read_excel("Ribosis.xlsx")[[1]]
common_genes <- intersect(ribo_genes, rownames(expData))
geneExp <- expData[common_genes, ]
write.table(geneExp, "GeneExp.txt", sep = "\t", quote = FALSE, col.names = NA)

clinical <- read.csv("merged_GSE_clinical.csv", row.names = "sample_id")
common_samples <- intersect(colnames(geneExp), rownames(clinical))
data <- t(geneExp[, common_samples])
cli <- clinical[common_samples, c("futime", "fustat")]
rt <- cbind(cli, data)
rt <- rt %>% 
  filter(futime > 0) %>% 
  mutate(futime = futime / 365)

outTab=data.frame()
km=c()
gene=c()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
  data=rt[,c("futime", "fustat", i)]
  colnames(data)=c("futime", "fustat", "gene")
  res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
  res.cat=surv_categorize(res.cut)
  fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
  diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
  pValue=1-pchisq(diff$chisq, df=1)
  km=c(km, pValue)
  if(pValue<0.05){
    gene=c(gene,i)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
  }
}

outTab=cbind(outTab, km)
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

Coxfile <- "uniCox.txt"
gene.group <- read.table("uniCox.txt",header=T,sep="\t")
gene.group=dplyr::filter(gene.group,pvalue < 0.05)
gene.group=as.data.frame(gene.group[,1]) 
colnames(gene.group) <- c('id')
gene.group$group<-"RRGs"
gene.exp <- read.csv("prognostic_gene_expression_matrix.csv", row.names = 1, check.names = FALSE)
gene.cox <- read.table(Coxfile,header=T,sep="\t")

genelist <- intersect(gene.group$id, gene.cox$id)
genelist <- intersect(genelist, rownames(gene.exp))
gene.group <- gene.group[match(genelist,gene.group$id),]
gene.group <- gene.group[order(gene.group$group),]
gene.exp <- gene.exp[match(gene.group$id,rownames(gene.exp)),]
gene.cox <- gene.cox[match(gene.group$id,gene.cox$id),]

gene.cor <- corr.test(t(gene.exp))
gene.cor.cor <- gene.cor$r
gene.cor.pvalue <- gene.cor$p
gene.cor.cor[upper.tri(gene.cor.cor)] = NA
gene.cor.pvalue[upper.tri(gene.cor.pvalue)] = NA
gene.cor.cor.melt <- melt(gene.cor.cor)
gene.cor.pvalue.melt <- melt(gene.cor.pvalue)
gene.melt <- data.frame(from = gene.cor.cor.melt$Var2,to=gene.cor.cor.melt$Var1,cor=gene.cor.cor.melt$value,pvalue=gene.cor.pvalue.melt$value)
gene.melt <- gene.melt[gene.melt$from!=gene.melt$to&!is.na(gene.melt$pvalue),,drop=F]
gene.edge <- gene.melt[gene.melt$pvalue<0.05,,drop=F]
gene.edge$color <- ifelse(gene.edge$cor>0,'pink','#6495ED')
gene.edge$weight <- abs(gene.edge$cor)*6

gene.node <- gene.group
group.color <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(gene.node$group)))
gene.node$color <- group.color[as.numeric(as.factor(gene.node$group))]
gene.node$shape <- "circle"
gene.node$frame <- ifelse(gene.cox$HR>1,'#FFFF00',"#4169E1")
gene.node$pvalue <- gene.cox$pvalue
pvalue.breaks <- c(0,0.0001,0.001,0.01,0.05,1)
pvalue.size <- c(16,14,12,10,8)
cutpvalue <- cut(gene.node$pvalue,breaks=pvalue.breaks)
gene.node$size <- pvalue.size[as.numeric(cutpvalue)]
nodefile <- "network.node.txt"
edgefile <- "network.edge.txt"
write.table(gene.node, nodefile, sep="\t", col.names=T, row.names=F, quote=F)
write.table(gene.edge, edgefile, sep="\t", col.names=T, row.names=F, quote=F)

node = read.table(nodefile, header=T, sep="\t", comment.char="")
edge = read.table(edgefile, header=T, sep="\t", comment.char="")

g = graph.data.frame(edge,directed = FALSE)
node = node[match(names(components(g)$membership),node$id),]

if(!is.na(match('color',colnames(node)))) V(g)$color = node$color
if(!is.na(match('size',colnames(node)))) V(g)$size = node$size
if(!is.na(match('shape',colnames(node)))) V(g)$shape = node$shape
if(!is.na(match('frame',colnames(node)))) V(g)$frame = node$frame

pdf(file="新network.pdf", width=15, height=12)
par(mar=c(0,0,0,0))
layout(matrix(c(1,1,4,2,3,4),nc=2),height=c(4,4,2),width=c(8,3))

coord = layout_in_circle(g)
degree.x = acos(coord[,1])
degree.y = asin(coord[,2])
degree.alpha = c()
for(i in 1:length(degree.x)){
  if(degree.y[i]<0) degree.alpha=c(degree.alpha,2*pi-degree.x[i]) else degree.alpha=c(degree.alpha,degree.x[i])
}
degree.cut.group = (0:8)/4*pi
degree.cut.group[1] = -0.0001
degree.cut = cut(degree.alpha,degree.cut.group)
degree.degree = c(-pi/4,-pi/4,-pi/2,-pi/2,pi/2,pi/2,pi/2,pi/4)
degree = degree.degree[as.numeric(degree.cut)]

values <- lapply(node$id,function(x)c(1,1))
V(g)$pie.color = lapply(1:nrow(node),function(x)c(node$color[x],node$frame[x]))
V(g)$frame = NA 

plot(g,layout=layout_in_circle,vertex.shape="pie",vertex.pie=values,
     edge.width = E(g)$weight,edge.arrow.size=0,
     vertex.label.color=V(g)$color,vertex.frame.color=V(g)$frame,edge.color=E(g)$color,
     vertex.label.cex=2,vertex.label.font=2,vertex.size=V(g)$size,edge.curved=0.4,
     vertex.color=V(g)$color,vertex.label.dist=1,vertex.label.degree=degree)

par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
groupinfo = unique(data.frame(group=node$group,color=node$color))
legend("left",legend=groupinfo$group,col=groupinfo$color,pch=16,bty="n",cex=3)

par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
legend("left",legend=c('Risk factors','Favorable factors'),col=c('#FFFF00',"#4169E1"),pch=16,bty="n",cex=2.5)

par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",axes=F,ylab="")
legend("top",legend=c('Postive correlation with P<0.05','Negative correlation with P<0.05'),lty=1,lwd=4,col=c('pink','#6495ED'),bty="n",cex=2.2)
legend('bottom',legend=c(0.0001,0.001,0.01,0.05,1),pch=16,pt.cex=c(1.6,1.4,1.2,1,0.8)*6,bty="n",ncol=5,cex=2.2,col="black",title="Cox test, pvalue")
dev.off()
##B&C
rm(list = ls())
library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(tidyverse)

gene_exp <- read.csv("prognostic_gene_expression_matrix.csv", row.names = 1, check.names = FALSE)
uni_cox <- read.table("uniCox.txt", header = TRUE, sep = "\t")

sig_genes <- uni_cox$id[uni_cox$pvalue < 0.05]
AAA <- intersect(sig_genes, rownames(gene_exp))
geneexpr <- gene_exp[AAA, ]

write.table(geneexpr, "geneexpr.txt", quote = FALSE, row.names = TRUE, sep = "\t")
save(geneexpr, file = "geneexpr.RDATA")

data <- as.matrix(geneexpr)
celltype <- "RRGs"
maxK <- 9

results <- ConsensusClusterPlus(
  data,
  maxK = maxK,
  reps = 50,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "pam",
  distance = "euclidean",
  seed = 123456,
  title = celltype,
  plot = "png"
)

clusterNum <- 2
cluster <- results[[clusterNum]][["consensusClass"]]
cluster <- as.data.frame(cluster)
colnames(cluster) <- c("cluster")
letter <- c("A", "B", "C", "D", "E", "F", "G")
uniqClu <- levels(factor(cluster$cluster))
cluster$cluster <- letter[match(cluster$cluster, uniqClu)]
clusterOut <- rbind(ID = rownames(cluster), cluster)
write.table(clusterOut, file = "Cluster.txt", sep = "\t", quote = FALSE, col.names = FALSE)

cli <- read.csv("merged_GSE_clinical.csv", row.names = "sample_id", check.names = FALSE)
cli <- cli[, c("futime", "fustat")]
cli$futime <- cli$futime / 365

sameSample <- intersect(rownames(cluster), rownames(cli))
rt <- cbind(cli[sameSample, , drop = FALSE], cluster[sameSample, , drop = FALSE])
rt <- dplyr::filter(rt, futime > 0)

length <- length(levels(factor(rt$cluster)))
diff <- survdiff(Surv(futime, fustat) ~ cluster, data = rt)
pValue <- 1 - pchisq(diff$chisq, df = length - 1)
if (pValue < 0.001) {
  pValue <- "p<0.001"
} else {
  pValue <- paste0("p=", sprintf("%.03f", pValue))
}
fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)

bioCol <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
bioCol <- bioCol[1:length]
surPlot <- ggsurvplot(
  fit,
  data = rt,
  conf.int = FALSE,
  pval = pValue,
  pval.size = 6,
  legend.title = "cluster",
  legend.labs = levels(factor(rt[,"cluster"])),
  legend = c(0.8, 0.8),
  font.legend = 10,
  xlab = "Time(years)",
  break.time.by = 2,
  palette = bioCol,
  surv.median.line = "hv",
  risk.table = TRUE,
  cumevents = FALSE,
  risk.table.height = .25
)

pdf(file = "fig5C_survival.pdf", onefile = FALSE, width = 7, height = 6)
print(surPlot)
dev.off()
##D
rm(list = ls())
library(limma)
library(reshape2)
library(ggpubr)
library(tidyverse)

expFile <- "geneexpr.txt"
geneCluFile <- "Cluster.txt"

rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
dimnames <- list(rownames(rt), colnames(rt))
data <- matrix(as.numeric(as.matrix(rt)), nrow = nrow(rt), dimnames = dimnames)
data <- avereps(data)
data <- t(data)

geneClu <- read.table(geneCluFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

sameSample <- intersect(row.names(data), row.names(geneClu))
expClu <- cbind(data[sameSample, , drop = FALSE], geneClu[sameSample, , drop = FALSE])

data_melt <- melt(expClu, id.vars = c("cluster"))
colnames(data_melt) <- c("Cluster", "Gene", "Expression")

bioCol <- c("#0069A7", "#FBA784", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
bioCol <- bioCol[1:length(levels(factor(data_melt[,"Cluster"])))]

data_melt$Gene <- str_replace_all(data_melt$Gene, "REACTOME_", "")

p <- ggboxplot(data_melt, 
               x = "Gene", 
               y = "Expression", 
               fill = "Cluster",
               color = "black",
               ylab = "Gene expression",
               xlab = "",
               legend.title = "Cluster",
               palette = bioCol,
               width = 1) 

p <- p + rotate_x_text(60)

p1 <- p + stat_compare_means(
  aes(group = Cluster),
  symnum.args = list(
    cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
    symbols = c("***", "**", "*", "ns")
  ),
  label = "p.signif",
  color = "black"
)

pdf(file = "fig5D_boxplot.pdf", width = 15, height = 6)
print(p1)
dev.off()
##E
```r
rm(list = ls())
library(pheatmap)
library(tidyverse)

expFile <- "geneexpr.txt"
clusterFile <- "Cluster.txt"
cliFile <- "merged_GSE_clinical.csv"

exp <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
exp <- t(exp)
cluster <- read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

sameSample <- intersect(row.names(exp), row.names(cluster))
exp <- exp[sameSample, , drop = FALSE]
cluster <- cluster[sameSample, , drop = FALSE]
expCluster <- cbind(exp, cluster)

Project <- "AML"
expCluster <- cbind(expCluster, Project)

cli <- read.csv(cliFile, row.names = "sample_id", check.names = FALSE)
sameSample <- intersect(row.names(expCluster), row.names(cli))
expCluster <- expCluster[sameSample, , drop = FALSE]
cli <- cli[sameSample, , drop = FALSE]
data <- cbind(expCluster, cli)

data$fustat[data$fustat == 1] <- "Dead"
data$fustat[data$fustat == 0] <- "Alive"

data <- data[order(data$cluster), ]
Type <- data[, ((ncol(exp) + 1):ncol(data))]
data <- t(data[, 1:ncol(exp)])

bioCol <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
ann_colors <- list()
CluCol <- bioCol[1:length(levels(factor(Type$cluster)))]
names(CluCol) <- levels(factor(Type$cluster))
ann_colors[["cluster"]] <- CluCol

pdf("fig5E_heatmap.pdf", height = 8, width = 10)
pheatmap(data,
         annotation = Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("purple",5), "white", rep("orange",5)))(100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "row",
         show_colnames = FALSE,
         fontsize = 6,
         fontsize_row = 6,
         fontsize_col = 6)
dev.off()

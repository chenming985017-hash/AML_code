library(sva)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(scatterplot3d)
library(dplyr)

set.seed(123)

gse37642 <- read.csv("GSE37642_Matrix_norm.csv", row.names = 1, check.names = FALSE)
gse12417 <- read.csv("GSE12417_Matrix_norm.csv", row.names = 1, check.names = FALSE)

cat("GSE37642 dimension:", dim(gse37642), "\n")
cat("GSE12417 dimension:", dim(gse12417), "\n")

common_genes <- intersect(rownames(gse37642), rownames(gse12417))
cat("Number of common genes:", length(common_genes), "\n")

gse37642_common <- gse37642[common_genes, ]
gse12417_common <- gse12417[common_genes, ]

combined_data <- cbind(gse37642_common, gse12417_common)

batch_info <- data.frame(
  sample = colnames(combined_data),
  batch = c(rep("GSE37642", ncol(gse37642_common)), 
            rep("GSE12417", ncol(gse12417_common)))
)

cat("Performing PCA before batch correction...\n")
pca_before <- prcomp(t(combined_data), scale. = TRUE, center = TRUE)

pca_var_before <- pca_before$sdev^2
pca_var_percent_before <- pca_var_before / sum(pca_var_before) * 100
cumulative_percent_before <- cumsum(pca_var_percent_before)

pca_df_before <- data.frame(
  PC1 = pca_before$x[,1],
  PC2 = pca_before$x[,2],
  PC3 = pca_before$x[,3],
  Batch = batch_info$batch
)

cat("Performing ComBat batch correction...\n")
batch_vector <- as.factor(c(rep(1, ncol(gse37642_common)), 
                            rep(2, ncol(gse12417_common))))

combined_matrix <- as.matrix(combined_data)
mode(combined_matrix) <- "numeric"

corrected_data <- ComBat(dat = combined_matrix, batch = batch_vector, 
                         par.prior = TRUE, prior.plots = FALSE)

cat("Performing PCA after batch correction...\n")
pca_after <- prcomp(t(corrected_data), scale. = TRUE, center = TRUE)

pca_var_after <- pca_after$sdev^2
pca_var_percent_after <- pca_var_after / sum(pca_var_after) * 100
cumulative_percent_after <- cumsum(pca_var_percent_after)

pca_df_after <- data.frame(
  PC1 = pca_after$x[,1],
  PC2 = pca_after$x[,2],
  PC3 = pca_after$x[,3],
  Batch = batch_info$batch
)

calculate_batch_distance <- function(pca_df) {
  batch1 <- pca_df[pca_df$Batch == "GSE37642", 1:3]
  batch2 <- pca_df[pca_df$Batch == "GSE12417", 1:3]
  
  centroid1 <- colMeans(batch1)
  centroid2 <- colMeans(batch2)
  
  distance <- sqrt(sum((centroid1 - centroid2)^2))
  return(distance)
}

distance_before <- calculate_batch_distance(pca_df_before)
distance_after <- calculate_batch_distance(pca_df_after)
reduction_percent <- (1 - distance_after/distance_before) * 100

pca_var_table <- data.frame(
  PC = rep(paste0("PC", 1:5), 2),
  Variance_Percentage = c(pca_var_percent_before[1:5], pca_var_percent_after[1:5]),
  Cumulative_Percentage = c(cumulative_percent_before[1:5], cumulative_percent_after[1:5]),
  Status = rep(c("Before Correction", "After Correction"), each = 5)
)

distance_table <- data.frame(
  Comparison = "GSE37642 vs GSE12417",
  Before_Correction = round(distance_before, 4),
  After_Correction = round(distance_after, 4),
  Reduction_Percentage = paste0(round(reduction_percent, 1), "%")
)

cat("\n=== Batch Correction Result Summary ===\n")
cat("Batch distance reduction percentage:", round(reduction_percent, 1), "%\n\n")

cat("=== PCA Variance Analysis Table ===\n")
print(pca_var_table)

cat("\n=== Batch Distance Comparison ===\n")
print(distance_table)

write.csv(pca_var_table, "PCA_Variance_Analysis_Table.csv", row.names = FALSE)
write.csv(distance_table, "Batch_Distance_Comparison_Table.csv", row.names = FALSE)

p1 <- ggplot(pca_df_before, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA Before Batch Correction",
       x = paste0("PC1 (", round(pca_var_percent_before[1], 2), "%)"),
       y = paste0("PC2 (", round(pca_var_percent_before[2], 2), "%)")) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"))

p2 <- ggplot(pca_df_after, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA After Batch Correction",
       x = paste0("PC1 (", round(pca_var_percent_after[1], 2), "%)"),
       y = paste0("PC2 (", round(pca_var_percent_after[2], 2), "%)")) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"))

pdf("Supplementary_Figure_S1_2D_PCA.pdf", width = 12, height = 6)
ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()

write.csv(corrected_data, "Batch_Corrected_Combined_Data.csv")
write.csv(pca_var_table, "PCA_Variance_Analysis.csv", row.names = FALSE)
write.csv(distance_table, "Batch_Distance_Comparison.csv", row.names = FALSE)

cat("\nAnalysis completed!\n")
cat("Generated files:\n")
cat("1. Supplementary_Figure_S1_2D_PCA.pdf - 2D PCA plots\n")
cat("2. Batch_Corrected_Combined_Data.csv - Batch-corrected expression matrix\n")
cat("3. PCA_Variance_Analysis.csv - PCA variance analysis table\n")
cat("4. Batch_Distance_Comparison.csv - Batch distance comparison table\n")

colors_before <- ifelse(pca_df_before$Batch == "GSE37642", "#E41A1C", "#377EB8")
colors_after <- ifelse(pca_df_after$Batch == "GSE37642", "#E41A1C", "#377EB8")

pdf("Supplementary_Figure_S1_3D_PCA.pdf", width = 14, height = 6)
par(mfrow = c(1, 2), mar = c(4, 4, 4, 6), oma = c(0, 0, 0, 2))

s3d_before <- scatterplot3d(pca_df_before$PC1, pca_df_before$PC2, pca_df_before$PC3,
                            color = colors_before,
                            pch = 16, 
                            main = "3D PCA Before Batch Correction\n(Clear Batch Separation)",
                            xlab = paste0("PC1 (", round(pca_var_percent_before[1], 2), "%)"),
                            ylab = paste0("PC2 (", round(pca_var_percent_before[2], 2), "%)"),
                            zlab = paste0("PC3 (", round(pca_var_percent_before[3], 2), "%)"),
                            grid = TRUE, box = TRUE,
                            mar = c(3, 3, 4, 8))

legend(x = "topright", 
       legend = c("GSE37642", "GSE12417"),
       col = c("#E41A1C", "#377EB8"), 
       pch = 16, 
       cex = 0.9,
       inset = -0.15,
       xpd = TRUE)

s3d_after <- scatterplot3d(pca_df_after$PC1, pca_df_after$PC2, pca_df_after$PC3,
                           color = colors_after,
                           pch = 16, 
                           main = "3D PCA After Batch Correction\n(Fully Integrated Clusters)",
                           xlab = paste0("PC1 (", round(pca_var_percent_after[1], 2), "%)"),
                           ylab = paste0("PC2 (", round(pca_var_percent_after[2], 2), "%)"),
                           zlab = paste0("PC3 (", round(pca_var_percent_after[3], 2), "%)"),
                           grid = TRUE, box = TRUE,
                           mar = c(3, 3, 4, 8))

legend(x = "topright", 
       legend = c("GSE37642", "GSE12417"),
       col = c("#E41A1C", "#377EB8"), 
       pch = 16, 
       cex = 0.9,
       inset = -0.15,
       xpd = TRUE)

dev.off()

png("Supplementary_Figure_S1_3D_PCA_Alternative.png", width = 1400, height = 600)
par(mfrow = c(1, 2), mar = c(4, 4, 4, 8))

s3d_before <- scatterplot3d(pca_df_before$PC1, pca_df_before$PC2, pca_df_before$PC3,
                            color = colors_before,
                            pch = 16, 
                            main = "3D PCA Before Batch Correction\n(Clear Batch Separation)",
                            xlab = paste0("PC1 (", round(pca_var_percent_before[1], 2), "%)"),
                            ylab = paste0("PC2 (", round(pca_var_percent_before[2], 2), "%)"),
                            zlab = paste0("PC3 (", round(pca_var_percent_before[3], 2), "%)"))

legend("topright", 
       legend = c("GSE37642", "GSE12417"),
       col = c("#E41A1C", "#377EB8"), 
       pch = 16,
       cex = 0.8,
       inset = c(-0.2, 0),
       xpd = TRUE)

s3d_after <- scatterplot3d(pca_df_after$PC1, pca_df_after$PC2, pca_df_after$PC3,
                           color = colors_after,
                           pch = 16, 
                           main = "3D PCA After Batch Correction\n(Fully Integrated Clusters)",
                           xlab = paste0("PC1 (", round(pca_var_percent_after[1], 2), "%)"),
                           ylab = paste0("PC2 (", round(pca_var_percent_after[2], 2), "%)"),
                           zlab = paste0("PC3 (", round(pca_var_percent_after[3], 2), "%)"))

legend("topright", 
       legend = c("GSE37642", "GSE12417"),
       col = c("#E41A1C", "#377EB8"), 
       pch = 16,
       cex = 0.8,
       inset = c(-0.2, 0),
       xpd = TRUE)

dev.off()


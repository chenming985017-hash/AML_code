library(survival)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)

# Load data
expression_matrix_TCGA <- read.csv("LAML_TPM.csv", row.names = 1)
clinical_data_TCGA <- read.csv("LAML_clinical.csv")
gene_list <- read_excel("Ribosis.xlsx")[[1]]
gene_list <- as.character(gene_list)

colnames(expression_matrix_TCGA) <- gsub("\\.", "-", colnames(expression_matrix_TCGA))

# Function to filter prognostic genes
filter_prognostic_genes <- function(expression_matrix, clinical_data, gene_list) {
  results <- data.frame(
    Gene = character(),
    HR = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  for (gene in gene_list) {
    if (gene %in% rownames(expression_matrix)) {
      gene_expression <- as.numeric(expression_matrix[gene, ])
      data <- clinical_data %>% mutate(GeneExpression = gene_expression)
      cox_model <- coxph(Surv(OS.time, OS) ~ GeneExpression, data = data)
      cox_summary <- summary(cox_model)
      results <- rbind(results, data.frame(
        Gene = gene,
        HR = cox_summary$coefficients[1, 2],
        CI_lower = cox_summary$conf.int[1, 3],
        CI_upper = cox_summary$conf.int[1, 4],
        P_value = cox_summary$coefficients[1, 5]
      ))
    }
  }
  significant_genes <- results %>% filter(P_value < 0.05)
  return(significant_genes)
}

# Analyze individual datasets
significant_genes_TCGA <- filter_prognostic_genes(expression_matrix_TCGA, clinical_data_TCGA, gene_list)

GSE37642_expression <- read.csv("GSE37642_Matrix_norm.csv", row.names = 1, check.names = FALSE)
GSE37642_clinical <- read.csv("GSE37642-GPL570_cli.csv", check.names = FALSE)
significant_genes_GSE37642 <- filter_prognostic_genes(GSE37642_expression, GSE37642_clinical, gene_list)

GSE12417_expression <- read.csv("GSE12417_Matrix_norm.csv", row.names = 1, check.names = FALSE)
GSE12417_clinical <- read.csv("GSE12417-GPL570_cli.csv", check.names = FALSE)
significant_genes_GSE12417 <- filter_prognostic_genes(GSE12417_expression, GSE12417_clinical, gene_list)

# Identify common prognostic genes across all datasets
common_genes <- intersect(significant_genes_TCGA$Gene, intersect(significant_genes_GSE37642$Gene, significant_genes_GSE12417$Gene))

# Prepare merged expression matrix (GSE37642 + GSE12417)
colnames(GSE37642_expression) <- paste0("GSE37642_", colnames(GSE37642_expression))
colnames(GSE12417_expression) <- paste0("GSE12417_", colnames(GSE12417_expression))
common_genes_expr <- intersect(rownames(GSE37642_expression), rownames(GSE12417_expression))
GSE37642_common <- GSE37642_expression[common_genes_expr, ]
GSE12417_common <- GSE12417_expression[common_genes_expr, ]
merged_expression <- cbind(GSE37642_common, GSE12417_common)

# Prepare merged clinical data
GSE37642_clinical$sample_id <- paste0("GSE37642_", GSE37642_clinical$sample_id)
GSE12417_clinical$sample_id <- paste0("GSE12417_", GSE12417_clinical$sample_id)
merged_clinical <- bind_rows(GSE37642_clinical, GSE12417_clinical)

# Analyze merged dataset
significant_genes_merged <- filter_prognostic_genes(merged_expression, merged_clinical, common_genes)
prognostic_expression_matrix <- merged_expression[significant_genes_merged$Gene, ]

# Save results
write.csv(significant_genes_TCGA, "significant_genes_TCGA.csv", row.names = FALSE)
write.csv(significant_genes_GSE37642, "significant_genes_GSE37642.csv", row.names = FALSE)
write.csv(significant_genes_GSE12417, "significant_genes_GSE12417.csv", row.names = FALSE)
write.csv(common_genes, "common_prognostic_genes.csv", row.names = FALSE)
write.csv(merged_expression, "merged_GSE_expression.csv")
write.csv(merged_clinical, "merged_GSE_clinical.csv", row.names = FALSE)
write.csv(significant_genes_merged, "significant_genes_merged.csv", row.names = FALSE)
write.csv(prognostic_expression_matrix, "prognostic_gene_expression_matrix.csv")
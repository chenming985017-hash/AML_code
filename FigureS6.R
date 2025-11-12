# =====================================================
# SHAP Analysis for AML Ribosome Biogenesis Signature
# Validation Sets - Supplementary Materials
# =====================================================

# Load required packages
library(randomForestSRC)
library(survival)
library(ggplot2)
library(kernelshap)
library(shapviz)
library(dplyr)
library(patchwork)
library(RColorBrewer)

# Set environment
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

# =====================================================
# 1. LOAD TRAINED MODEL AND BACKGROUND DATA
# =====================================================

# Load trained model
fit <- readRDS("rsf_model_GSE37642.rds")

# Load pre-processed data and prepare background data
load("filtered_gene_expression_clinical.RDATA")
OSdata <- dplyr::filter(OSdata, OS.time > 0)
train_data <- OSdata %>%
  dplyr::filter(cohort == "GSE37642") %>%
  dplyr::select("ID", "OS.time", "OS", badgene, goodgene) %>%
  na.omit()
est_dd <- train_data[, c('OS.time', 'OS', colnames(train_data)[-c(1:3)])]

set.seed(123)
bg_sample_size <- min(50, nrow(est_dd))
bg_indices <- sample(nrow(est_dd), bg_sample_size)
bg_data <- est_dd[bg_indices, ] %>% select(-OS.time, -OS)

# Prediction function for RSF model
pred_fun <- function(object, newdata) {
  pred_result <- predict(object, newdata = newdata)
  risk_scores <- as.numeric(pred_result$predicted)
  return(risk_scores)
}

# =====================================================
# 2. PREPARE VALIDATION DATASETS
# =====================================================

# Function to standardize and prepare validation data
prepare_val_data <- function(data, features) {
  data <- data[, c('OS.time', 'OS', features)]
  data[, -c(1:3)] <- scale(data[, -c(1:3)])
  features <- data %>% select(-OS.time, -OS)
  return(features)
}

# Validation datasets list
val_datasets <- list(
  GSE12417 = {
    OSdata %>%
      dplyr::filter(cohort == "GSE12417") %>%
      dplyr::select("ID", "OS.time", "OS", badgene, goodgene) %>%
      na.omit()
  },
  Merge = {
    OSdata %>%
      dplyr::select("ID", "OS.time", "OS", badgene, goodgene) %>%
      na.omit()
  },
  BeatAML = {
    # Load and process BeatAML data (adjust paths and columns as needed)
    beataml_exp <- read.csv("BeatAML.csv", row.names = 1)
    beataml_cli <- read.csv("BeatAML_clean.csv")
    target_genes <- colnames(bg_data)
    beataml_exp_filtered <- beataml_exp[intersect(target_genes, rownames(beataml_exp)), ]
    beataml_exp_t <- as.data.frame(t(beataml_exp_filtered))
    beataml_exp_t$SampleID <- rownames(beataml_exp_t)
    beataml_data <- beataml_cli %>%
      inner_join(beataml_exp_t, by = "SampleID") %>%
      select(-SampleID)
    beataml_data$OS.time <- beataml_data$OS.time / 365
    beataml_data
  },
  TCGA = {
    # Load and process TCGA data (adjust paths and columns as needed)
    tcga_exp <- read.csv("sorted_expression_matrix.csv", row.names = 1, check.names = FALSE)
    tcga_cli <- read.csv("TCGA.csv")
    if (any(duplicated(rownames(tcga_exp)))) {
      tcga_exp_agg <- aggregate(. ~ rownames(tcga_exp), data = tcga_exp, mean)
      rownames(tcga_exp_agg) <- tcga_exp_agg[,1]
      tcga_exp <- tcga_exp_agg[-1]
    }
    colnames(tcga_exp) <- substr(colnames(tcga_exp), 1, 12)
    target_genes <- colnames(bg_data)
    tcga_exp_filtered <- tcga_exp[intersect(target_genes, rownames(tcga_exp)), ]
    tcga_exp_t <- as.data.frame(t(tcga_exp_filtered))
    tcga_exp_t$sample_id <- rownames(tcga_exp_t)
    tcga_data <- tcga_cli %>%
      inner_join(tcga_exp_t, by = "sample_id") %>%
      select(-sample_id)
    tcga_data
  },
  OHSU = {
    # Load and process OHSU data (adjust paths and columns as needed)
    ohsu_exp <- read.csv("OSHU_final_expression_matrix.csv")
    ohsu_cli <- read.csv("OSHU_clean.csv")
    target_genes <- colnames(bg_data)
    ohsu_exp_filtered <- ohsu_exp[, c("sample_id", intersect(target_genes, colnames(ohsu_exp)))]
    ohsu_data <- ohsu_cli %>%
      inner_join(ohsu_exp_filtered, by = "sample_id") %>%
      select(-sample_id)
    ohsu_data
  }
)

# Standardize and extract features for each validation dataset
val_features_list <- lapply(val_datasets, function(data) {
  prepare_val_data(data, colnames(bg_data))
})

# =====================================================
# 3. SHAP ANALYSIS - VALIDATION SETS
# =====================================================

# Function to compute SHAP for validation sets
compute_val_shap <- function(model, features, background_data, cohort_name) {
  shap_values <- kernelshap(
    model,
    X = features,
    bg_X = background_data,
    pred_fun = pred_fun
  )
  shp <- shapviz(shap_values)
  p_bar <- sv_importance(
    shp,
    kind = "bar",
    max_display = 9,
    fill = brewer.pal(9, "Reds"),
    bar_width = 0.8,
    show_numbers = TRUE,
    number_size = 3.5
  ) + 
    labs(title = paste("SHAP Feature Importance -", cohort_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(p_bar, file = paste0("shap_", cohort_name, "_supplementary.pdf"), width = 9, height = 6)
  return(p_bar)
}

# Run SHAP analysis for all validation sets
val_shap_plots <- lapply(names(val_features_list), function(cohort) {
  compute_val_shap(fit, val_features_list[[cohort]], bg_data, cohort)
})

# Combine all validation plots
combined_val_plots <- wrap_plots(val_shap_plots, ncol = 2) +
  plot_annotation(title = "SHAP Feature Importance - Validation Cohorts (Supplementary)") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

ggsave(combined_val_plots, file = "shap_ALL_validation_supplementary.pdf", width = 16, height = 12)
cat("Validation sets SHAP analysis completed. Plots saved for Supplementary Materials.\n")
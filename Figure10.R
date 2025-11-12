# =====================================================
# SHAP Analysis for AML Ribosome Biogenesis Signature
# Training Set (GSE37642) - Figure 10
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
# 1. DATA PREPARATION
# =====================================================

# Load pre-processed data
load("filtered_gene_expression_clinical.RDATA")

# Filter and prepare data
OSdata <- dplyr::filter(OSdata, OS.time > 0)
train_data <- OSdata %>%
  dplyr::filter(cohort == "GSE37642") %>%
  dplyr::select("ID", "OS.time", "OS", badgene, goodgene) %>%
  na.omit()

# Standardize data
train_data[, -c(1:3)] <- scale(train_data[, -c(1:3)])

# =====================================================
# 2. MODEL TRAINING
# =====================================================

pre_var <- colnames(train_data)[-c(1:3)]
est_dd <- train_data[, c('OS.time', 'OS', pre_var)]

# Train Random Survival Forest model
set.seed(123)
fit <- rfsrc(Surv(OS.time, OS) ~ ., 
             data = est_dd, 
             ntree = 1000, 
             nodesize = 5,
             splitrule = 'logrank', 
             importance = TRUE, 
             proximity = TRUE, 
             forest = TRUE)

# Save model
saveRDS(fit, "rsf_model_GSE37642.rds")

# =====================================================
# 3. SHAP ANALYSIS FUNCTIONS
# =====================================================

# Prediction function for RSF model
pred_fun <- function(object, newdata) {
  pred_result <- predict(object, newdata = newdata)
  risk_scores <- as.numeric(pred_result$predicted)
  return(risk_scores)
}

# Function to compute and visualize SHAP values for training set
compute_train_shap <- function(model, features, background_data, save_plot = TRUE) {
  
  # Compute SHAP values
  shap_values <- kernelshap(
    model,
    X = features,
    bg_X = background_data,
    pred_fun = pred_fun
  )
  
  # Create shapviz object
  shp <- shapviz(shap_values)
  
  # Feature importance bar plot
  p_bar <- sv_importance(
    shp,
    kind = "bar",
    max_display = 9,
    fill = brewer.pal(9, "Blues"),
    bar_width = 0.8,
    show_numbers = TRUE,
    number_size = 3.5
  ) + 
    labs(title = "SHAP Feature Importance - GSE37642 (Training Set)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Beeswarm plot
  gradient_colors <- c("darkblue", "dodgerblue", "white", "orangered", "#F60")
  p_beeswarm <- sv_importance(
    shp,
    kind = "beeswarm",
    max_display = 9,
    bee_width = 0.3,
    bee_adjust = 0.6
  ) + 
    scale_color_gradientn(
      colors = gradient_colors,
      values = scales::rescale(c(0, 0.3, 0.5, 0.7, 1)),
      name = "Expression Level",
      breaks = c(0, 1),
      labels = c("Low", "High")
    ) +
    labs(title = "SHAP Beeswarm Plot - GSE37642 (Training Set)") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  
  # Combined bar and beeswarm plot
  combined_bar_beeswarm <- p_bar + plot_spacer() + p_beeswarm + 
    plot_layout(widths = c(1.5, 0.3, 1.5)) + 
    plot_annotation(title = "SHAP Analysis of Gene Features - GSE37642 (Training Set)") & 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95")
    )
  
  # Save plots
  if (save_plot) {
    ggsave(p_bar, file = "shap_training_bar.pdf", width = 9, height = 6)
    ggsave(p_beeswarm, file = "shap_training_beeswarm.pdf", width = 9, height = 6)
    ggsave(combined_bar_beeswarm, file = "shap_training_combined.pdf", width = 16, height = 8)
  }
  
  # Waterfall and force plots for representative samples
  risk_scores <- pred_fun(fit, features)
  sample_high_risk <- which.max(risk_scores)
  sample_medium_risk <- which.min(abs(risk_scores - median(risk_scores)))
  sample_low_risk <- which.min(risk_scores)
  
  create_waterfall_force_plot <- function(shap_obj, row_id, sample_type) {
    p_waterfall <- sv_waterfall(shap_obj, row_id = row_id) + 
      ggtitle(paste("Waterfall Plot -", sample_type, "Risk Sample")) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    p_force <- sv_force(shap_obj, row_id = row_id) + 
      ggtitle(paste("Force Plot -", sample_type, "Risk Sample")) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    combined_plot <- p_waterfall + p_force + 
      plot_layout(widths = c(1.3, 1.5)) + 
      plot_annotation(title = paste("SHAP Explanation for", sample_type, "Risk Sample (ID:", row_id, ")")) & 
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    
    return(combined_plot)
  }
  
  high_risk_plot <- create_waterfall_force_plot(shp, sample_high_risk, "High")
  medium_risk_plot <- create_waterfall_force_plot(shp, sample_medium_risk, "Medium")
  low_risk_plot <- create_waterfall_force_plot(shp, sample_low_risk, "Low")
  
  all_samples_combined <- high_risk_plot / medium_risk_plot / low_risk_plot +
    plot_annotation(title = "SHAP Explanations Across Different Risk Levels (GSE37642)") &
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  
  if (save_plot) {
    ggsave(high_risk_plot, file = "shap_high_risk_sample.pdf", width = 16, height = 8)
    ggsave(medium_risk_plot, file = "shap_medium_risk_sample.pdf", width = 16, height = 8)
    ggsave(low_risk_plot, file = "shap_low_risk_sample.pdf", width = 16, height = 8)
    ggsave(all_samples_combined, file = "shap_all_risk_samples.pdf", width = 16, height = 20)
  }
  
  return(list(shap_object = shp, plots = list(bar = p_bar, beeswarm = p_beeswarm, combined = combined_bar_beeswarm, risk_samples = all_samples_combined)))
}

# =====================================================
# 4. PREPARE BACKGROUND DATA
# =====================================================

set.seed(123)
bg_sample_size <- min(50, nrow(est_dd))
bg_indices <- sample(nrow(est_dd), bg_sample_size)
bg_data <- est_dd[bg_indices, ] %>% select(-OS.time, -OS)

# =====================================================
# 5. SHAP ANALYSIS - TRAINING SET (GSE37642)
# =====================================================

train_features <- est_dd %>% select(-OS.time, -OS)
train_shap <- compute_train_shap(fit, train_features, bg_data)

cat("Training set SHAP analysis completed. Plots saved for Figure 10.\n")
library(survival)
library(survminer)
library(timeROC)
library(ranger)
library(dplyr)
library(forestplot)

set.seed(123)

load("RS.RDATA")

rsf_model_GSE37642 <- rs[["GSE37642"]]

cat("GSE37642 model class:\n")
print(class(rsf_model_GSE37642))

predict_risk_scores <- function(model, new_data, train_feature_names) {
  valid_features <- new_data[, train_feature_names, drop = FALSE]
  valid_features <- valid_features[, train_feature_names]
  valid_features <- as.data.frame(sapply(valid_features, as.numeric))
  
  predictions <- predict(model, newdata = valid_features)
  risk_scores <- as.vector(t(predictions$predicted))
  
  return(risk_scores)
}

perform_survival_analysis <- function(data, risk_scores, dataset_name) {
  data$risk_score <- risk_scores
  
  score_data <- data.frame(
    OS.time = data$OS.time / 365,
    OS = data$OS,
    Score = data$risk_score
  )
  
  cutpoint <- surv_cutpoint(
    data = score_data,
    time = "OS.time",
    event = "OS",
    variables = "Score"
  )
  
  cat_data <- surv_categorize(cutpoint)
  data$risk_group <- cat_data$Score
  
  fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = data)
  
  km_plot <- ggsurvplot(
    fit,
    data = data,
    palette = "jco",
    size = 1.3,
    pval = TRUE,
    legend.labs = c("High", "Low"),
    legend.title = "Risk Score",
    xlab = "Time (years)",
    ylab = "Overall survival",
    ggtheme = theme_survminer(
      font.legend = c(14, "plain", "black"),
      font.x = c(14, "plain", "black"),
      font.y = c(14, "plain", "black")
    ),
    break.time.by = 2,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.title = "",
    risk.table.height = 0.25
  )
  
  roc_data <- data.frame(
    OS.time = data$OS.time / 365,
    OS = data$OS,
    Score = data$risk_score
  )
  
  tt <- timeROC(
    T = roc_data$OS.time,
    delta = roc_data$OS,
    marker = roc_data$Score,
    cause = 1,
    weighting = "marginal",
    times = c(1, 3, 5),
    ROC = TRUE,
    iid = TRUE
  )
  
  col <- c("#0073C2FF", "firebrick1", "orange")
  roc_plot <- function() {
    plot(tt, time = 1, title = FALSE, lwd = 1.5, col = col[1])
    plot(tt, time = 3, col = col[2], add = TRUE, title = FALSE, lwd = 1.5)
    plot(tt, time = 5, col = col[3], add = TRUE, title = FALSE, lwd = 1.5)
    legend(
      "bottomright",
      legend = c(
        paste0("1-year AUC = ", round(tt$AUC[1], 3)),
        paste0("3-year AUC = ", round(tt$AUC[2], 3)),
        paste0("5-year AUC = ", round(tt$AUC[3], 3))
      ),
      fill = col[1:3],
      bty = "o",
      cex = 1,
      border = NA
    )
    abline(0, 1, lty = 2, lwd = 0.5)
  }
  
  return(list(
    data = data,
    km_plot = km_plot,
    roc_object = tt,
    roc_plot = roc_plot,
    auc_results = data.frame(
      Time = c("1-year", "3-year", "5-year"),
      AUC = round(tt$AUC, 3)
    )
  ))
}

create_forest_plot <- function(cox_results, title) {
  results_df <- cox_results
  
  if ("Variable" %in% colnames(results_df)) {
    labels <- results_df$Variable
    hr <- results_df$HR
    lower <- results_df$CI_lower
    upper <- results_df$CI_upper
    p_values <- results_df$P_value
    
    tabletext <- cbind(
      c("Variable", as.character(labels)),
      c("HR", sprintf("%.2f", hr)),
      c("95% CI", sprintf("(%.2f-%.2f)", lower, upper)),
      c("P-value", sprintf("%.3f", p_values))
    )
    
    mean <- c(NA, hr)
    lower_ci <- c(NA, lower)
    upper_ci <- c(NA, upper)
    
    pdf(paste0(gsub(" ", "_", title), "_forest_plot.pdf"), width = 10, height = 6)
    forestplot(
      labeltext = tabletext,
      mean = mean,
      lower = lower_ci,
      upper = upper_ci,
      zero = 1,
      boxsize = 0.2,
      graph.pos = 2,
      xticks = c(0.5, 1, 2, 4),
      title = title,
      xlab = "Hazard Ratio"
    )
    dev.off()
  }
}

perform_cox_analysis <- function(data, variables, dataset_name) {
  cox_data <- data
  
  factor_vars <- variables[sapply(cox_data[, variables], function(x) !is.numeric(x))]
  for(var in factor_vars) {
    cox_data[[var]] <- as.factor(cox_data[[var]])
  }
  
  univ_results <- data.frame()
  for(var in variables) {
    formula <- as.formula(paste("Surv(OS.time, OS) ~", var))
    model <- coxph(formula, data = cox_data)
    summary_model <- summary(model)
    
    result_row <- data.frame(
      Variable = var,
      HR = round(exp(coef(model)), 3),
      CI_lower = round(exp(confint(model))[1], 3),
      CI_upper = round(exp(confint(model))[2], 3),
      P_value = round(summary_model$coefficients[5], 4),
      n = model$n,
      Events = model$nevent
    )
    univ_results <- rbind(univ_results, result_row)
  }
  
  multiv_formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(variables, collapse = " + ")))
  multiv_model <- coxph(multiv_formula, data = cox_data)
  summary_multiv <- summary(multiv_model)
  
  multiv_results <- data.frame()
  for(i in 1:length(coef(multiv_model))) {
    var_name <- names(coef(multiv_model))[i]
    result_row <- data.frame(
      Variable = var_name,
      HR = round(exp(coef(multiv_model)[i]), 3),
      CI_lower = round(exp(confint(multiv_model)[i, 1]), 3),
      CI_upper = round(exp(confint(multiv_model)[i, 2]), 3),
      P_value = round(summary_multiv$coefficients[i, 5], 4)
    )
    multiv_results <- rbind(multiv_results, result_row)
  }
  
  create_forest_plot(univ_results, paste(dataset_name, "Univariate Analysis"))
  create_forest_plot(multiv_results, paste(dataset_name, "Multivariate Analysis"))
  
  return(list(
    univariate = univ_results,
    multivariate = multiv_results
  ))
}

analyze_dataset <- function(dataset_name, clinical_file, model, train_features, 
                            cox_variables, output_prefix) {
  
  data <- read.csv(clinical_file)
  
  risk_scores <- predict_risk_scores(model, data, train_features)
  
  survival_results <- perform_survival_analysis(data, risk_scores, dataset_name)
  
  cox_results <- perform_cox_analysis(survival_results$data, cox_variables, dataset_name)
  
  write.csv(survival_results$data, paste0(output_prefix, "_with_risk_scores.csv"), row.names = FALSE)
  write.csv(cox_results$univariate, paste0(output_prefix, "_univariate_cox.csv"), row.names = FALSE)
  write.csv(cox_results$multivariate, paste0(output_prefix, "_multivariate_cox.csv"), row.names = FALSE)
  write.csv(survival_results$auc_results, paste0(output_prefix, "_auc_results.csv"), row.names = FALSE)
  
  pdf(paste0(output_prefix, "_KM_curve.pdf"), width = 8, height = 8)
  print(survival_results$km_plot)
  dev.off()
  
  pdf(paste0(output_prefix, "_ROC_curve.pdf"), width = 8, height = 8)
  survival_results$roc_plot()
  dev.off()
  
  return(list(
    survival = survival_results,
    cox = cox_results
  ))
}

tcga_variables <- c("Age", "ELN_simple", "FLT3_mut", "RUNX1_mut", "TP53_mut", "risk_score")
ohsu_variables <- c("FAB", "Age", "ELN_simple", "FLT3_ITD", "NPM1", "DNMT3A_MUTATION", 
                    "FLT3_MUTATION", "RUNX1_MUTATION", "TP53_MUTATION", "risk_score")
beataml_variables <- c("FAB", "Age", "ELN_simple", "TP53_binary", "FLT3.ITD", "NPM1", "risk_score")

tcga_results <- analyze_dataset(
  dataset_name = "TCGA",
  clinical_file = "TCGA_clinical_data.csv",
  model = rsf_model_GSE37642,
  train_features = train_feature_names,
  cox_variables = tcga_variables,
  output_prefix = "TCGA_validation"
)

ohsu_results <- analyze_dataset(
  dataset_name = "OHSU", 
  clinical_file = "OHSU_clinical_data.csv",
  model = rsf_model_GSE37642,
  train_features = train_feature_names,
  cox_variables = ohsu_variables,
  output_prefix = "OHSU_validation"
)

beataml_results <- analyze_dataset(
  dataset_name = "BeatAML",
  clinical_file = "BeatAML_clinical_data.csv", 
  model = rsf_model_GSE37642,
  train_features = train_feature_names,
  cox_variables = beataml_variables,
  output_prefix = "BeatAML_validation"
)

summary_report <- data.frame(
  Dataset = c("TCGA", "OHSU", "BeatAML"),
  Samples = c(
    nrow(tcga_results$survival$data),
    nrow(ohsu_results$survival$data), 
    nrow(beataml_results$survival$data)
  ),
  Events = c(
    sum(tcga_results$survival$data$OS),
    sum(ohsu_results$survival$data$OS),
    sum(beataml_results$survival$data$OS)
  ),
  AUC_1year = c(
    tcga_results$survival$auc_results$AUC[1],
    ohsu_results$survival$auc_results$AUC[1],
    beataml_results$survival$auc_results$AUC[1]
  ),
  AUC_3year = c(
    tcga_results$survival$auc_results$AUC[2],
    ohsu_results$survival$auc_results$AUC[2],
    beataml_results$survival$auc_results$AUC[2]
  ),
  AUC_5year = c(
    tcga_results$survival$auc_results$AUC[3],
    ohsu_results$survival$auc_results$AUC[3],
    beataml_results$survival$auc_results$AUC[3]
  )
)

write.csv(summary_report, "GSE37642_Model_Validation_Summary.csv", row.names = FALSE)

save(tcga_results, ohsu_results, beataml_results, summary_report, rsf_model_GSE37642,
     file = "Supplementary_Figure_5_Results.RData")

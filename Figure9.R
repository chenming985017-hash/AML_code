# A&B
library(survival)
library(forestplot)
library(dplyr)

# 
data <- read.csv("GSE37642RS.csv", stringsAsFactors = TRUE)

# 
data$Fab <- factor(data$Fab)
data$runx1.runx1t1_fusion <- factor(data$runx1.runx1t1_fusion, levels = c("No", "Yes"))
data$runx1_mutation <- factor(data$runx1_mutation, levels = c("No", "Yes"))

#### 
univ_results <- list()

# 
fab_levels <- levels(data$Fab)
for (level in fab_levels) {
  data$fab_level <- ifelse(data$Fab == level, 1, 0)
  fit <- coxph(Surv(OS.time, OS) ~ fab_level, data = data)
  sfit <- summary(fit)
  univ_results[[paste0("Fab_", level)]] <- data.frame(
    variable = paste0("Fab (", level, ")"),
    hr = sfit$coefficients[1, "exp(coef)"],
    lower = sfit$conf.int[1, "lower .95"],
    upper = sfit$conf.int[1, "upper .95"],
    pval = sfit$coefficients[1, "Pr(>|z|)"]
  )
}

# 
fit <- coxph(Surv(OS.time, OS) ~ age, data = data)
sfit <- summary(fit)
univ_results[["age"]] <- data.frame(
  variable = "age",
  hr = sfit$coefficients[1, "exp(coef)"],
  lower = sfit$conf.int[1, "lower .95"],
  upper = sfit$conf.int[1, "upper .95"],
  pval = sfit$coefficients[1, "Pr(>|z|)"]
)

# 
fit <- coxph(Surv(OS.time, OS) ~ runx1.runx1t1_fusion, data = data)
sfit <- summary(fit)
univ_results[["runx1.runx1t1_fusion"]] <- data.frame(
  variable = "runx1-runx1t1_fusion",
  hr = sfit$coefficients[1, "exp(coef)"],
  lower = sfit$conf.int[1, "lower .95"],
  upper = sfit$conf.int[1, "upper .95"],
  pval = sfit$coefficients[1, "Pr(>|z|)"]
)

# 
fit <- coxph(Surv(OS.time, OS) ~ runx1_mutation, data = data)
sfit <- summary(fit)
univ_results[["runx1_mutation"]] <- data.frame(
  variable = "runx1_mutation",
  hr = sfit$coefficients[1, "exp(coef)"],
  lower = sfit$conf.int[1, "lower .95"],
  upper = sfit$conf.int[1, "upper .95"],
  pval = sfit$coefficients[1, "Pr(>|z|)"]
)

# 
fit <- coxph(Surv(OS.time, OS) ~ RS, data = data)
sfit <- summary(fit)
univ_results[["RS"]] <- data.frame(
  variable = "riskScore",
  hr = sfit$coefficients[1, "exp(coef)"],
  lower = sfit$conf.int[1, "lower .95"],
  upper = sfit$conf.int[1, "upper .95"],
  pval = sfit$coefficients[1, "Pr(>|z|)"]
)

# 
univ_df <- bind_rows(univ_results)

##
fit_multi <- coxph(Surv(OS.time, OS) ~ Fab + age + runx1.runx1t1_fusion + runx1_mutation + RS, data = data)
sfit_multi <- summary(fit_multi)

multi_results <- list()

# 
for (level in fab_levels) {
  coef_name <- paste0("Fab", level)
  if (coef_name %in% rownames(sfit_multi$coefficients)) {
    multi_results[[paste0("Fab_", level)]] <- data.frame(
      variable = paste0("Fab (", level, ")"),
      hr = sfit_multi$coefficients[coef_name, "exp(coef)"],
      lower = sfit_multi$conf.int[coef_name, "lower .95"],
      upper = sfit_multi$conf.int[coef_name, "upper .95"],
      pval = sfit_multi$coefficients[coef_name, "Pr(>|z|)"]
    )
  }
}

# 
multi_results[["age"]] <- data.frame(
  variable = "age",
  hr = sfit_multi$coefficients["age", "exp(coef)"],
  lower = sfit_multi$conf.int["age", "lower .95"],
  upper = sfit_multi$conf.int["age", "upper .95"],
  pval = sfit_multi$coefficients["age", "Pr(>|z|)"]
)

# 
multi_results[["runx1.runx1t1_fusion"]] <- data.frame(
  variable = "runx1-runx1t1_fusion",
  hr = sfit_multi$coefficients["runx1.runx1t1_fusionYes", "exp(coef)"],
  lower = sfit_multi$conf.int["runx1.runx1t1_fusionYes", "lower .95"],
  upper = sfit_multi$conf.int["runx1.runx1t1_fusionYes", "upper .95"],
  pval = sfit_multi$coefficients["runx1.runx1t1_fusionYes", "Pr(>|z|)"]
)

# 
multi_results[["runx1_mutation"]] <- data.frame(
  variable = "runx1_mutation",
  hr = sfit_multi$coefficients["runx1_mutationYes", "exp(coef)"],
  lower = sfit_multi$conf.int["runx1_mutationYes", "lower .95"],
  upper = sfit_multi$conf.int["runx1_mutationYes", "upper .95"],
  pval = sfit_multi$coefficients["runx1_mutationYes", "Pr(>|z|)"]
)

# 
multi_results[["RS"]] <- data.frame(
  variable = "riskScore",
  hr = sfit_multi$coefficients["RS", "exp(coef)"],
  lower = sfit_multi$conf.int["RS", "lower .95"],
  upper = sfit_multi$conf.int["RS", "upper .95"],
  pval = sfit_multi$coefficients["RS", "Pr(>|z|)"]
)

# 
multi_df <- bind_rows(multi_results)

##
univ_forest <- rbind(
  data.frame(variable = "Characteristics", hr = NA, lower = NA, upper = NA, pval = NA),
  univ_df
)

univ_forest$pval_label <- ifelse(
  univ_forest$pval < 0.001, "***",
  ifelse(univ_forest$pval < 0.01, "**",
         ifelse(univ_forest$pval < 0.05, "*", sprintf("%.3f", univ_forest$pval)))
)

univ_tabletext <- cbind(
  univ_forest$variable,
  paste0(
    sprintf("%.3f", univ_forest$hr),
    " (", sprintf("%.3f", univ_forest$lower),
    " - ", sprintf("%.3f", univ_forest$upper), ")"
  ),
  univ_forest$pval_label
)

univ_mean <- c(NA, univ_forest$hr)
univ_lower <- c(NA, univ_forest$lower)
univ_upper <- c(NA, univ_forest$upper)

pdf("univariate_forestplot.pdf", width = 10, height = 8)
forestplot(
  labeltext = univ_tabletext,
  mean = univ_mean,
  lower = univ_lower,
  upper = univ_upper,
  new_page = TRUE,
  is.summary = c(TRUE, rep(FALSE, nrow(univ_forest)-1)),
  xlog = FALSE,
  col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
  xlab = "Hazard Ratio (HR)",
  title = "Univariate Cox Regression Analysis",
  csvmx = 1,
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 10),
    ticks = gpar(fontsize = 10),
    xlab = gpar(fontsize = 12),
    title = gpar(fontsize = 14, fontface = "bold")
  ),
  align = c("left", "center", "center")
)
dev.off()

##
multi_forest <- rbind(
  data.frame(variable = "Characteristics", hr = NA, lower = NA, upper = NA, pval = NA),
  multi_df
)

multi_forest$pval_label <- ifelse(
  multi_forest$pval < 0.001, "***",
  ifelse(multi_forest$pval < 0.01, "**",
         ifelse(multi_forest$pval < 0.05, "*", sprintf("%.3f", multi_forest$pval)))
)

multi_tabletext <- cbind(
  multi_forest$variable,
  paste0(
    sprintf("%.3f", multi_forest$hr),
    " (", sprintf("%.3f", multi_forest$lower),
    " - ", sprintf("%.3f", multi_forest$upper), ")"
  ),
  multi_forest$pval_label
)

multi_mean <- c(NA, multi_forest$hr)
multi_lower <- c(NA, multi_forest$lower)
multi_upper <- c(NA, multi_forest$upper)

pdf("multivariate_forestplot.pdf", width = 10, height = 8)
forestplot(
  labeltext = multi_tabletext,
  mean = multi_mean,
  lower = multi_lower,
  upper = multi_upper,
  new_page = TRUE,
  is.summary = c(TRUE, rep(FALSE, nrow(multi_forest)-1)),
  xlog = FALSE,
  col = fpColors(box = "seagreen", line = "darkgreen", summary = "seagreen"),
  xlab = "Hazard Ratio (HR)",
  title = "Multivariate Cox Regression Analysis",
  csvmx = 1,
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 10),
    ticks = gpar(fontsize = 10),
    xlab = gpar(fontsize = 12),
    title = gpar(fontsize = 14, fontface = "bold")
  ),
  align = c("left", "center", "center")
)
dev.off()
##C&D&E
install.packages(c("rms", "rmda", "survival"))
library(rms)
library(rmda)
library(survival)

data <- read.csv("GSE37642RS.csv", stringsAsFactors = TRUE)
data$runx1_mutation <- factor(data$runx1_mutation, levels = c("No", "Yes"))
data$OS <- as.integer(data$OS)

dd <- datadist(data)
options(datadist = "dd")
fit_nomo <- cph(Surv(OS.time, OS) ~ age + runx1_mutation + RS, data = data, x = TRUE, y = TRUE, surv = TRUE)

pdf("nomogram.pdf", width = 8, height = 10)
nomogram(fit_nomo, 
         fun = list(function(x) survprob(fit_nomo, x, times = 1),
                    function(x) survprob(fit_nomo, x, times = 3),
                    function(x) survprob(fit_nomo, x, times = 5)),
         funlabel = c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
         lp = FALSE,
         xfrac = 0.4)
dev.off()

cal <- calibrate(fit_nomo, u = c(1, 3, 5), m = 30, B = 1000)
pdf("calibration_curve.pdf", width = 8, height = 8)
plot(cal, xlim = c(0, 1), ylim = c(0, 1), xlab = "Predicted Survival Probability", ylab = "Observed fraction survival probability")
abline(0, 1, lty = 2, col = "gray")
dev.off()

dc_age <- decision_curve(Surv(OS.time, OS) ~ age, data = data, family = "survival", times = 1, thresholds = seq(0, 1, 0.01), model = "logistic")
dc_runx1 <- decision_curve(Surv(OS.time, OS) ~ runx1_mutation, data = data, family = "survival", times = 1, thresholds = seq(0, 1, 0.01), model = "logistic")
dc_risk <- decision_curve(Surv(OS.time, OS) ~ RS, data = data, family = "survival", times = 1, thresholds = seq(0, 1, 0.01), model = "logistic")
dc_all <- decision_curve(Surv(OS.time, OS) ~ age + runx1_mutation + RS, data = data, family = "survival", times = 1, thresholds = seq(0, 1, 0.01), model = "logistic")
dc_all_positive <- decision_curve(Surv(OS.time, OS) ~ I(age >= 60) + runx1_mutation + I(RS >= median(RS)), data = data, family = "survival", times = 1, thresholds = seq(0, 1, 0.01), model = "logistic")
dc_all_negative <- decision_curve(Surv(OS.time, OS) ~ I(age < 60) + runx1_mutation + I(RS < median(RS)), data = data, family = "survival", times = 1, thresholds = seq(0, 1, 0.01), model = "logistic")

pdf("decision_curve.pdf", width = 8, height = 8)
plot_decision_curve(
  dc_age, dc_runx1, dc_risk, dc_all, dc_all_positive, dc_all_negative,
  curve.names = c("Age", "runx1_mutation", "riskScore", "All", "All positive", "All negative"),
  xlim = c(0, 1), ylim = c(0, 0.5),
  xlab = "Threshold Probability", ylab = "Net Benefit",
  legend.position = "topright"
)
dev.off()
##F&G&H
library(survival)
library(survminer)
library(timeROC)
library(RColorBrewer)

# Define datasets to analyze
datasets <- c("GSE37642", "GSE12417", "merge")

# Loop through each dataset
for (sel in datasets) {
  # Prepare data for current dataset
  score_t <- data.frame(OS.time = get(sel)$OS.time, OS = get(sel)$OS, Score = rs[[sel]]$RS)
  score_t <- na.omit(score_t)
  
  # Determine optimal cutpoint for risk score
  cut <- surv_cutpoint(score_t, time = "OS.time", event = "OS", variable = "Score")
  cat <- surv_categorize(cut)
  
  # Survival analysis
  fit1 <- survfit(Surv(OS.time, OS) ~ Score, data = cat)
  
  # Custom theme for plots
  mytheme <- theme_survminer(font.legend = c(14, "plain", "black"),
                             font.x = c(14, "plain", "black"),
                             font.y = c(14, "plain", "black"))
  
  # Save survival curve with risk table
  pdf(paste0(sel, "_survival_curve.pdf"), width = 8, height = 8)
  ggsurvplot(fit1, cat,
             palette = "jco",
             size = 1.3,
             pval = TRUE,
             legend.labs = c("High", "Low"),
             legend.title = "Score",
             xlab = "Time(years)",
             ylab = "Overall survival",
             ggtheme = mytheme,
             break.time.by = 2,
             conf.int = TRUE,
             risk.table = TRUE,
             risk.table.title = "",
             risk.table.height = .25)
  dev.off()
  
  # Time-dependent ROC analysis
  tt <- timeROC(T = score_t$OS.time, delta = score_t$OS, marker = score_t$Score,
                cause = 1, weighting = 'marginal',
                times = seq(1, 5, 1), ROC = TRUE, iid = TRUE)
  
  # Save time ROC curve plot
  pdf(paste0(sel, "_time_ROC_curve.pdf"), width = 8, height = 8)
  plot(tt, time = 1, title = FALSE, lwd = 1.5, col = brewer.pal(3, "Set1")[1])
  plot(tt, time = 3, col = brewer.pal(3, "Set1")[2], add = TRUE, title = FALSE, lwd = 1.5)
  plot(tt, time = 5, col = brewer.pal(3, "Set1")[3], add = TRUE, title = FALSE, lwd = 1.5)
  legend("bottomright", 
         legend = c(paste0("1-year AUC = ", round(tt$AUC[1], 3)),
                    paste0("3-year AUC = ", round(tt$AUC[3], 3)),
                    paste0("5-year AUC = ", round(tt$AUC[5], 3))),
         fill = brewer.pal(3, "Set1")[1:3],
         bty = "o", cex = 1,
         border = NA)
  abline(0, 1, lty = 2, lwd = 0.5)
  dev.off()
  
  # Plot AUC curve
  pdf(paste0(sel, "_AUC_curve.pdf"), width = 8, height = 8)
  plotAUCcurve(tt, conf.int = TRUE, col = "tomato")
  dev.off()
}
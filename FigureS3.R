rm(list = ls())
library(Seurat)
library(copykat)
library(tidyr)
library(gtools)
library(ggplot2)
library(cowplot)

setwd("/home/cm985017/zenodo/zenodo/zenodo/")

sce.all <- readRDS("data/zenodo_final.rds")

table(sce.all$celltype)

cells_to_analyze <- subset(sce.all, subset = celltype %in% c("HSC", "GMP", "MEP", "T/NK"))

set.seed(123)

if (ncol(cells_to_analyze) > 10000) {
  sampled_cells <- sample(colnames(cells_to_analyze), size = 10000)
  cells_to_analyze <- subset(cells_to_analyze, cells = sampled_cells)
}

counts_matrix <- as.matrix(GetAssayData(
  object = cells_to_analyze, 
  assay = "RNA",
  layer = "counts"
))

ref_cells <- colnames(cells_to_analyze)[cells_to_analyze$celltype == "T/NK"]

outputdir <- "/home/cm985017/zenodo/zenodo/zenodo/copykat_analysis/"
dir.create(outputdir, showWarnings = FALSE)
dir.create(file.path(outputdir, "01_copykat"), showWarnings = FALSE)

genome_coords <- read.table("chromosome_coordinates_hg19.txt", header = TRUE)

copykat_res <- copykat(
  rawmat = counts_matrix,
  id.type = "S",
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.1,
  sam.name = file.path(outputdir, "01_copykat", "AML_copykat"),
  norm.cell.names = ref_cells,
  distance = "euclidean",
  n.cores = 8
)

saveRDS(copykat_res, file.path(outputdir, "01_copykat", "copykat_fullresult.rds"))

predictions <- data.frame(copykat_res$prediction)
CNA_data <- data.frame(copykat_res$CNAmat)

convert_cna_colnames <- function(colnames) {
  new_names <- gsub("^X", "", colnames)
  new_names <- gsub("\\.", "-", new_names)
  return(new_names)
}

original_colnames <- colnames(CNA_data)
new_colnames <- c(original_colnames[1:3], convert_cna_colnames(original_colnames[4:length(original_colnames)]))
colnames(CNA_data) <- new_colnames

common_cells <- intersect(colnames(CNA_data), predictions$cell.names)

assess_chromarm <- function(chrom, arm, CNA_data, genome_coords, cell_labels, 
                            calling_threshold = 0.8, pvalue_signif_cutoff = 0.01) {
  
  cent_start = genome_coords$cen.left.base[genome_coords$chr == chrom]
  cent_end = genome_coords$cen.right.base[genome_coords$chr == chrom]
  
  aneuploid_cells <- intersect(
    cell_labels$cell.names[cell_labels$copykat.pred == "aneuploid"], 
    colnames(CNA_data)
  )
  diploid_cells <- intersect(
    cell_labels$cell.names[cell_labels$copykat.pred == "diploid"], 
    colnames(CNA_data)
  )
  
  if (length(aneuploid_cells) == 0 | length(diploid_cells) == 0) {
    return(list(
      chromname = paste0(chrom, arm),
      cells_with_alteration = character(0)
    ))
  }
  
  if (arm == "p") {
    chromname = paste0(chrom, "p")
    region_filter <- CNA_data$chrom == chrom & CNA_data$chrompos < cent_start
    chrom_aneuploid = CNA_data[region_filter, aneuploid_cells, drop = FALSE]
    chrom_diploid = CNA_data[region_filter, diploid_cells, drop = FALSE]
  } else {
    chromname = paste0(chrom, "q")
    region_filter <- CNA_data$chrom == chrom & CNA_data$chrompos > cent_end
    chrom_aneuploid = CNA_data[region_filter, aneuploid_cells, drop = FALSE]
    chrom_diploid = CNA_data[region_filter, diploid_cells, drop = FALSE]
  }
  
  if (nrow(chrom_aneuploid) == 0 | nrow(chrom_diploid) == 0) {
    return(list(
      chromname = chromname,
      cells_with_alteration = character(0)
    ))
  }
  
  chrom_aneuploid_m = tidyr::pivot_longer(
    chrom_aneuploid, 
    cols = colnames(chrom_aneuploid),
    names_to = "cell",
    values_to = "value"
  )
  chrom_diploid_m = tidyr::pivot_longer(
    chrom_diploid, 
    cols = colnames(chrom_diploid),
    names_to = "cell",
    values_to = "value"
  )
  
  chrom_diploid_mean = mean(chrom_diploid_m$value, na.rm = TRUE)
  chrom_diploid_sd = sd(chrom_diploid_m$value, na.rm = TRUE)
  
  chrom_aneuploid_means <- colMeans(chrom_aneuploid, na.rm = TRUE)
  pvalues <- sapply(chrom_aneuploid_means, function(x) {
    z <- (x - chrom_diploid_mean) / chrom_diploid_sd
    2 * pnorm(-abs(z))
  })
  
  is_signif <- pvalues < pvalue_signif_cutoff
  cells_with_alteration <- names(chrom_aneuploid_means)[is_signif]
  
  if (nrow(chrom_aneuploid_m) > 0) {
    chrom_aneuploid_m$altered = chrom_aneuploid_m$cell %in% cells_with_alteration
  }
  
  return(list(
    chromname = chromname,
    chrom_aneuploid_m = chrom_aneuploid_m,
    chrom_diploid_m = chrom_diploid_m,
    chrom_diploid_mean = chrom_diploid_mean,
    chrom_diploid_sd = chrom_diploid_sd,
    cells_with_alteration = cells_with_alteration
  ))
}

aml_chromosomes <- c(5, 7, 17)
aml_results <- list()

for (chrom in aml_chromosomes) {
  res_p <- assess_chromarm(
    chrom = chrom, arm = "p", 
    CNA_data = CNA_data, genome_coords = genome_coords, 
    cell_labels = predictions, calling_threshold = 0.8
  )
  
  res_q <- assess_chromarm(
    chrom = chrom, arm = "q", 
    CNA_data = CNA_data, genome_coords = genome_coords, 
    cell_labels = predictions, calling_threshold = 0.8
  )
  
  aml_results[[as.character(chrom)]] <- list(p_arm = res_p, q_arm = res_q)
}

saveRDS(aml_results, file.path(outputdir, "AML_specific_CNV.rds"))

generate_aml_report <- function(aml_results, predictions) {
  report <- data.frame()
  total_aneuploid <- sum(predictions$copykat.pred == "aneuploid")
  
  for (chrom in names(aml_results)) {
    p_data <- aml_results[[chrom]]$p_arm
    q_data <- aml_results[[chrom]]$q_arm
    
    p_del_cells <- length(p_data$cells_with_alteration)
    q_del_cells <- length(q_data$cells_with_alteration)
    
    chrom_report <- data.frame(
      Chromosome = chrom,
      Arm = c("p", "q"),
      Altered_Cells = c(p_del_cells, q_del_cells),
      Percentage_Aneuploid = c(
        ifelse(total_aneuploid > 0, round(p_del_cells/total_aneuploid * 100, 2), 0),
        ifelse(total_aneuploid > 0, round(q_del_cells/total_aneuploid * 100, 2), 0)
      )
    )
    
    report <- rbind(report, chrom_report)
  }
  
  return(report)
}

aml_report <- generate_aml_report(aml_results, predictions)
write.csv(aml_report, file.path(outputdir, "AML_CNV_Report.csv"), row.names = FALSE)

plot_aml_cnv <- function(aml_results) {
  plots <- list()
  
  for (chrom in names(aml_results)) {
    p_data <- aml_results[[chrom]]$p_arm
    q_data <- aml_results[[chrom]]$q_arm
    
    if (!is.null(p_data$chrom_diploid_m) && !is.null(p_data$chrom_aneuploid_m)) {
      p_plot <- ggplot() +
        geom_density(data = p_data$chrom_diploid_m,
                     aes(x = value), 
                     fill = "lightblue", alpha = 0.7, color = "blue") +
        geom_density(data = p_data$chrom_aneuploid_m,
                     aes(x = value), 
                     fill = "red", alpha = 0.7, color = "darkred") +
        ggtitle(paste("Chromosome", chrom, "p arm")) +
        xlab("CNA value") +
        ylab("Density") +
        theme_cowplot()
      
      plots[[paste0("chr", chrom, "p")]] <- p_plot
    }
    
    if (!is.null(q_data$chrom_diploid_m) && !is.null(q_data$chrom_aneuploid_m)) {
      q_plot <- ggplot() +
        geom_density(data = q_data$chrom_diploid_m,
                     aes(x = value), 
                     fill = "lightblue", alpha = 0.7, color = "blue") +
        geom_density(data = q_data$chrom_aneuploid_m,
                     aes(x = value), 
                     fill = "red", alpha = 0.7, color = "darkred") +
        ggtitle(paste("Chromosome", chrom, "q arm")) +
        xlab("CNA value") +
        ylab("Density") +
        theme_cowplot()
      
      plots[[paste0("chr", chrom, "q")]] <- q_plot
    }
  }
  
  return(plots)
}

aml_plots <- plot_aml_cnv(aml_results)
for (plot_name in names(aml_plots)) {
  ggsave(
    filename = file.path(outputdir, paste0(plot_name, "_CNV.png")),
    plot = aml_plots[[plot_name]],
    width = 8, height = 6
  )
}

simple_summary <- data.frame(
  CNV_Type = c("Chromosome 7 Deletion (-7)", "Chromosome 5q Deletion (del(5q))", "Chromosome 17p Aberration (abn(17p))"),
  Affected_Cells = c(
    ifelse("7" %in% names(aml_results), length(aml_results[["7"]]$p_arm$cells_with_alteration), 0),
    ifelse("5" %in% names(aml_results), length(aml_results[["5"]]$q_arm$cells_with_alteration), 0),
    ifelse("17" %in% names(aml_results), length(aml_results[["17"]]$p_arm$cells_with_alteration), 0)
  ),
  Percentage_of_Aneuploid = c(
    ifelse("7" %in% names(aml_results), 
           round(length(aml_results[["7"]]$p_arm$cells_with_alteration)/sum(predictions$copykat.pred == "aneuploid") * 100, 2), 0),
    ifelse("5" %in% names(aml_results), 
           round(length(aml_results[["5"]]$q_arm$cells_with_alteration)/sum(predictions$copykat.pred == "aneuploid") * 100, 2), 0),
    ifelse("17" %in% names(aml_results), 
           round(length(aml_results[["17"]]$p_arm$cells_with_alteration)/sum(predictions$copykat.pred == "aneuploid") * 100, 2), 0)
  ),
  Clinical_Significance = c(
    "High-risk genetic abnormality, associated with poor prognosis",
    "Poor prognosis marker, common in therapy-related AML", 
    "TP53-associated, linked to treatment resistance"
  )
)

write.csv(simple_summary, file.path(outputdir, "AML_Key_CNV_Summary.csv"), row.names = FALSE)
rm(list = ls())
load("R_Datasets/abundance_by_group.rdata")

datasets_list <- list(
  fish = fish[, -(1:8)],
  benthic_invertebrate = benthic_invertebrate[, -(1:8)],
  insect = insect[, -(1:8)],
  zooplankton = zooplankton[, -(1:8)],
  fungi = fungi[, -(1:8)],
  bacteria = bacteria[, -(1:8)],
  phytoplankton = phytoplankton[, -(1:8)]
)

perform_procrustes <- function(name1, name2, data1, data2) {
  transpose <- function(y) {
    current_df <- as.data.frame(y)
    species <- current_df[, 1]
    mat <- t(current_df[, -1])
    colnames(mat) <- species
    return(as.data.frame(mat))
  }
  
  pcoa1 <- cmdscale(vegdist(decostand(transpose(data1), method = "hellinger"), method = "bray"))
  pcoa2 <- cmdscale(vegdist(decostand(transpose(data2), method = "hellinger"), method = "bray"))
  
  set.seed(123)
  pro_test <- protest(pcoa1, pcoa2, permutations = 999)
  Pro_coords <- cbind(data.frame(pro_test$Yrot), data.frame(pro_test$X))
  Pro_rotation <- data.frame(pro_test$rotation)
  colnames(Pro_coords) <- c("X1", "X2", "Dim1", "Dim2")
  
  plotdf_group <- data.frame(group = c(name1, name2))
  legend_data <- data.frame(
    group = c(name1, name2),
    color = c("#E41A1C", "#377EB8")
  )
  
  p <- ggplot(data = Pro_coords) +
    geom_point(aes(X1, X2, color = name1), size = 7, shape = 16) +
    geom_point(aes(Dim1, Dim2, color = name2), size = 7, shape = 16) +
    geom_segment(
      aes(x = X1, y = X2, xend = (X1 + Dim1) / 2, yend = (X2 + Dim2) / 2),
      arrow = arrow(length = unit(0, "cm")),
      color = "#E41A1C",
      linewidth = 0.9
    ) +
    geom_segment(
      aes(x = (X1 + Dim1) / 2, y = (X2 + Dim2) / 2, xend = Dim1, yend = Dim2),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "#377EB8",
      linewidth = 0.9
    ) +
    scale_color_manual(
      values = setNames(c("#E41A1C", "#377EB8"), c(name1, name2)),
      name = "Taxonomic Group",
      limits = c(name1, name2)
    ) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x)) +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
    geom_vline(xintercept = 0, color = "gray", linetype = 2, linewidth = 0.95) +
    geom_hline(yintercept = 0, color = "gray", linetype = 2, linewidth = 0.95) +
    geom_abline(intercept = 0, slope = Pro_rotation[1, 2] / Pro_rotation[1, 1], linewidth = 0.85) +
    geom_abline(intercept = 0, slope = Pro_rotation[2, 2] / Pro_rotation[2, 1], linewidth = 0.85) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      axis.ticks.length = unit(0.4, "lines"),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_text(colour = "black", size = 37),
      axis.title.y = element_text(colour = "black", size = 37),
      axis.text = element_text(colour = "black", size = 26),
      legend.position = c(0.235, 0.15),
      legend.key.size = unit(1.5, "cm"),
      legend.title = element_text(size = 24.5, face = "bold"),
      legend.text = element_text(size = 24.5),
      legend.background = element_rect(color = "transparent", fill = alpha("white", 0.1)),
      legend.key = element_rect(fill = "transparent", color = "transparent"),
      plot.margin = unit(c(0.55, 0.55, 0.55, 0.55), "cm")
    ) +
    annotate(
      "text",
      label = sprintf("M2 = %.3f\np = %.3f", pro_test$ss, pro_test$signif),
      x = Inf,
      y = Inf,
      hjust = 1.2,
      vjust = 1.35,
      size = 27.5
    )
  
  ggsave(
    sprintf("Result_Figures/Procrustes/%s-%s.png", name1, name2),
    plot = p,
    width = 11,
    height = 11,
    dpi = 200
  )
  
  return(list(
    pair = paste(name1, "-", name2),
    M2 = pro_test$ss,
    p_value = pro_test$signif,
    residuals = residuals(procrustes(pcoa1, pcoa2, symmetric = TRUE))
  ))
}

results <- list()
dataset_names <- names(datasets_list)

for (i in 1:(length(dataset_names) - 1)) {
  name1 <- dataset_names[i]
  for (j in (i + 1):length(dataset_names)) {
    name2 <- dataset_names[j]
    result <- perform_procrustes(name1, name2, datasets_list[[name1]], datasets_list[[name2]])
    results[[length(results) + 1]] <- result
  }
}

for (i in 1:length(results)) {
  result <- results[[i]]
  cat(sprintf("%s: M2 = %.3f, p = %.3f\n", result$pair, result$M2, result$p_value))
}
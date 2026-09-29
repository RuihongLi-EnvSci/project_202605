rm(list = ls())
set.seed(123)
load("R_Datasets/abundance_by_group.rdata")
load("R_Datasets/metrics.rdata")

multigroups <- rbind(fish, benthic_invertebrate, insect, zooplankton, fungi, bacteria, phytoplankton)

datasets_list <- list(
  fish = fish[, -(1:8)],
  benthic_invertebrate = benthic_invertebrate[, -(1:8)],
  insect = insect[, -(1:8)],
  zooplankton = zooplankton[, -(1:8)],
  fungi = fungi[, -(1:8)],
  bacteria = bacteria[, -(1:8)],
  phytoplankton = phytoplankton[, -(1:8)],
  multigroups = multigroups[, -(1:8)]
)

data_original <- data_original[order(data_original[, 1]), ]
env_vars <- c("T", "pH", "DO", "COD", "NH3N", "TP", "TN", "EC", "TUB")
env_data <- data_original[, env_vars]

results_list <- list()

for (group_name in names(datasets_list)) {
  df <- datasets_list[[group_name]]
  
  abundance <- df[, 2:73]
  abundance <- t(abundance)
  abundance <- abundance[order(rownames(abundance)), ]
  
  rda_model <- rda(abundance ~ ., data = env_data)
  rsq <- RsquareAdj(rda_model)
  r2 <- rsq$r.squared
  r2_adj <- rsq$adj.r.squared
  
  anova_model <- anova(rda_model, permutations = 999)
  model_p <- anova_model[1, "Pr(>F)"]
  
  anova_axis <- anova(rda_model, by = "axis", permutations = 999)
  rda1_p <- anova_axis[1, "Pr(>F)"]
  rda2_p <- ifelse(nrow(anova_axis) >= 2, anova_axis[2, "Pr(>F)"], NA)
  
  anova_env <- anova(rda_model, by = "margin", permutations = 999)
  env_p_values <- anova_env[, "Pr(>F)", drop = FALSE]
  colnames(env_p_values) <- "p_value"
  
  env_scores <- vegan::scores(rda_model, display = "bp", choices = c(1, 2), scaling = 2)
  env_scores_df <- data.frame(
    Environmental_factor = rownames(env_scores),
    RDA1 = env_scores[, 1],
    RDA2 = env_scores[, 2],
    row.names = NULL
  )
  
  env_result_df <- merge(
    env_scores_df,
    data.frame(
      Environmental_factor = rownames(env_p_values),
      p_value = env_p_values[, 1],
      row.names = NULL
    ),
    by = "Environmental_factor",
    all.x = TRUE
  )
  
  env_result_df <- env_result_df[match(env_vars, env_result_df$Environmental_factor), ]
  
  group_result <- list(
    Group = group_name,
    R2 = r2,
    Adj_R2 = r2_adj,
    Model_p = model_p,
    RDA1_p = rda1_p,
    RDA2_p = rda2_p,
    Env_p_values = env_p_values,
    Env_scores = env_scores_df,
    Env_result = env_result_df
  )
  
  results_list[[group_name]] <- group_result
  
  cat("======", group_name, "======", "\n")
  cat("R² =", round(r2, 4), "\n")
  cat("Adjusted R² =", round(r2_adj, 4), "\n")
  cat("Model p-value =", model_p, "\n")
  cat("RDA1 p-value =", rda1_p, "\n")
  cat("RDA2 p-value =", rda2_p, "\n")
  cat("\n--- comprehensive results---\n")
  print(env_result_df)
  cat("\n")
}

model_tables <- list()

for (group_name in names(results_list)) {
  group_result <- results_list[[group_name]]
  
  result_df <- group_result$Env_result
  result_df$Group <- group_name
  result_df$R2 <- group_result$R2
  result_df$Adj_R2 <- group_result$Adj_R2
  result_df$Model_p <- group_result$Model_p
  result_df$RDA1_p <- group_result$RDA1_p
  result_df$RDA2_p <- group_result$RDA2_p
  
  model_tables[[group_name]] <- result_df
}

write.xlsx(model_tables, "Result_Tables/RDA Results.xlsx")
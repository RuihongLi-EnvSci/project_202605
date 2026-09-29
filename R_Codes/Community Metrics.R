rm(list = ls())

load("R_Datasets/abundance_by_group.rdata")

all_groups <- list(
  bacteria = bacteria[, -(1:8)],
  benthic_invertebrate = benthic_invertebrate[, -(1:8)],
  fish = fish[, -(1:8)],
  fungi = fungi[, -(1:8)],
  insect = insect[, -(1:8)],
  phytoplankton = phytoplankton[, -(1:8)],
  zooplankton = zooplankton[, -(1:8)]
)

hellinger_transform <- function(df) {
  species_col <- df[, 1, drop = FALSE]
  abundance <- df[, -1, drop = FALSE]
  
  hellinger_abundance <- apply(abundance, 2, function(col) {
    col_sum <- sum(col, na.rm = TRUE)
    if (col_sum == 0) {
      return(rep(0, length(col)))
    } else {
      return(sqrt(col / col_sum))
    }
  })
  
  hellinger_abundance <- as.data.frame(hellinger_abundance)
  colnames(hellinger_abundance) <- colnames(abundance)
  
  result <- cbind(species_col, hellinger_abundance)
  
  return(result)
}

all_groups_hellinger <- lapply(all_groups, hellinger_transform)

rm(list = setdiff(ls(), c("all_groups", "all_groups_hellinger")))

transpose <- function(x) {
  current_df <- as.data.frame(x)
  col_names <- current_df[, 1]
  current_df <- as.data.frame(t(current_df[, -1]))
  current_df[] <- lapply(current_df, function(x) as.numeric(as.character(x)))
  colnames(current_df) <- col_names
  return(current_df)
}

all_groups_transposed <- lapply(all_groups, transpose)

richness <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]]))
  
  for (y in names(x)) {
    richness_values <- vegan::specnumber(x[[y]])
    temp_df <- data.frame(richness = as.numeric(richness_values))
    colnames(temp_df)[1] <- y
    result_df <- cbind(result_df, temp_df)
  }
  
  return(result_df)
}

richness_single_groups <- richness(all_groups_transposed)

diversity <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]]))
  
  for (y in names(x)) {
    diversity_values <- vegan::diversity(x[[y]])
    temp_df <- data.frame(diversity = as.numeric(diversity_values))
    colnames(temp_df)[1] <- y
    result_df <- cbind(result_df, temp_df)
  }
  
  return(result_df)
}

diversity_single_groups <- diversity(all_groups_transposed)

convert_to_longdata <- function(df) {
  colnames(df)[1] <- "species"
  
  df[, -1] <- lapply(df[, -1], function(x) as.numeric(as.character(x)))
  
  long_df <- df %>%
    pivot_longer(
      cols = -species,
      names_to = c("sampling_site", "sampling_order"),
      names_sep = "_",
      values_to = "abundance"
    ) %>%
    mutate(sampling_order = as.integer(sampling_order))
  
  return(long_df)
}

all_groups_longdata <- lapply(all_groups_hellinger, convert_to_longdata)

synchrony <- function(x) {
  result_df <- data.frame(sampling_site = x[[1]][(1:18), 2])
  
  for (y in names(x)) {
    synchrony_values <- codyn::synchrony(
      x[[y]],
      species.var = "species",
      time.var = "sampling_order",
      abundance.var = "abundance",
      replicate.var = "sampling_site"
    )
    colnames(synchrony_values)[2] <- y
    result_df <- cbind(result_df, synchrony_values[2])
  }
  
  return(result_df)
}

synchrony_single_groups <- synchrony(all_groups_longdata)

stability <- function(x) {
  result_df <- data.frame(sampling_site = x[[1]][(1:18), 2])
  
  for (y in names(x)) {
    stability_values <- codyn::community_stability(
      x[[y]],
      time.var = "sampling_order",
      abundance.var = "abundance",
      replicate.var = "sampling_site"
    )
    colnames(stability_values)[2] <- y
    result_df <- cbind(result_df, stability_values[2])
  }
  
  return(result_df)
}

stability_single_groups <- stability(all_groups_longdata)

multi_groups_transposed <- do.call(cbind, all_groups_transposed)
multi_groups_longdata <- do.call(rbind, all_groups_longdata)
richness_multi_groups <- data.frame(multi_groups = vegan::specnumber(multi_groups_transposed))
diversity_multi_groups <- data.frame(multi_groups = vegan::diversity(multi_groups_transposed))
synchrony_multi_groups <- synchrony(list(multi_groups = multi_groups_longdata))
stability_multi_groups <- stability(list(multi_groups = multi_groups_longdata))

community_metrics_collection <- list(
  richness = cbind(richness_single_groups, richness_multi_groups),
  diversity = cbind(diversity_single_groups, diversity_multi_groups),
  synchrony = cbind(synchrony_single_groups, multi_groups = synchrony_multi_groups[, 2]),
  stability = cbind(stability_single_groups, multi_groups = stability_multi_groups[, 2])
)

write.xlsx(community_metrics_collection, "Result_Tables/Community Metrics.xlsx")
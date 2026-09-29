### Box Plot

rm(list = ls())

load("R_Datasets/metrics.rdata")

data_processed <- data_original %>%
  separate(sampling_site_order, into = c("site", "year_code"), sep = "_", remove = FALSE) %>%
  mutate(
    year = case_when(
      year_code == "1" ~ "1st",
      year_code == "2" ~ "2nd",
      year_code == "3" ~ "3rd",
      year_code == "4" ~ "4th",
      TRUE ~ year_code
    ),
    year = factor(year, levels = c("1st", "2nd", "3rd", "4th")),
    year_numeric = as.numeric(year)
  )

data_long <- data_processed %>%
  pivot_longer(
    cols = c(connectance, mean_path_length, modularity, nestedness, vulnerability, robustness),
    names_to = "index",
    values_to = "value"
  ) %>%
  mutate(
    index = factor(index, levels = c("connectance", "mean_path_length", "modularity", "nestedness", "vulnerability", "robustness"))
  )

year_colors <- c(
  "1st" = "#647ADD",
  "2nd" = "#C0C0C0",
  "3rd" = "#FFA500",
  "4th" = "#F5664D"
)

index_names <- c(
  "connectance" = "Connectance",
  "mean_path_length" = "Mean Path Length",
  "modularity" = "Modularity",
  "nestedness" = "Nestedness",
  "vulnerability" = "Vulnerability",
  "robustness" = "Robustness"
)

y_labels <- c(
  "connectance" = "Connectance",
  "mean_path_length" = "Mean Path Length",
  "modularity" = "Modularity",
  "nestedness" = "Nestedness",
  "vulnerability" = "Vulnerability",
  "robustness" = "Robustness"
)

y_params <- list(
  connectance = list(limits = c(0, 0.15)),
  mean_path_length = list(limits = c(1, 2)),
  modularity = list(limits = c(0, 0.18)),
  nestedness = list(limits = c(0, 6.5)),
  robustness = list(limits = c(0, 0.35)),
  vulnerability = list(limits = c(10, 75))
)

for (current_index in levels(data_long$index)) {
  index_data <- data_long %>%
    filter(index == current_index) %>%
    drop_na(value)
  
  index_data$year <- factor(index_data$year, levels = c("1st", "2nd", "3rd", "4th"))
  
  lm_model <- lm(value ~ year, data = index_data)
  anova_res <- anova(lm_model)
  p_value <- anova_res$`Pr(>F)`[1]
  f_value <- anova_res$`F value`[1]
  f_value_text <- sprintf("F = %.2f", f_value)
  
  p_label <- if (is.na(p_value)) {
    "p = NA"
  } else if (p_value < 0.001) {
    "P < 0.001"
  } else if (p_value < 0.01) {
    "P < 0.01"
  } else if (p_value < 0.05) {
    "P < 0.05"
  } else {
    paste0("P = ", sprintf("%.3f", p_value))
  }
  
  anova_text <- paste(f_value_text, "", p_label)
  
  anova_model <- aov(value ~ year, data = index_data)
  tukey_result <- HSD.test(anova_model, "year", group = TRUE, console = FALSE)
  tukey_groups <- tukey_result$groups
  tukey_groups <- tukey_groups[levels(index_data$year), , drop = FALSE]
  
  boxplot_stats <- index_data %>%
    group_by(year) %>%
    dplyr::summarize(
      Q1 = quantile(value, 0.25, na.rm = TRUE),
      Q3 = quantile(value, 0.75, na.rm = TRUE),
      Median = median(value, na.rm = TRUE),
      IQR = IQR(value, na.rm = TRUE),
      Upper_whisker = min(max(value, na.rm = TRUE), Q3 + 1.5 * IQR),
      Lower_whisker = max(min(value, na.rm = TRUE), Q1 - 1.5 * IQR),
      Max_data = max(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  letter_df <- data.frame(
    year = rownames(tukey_groups),
    letter = tukey_groups$groups,
    stringsAsFactors = FALSE
  ) %>%
    mutate(year = factor(year, levels = levels(index_data$year))) %>%
    left_join(boxplot_stats, by = "year")
  
  y_range <- diff(y_params[[current_index]]$limits)
  
  letter_df <- letter_df %>%
    mutate(
      base_y = pmax(Upper_whisker, Q3),
      letter_y = base_y + y_range * 0.085
    )
  
  if (current_index == "nestedness") {
    letter_df <- letter_df %>%
      mutate(letter_y = ifelse(letter_y > 1.5 & letter_y < 2.5, 2.55, letter_y))
  } else if (current_index == "connectance") {
    letter_df <- letter_df %>%
      mutate(letter_y = ifelse(letter_y > 0.3 & letter_y < 0.4, 0.405, letter_y))
  }
  
  y_limit_max <- max(letter_df$letter_y, na.rm = TRUE)
  final_y_limit_max <- max(y_params[[current_index]]$limits[2], y_limit_max) + y_range * 0.09
  y_limits <- c(y_params[[current_index]]$limits[1], final_y_limit_max)
  
  anova_x <- 0.55
  anova_y <- y_params[[current_index]]$limits[1] + y_range * 0.04
  
  p_base <- ggplot(index_data, aes(x = year, y = value, fill = year, color = year)) +
    geom_boxplot(
      width = 0.575,
      alpha = 0,
      color = "black",
      linewidth = 1.1,
      outlier.shape = NA
    ) +
    geom_point(
      position = position_jitter(width = 0.2, height = 0),
      size = 7,
      alpha = 0.7,
      aes(color = year)
    ) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.4,
      linewidth = 0.6,
      color = "black",
      alpha = 0.8
    ) +
    scale_fill_manual(values = year_colors) +
    scale_color_manual(values = year_colors) +
    scale_y_continuous(
      name = y_labels[current_index],
      labels = function(x) {
        ifelse(round(x, 0) == x, as.character(round(x)), sprintf("%.2f", x))
      },
      limits = y_limits,
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_discrete(name = "Year") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 2),
      axis.line.x = element_line(color = "black", linewidth = 1.25),
      axis.line.y = element_line(color = "black", linewidth = 1.25),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 50, color = "black"),
      axis.text.y = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.ticks.length = unit(0.3, "cm"),
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  output_file <- paste0(
    "Result_Figures/Box Plots/",
    current_index,
    ".png"
  )
  
  ggsave(
    filename = output_file,
    plot = p_base,
    width = 8,
    height = 6.75,
    dpi = 600
  )
  
  cat("\n", index_names[current_index])
  cat("\n\n", "ANOVA Result", "\n")
  cat(anova_text, "\n\n")
  cat("Tukey HSD Result", "\n\n")
  print(tukey_groups)
}

### Scatter Plot

rm(list = ls())

load("R_Datasets/metrics.rdata")

richness_variables <- c(
  "sampling_site_order",
  "phytoplankton_richness",
  "bacteria_richness",
  "fungi_richness",
  "zooplankton_richness",
  "insect_richness",
  "benthic_invertebrate_richness",
  "fish_richness",
  "multigroups_richness"
)

diversity_variables <- c(
  "sampling_site_order",
  "phytoplankton_diversity",
  "bacteria_diversity",
  "fungi_diversity",
  "zooplankton_diversity",
  "insect_diversity",
  "benthic_invertebrate_diversity",
  "fish_diversity",
  "multigroups_diversity"
)

data_richness <- data_original[, richness_variables]
data_diversity <- data_original[, diversity_variables]

all_years <- c("1st", "2nd", "3rd", "4th")

species_base_colors_diversity <- c(
  phytoplankton_diversity        = "#40904F",
  bacteria_diversity             = "#9F7C66",
  benthic_invertebrate_diversity = "#CCB254",
  fish_diversity                 = "#A44A4B",
  fungi_diversity                = "#A484A7",
  insect_diversity               = "#CC8725",
  zooplankton_diversity          = "#28466B",
  multigroups_diversity          = "#656565"
)

species_base_colors_richness <- c(
  phytoplankton_richness        = "#40904F",
  bacteria_richness             = "#9F7C66",
  benthic_invertebrate_richness = "#CCB254",
  fish_richness                 = "#A44A4B",
  fungi_richness                = "#A484A7",
  insect_richness               = "#CC8725",
  zooplankton_richness          = "#28466B",
  multigroups_richness          = "#656565"
)

generate_year_gradient <- function(base_color, years = all_years) {
  pal <- colorRampPalette(c("white", base_color))(10)
  year_cols <- pal[c(3.75, 5.75, 8.5, 10)]
  names(year_cols) <- years
  return(year_cols)
}

data_richness_processed <- data_richness %>%
  separate(sampling_site_order, into = c("site", "year_code"), sep = "_", remove = FALSE) %>%
  mutate(
    year = case_when(
      year_code == "1" ~ "1st",
      year_code == "2" ~ "2nd",
      year_code == "3" ~ "3rd",
      year_code == "4" ~ "4th",
      TRUE ~ year_code
    ),
    year = factor(year, levels = all_years),
    year_numeric = as.numeric(year)
  )

data_richness_long <- data_richness_processed %>%
  pivot_longer(
    cols = phytoplankton_richness:multigroups_richness,
    names_to = "species",
    values_to = "value"
  ) %>%
  mutate(
    species = factor(
      species,
      levels = c(
        "phytoplankton_richness",
        "bacteria_richness",
        "benthic_invertebrate_richness",
        "fish_richness",
        "fungi_richness",
        "insect_richness",
        "zooplankton_richness",
        "multigroups_richness"
      )
    )
  )

richness_breaks_list <- list(
  phytoplankton_richness        = c(30, 60, 90),
  bacteria_richness             = c(0, 250, 500),
  fungi_richness                = c(0, 40, 80),
  fish_richness                 = c(20, 40, 60),
  benthic_invertebrate_richness = c(0, 10, 20, 30),
  insect_richness               = c(0, 20, 40),
  zooplankton_richness          = c(0, 10, 20, 30),
  multigroups_richness          = c(200, 400, 600, 800)
)

richness_limits_list <- list(
  phytoplankton_richness        = c(25, 105),
  bacteria_richness             = c(0, 520),
  fungi_richness                = c(0, 80),
  fish_richness                 = c(10, 60),
  benthic_invertebrate_richness = c(0, 35),
  insect_richness               = c(0, 40),
  zooplankton_richness          = c(0, 35),
  multigroups_richness          = c(150, 800)
)

data_diversity_processed <- data_diversity %>%
  separate(sampling_site_order, into = c("site", "year_code"), sep = "_", remove = FALSE) %>%
  mutate(
    year = case_when(
      year_code == "1" ~ "1st",
      year_code == "2" ~ "2nd",
      year_code == "3" ~ "3rd",
      year_code == "4" ~ "4th",
      TRUE ~ year_code
    ),
    year = factor(year, levels = all_years),
    year_numeric = as.numeric(year)
  )

data_diversity_long <- data_diversity_processed %>%
  pivot_longer(
    cols = phytoplankton_diversity:multigroups_diversity,
    names_to = "species",
    values_to = "value"
  ) %>%
  mutate(
    species = factor(
      species,
      levels = c(
        "phytoplankton_diversity",
        "bacteria_diversity",
        "benthic_invertebrate_diversity",
        "fish_diversity",
        "fungi_diversity",
        "insect_diversity",
        "zooplankton_diversity",
        "multigroups_diversity"
      )
    )
  )

diversity_breaks_list <- list(
  phytoplankton_diversity        = c(1, 2, 3),
  bacteria_diversity             = c(1, 2, 3, 4),
  fungi_diversity                = c(0, 1, 2, 3),
  fish_diversity                 = c(1, 2, 3),
  benthic_invertebrate_diversity = c(0, 1, 2, 3),
  insect_diversity               = c(0, 1, 2),
  zooplankton_diversity          = c(0, 1, 2, 3),
  multigroups_diversity          = c(2, 3, 4, 5)
)

diversity_limits_list <- list(
  phytoplankton_diversity        = c(1, 3.2),
  bacteria_diversity             = c(1, 4.5),
  fungi_diversity                = c(0, 3.5),
  fish_diversity                 = c(0.5, 3.1),
  benthic_invertebrate_diversity = c(0, 3),
  insect_diversity               = c(0, 2.5),
  zooplankton_diversity          = c(0, 3),
  multigroups_diversity          = c(2, 5)
)

plot_point_style <- function(data_long,
                             current_species,
                             y_breaks_list,
                             y_limits_list = NULL,
                             species_base_colors,
                             file_prefix,
                             suffix_to_remove) {
  species_data <- data_long %>%
    filter(species == current_species) %>%
    drop_na(value) %>%
    mutate(year = factor(as.character(year), levels = all_years))
  
  if (nrow(species_data) == 0) {
    return(NULL)
  }
  
  base_color <- species_base_colors[[as.character(current_species)]]
  if (is.null(base_color) || is.na(base_color)) {
    base_color <- "#4C78A8"
  }
  
  year_colors <- generate_year_gradient(base_color, years = all_years)
  
  years_present <- species_data %>%
    filter(!is.na(value)) %>%
    pull(year) %>%
    as.character() %>%
    unique()
  
  n_years_present <- length(years_present)
  
  y_breaks <- y_breaks_list[[as.character(current_species)]]
  
  data_min <- min(species_data$value, na.rm = TRUE)
  data_max <- max(species_data$value, na.rm = TRUE)
  
  manual_limits <- NULL
  if (!is.null(y_limits_list)) {
    manual_limits <- y_limits_list[[as.character(current_species)]]
  }
  
  if (!is.null(manual_limits) && length(manual_limits) == 2 && all(!is.na(manual_limits))) {
    base_y_min <- manual_limits[1]
    base_y_max <- manual_limits[2]
  } else {
    y_min_break <- min(y_breaks, na.rm = TRUE)
    y_max_break <- max(y_breaks, na.rm = TRUE)
    base_y_min <- min(data_min, y_min_break, na.rm = TRUE)
    base_y_max <- max(data_max, y_max_break, na.rm = TRUE)
  }
  
  y_range <- base_y_max - base_y_min
  if (is.na(y_range) || y_range == 0) {
    y_range <- 1
  }
  
  if (!is.null(manual_limits) && length(manual_limits) == 2 && all(!is.na(manual_limits))) {
    final_y_max <- manual_limits[2]
  } else {
    final_y_max <- base_y_max + y_range * 0.08
  }
  
  y_limits <- c(base_y_min, final_y_max)
  
  center_points <- species_data %>%
    group_by(year) %>%
    summarise(
      center_value = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(year = factor(as.character(year), levels = all_years))
  
  species_label <- gsub(suffix_to_remove, "", as.character(current_species))
  
  if (n_years_present >= 2) {
    anova_ok <- FALSE
    
    try({
      lm_model <- lm(value ~ year, data = species_data)
      anova_res <- anova(lm_model)
      
      p_value <- anova_res$`Pr(>F)`[1]
      f_value <- anova_res$`F value`[1]
      f_value_text <- sprintf("F = %.2f", f_value)
      
      p_label <- if (is.na(p_value)) {
        "p = NA"
      } else if (p_value < 0.001) {
        "P < 0.001"
      } else if (p_value < 0.01) {
        "P < 0.01"
      } else if (p_value < 0.05) {
        "P < 0.05"
      } else {
        paste0("P = ", sprintf("%.3f", p_value))
      }
      
      anova_text <- paste(f_value_text, "", p_label)
      
      lsd_res <- agricolae::LSD.test(
        y = lm_model,
        trt = "year",
        alpha = 0.05,
        group = TRUE,
        console = FALSE
      )
      
      letter_df <- lsd_res$groups
      letter_df$year <- rownames(letter_df)
      rownames(letter_df) <- NULL
      
      letter_df <- letter_df %>%
        mutate(year = factor(year, levels = all_years)) %>%
        arrange(year)
      
      cat("\n", species_label)
      cat("\n\n", "ANOVA Result", "\n")
      cat(anova_text, "\n\n")
      cat("LSD Result", "\n\n")
      print(letter_df[, c("year", "groups")])
      
      anova_ok <- TRUE
    }, silent = TRUE)
    
    if (!anova_ok) {
      cat("\n", species_label)
      cat("\n\n", "ANOVA Result", "\n")
      cat("ANOVA unavailable for:", species_label, "\n")
    }
  } else {
    cat("\n", species_label)
    cat("\n\n", "ANOVA Result", "\n")
    cat("Insufficient years for ANOVA\n")
  }
  
  p <- ggplot(species_data, aes(x = year, y = value, color = year)) +
    geom_point(
      position = position_jitter(width = 0.125, height = 0),
      size = 7,
      alpha = 0.7,
      stroke = 0
    ) +
    geom_point(
      data = center_points,
      aes(x = year, y = center_value, fill = year),
      inherit.aes = FALSE,
      shape = 21,
      size = 10,
      stroke = 1.2,
      color = "black"
    ) +
    scale_color_manual(values = year_colors, drop = FALSE) +
    scale_fill_manual(values = year_colors, drop = FALSE) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = y_limits,
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 2),
      axis.line.x = element_line(color = "black", linewidth = 1.35),
      axis.line.y = element_line(color = "black", linewidth = 1.35),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 40, color = "black"),
      axis.text.y = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.ticks.length = unit(0.45, "cm"),
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  species_simple_name <- gsub(suffix_to_remove, "", as.character(current_species))
  output_file <- paste0("Result_Figures/Scatter Plots/", file_prefix, species_simple_name, ".png")
  
  ggsave(
    filename = output_file,
    plot = p,
    width = 5.5,
    height = 5.5,
    dpi = 600
  )
}

species_list_richness <- levels(data_richness_long$species)
for (current_species in species_list_richness) {
  plot_point_style(
    data_long = data_richness_long,
    current_species = current_species,
    y_breaks_list = richness_breaks_list,
    y_limits_list = richness_limits_list,
    species_base_colors = species_base_colors_richness,
    file_prefix = "richness_",
    suffix_to_remove = "_richness"
  )
}

species_list_diversity <- levels(data_diversity_long$species)
for (current_species in species_list_diversity) {
  plot_point_style(
    data_long = data_diversity_long,
    current_species = current_species,
    y_breaks_list = diversity_breaks_list,
    y_limits_list = diversity_limits_list,
    species_base_colors = species_base_colors_diversity,
    file_prefix = "diversity_",
    suffix_to_remove = "_diversity"
  )
}
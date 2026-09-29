rm(list = ls())

load("R_Datasets/abundance_by_time.rdata")

data_all <- list(data_2022, data_2023, data_2024, data_2025) %>%
  reduce(full_join, by = c("group", "species")) %>%
  mutate(across(-c(group, species), ~ replace_na(., 0))) %>%
  arrange(group, species)

rm(list = setdiff(ls(), "data_all"))

for (colname in names(data_all)[-1]) {
  tmp <- data_all[data_all[[colname]] != 0, "species", drop = FALSE]
  assign(colname, tmp, envir = .GlobalEnv)
}

rm(data_all, tmp, colname, species)

all_objects <- ls()
occurrences <- list()

for (x in all_objects) {
  if (is.data.frame(get(x))) {
    occurrences[[x]] <- get(x)
  }
}

occurrences <- lapply(occurrences, as.data.frame)

rm(list = setdiff(ls(), "occurrences"))

load("R_Datasets/adjacent_matrix.rdata")
original_network <- original_network[, -1]

diag(original_network) <- 0
original_network[is.na(original_network)] <- 0
rownames(original_network) <- colnames(original_network)

dfs <- list()
matrices <- list()
graphs <- list()

get_df <- function(x) {
  occurrence_list <- x[[ncol(x)]]
  valid_species <- intersect(occurrence_list, colnames(original_network))
  if (length(valid_species) > 0) {
    current_df <- original_network[valid_species, valid_species, drop = FALSE]
    dfs <<- c(dfs, list(current_df))
  }
}

lapply(occurrences, get_df)

get_matrix_graph <- function(x) {
  mat <- as.matrix(x)
  graph <- igraph::graph_from_adjacency_matrix(
    mat,
    mode = "directed",
    weighted = TRUE,
    diag = FALSE
  )
  none_isolated_nodes <- igraph::V(graph)[igraph::degree(graph) > 0]
  graph2 <- igraph::induced_subgraph(graph, vids = none_isolated_nodes)
  matrices <<- c(matrices, list(mat))
  graphs <<- c(graphs, list(graph2))
}

lapply(dfs, get_matrix_graph)

names(dfs) <- names(occurrences)
names(matrices) <- names(occurrences)
names(graphs) <- names(occurrences)

cal_foodweb_1 <- function(x) {
  result_df <- data.frame(
    sampling_site_order = character(0),
    connectance = numeric(0),
    max_path_length = numeric(0),
    mean_path_length = numeric(0),
    no_nodes = numeric(0)
  )
  for (y in names(x)) {
    current_graph <- x[[y]]
    metrics <- data.frame(
      sampling_site_order = y,
      connectance = edge_density(current_graph),
      max_path_length = diameter(
        current_graph,
        directed = TRUE,
        weights = E(current_graph)$weight
      ),
      mean_path_length = mean_distance(
        current_graph,
        directed = TRUE,
        weights = E(current_graph)$weight
      ),
      no_nodes = sum(igraph::degree(current_graph) > 0)
    )
    result_df <- rbind(result_df, metrics)
  }
  return(result_df)
}

network_metrics_1 <- cal_foodweb_1(graphs)

cal_foodweb_2 <- function(x, y) {
  result_df1 <- data.frame(
    sampling_site_order = character(0),
    modularity = numeric(0)
  )
  for (z in names(y)) {
    current_community <- cluster_walktrap(y[[z]])
    mod_df <- data.frame(
      sampling_site_order = z,
      modularity = modularity(current_community)
    )
    result_df1 <- rbind(result_df1, mod_df)
  }
  result_df2 <- data.frame(
    nestedness = numeric(0),
    vulnerability = numeric(0),
    robustness = numeric(0)
  )
  get_vulnerability <- function(n) {
    gen_vul <- networklevel(n, index = "vulnerability", weighted = TRUE)
    gen_vul[!grepl("generality", names(gen_vul))]
  }
  for (w in names(x)) {
    current_matrix <- x[[w]]
    net_metrics <- data.frame(
      nestedness = bipartite::nested(current_matrix),
      vulnerability = get_vulnerability(current_matrix),
      robustness = robustness(
        second.extinct(
          current_matrix,
          participant = "higher",
          nrep = 100,
          method = "random",
          details = FALSE
        )
      )
    )
    result_df2 <- rbind(result_df2, net_metrics)
  }
  return(cbind(result_df1, result_df2))
}

network_metrics_2 <- cal_foodweb_2(matrices, graphs)

network_metrics_collection <- cbind(network_metrics_1, network_metrics_2[, -1])

write.xlsx(
  network_metrics_collection,
  "Result_Tables/Foodweb Metrics.xlsx"
)
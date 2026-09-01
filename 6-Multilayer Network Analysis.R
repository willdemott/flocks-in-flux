# Tree Swallow Multilayer Network made with Infomap Ecology.
# This file constructs our empirical multiplex network, along with 100 randomized networks using the same set of edge weights, only shuffled across existing node-to-node connections. The randomized networks are then compared. Finally, network nodes and edges are exported for visualization construction in ArcPro and temporal analysis in "Module Analysis 1.R".

# Libraries 
library(tidyverse)
library(igraph)
library(infomapecology)
library(ggplot2)
library(sf)
library(scales)
library(ineq)
library(dplyr)
library(showtext)

# Add font for random histogram figure
font_add("Kigela",
         "C:/USERS/WILLD/APPDATA/LOCAL/MICROSOFT/WINDOWS/FONTS/FONNTS.COM-KIGELIA_VAI_REGULAR.OTF")
showtext_auto()

# Load data
setwd("C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/Radar Data")
nodes <- read.csv("nodes_df.csv")
edges <- read.csv("edges_df.csv") %>% rename(year_from = year) %>%
  mutate(roost_id      = str_remove(from, "_\\d{4}$"),
         roost_id_next = str_remove(to,   "_\\d{4}$"),
         year_to       = str_extract(to, "\\d{4}$") %>% as.integer())

# Build spatial and temporal tables
spatial_tbl <- edges %>%
  filter(type=="spatial") %>%
  transmute(
    year      = year_from,
    node_from = roost_id,
    node_to  = roost_id_next,
    weight = weight,
    degree = annual_degree)

temporal_tbl <- edges %>%
  filter(type=="temporal") %>%
  transmute(
    layer_from = year_from,
    node_from       = roost_id,
    layer_to    = year_to,
    node_to  = roost_id_next,
    weight = weight)

# Make list of years
years <- c(2005:2023)

# Build intralayer adjacency matrices (•᷄- •᷅ ;)
year_layers <- lapply(years, function(y) {
  df <- spatial_tbl %>% filter(year==y)
  nodes <- unique(c(
    df$node_from, 
    df$node_to,
    temporal_tbl %>% filter(layer_from==y) %>% pull(node_from),
    temporal_tbl %>% filter(layer_to==y)   %>% pull(node_to)))
  mat <- matrix(0, nrow=length(nodes), ncol=length(nodes), dimnames=list(nodes,nodes))
  for(i in seq_len(nrow(df))) {
    f <- df$node_from[i]; 
    t <- df$node_to[i]; 
    w <- df$weight[i];
    d <- df$degree
    mat[f,t] <- w
    }
  mat
  })

names(year_layers) <- as.character(years)

# Build "layer attributes"
layer_attributes <- tibble(
  layer_id   = years,
  layer_name = as.character(years))

# Create network!
ml_net <- create_multilayer_network(
  list_of_layers   = year_layers,
  interlayer_links = temporal_tbl,
  layer_attributes = layer_attributes,
  bipartite        = FALSE,
  directed         = TRUE)

# Run community-detection algorithm
infomap_res <- run_infomap_multilayer(
  M                        = ml_net,
  flow_model               = "directed",
  temporal_network         = T,
  two_level                = T,
  relax                    = T,
  multilayer_relax_rate    = 1,
  multilayer_relax_limit_up   = -1,
  multilayer_relax_limit_down = 0,
  trials                   = 50,
  silent                   = T,
  seed                     = 42)

# Extract node "strengths" for calculating gini coefficients of communities
node_strengths <- nodes %>%
  group_by(year, CLUSTER_ID) %>%
  summarise(in_strength = first(in_strength), out_strength = first(out_strength), strength = first(strength), mean_size = first(mean_size), degree = first(annual_degree), .groups = "drop") %>% 
  mutate(CLUSTER_ID = as.character(CLUSTER_ID)) %>%
  rename(layer_id = year, node_name = CLUSTER_ID)

module_strengths <- node_strengths %>%
  left_join(infomap_res$modules, by = c("layer_id", "node_name")) %>% mutate(module = ifelse(is.na(module), 0, module))

gini_by_module_year <- module_strengths %>% filter(module != 0) %>%
  group_by(layer_id, module) %>%
  filter(n_distinct(node_name) > 2) %>%
  summarise(
    gini               = ineq::Gini(strength),
    roosts_in_comm = n_distinct(node_name),
    .groups            = "drop")

# Create randomized networks!
run_null <- function(iter = 1) {
  print(c("ITERATION:", iter))
  set.seed(42 + iter)
  
  random_layers <- lapply(year_layers, function(mat) {
    non_zero_idx <- which(mat > 0, arr.ind = TRUE)
    mat[non_zero_idx] <- sample(mat[non_zero_idx])
    mat})
  
  shuffled_temporal <- temporal_tbl %>% mutate(weight = sample(weight))
  
  null_net <- create_multilayer_network(
    list_of_layers   = random_layers,
    interlayer_links = shuffled_temporal,
    layer_attributes = layer_attributes,
    bipartite        = FALSE,
    directed         = TRUE)
  
  res <- run_infomap_multilayer(
    M                           = null_net,
    flow_model                  = "directed",
    temporal_network            = TRUE,
    two_level                   = TRUE,
    relax                       = TRUE,
    multilayer_relax_rate       = 1,
    multilayer_relax_limit_up   = -1,
    multilayer_relax_limit_down = 0,
    silent                      = TRUE,
    trials                      = 50,
    seed                        = 42)
  
  node_metrics <- purrr::imap_dfr(random_layers, ~{
    mat <- .x
    tibble(
      layer_id     = as.integer(.y),
      node_name    = rownames(mat),
      out_degree   = rowSums(mat > 0),
      in_degree    = colSums(mat > 0),
      degree       = rowSums(mat > 0) + colSums(mat > 0),
      out_strength = rowSums(mat),
      in_strength  = colSums(mat),
      strength     = rowSums(mat) + colSums(mat))})
  
  module_metrics <- res$modules %>%
    left_join(node_metrics, by = c("layer_id", "node_name"))
  
  gini_by_module_year <- module_metrics %>%
    filter(module != 0) %>%
    group_by(layer_id, module) %>%
    filter(n_distinct(node_name) > 2) %>%
    summarise(
      gini               = ineq::Gini(strength),
      avg_roosts_per_com = n_distinct(node_name),
      .groups            = "drop")
  
  avg_gini_coeff <- mean(gini_by_module_year$gini)
  
  metrics <- tibble(
    iter       = iter,
    gini_coeff = avg_gini_coeff,
    codelength = res$L)
  
  list(metrics = metrics)
}


# Run randomizer function
all_results <- purrr::map(1:100, run_null)
metrics_results      <- purrr::map_dfr(all_results, "metrics") %>% dplyr::select(-iter)

# Write metrics to CSV
write.csv(metrics_results, "null_metrics.csv", row.names = FALSE)
metrics_results <- read.csv("null_metrics.csv")
summary(metrics_results)

# Empirical table
empirical_summary <- tibble(
  codelength = infomap_res$L,
  gini_coeff = mean(gini_by_module_year$gini))
empirical_summary

# Are differences between random and true significantly different (<= 2.5% of dist or >= 97.5% of dist)?
pvals <- map2_dfr(metrics_results, empirical_summary, ~{
  null_vals <- .x
  obs_val <- .y
  lower <- quantile(null_vals, 0.025)
  upper <- quantile(null_vals, 0.975)
  p_lower <- mean(null_vals <= obs_val)
  p_upper <- mean(null_vals >= obs_val)
  tibble(lower_2.5 = lower,
         upper_97.5 = upper,
         p_lower = round(p_lower, 3),
         p_upper = round(p_upper, 3),
         significant = (obs_val < lower) | (obs_val > upper))}, .id = "metric")

null_long <- metrics_results %>% pivot_longer(cols = c(codelength, gini_coeff), names_to = "metric", values_to = "value")
emp_long <- empirical_summary %>% pivot_longer(cols = c(codelength, gini_coeff), names_to = "metric", values_to = "obs_value")

label_names <- c(codelength  = "Codelength", gini_coeff  = "Gini Coefficient")

# Define a Nature-friendly color palette
palette <- c("#6b95c8", "#df9b56")

p <- ggplot(null_long, aes(x = value, fill = metric)) +
  geom_histogram(bins = 30, alpha = 0.9) +
  geom_vline(
    data = emp_long,
    aes(xintercept = obs_value),
    color = "#c56d68",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_x",
    ncol = 2,
    labeller = as_labeller(label_names)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = palette) +
  labs(x = "Metric value", y = "Frequency") +
  theme_classic(base_size = 10, base_family = "Kigela") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = rel(1.2),
      hjust = .5,
      margin = margin(b = 6, t = 6)),
    axis.line  = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.text  = element_text(color = "black", size = rel(0.9)),
    axis.title = element_text(face = "bold", size = rel(1)),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4)),
    panel.spacing.x = unit(0.6, "lines"),
    plot.margin = margin(10, 12, 8, 8)
  )
p
ggsave("random_histogram.png", p, width = 180, height = 100, units = "mm", dpi = 600)

# Make modules_df
modules_df <- infomap_res$modules %>% 
  as_tibble() %>% rename(year = layer_id, CLUSTER_ID = node_name)

annual_density <- spatial_tbl %>%
  group_by(year) %>%
  summarise(
    n_nodes = n_distinct(c(node_from, node_to)),
    n_edges = sum(node_from != node_to),
    n_potential_edges = n_nodes * (n_nodes - 1) / 2,
    density = n_edges / n_potential_edges,
    .groups = "drop")

plot(annual_density$year, annual_density$density)

# Calculate Clustering Coefficient for each year
annual_clustering <- lapply(names(year_layers), function(y) {
  g <- graph_from_adjacency_matrix(year_layers[[y]], mode = "directed", weighted = TRUE)
  # Calculate Transitivity (Global Clustering Coefficient)
  cc <- transitivity(g, type = "global")
  tibble(year = as.integer(y), clustering_coeff = cc)
}) %>% bind_rows()

# Final nodes with modules & flow
connected_nodes <- edges %>%
  filter(type == "spatial")

nodes_final <- nodes %>%
  mutate(CLUSTER_ID = as.character(CLUSTER_ID)) %>%
  left_join(modules_df, by = c("CLUSTER_ID", "year")) %>%
  mutate(across(c(annual_degree, strength, in_strength, 
                  out_strength, flow, module), ~replace_na(.x,0))) %>%
  group_by(year, module) %>%
  mutate(n_module_roosts = ifelse(module != 0, n_distinct(CLUSTER_ID), NA), 
         gini_coeff = ifelse(module == 0 | n_module_roosts <= 2, NA, ineq::Gini(annual_degree))) %>%
  ungroup() %>% left_join(annual_density, by = "year") %>% left_join(annual_clustering, by = "year")

######### EXPORTING ##############
Mode <- function(x) {
  x <- x[x != 0]
  if (length(x) == 0) return(0)
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}
nodes_final <- nodes_final %>% group_by(CLUSTER_ID) %>% mutate(mode_module = Mode(module)) %>% ungroup()

write.csv(nodes_final, "nodes_final.csv")

# turn into sf (keep all columns as attributes)
nodes_sf <- st_as_sf(
  nodes_final,
  coords    = c("XCoord", "YCoord"),
  crs       = 4326,
  remove    = F)

# Spatial edges with all attributes
all_edges_sf <- edges %>%
  filter(type == "spatial") %>%
  left_join(nodes_final %>% dplyr::select(name, year, XCoord, YCoord, module),
            by = c("from" = "name", "year_from" = "year")) %>%
  rename(from_x = XCoord, from_y = YCoord, from_module = module) %>%
  left_join(nodes_final %>% dplyr::select(name, year, XCoord, YCoord, module),
            by = c("to" = "name", "year_from" = "year")) %>%
  rename(to_x   = XCoord, to_y = YCoord, to_module = module, year = year_from) %>%
  group_by(CLUSTER_ID_from, CLUSTER_ID_to, year) %>% 
  summarize(dist_km = first(dist_km),
            from_x = first(from_x),
            from_y = first(from_y),
            to_x = first(to_x),
            to_y = first(to_y),
            from_module = first(from_module),
            to_module = first(to_module)) %>%
  ungroup() #%>% mutate(geometry = st_sfc(lapply(seq_len(n()), function(i) {st_linestring(matrix(c(from_x[i], from_y[i], to_x[i], to_y[i]),ncol = 2, byrow = T))}),crs = 4326)) %>% st_as_sf()

write.csv(all_edges_sf,"final_weightless_edges.csv")

summary_edges <- edges %>%
  filter(type == "spatial") %>%
  left_join(nodes_final %>% dplyr::select(name, year, XCoord, YCoord, mode_module),
            by = c("from" = "name", "year_from" = "year")) %>%
  rename(from_x = XCoord, from_y = YCoord, from_module = mode_module) %>%
  left_join(nodes_final %>% dplyr::select(name, year, XCoord, YCoord, mode_module),
            by = c("to" = "name", "year_from" = "year")) %>%
  rename(to_x   = XCoord, to_y = YCoord, to_module = mode_module, year = year_from) %>%
  group_by(CLUSTER_ID_from, CLUSTER_ID_to) %>% 
  summarize(dist_km = first(dist_km),
            from_x = first(from_x),
            from_y = first(from_y),
            to_x = first(to_x),
            to_y = first(to_y),
            from_module = first(from_module),
            to_module = first(to_module)) %>%
  mutate(
    geometry = st_sfc(
      lapply(seq_len(n()), function(i) {
        st_linestring(matrix(
          c(from_x[i], from_y[i], to_x[i], to_y[i]),
          ncol = 2, byrow = T))}),
      crs = 4326)) %>%
  st_as_sf()

st_write(all_edges_sf,"C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/NE ROOSTS GIS/all_edges.gpkg",layer = "edges", delete_layer = T)
st_write(summary_edges,"C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/NE ROOSTS GIS/mode_edges.gpkg",layer = "summaryedges", delete_layer=T)
st_write(nodes_sf,"C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/NE ROOSTS GIS/nodes.gpkg",  layer = "nodes", delete_layer = T)




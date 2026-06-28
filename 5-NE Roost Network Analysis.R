# Here, we construct our nodes and edges for the multiplex network. Note the network is not made here, only the parameters for it. This file also produces weight-distribution histograms for varying roost-switching rates.

# Libraries
library(dplyr)
library(purrr)
library(tidyverse)
library(sf)
library(units)

# Load & prep
setwd("C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/Radar Data")

# Years prior to 2005 had missing radar-station data, so we limited our sample period to >2005.
roosts <- read.csv("NE_CLEANED_CLUSTERS.csv") %>% filter(tres_ratio >= .9 ) %>% filter(year > 2004) %>%
  mutate(tres_count = max_scan_sum * tres_ratio) 

# Summarize per‐roost/year (size + season window)
#sizes_by_year
nodes_multilayer <- roosts %>%
  group_by(CLUSTER_ID, year) %>%
  summarise(
    size      = mean(tres_count, na.rm = TRUE),
    mean_size = mean(tres_count),
    max_size = max(tres_count),
    annual_sum = sum(tres_count),
    station = first(station),
    XCoord = first(XCoord),
    YCoord = first(YCoord),
    .groups   = "drop"
  ) %>% mutate(name = paste0(CLUSTER_ID, "_", year)) %>% ungroup()
plot(nodes_multilayer$annual_sum, nodes_multilayer$mean_size)
#sf list to append weights
annual_nodes_sf <- nodes_multilayer %>% filter(!is.na(CLUSTER_ID)) %>%
  st_as_sf(coords = c("XCoord","YCoord"), crs = 4326) %>%
  split(.$year)

# Create spatial edges and weights!!
switch_rate <- 0.22
make_spatial_edges <- function(nodes_sf, sizes_df, switch_rate) {
  
  all_ids <- as.integer(nodes_sf$CLUSTER_ID)
  year <- nodes_sf$year[1]
  
  D   <- st_distance(nodes_sf) %>% set_units("km") %>% drop_units()
  idx <- which(D <= 60 & row(D) != col(D), arr.ind = TRUE)
  
  edge_base <- tibble(
    CLUSTER_ID_from = nodes_sf$CLUSTER_ID[idx[, 1]],
    CLUSTER_ID_to   = nodes_sf$CLUSTER_ID[idx[, 2]],
    year            = nodes_sf$year[1],
    dist_km         = D[idx],
    base_p          = exp(-0.045 * D[idx]))
  
  src_days <- sizes_df %>%
    dplyr::select(CLUSTER_ID, doy, tres_count) %>%
    rename(CLUSTER_ID_from = CLUSTER_ID)
  
  dst_days <- sizes_df %>%
    dplyr::select(CLUSTER_ID, doy) %>%
    rename(CLUSTER_ID_to = CLUSTER_ID)
  
  flows <- edge_base %>%
    inner_join(src_days, by = "CLUSTER_ID_from", relationship = "many-to-many") %>%
    semi_join(dst_days, by = c("CLUSTER_ID_to", "doy")) %>%
    group_by(CLUSTER_ID_from, year, doy) %>%
    mutate(denom = sum(base_p, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      day_out = switch_rate * tres_count,
      alloc   = ifelse(denom > 0, (day_out * base_p) / denom, 0)
    ) %>%
    group_by(CLUSTER_ID_from, CLUSTER_ID_to, year) %>%
    summarise(
      dist_km = first(dist_km),
      weight  = sum(alloc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      type = "spatial",
      from = paste0(CLUSTER_ID_from, "_", year),
      to   = paste0(CLUSTER_ID_to,   "_", year)
    ) %>%
    filter(weight > 0) %>%
    dplyr::select(from, to, CLUSTER_ID_from, CLUSTER_ID_to, year, dist_km, weight, type)
  
  # Identify isolates: nodes with no incident edges in 'flows'
  incident_ids <- unique(c(flows$CLUSTER_ID_from, flows$CLUSTER_ID_to))
  isolate_ids  <- setdiff(all_ids, incident_ids)
  
  isolates <- tibble(
    from            = paste0(isolate_ids, "_", year),
    to              = paste0(isolate_ids, "_", year),
    CLUSTER_ID_from = isolate_ids,
    CLUSTER_ID_to   = isolate_ids,
    year            = year,
    dist_km         = NA_real_,
    weight          = 0,
    type            = "spatial"
  )
  bind_rows(flows, isolates)
}

# Run funtion that makes spatial edges
intralayer_edges <- purrr::imap_dfr(
  annual_nodes_sf,
  ~ make_spatial_edges(.x, sizes_df = roosts %>% filter(year == .y), switch_rate = 0.22)
) %>% mutate(CLUSTER_ID = CLUSTER_ID_from) %>% group_by(CLUSTER_ID, year) %>% mutate(annual_degree = n()) %>% ungroup()

# Compute edges for different switch rates
switch_rates <- c(0.10, 0.20, 0.30, 0.40)
rate_labels  <- c("10%", "20%", "30%", "40%")

weight_combined <- purrr::map2_dfr(switch_rates, rate_labels, function(sr, label) {
  edges <- purrr::imap_dfr(
    annual_nodes_sf,
    ~ make_spatial_edges(.x, sizes_df = roosts %>% filter(year == .y), switch_rate = sr)
  )
  edges %>%
    filter(weight > 0) %>%
    transmute(
      log_weight = log(weight),
      group      = label
    )
}) %>%
  mutate(group = factor(group, levels = rate_labels))

# --- Faceted histogram ---
rate_label_names <- setNames(
  paste0("Switch Rate: ", rate_labels),
  rate_labels
)

ggplot(weight_combined, aes(x = log_weight, fill = group)) +
  geom_histogram(bins = 500, alpha = 0.9) +
  facet_wrap(
    ~ group,
    scales = "free_y",
    ncol   = 2,
    labeller = as_labeller(rate_label_names)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("#6b95c8", "#6b95c8", "#6b95c8", "#6b95c8")) +
  labs(x = "log(weight)", y = "Count") +
  theme_classic(base_size = 10, base_family = "Kigela") +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(
      face = "bold",
      size = rel(1.2),
      hjust = 0.5,
      margin = margin(b = 6, t = 6)
    ),
    axis.line  = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.text  = element_text(color = "black", size = rel(0.9)),
    axis.title = element_text(face = "bold", size = rel(1)),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4)),
    panel.spacing.x = unit(0.6, "lines"),
    plot.margin = margin(10, 12, 8, 8)
  )
##################

# Inter‐year edges (year‐to‐year flux = next year’s size)
interlayer_edges <- nodes_multilayer %>%
  arrange(CLUSTER_ID, year) %>%
  group_by(CLUSTER_ID) %>%
  mutate(next_year = lead(year), 
         sum_next = lead(annual_sum),
         n_years = n_distinct(year)) %>%
  filter(!is.na(next_year)) %>%
  transmute(
    from              = paste(CLUSTER_ID, year,      sep = "_"),
    to                = paste(CLUSTER_ID, next_year, sep = "_"),
    CLUSTER_ID        = CLUSTER_ID,
    CLUSTER_ID_from   = CLUSTER_ID,
    CLUSTER_ID_to     = CLUSTER_ID,
    year              = year,
    dist_km           = NA_real_,
    weight            = sum_next * .1,
    type              = "temporal"
  )  %>% ungroup()

# Final tables for multilayer analysis:
edges_df <- bind_rows(intralayer_edges, interlayer_edges)

# nodes_df
node_strengths <- intralayer_edges %>%
  group_by(CLUSTER_ID_from, year) %>%
  mutate(out_strength = sum(weight)) %>% 
  group_by(CLUSTER_ID_to, year) %>%
  mutate(in_strength = sum(weight)) %>% 
  group_by(CLUSTER_ID, year) %>%
  summarise(in_strength = first(in_strength),
            out_strength = first(out_strength),
            strength = first(in_strength) + first(out_strength),
            annual_degree = n(),
            CLUSTER_ID = first(CLUSTER_ID)) %>% ungroup()
  
nodes_df <- left_join(nodes_multilayer, node_strengths, by = c("CLUSTER_ID", "year")) %>% ungroup()

# Export final node & edge tables
write.csv(nodes_df, "nodes_df.csv", row.names = FALSE)
write.csv(edges_df, "edges_df.csv", row.names = FALSE)

# Now we can move onto "Multilayer Network Analysis.R" to construct the multiplex network using our nodes and edges.
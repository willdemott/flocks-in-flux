# NE Roost Cluster Cleaning
# Here we used 2-dimensional kernel-density estimation to remove erroneous roost-detections that fall outside of the standard DOY x size distribution for each roost site.

#set wd
setwd("C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/Radar Data")

#read clusters
raw_clusters <- read.csv("NE_ROOSTS_RAW_CLUSTERS.csv")

library(dplyr)
library(ggplot2)
library(MASS)
library(fields)
library(lubridate)

# Remove roosts outside of Tree Swallow migratory season
valid_clusters <- raw_clusters %>%
  filter(CLUSTER_ID != -1, doy >= 186 & doy <= 333)
plot(valid_clusters$adjusted_bird_count)
summary(valid_clusters$adjusted_bird_count)

# Filter impossible outliers by keeping lower 98% of roosts by size
q98 <- quantile(valid_clusters$adjusted_bird_count, 0.98, na.rm = TRUE)
valid_clusters <- valid_clusters %>%
  filter(birds < q98)
plot(valid_clusters$adjusted_bird_count)

cluster_sizes <- valid_clusters %>% count(CLUSTER_ID, name = "n_points")

ggplot(cluster_sizes, aes(x = n_points)) +
  geom_histogram(binwidth = 10, fill = "steelblue", colour = "black") +
  labs(title = "Distribution of Cluster Sizes",
    x     = "Number of points in cluster",
    y     = "Number of clusters") +
  theme_minimal()

# Filter out roosts that are less than 90% Tree Swallows, determined via 
raw_year_data <- valid_clusters %>% filter(tres_ratio > 0.9) %>% group_by(year) %>%summarise(n_birds = sum(max_scan_sum))
ggplot(raw_year_data, aes(x = year, y = n_birds)) +
  geom_point() +
  geom_smooth()

### ELIMINATE OUTLIERS USING 2D-Kernel density analysis ##
## Helper to clamp values
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# KDE parameter function: uses log10(adjusted_bird_count) range for grid. get_kde_params values were set by manually testing different values until sample outputs (roost sites with varying sizes and occurrences) were appropriately filtered
get_kde_params <- function(cluster) {
  n_points <- nrow(cluster)
  anchor_n      <- c(5, 50, 100, 250, 500, 1000)
  anchor_outpct <- c(0.1, 0.05, 0.05, 0.05, 0.05, 0.01)
  outlier_raw   <- approx(anchor_n, anchor_outpct, xout = n_points, rule = 2)$y
  outlier_pct   <- clamp(outlier_raw, 0.01, 0.40)
  cluster       <- cluster %>% mutate(log_count = log10(adjusted_bird_count + 1))
  range_log     <- diff(range(cluster$log_count, na.rm = TRUE))
  n_grid        <- clamp(round(range_log * 50), 25, 200)
  list(
    n_grid      = n_grid,
    outlier_pct = outlier_pct)}

# Example for manual interpretation. The for loop for printing figures of all roost sites is at the end of the script, which was also used to tune parameters
cluster <- valid_clusters %>%
  filter(CLUSTER_ID == 208) %>%
  mutate(log_count = log10(adjusted_bird_count + 1))
params <- get_kde_params(cluster)

# Compute 2D KDE on (doy, log_count)
dens <- kde2d(
  x = cluster$doy,
  y = cluster$log_count,
  n = params$n_grid)

# Interpolate densities back to points
interp_dens <- interp.surface(
  obj = dens,
  loc = cbind(cluster$doy, cluster$log_count))
cluster$density <- interp_dens

# Flag lowest-density outliers
threshold        <- quantile(cluster$density, params$outlier_pct, na.rm = TRUE)
cluster$outlier <- cluster$density < threshold

# Plot with outliers highlighted
ggplot(cluster, aes(x = doy, y = adjusted_bird_count, color = outlier)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(
    title = "Cluster 50 — 2D KDE Outlier Detection",
    x = "Day of Year",
    y = "Adjusted Bird Count",
    color = "Outlier"
  ) +
  theme_minimal()

# When satisfied with parameter tuning, clean and filter
cleaned_clusters <- valid_clusters %>%
  mutate(log_count = log10(adjusted_bird_count + 1)) %>%
  group_split(CLUSTER_ID) %>%
  purrr::map_dfr(function(cluster) {
    params <- get_kde_params(cluster)
    kde <- kde2d(x = cluster$doy, y = cluster$log_count, n = params$n_grid)
    interp <- interp.surface(kde, cbind(cluster$doy, cluster$log_count))
    cluster$density <- interp
    threshold <- quantile(interp, params$outlier_pct, na.rm = TRUE)
    cluster %>% filter(density >= threshold)
  })

# Filter clusters with >= 10 points
cleaned_filtered_clusters <- cleaned_clusters %>%
  group_by(CLUSTER_ID) %>%
  filter(n() >= 10)

cleaned_cluster_list <- sort(unique(cleaned_clusters$CLUSTER_ID))
cleaned_filtered_cluster_list <- sort(unique(cleaned_filtered_clusters$CLUSTER_ID))
tossed_clusters <- sort(setdiff(cleaned_cluster_list, cleaned_filtered_cluster_list))
tossed_clusters

# Create output directory
#dir.create("cluster_doy_bird_plots", showWarnings = FALSE)
# Manually screen the rest for abnormally large/unevenly distributed roosts
for (cid in cleaned_filtered_cluster_list) {
  cluster <- cleaned_filtered_clusters %>% filter(CLUSTER_ID == cid)
  p <- ggplot(cluster, aes(x = doy, y = adjusted_bird_count, colour = factor(year))) +
    geom_point(alpha = 0.7, size = 1.8) +
    labs(
      title = paste("Cluster", cid, "– Bird Counts over DOY"),
      x     = "Day of Year",
      y     = "Bird Count",
      colour = "Year") +
    scale_colour_viridis_d(option = "mako") +
    theme_bw(base_size = 25)
  ggsave(
    filename = paste0("cluster_doy_bird_plots/supfig4_", cid, ".png"),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300)}

newlist <- c(95, 162, 58, 122)
panel_labels <- setNames(LETTERS[1:4], as.character(newlist))
plot_data <- cleaned_filtered_clusters %>%
  filter(CLUSTER_ID %in% newlist) %>%
  mutate(CLUSTER_ID = factor(as.character(CLUSTER_ID), 
                             levels = as.character(newlist)))

p <- ggplot(plot_data, aes(x = doy, y = adjusted_bird_count, colour = "#6b95c8")) +
  geom_point(alpha = 0.7, size = 1.5, color = "#6b95c8") +
  facet_wrap(~ CLUSTER_ID, scales = "free_y",
             labeller = labeller(CLUSTER_ID = panel_labels)) +
  labs(
    x      = "Day of Year",
    y      = "Bird Count") +
  theme_minimal(base_size = 50) +
  theme(legend.position = "right")

ggsave(
  filename = "cluster_doy_bird_plots/supfig4_faceted.png",
  plot     = p,
  width    = 16,
  height   = 10,
  dpi      = 300)

# Filter out manually-determined non-roosts
cleaned_filtered_clusters <- cleaned_filtered_clusters %>%
  filter(!CLUSTER_ID %in% c(1, 2, 5, 6, 11, 20, 
                            30, 45, 60, 69, 70, 93, 94, 135, 166))

write.csv(cleaned_filtered_clusters, "NE_CLEANED_CLUSTERS.csv")

# Now we can finally start to make our network in "NE Roost Network Analysis.R" with "NE_CLEANED_CLUSTERS.csv".
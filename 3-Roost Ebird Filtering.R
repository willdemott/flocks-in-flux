setwd("C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/Radar Data")

roosts <- read.csv("NE_ROOSTS.csv")
head(roosts)

library(terra)
library(dplyr)
library(ebirdst)
key <- "41ksdggr9mie" 
set_ebirdst_access_key(key)
# Define WGS84 projection
wgs84_crs <- "EPSG:4326"

####### CREATE POISSON BUFFERS IN RASTERS ###########

# Download the data
treswa_data <- ebirdst_download_status(species = "treswa")
purmar_data <- ebirdst_download_status(species = "purmar")

# Load the weekly abundance raster
purmar_abundance_raster <- load_raster("purmar", product = "abundance", period = "weekly", resolution = "27km")
treswa_abundance_raster <- load_raster("treswa", product = "abundance", period = "weekly", resolution = "27km")

# Reproject the rasters (switch to _smoothed if you want buffered rasters)
purmar_wgs84 <- project(purmar_abundance_raster, wgs84_crs)
treswa_wgs84 <- project(treswa_abundance_raster, wgs84_crs)

raster_dates <- names(purmar_wgs84)

# Function to find the most recent Tuesday in 2022 for a given roost date
find_raster_week <- function(roost_date, raster_dates) {
  roost_date <- as.Date(roost_date)  # Convert to Date format
  
  # Convert roost date to 2022 (ignoring the actual year)
  roost_date_2022 <- as.Date(paste0("2022-", format(roost_date, "%m-%d")))
  
  # Find most recent Tuesday
  offset <- as.integer(format(roost_date_2022, "%w")) - 2  # Days since last Tuesday
  if (offset < 0) offset <- offset + 7  # Move back a full week if before Tuesday
  raster_week <- as.character(roost_date_2022 - offset)  # Compute correct raster week
}
roosts$date <- as.Date(as.character(roosts$date), format = "%Y%m%d")
roosts$closest_week <- sapply(roosts$date, function(d) find_raster_week(d, raster_dates))

# Convert roost data to spatial points (coordinates)
roosts_sp <- vect(roosts[, c("lon", "lat")], crs = "EPSG:4326")

# Run for loop appending tres & puma abundances with matching dates and locations
roosts$puma_abundance <- NA
roosts$tres_abundance <- NA

for (i in 1:nrow(roosts)){
  week_date <- roosts$closest_week[i]
  
  purmar_raster_week <- purmar_wgs84[[which(names(purmar_wgs84) == week_date)]]
  treswa_raster_week <- treswa_wgs84[[which(names(treswa_wgs84) == week_date)]]
  
  roost_coords <- c(roosts$lon[i], roosts$lat[i])
  roost_coords_spat <- vect(data.frame(lon = roost_coords[1], lat = roost_coords[2]), 
                            geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # Now try extracting values from the raster
  puma_abundance_value <- terra::extract(purmar_raster_week, roost_coords_spat)[ , 2]
  treswa_abundance_value <- terra::extract(treswa_raster_week, roost_coords_spat)[ , 2]
  
  roosts$puma_abundance[i] <- puma_abundance_value
  roosts$tres_abundance[i] <- treswa_abundance_value
}

# Calculate the total abundance for each roost
roosts$total_abundance <- roosts$puma_abundance + roosts$tres_abundance

# Calculate the species ratios
roosts$puma_ratio <- roosts$puma_abundance / roosts$total_abundance
roosts$tres_ratio <- roosts$tres_abundance / roosts$total_abundance

# Constants for bird weights (in grams)
tres_weight <- 21.2  # g (Tree Swallow)
puma_weight <- 53.8  # g (Purple Martin)

# Calculate the adjusted bird count
roosts$adjusted_bird_count <- roosts$birds * (
  roosts$tres_ratio + (roosts$puma_ratio * tres_weight / puma_weight))

roosts$tres_count <- roosts$adjusted_bird_count * roosts$tres_ratio
roosts$puma_count <- roosts$adjusted_bird_count * roosts$puma_ratio

roosts$adjusted_max_count <- roosts$max_scan_sum * (
  roosts$tres_ratio + (roosts$puma_ratio * tres_weight / puma_weight))
roosts$max_tres_count <- roosts$adjusted_max_count * roosts$tres_ratio
# View the updated dataframe
head(roosts)

# Make TRES/PUMA ratio columns
# Filter out rows where both puma_abundance and tres_abundance are 0
roosts_filtered <- roosts %>%
  filter(!(puma_abundance == 0 & tres_abundance == 0))

roosts_filtered <- read.csv("NE_ROOSTS_EBIRD_FILTERED_27km.csv")
ggplot(roosts_filtered, aes(tres_ratio, fill = factor(round(tres_ratio, digits = 1)))) + 
  geom_histogram(color = "white", binwidth = 0.1, show.legend = F) + 
  scale_fill_discrete(h = c(290, 250), c = 150, l = 50) +
  labs(x = "Ratio of Tree Swallows to Purple Martins", y = "Roost Detections") +
  theme_minimal() + ggtitle("Roost Species Abundance Ratio Counts")

roosts_filtered %>%
  summarise(
    total_roosts = n(),
    below_thresh = sum(tres_ratio < 0.9, na.rm = TRUE),
    percentage = (below_thresh / total_roosts) * 100
  )
write.csv(roosts, "NE_ROOSTS_EBIRD_UNFILTERED_27km.csv")
write.csv(roosts_filtered, "NE_ROOSTS_EBIRD_FILTERED_27km.csv")

# Move "NE_ROOSTS_EBIRD_UNFILTERED_27km.csv" to ArcPro to create cluster ID's and coords using Heirarchical Density-based Clustering with minimum features per cluster set to 10.
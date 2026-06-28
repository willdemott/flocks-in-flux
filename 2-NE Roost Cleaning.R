# Now that we have raw roosts with estimates of bird counts, we need to combine all station-year .csv's and do initial steps for cleaning the data. 
# Here we will be filtering out jumpy-track frames, removing duplicate overlapping tracks within stations, combining tracks of the same roost with different track id's, and removing duplicate roosts between radar stations (keeping the detection from the radar thats closer)

# Load libraries
library(dplyr)         
library(raster)
library(sf)
library(scales)
library(lubridate)
library(stringr)
library(geosphere)

# Set home directory and load data
setwd("C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/Radar Data/")

#Load data
roost_files <- list.files(path = "./Counted Roosts", pattern = ".csv", full.names = T)

# See column differences in old KOKX data that was manually screened, and exemplary new file
KOKX2017 <- read.csv(roost_files[415])
KAKQ1995 <- read.csv(roost_files[1])
str(KOKX2017)
str(KAKQ1995)

# Create and apply function that removes OL_ columns (a dead-end product I created that accounted for overlapping roosts) and updates label columns to match
roost_list <- lapply(roost_files, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # If KOKX, drop OL_ columns
  if (grepl("KOKX", file)) {
    df <- df[ , !grepl("^OL_", names(df))]
  }
  
  # Coerce problem columns to character
  if ("viewed" %in% names(df)) df$viewed <- as.character(df$viewed)
  if ("user_labeled" %in% names(df)) df$user_labeled <- as.character(df$user_labeled)
  
  return(df)
})

# Get unified column set
all_cols <- unique(unlist(lapply(roost_list, names)))

# Fill in missing columns & align
roost_list <- lapply(roost_list, function(df) {
  df[setdiff(all_cols, names(df))] <- NA
  df <- df[all_cols]
  df
})

# Bind all data
raw_roosts <- bind_rows(roost_list)
n_distinct(raw_roosts$track_id)
# Test for the missing point
track_test <- "KOKX20170826-5"
cat("Check after STEP_NAME:\n")
print(any(raw_roosts$track_id == track_test))

################ FINALLY move on to processing!################
################

# Filter out KOKX rows where the label is "non-roost" (all other station-years did not count non-roosts)
roosts <- raw_roosts %>%  filter(!(label %in% c("non-roost")))

# We first need to create a function that calculates the distance between the lat and lon points of each roost

# Define function to calculate distance in km between two lat/lon sets
calc_dist_km <- function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # meters to km
}

############# REMOVE JUMPY FRAMES

# Set max allowable movement between consecutive frames
max_jump_km <- 40
# Remove track frames that jump
roosts <- roosts %>%
  group_by(track_id) %>%
  arrange(time) %>%
  group_split() %>%
  lapply(function(df) {
    if (nrow(df) < 2) return(df)
    keep_idx <- 1  # always keep the first row
    for (i in 1:(nrow(df) - 1)) {
      dist_km <- calc_dist_km(df$lon[i], df$lat[i], df$lon[i + 1], df$lat[i + 1])
      
      if (!is.na(dist_km) && dist_km <= max_jump_km) {
        keep_idx <- c(keep_idx, i + 1)
      } else {
        cat("Removed row in", df$track_id[i + 1], "- jump:", round(dist_km, 2), "km\n")
      }
    }
    df[sort(unique(keep_idx)), ]
  }) %>%
  bind_rows() %>%
  ungroup()


################### Remove duplicate tracks (true spatial/temporal overlap)
roosts <- roosts %>%
  group_by(local_time, station) %>%
  arrange(desc(length)) %>%
  group_split() %>%
  lapply(function(df) {
    n <- nrow(df)
    remove_indices <- c()
    for (i in seq_len(n - 1)) {
      for (j in seq((i + 1), n)) {
        distance <- calc_dist_km(df$lon[i], df$lat[i], df$lon[j], df$lat[j])
        if (distance <= 2) {
          if (df$length[i] < df$length[j]) {
            cat("Removing", df$track_id[i], "\n")
            remove_indices <- c(remove_indices, i)
          } else {
            cat("Removing", df$track_id[j], "\n")
            remove_indices <- c(remove_indices, j)
          }
        }
      }
    }
    if (length(remove_indices) > 0) {
      df <- df[-unique(remove_indices), ]
    }
    return(df)
  }) %>%
  bind_rows() %>%
  ungroup()


################### Synchronize double-tracks (two unique track IDs covering the same roost at different time frames)

# Helper to split track_id into base and suffix
split_track_id <- function(id) {
  parts <- stringr::str_match(id, "^(.*)-(\\d+)$")
  list(base = parts[, 2], suffix = parts[, 3])
}

# Synchronize track IDs for overlapping tracks from the same station/date
roosts <- roosts %>%
  group_by(date, station) %>%
  arrange(time, desc(length)) %>%
  group_split() %>%
  lapply(function(df) {
    n <- nrow(df)
    for (i in seq_len(n - 1)) {
      for (j in seq(i + 1, n)) {
        
        # Skip if already same track_id and/or time
        if (df$track_id[i] == df$track_id[j]) next
        if (df$local_time[i] == df$local_time[j]) next
        
        # Calculate distance
        distance <- calc_dist_km(df$lon[i], df$lat[i], df$lon[j], df$lat[j])
        
        if (!is.na(distance) && distance <= 3) {
           # Always keep i, merge j into i
          base_id    <- split_track_id(df$track_id[i])$base
          new_suffix <- split_track_id(df$track_id[i])$suffix
          old_id     <- df$track_id[j]
          new_id     <- paste0(base_id, "-", new_suffix)
          cat("Merged", old_id, "into", new_id, "\n")
          df$track_id[df$track_id == old_id] <- new_id
        }
      }
    }
    return(df)
  }) %>%
  bind_rows() %>%
  ungroup()

###### SUMMING AND CONDENSING ROOSTS ##########

# Remove additional track frames on tracks with >5 frames
roosts_sliced <- roosts %>%
  group_by(track_id) %>%
  arrange(time, .by_group = TRUE) %>%
  slice_min(order_by = row_number(), n = 5) %>%
  ungroup()

# Identify the "sum" columns and ensure they're numeric
sum_cols <- grep("^sum", names(roosts_sliced), value = TRUE)
numeric_sum_cols <- sum_cols[sapply(roosts_sliced[, sum_cols, drop = FALSE], is.numeric)]

# Compute the row-wise sum for the selected columns and create 'scan_sum'
roosts_sliced$scan_sum <- rowSums(roosts_sliced[, numeric_sum_cols, drop = FALSE], na.rm = TRUE)


head(roosts_sliced)

# Add scan order within each track_id
roosts_ordered <- roosts_sliced %>%
  arrange(track_id, from_sunrise) %>%
  group_by(track_id) %>%
  mutate(scan_num = row_number()) %>%
  ungroup()

# Reshape to long format and compute proportions
long_props <- roosts_ordered %>%
  dplyr::select(track_id, scan_num, starts_with("sum_"), scan_sum) %>%
  pivot_longer(
    cols = starts_with("sum_"),
    names_to = "bin",
    values_to = "bin_sum"
  ) %>%
  mutate(proportion = bin_sum / scan_sum)

# Calculate average proportion per bin per scan number
avg_props_by_scan <- long_props %>%
  group_by(scan_num, bin) %>%
  summarise(avg_proportion = mean(proportion, na.rm = TRUE), .groups = "drop")

# Plot
ggplot(avg_props_by_scan, aes(x = bin, y = avg_proportion, fill = as.factor(scan_num))) +
  geom_col(position = "dodge", color = "white") +
  labs(
    title = "Change in Vertical Distribution Across Scans",
    x = "Altitude Bin",
    y = "Average Proportion of scan_sum",
    fill = "Track Number"
  ) +
  theme_minimal()

# Filter for scans 1 and 5
scan_change <- avg_props_by_scan %>%
  filter(scan_num %in% c(1, 5)) %>%
  pivot_wider(names_from = scan_num, values_from = avg_proportion, names_prefix = "scan_") %>%
  mutate(change = scan_5 - scan_1)

# View the change
print(scan_change)


roosts_high_alt <- roosts_sliced %>%
  mutate(
    high_alt_sum = sum_2.5 + sum_3.5 + sum_4.5,
    high_alt_prop = high_alt_sum / scan_sum
  )

# Summary shows us the mean high altitude proportion is 5%, which is likely biased by weather and noise
# We will use the median (0.001%) to justify not counting higher elevation angles
summary(roosts_high_alt$high_alt_prop)

# Compute row-wise sum for the first two elevation angles
roosts_sliced$true_scan_sum <- rowSums(
  roosts_sliced[, c("sum_0.5", "sum_1.5")],
  na.rm = TRUE
)

# Remove all but the first row of each unique track_id, keeping all columns
roosts_condensed <- roosts_sliced %>%
  group_by(track_id) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(birds = sum(true_scan_sum, na.rm = TRUE),
         max_scan_sum = max(true_scan_sum)) %>%
  slice(1) %>%
  ungroup()

plot(roosts_condensed$birds)

############################## REMOVING DUPLICATE ROOSTS BETWEEN STATIONS #

# Create Thresholds (tracks will appear within 10 minutes of each other)
time_threshold <- 600   # 10 minutes in seconds
dist_threshold <- 5     # 5 km

# Make function to convert hhmmss to total seconds from midnight
hms_to_seconds <- function(hhmmss) {
  hms_str <- sprintf("%06d", hhmmss)  # ensure 6 digits
  hh <- as.numeric(substr(hms_str, 1, 2))
  mm <- as.numeric(substr(hms_str, 3, 4))
  ss <- as.numeric(substr(hms_str, 5, 6))
  hh * 3600 + mm * 60 + ss
}

# Remove near-duplicate roosts across stations
roosts_condensed <- roosts_condensed %>%
  mutate(time_sec = hms_to_seconds(time)) %>%
  group_by(date) %>%
  arrange(time_sec) %>% 
  group_split() %>%
  lapply(function(df) {
    remove_idx <- integer(0)
    n <- nrow(df)
    for (i in seq_len(n - 1)) {
      if (i %in% remove_idx) next
      for (j in seq(i + 1, n)) {
        if (j %in% remove_idx) next
        if (df$station[i] != df$station[j]) {
          if (abs(df$time_sec[i] - df$time_sec[j]) <= time_threshold) {
            dist_ij <- calc_dist_km(df$lon[i], df$lat[i], df$lon[j], df$lat[j])
            if (dist_ij <= dist_threshold) {
              if (df$geo_dist[i] > df$geo_dist[j]) {
                cat("Removing roost:", df$track_id[i], "keeping", df$track_id[j], "\n")
                remove_idx <- c(remove_idx, i)
              } else {
                cat("Removing roost:", df$track_id[j], "keeping", df$track_id[i], "\n")
                remove_idx <- c(remove_idx, j)
              }
            }
          }
        }
      }
    }
    if (length(remove_idx) > 0) {
      df <- df[-unique(remove_idx), ]
    }
    return(df)
  }) %>%
  bind_rows() %>%
  ungroup() %>%
  dplyr::select(-time_sec)

##### FINAL CLEANING STEPS ##########
roosts_condensed <- roosts_condensed %>%
  dplyr::select(
    -starts_with("sum_"),
    -starts_with("OL_sum_"),
    -starts_with("n_"), 
    -c("avg_score", "tot_score", "det_score", "notes", "day_notes", "label", "original_label", "viewed", "user_labeled", "scan_sum", "true_scan_sum")
  )

roosts_condensed <- roosts_condensed %>% 
  mutate(year = as.integer(substr(as.character(date), 1, 4)))

roosts_condensed$doy <- yday(as.Date(as.character(roosts_condensed$date), format = "%Y%m%d"))

# Write processed roosts as csv
write.csv(roosts_condensed, "NE_ROOSTS.csv")

###########################
# Move to Roost Ebird Filtering.R


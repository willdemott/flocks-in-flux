# This is the first script used for the Flocks in Flux project. It pulls .csv files provided by UMass researchers Wenlong Zhao, Subhransu Maji, Gustavo Perez, and Daniel Sheldon to turn their machine-learning-system's identified roost rings into roost bird count estimates. Initial .csv files contain one row for every radar scan with machine-learning-detected roost rings, along with track information, geometry, and coordinates. This file processes just one radar-year .csv file, and was used within bash to run iteratively through all radar-year .csv files, producing unique .csv files for each radar-year combination.

# process_one_csv.R
args <- commandArgs(trailingOnly = TRUE)
ui_idx <- args[1]  # this will be the CSV file passed in

# Prepare Workspace
library(bioRad)
library(dplyr)
library(lubridate)
library(raster)
library(tidyverse)
library(sf)
library(ggspatial)

HOME <- "C:/Users/willd/Desktop/Everything/SWALLOW QUANTITY PROJECT/Radar Data"
setwd("/home/radar_path/")
dir.create("./data_pvol")
csv_files <- list.files(path = "/home/radar_path/roosts", pattern = ".csv", full.names = T)

# Silence messages to lower memory load
sink(file = "/dev/null")
options(warn = -1)  # Suppress warnings
options(max.print = 1000)
options(device = function(...) pdf(NULL))


#########################################################################
# Load & filter data 
  
  # Load CSV file for sample year
  year <- read.csv(ui_idx)
  
  # Remove all non-roost tracks
  year <- year %>%
    filter(label != "non-roost")
  
  # Get radar from csv file
  radar <- year$station[1]
  
  # or just change year to season if using unfiltered data
  season <- year
  rm(year)
  # Add empty bird count columns (testing OL feature)
  sum_columns <- c("sum_0.5", "sum_1.5", "sum_2.5", "sum_3.5", "sum_4.5")
  season[sum_columns] <- NA
  
  ###########################################################################
  #Load all raw pvol files to hard-drive using a "for" loop referencing CSV file
  # Create temporary dataframe with reformatted dates from the filename column
  reformatted_dates_df <- season %>%
    mutate(formatted_date = as.POSIXct(strptime(substr(filename, 5, 19),format = "%Y%m%d_%H%M%S"))) %>%
    mutate(date_day = as.Date(formatted_date)) %>%
    group_by(date_day) %>%
    summarise(formatted_min_date = min(formatted_date), formatted_max_date = max(formatted_date)) %>%
    ungroup()
  
  # Iterate over reformatted_dates_df to download all radar data into the directory
  for (i in seq_len(nrow(reformatted_dates_df))) {
    # Get min and max dates for current iteration
    date_min <- reformatted_dates_df$formatted_min_date[i]
    date_max <- reformatted_dates_df$formatted_max_date[i]
    radar <- radar
    directory <- "/home/radar_path/data_pvol"
    
    # Download radar data (polar-volume files) via BioRad package
    suppressMessages(download_pvolfiles(date_min, date_max, radar, directory)
    )

    gc()
  }
  
  # Extract the year from the ui_idx for folder specification
  year_range <- substr(reformatted_dates_df$date_day[1], 1, 4)
  
  # Specify the folder in data_pvol using the extracted year
  pvol_path <- file.path("/home/radar_path/data_pvol", year_range)
  
  # Create a list of downloaded pvolfiles to iterate through
  pvolfiles <- list.files(pvol_path, recursive = TRUE, full.names = TRUE, pattern = radar)
  pvolfiles <- pvolfiles[!grepl("MDM$", basename(pvolfiles))]
  
  # Remove corrupting files (they are always small files, <50 kb)
  pvolfiles <- sapply(pvolfiles, function(file) {
    # Get the file size
    file_size <- file.info(file)$size
    # Return the file if its size is greater than or equal to 50 KB
    if (file_size >= 50 * 1024) {
      return(file)
    } else {
      return(NULL) # Remove the file if it doesn't meet the size condition
    }
  })
  
  # Remove NULL values from the list
  pvolfiles <- pvolfiles[!sapply(pvolfiles, is.null)]
  pvolfiles <- as.character(pvolfiles)
  
##################################
  # Pull updated pvolfiles, count birds!
  for (pvol_idx in seq_along(pvolfiles)) {
    
    # Extract the filename and take the first 19 characters
    current_filename <- substr(basename(pvolfiles[pvol_idx]), 1, 19)
    
    # Filter based on the first 19 characters of the filename
    current_pvol_circles <- season %>%
      filter(substr(filename, 1, 19) == current_filename)
    #message("Number of circles detected: ", nrow(current_pvol_circles))
    
    # Check if there are any circles for the current pvol (biorad downloads pvols through a date min/date max method, and roosts are not always detected between the first and last detection on a given day)
    if (nrow(current_pvol_circles) == 0) {
      next  # Skip this iteration if there are no circles
    }
    
    # Load one of our downloaded files. Occasionally pvol files cannot be read due to missing parameters or warnings.
    pvol <- tryCatch({
      read_pvolfile(pvolfiles[pvol_idx], verbose = FALSE)
    }, error = function(e) {
      message("Error: ", e$message)
      gc()
      return(NULL)
    })
    # If pvol is NULL due to an error, skip to the next iteration
    if (is.null(pvol)) {
      next
    }
    
    # Calculate the parameters used when calculating # of birds
    radar_wavelength <- pvol$attributes$how$wavelength
    
    # Radar cross section of the bird, based on weight (21.2g for TRES, 53.8 for PUMA)
    RCS <- 10^((0.699 * log10(21.2)))
    
    # eta is "true" reflectivity
    pvol <- calculate_param(pvol, eta = dbz_to_eta(DBZH, radar_wavelength, K = sqrt(0.93)))
    
    # N is number of birds per cubic kilometer
    pvol <- calculate_param(pvol, N = eta / RCS)
    
    # Create an object for pvol scans. Sometimes the pvol is missing a scan, if so then we skip that scan
    scan_list <- list()
    for (elev_angle in seq(0.5, 4.5, by = 1)) {
      scan <- tryCatch({
        get_scan(pvol, elev_angle, all = FALSE)
      }, error = function(e) {
        return(NULL)
      })
      if (!is.null(scan)) {
        scan_list[[length(scan_list) + 1]] <- scan 
      } else {
        next
      }
    }
    
    # Prepare spatial-point dataframe list to manipulate data and add layer for bird count
    spatial_scans <- list()
    
    # For each elevation angle: convert to spdf, calculate Vrad, use Vrad to calculate BIRDS
    for (scan_idx in seq_along(scan_list)) {
      
      if (scan_idx <= length(pvol$scans)) {
        
        # Create spatial-point dataframe of the scan
        spatial_scans[[scan_idx]] <- scan_to_spatial(scan_list[[scan_idx]])
        
        # Prepare for Vrad equation, set vertical beam width based on resolution which changes pre/post 2013.
        pi <- 3.141592653589793
        if (pvol$scans[[scan_idx]]$attributes$where$nrays == 720) {
          beamwH <- (pi / 360) * pvol$attributes$how$beamwH
        } else if (pvol$scans[[scan_idx]]$attributes$where$nrays == 360) {
          beamwH <- (pi / 180) * pvol$attributes$how$beamwH
        } else {
          print("nrays does not equal 720 or 360")
        }
        beamwV <- (pi / 180) * pvol$attributes$how$beamwV 
        range_gate_spacing <- pvol$scans[[scan_idx]]$attributes$where$rscale
        range <- spatial_scans[[scan_idx]]@data$range
        
        # Create Vrad
        Vrad <- ((0.35 * sqrt(2 * pi)) / (2 * log(2)) * (((pi * range^2 * beamwH * beamwV * range_gate_spacing) / 4))) / 1e9
        
        # Append Vrad to scan
        spatial_scans[[scan_idx]]@data$Vrad <- Vrad
        
        # Create and append BIRDS for the full scan
        spatial_scans[[scan_idx]]@data$birds <- Vrad * spatial_scans[[scan_idx]]@data$N
        
      } else {
        message("Skipping scan")
        next
      }
    }
    
    # Now that we created a new parameter indicating how many birds are in each point in the radar scans, we can filter out all the data that aren't roosts using masks... We used circular masks, and collected data on the full # of birds in each bounding circle
    #############################
    # CREATING THE MASKS
    #########
    
    ########################
    ## Create overlap coordinates ##
    circle_coords <- list()
    
    # Create sf object to hold all circle polygons for that frame
    polygons <- vector("list", nrow(current_pvol_circles))
    
    # Run for loop to create circles
    for (circle_idx in 1:nrow(current_pvol_circles)) {
      
      # Convert pixel coordinates to meters on Cartesian azimuth grid
      x_cartesian <- (current_pvol_circles$x[circle_idx] - 300) * 500
      y_cartesian <- 500 * (300 - current_pvol_circles$y[circle_idx])
      r_cartesian <- current_pvol_circles$r[circle_idx] * 1.2 * 500
      
      # Define the coordinates for the circular polygon (100 point resolution)
      coords <- matrix(0, ncol = 2, nrow = 101)
      theta <- seq(0, 2 * pi, length.out = 100)
      
      coords[1:100, 1] <- x_cartesian + r_cartesian * cos(theta)
      coords[1:100, 2] <- y_cartesian + r_cartesian * sin(theta)
      coords[101, ] <- coords[1, ]
      
      circle_coords[[circle_idx]] <- list(
        coordinates = coords,
        filename = current_pvol_circles$filename[circle_idx],
        track_id = current_pvol_circles$track_id[circle_idx],
        notes = current_pvol_circles$notes[circle_idx]
      )
      
      # Create a polygon and add it to the polygon list
      polygons[[circle_idx]] <- st_polygon(list(coords))
    }
    
    #################
    # COUNTING AND APPENDING
    #################
    
    # Loop through each circle
    for (coords_idx in seq_along(circle_coords)) {
      # Extract circle metadata
      current_track_id <- circle_coords[[coords_idx]]$track_id
      current_filename <- circle_coords[[coords_idx]]$filename
      
      # Process each scan
      for (scan_idx in seq_along(spatial_scans)) {
        # Initialize and apply masks
        spatial_scans[[scan_idx]]@data$circle_mask <- ifelse(
          sqrt((spatial_scans[[scan_idx]]@coords[, 1] - mean(circle_coords[[coords_idx]]$coordinates[, 1]))^2 + 
                 (spatial_scans[[scan_idx]]@coords[, 2] - mean(circle_coords[[coords_idx]]$coordinates[, 2]))^2) <= 
            (current_pvol_circles$r[coords_idx] * 1.2 * 500), 1, 0
        )
        
        # Calculate masked birds and sum for the scan
        masked_birds <- spatial_scans[[scan_idx]]@data$circle_mask * 
          spatial_scans[[scan_idx]]@data$birds
        bird_sum <- sum(masked_birds, na.rm = TRUE)  
        # Append the result to the correct column in the correct row in season
        row_index <- which(season$track_id == current_track_id & season$filename == current_filename)
        season[row_index, sum_columns[scan_idx]] <- bird_sum
      }
    }
    #message("Memory usage before clearing: ", gc()[, "used"][2], " MB")
    rm(list = setdiff(ls(), c("HOME", "directory", "csv_files", "directory", "reformatted_dates_df", "year", "radar", "season", "sum_columns", "pvolfiles", "pvol_path", "year_range")))
    gc()
  }
  # Initialize a column for the beam heights
  season$height <- NA  # Create a new column for beam heights
  # Loop through each roost in the season dataframe
  for (i in 1:nrow(season)) {
    h0 <- pvol$attributes$where$height
    distance <- season$geo_dist[i]
    elangle <- pvol$scans[[1]]$geo$elangle
    # Use the beam_height function to calculate the height for the given distance
    season$height[i] <- beam_height(distance, elangle) + h0
  }
  
  unlink(pvolfiles, force = TRUE)
  # Write CSV for season
  filename <- paste0(radar, "_", year_range, "_roosts.csv")
  write.csv(season, file = filename, row.names = FALSE)
  # Clean up the environment, keeping only the necessary variables
  rm(list = setdiff(ls(), c("HOME", "csv_files")))
  gc()
  sink()

# Next script used was "NE Roost Cleaning.R"
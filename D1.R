# ==============================================================================
# Day 1: Create your first SDM in R
# Target Species: Grus japonensis (Red-crowned crane)
# ==============================================================================

# Install necessary packages if not already installed
# install.packages(c("geodata", "terra", "sf", "dismo", "data.table"))

library(geodata)      # For downloading species and climate data
library(terra)        # Core package for modern raster data processing
library(sf)           # Core package for modern vector data processing
library(dismo)        # Classic SDM package
library(data.table)   # High-performance data manipulation
library(here)

project_root <- here()

data_dir <- file.path(project_root, "Data")
occ_dir <- file.path(data_dir, "Occurrences")
bioclim_dir <- file.path(data_dir, "Bioclim")

#setwd("~/GIT/ENM_curriculum/ENM_curriculum")
setwd(project_root)

# ------------------------------------------------------------------------------
# Step 1: Downloading the occurrences
# ------------------------------------------------------------------------------

# Fetch GBIF data using geodata package, saving to temporary directory
occ_file <- file.path(occ_dir, "crane_gbif.rda")

if (file.exists(occ_file)){
  cat("Reading GBIF data for Grus japonensis...\n")
  crane_dt <- readRDS(occ_file)
} else{
  # 创建目录
  dir.create(occ_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("Downloading GBIF data for Grus japonensis...\n")
  crane_gbif <- tryCatch({
    sp_occurrence(genus = "Grus", species = "japonensis", path = tempdir())
  }, error = function(e) {
    stop("Failed to download GBIF data: ", e$message, "\n")
  })
  
  # Convert to data.table for efficient processing
  crane_dt <- as.data.table(crane_gbif)
  saveRDS(crane_dt, occ_file)
}

# ------------------------------------------------------------------------------
# Step 2: Getting environmental predictors
# ------------------------------------------------------------------------------

# Download global Bioclim data (10 arc-minutes resolution, 19 variables)
bioclim_10m_dir <- file.path(bioclim_dir, "wc2.1_10m_bio")

if (dir.exists(bioclim_10m_dir)){
  cat("Read Bioclim data (10 arc-minutes resolution data\n")
  files <- list.files(bioclim_10m_dir, pattern = "*.tif", full.names = T)
  clim_global <- rast(files)
} else{
    # 创建目录
    dir.create(bioclim_10m_dir, recursive = TRUE, showWarnings = FALSE)

    cat("Downloading WorldClim Bioclimatic variables...\n")
    clim_global <- tryCatch({
      worldclim_global(var = "bio", res = 10, path = tempdir())
    }, error = function(e) {
      stop("Failed to download BioClim data: ", e$message, "\n")
    })
    
    cat("Saving BioClim data as .tif files...\n")
    for (i in 1:nlyr(clim_global)) {
      writeRaster(clim_global[[i]], 
                  filename = file.path(bioclim_10m_dir, paste0("bio", i, ".tif")),
                  overwrite = TRUE)
    }
    cat("Saved", nlyr(clim_global), "layers to", bioclim_10m_dir, "\n")
    cat("\n")
  }
  #https://www.worldclim.org/data/bioclim.html


# 创建空间范围对象
# 坐标范围：110°E - 150°E，30°N - 55°N
study_extent <- ext(110, 150, 30, 55) 
clim_study_area <- crop(clim_global, study_extent)

# ------------------------------------------------------------------------------
# Step 3: Advanced Spatial Cleaning (NA Removal & Spatial Thinning)
# ------------------------------------------------------------------------------
cat("Performing spatial cleaning based on raster resolution...\n")
occurrences_raw <- crane_dt[!is.na(lon) & !is.na(lat), .(lon, lat)]

# 3.1 Extract the raster cell ID for each occurrence coordinate
# terra::cellFromXY returns the exact cell number that a coordinate falls into
pts_matrix <- as.matrix(occurrences_raw[, .(lon, lat)])
occurrences_raw[, cell_id := cellFromXY(clim_study_area, pts_matrix)]

# 3.2 Extract environmental values to check for NAs (e.g., points in the ocean)
# We use the first layer (Bio1) as the mask. The returned object is a matrix/data.frame,
# where the second column contains the actual extracted values.
env_vals <- extract(clim_study_area[[1]], pts_matrix)
occurrences_raw[, env_val := env_vals[, 1]]

# 3.3 Remove points that fall in NA areas (ocean or outside raster bounds)
occ_no_na <- occurrences_raw[!is.na(env_val)]
cat("Occurrences after removing NAs:", nrow(occ_no_na), "\n")

# 3.4 Spatial Thinning: Keep only ONE point per raster cell
# data.table's 'unique' function with 'by' argument does this instantly
occurrences_clean <- unique(occ_no_na, by = "cell_id")
cat("Occurrences after spatial thinning (1 per cell):", nrow(occurrences_clean), "\n")

# ------------------------------------------------------------------------------
# Step 3.5: Visualizing the Cleaning Process
# ------------------------------------------------------------------------------
cat("Plotting before-and-after spatial cleaning...\n")

# Set up a side-by-side plotting area
par(mfrow = c(1, 2))

# Plot 1: Raw Data (Includes overlapping points and points in the ocean)
occurrences_raw[, cell_count := .N, by = cell_id]                      # calculate how many occ in each cell
occurrences_density <- unique(occurrences_raw, by = "cell_id")
occurrences_density[, point_size := 0.5 + 0.1443 * log1p(cell_count)]  # smooth point size with log function 
                                                                       # make coefficient = 0.1443, as when cell_count = 1, point_size = 0.5 + 0.1443*log1p(1) = 0.6, same as the default point size 
plot(clim_study_area[[1]], main = "Raw Occurrences", legend = FALSE)

points(
  occurrences_density$lon, 
  occurrences_density$lat, 
  col = adjustcolor("black"),
  pch = 16, 
  cex = occurrences_density$point_size
)

# Plot 2: Cleaned Data (No ocean NAs, thinned to 1 point per pixel)
plot(clim_study_area[[1]], main = "Cleaned & Thinned", legend = FALSE)
points(occurrences_clean$lon, occurrences_clean$lat, col = "red", pch = 16, cex = 0.6) # default point size = 0.6

# Reset plot parameters
par(mfrow = c(1, 1))

# ------------------------------------------------------------------------------
# Step 4: Running your first ENM (Bioclim approach)
# ------------------------------------------------------------------------------
cat("Fitting the Bioclim model...\n")

# Fit the Bioclim model using only our strictly cleaned occurrence data.frame
# dismo's bioclim accepts a Raster object and a data.frame of coordinates
bc_model <- bioclim(as(clim_study_area, "Raster"), as.data.frame(occurrences_clean[, .(lon, lat)]))

# ------------------------------------------------------------------------------
# Step 5: Visualizing the results
# ------------------------------------------------------------------------------
cat("Predicting spatial distribution...\n")

# Predict and convert back to terra's SpatRaster
suitability_map <- predict(as(clim_study_area, "Raster"), bc_model)
suitability_map <- rast(suitability_map)

# Plot the final suitability map
plot(suitability_map, main = "Habitat Suitability for Grus japonensis (Bioclim)", 
     col = terrain.colors(100, rev = TRUE))
points(occurrences_clean$lon, occurrences_clean$lat, col = adjustcolor("blue", alpha.f = 0.2), pch = 16, cex = 0.5)

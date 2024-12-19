library(terra)
library(dplyr)
library(ggplot2)

rm(list = ls())
select = dplyr::select

# df <- readRDS("data/eds_grid_100m.rds") %>% filter(unit == "Tutuila")
df <- readRDS("data/eds_grid_500m.rds") %>% filter(unit == "Tutuila")
# df <- readRDS("data/eds_grid_1km.rds") %>% filter(unit == "Tutuila")

names(df) <- gsub("Daily", "daily", names(df)); names(df)
names(df) <- gsub("Weekly", "weekly", names(df)); names(df)
names(df) <- gsub("_8Day_", "_8days_", names(df)); names(df)
names(df) <- gsub("Monthly", "monthly", names(df)); names(df)

names(df) <- gsub("DHW.", "dhw.", names(df)); names(df)
names(df) <- gsub("Major_", "major_", names(df)); names(df)
names(df) <- gsub("Np10y_", "np10y_", names(df)); names(df)
names(df) <- gsub("MeanMax_", "mean_max_", names(df)); names(df)
names(df) <- gsub("CI95Max_", "ci95_max_", names(df)); names(df)
names(df) <- gsub("MeanDur_", "mean_durnal_", names(df)); names(df)
names(df) <- gsub("MaxMax_", "max_max_", names(df)); names(df)

names(df) <- gsub("_DY01", "_01dy", names(df)); names(df)
names(df) <- gsub("_WK01", "_01wk", names(df)); names(df)
names(df) <- gsub("_MO01", "_01mo", names(df)); names(df)
names(df) <- gsub("_MO03", "_03mo", names(df)); names(df)
names(df) <- gsub("_MO06", "_06mo", names(df)); names(df)
names(df) <- gsub("_YR10YR01", "_10yr_01yr", names(df)); names(df)
names(df) <- gsub("_YR01", "_01yr", names(df)); names(df)
names(df) <- gsub("_YR03", "_03yr", names(df)); names(df)
names(df) <- gsub("_YR05", "_05yr", names(df)); names(df)
names(df) <- gsub("_YR10", "_10yr", names(df)); names(df)
names(df) <- gsub("_ALLB4", "_all_before", names(df)); names(df)

# df <- df %>% select(-contains("kd490")); names(df)
# df <- df %>% select(-contains("chlorophyll_a_npp_viirs")); names(df)
df <- df %>% select(-contains("hi_otp")); names(df)
df <- df %>% select(-contains("mhi")); names(df)
df <- df %>% select(-contains("biweekly_range")); names(df)
df <- df %>% select(-contains("monthly_range")); names(df)
df <- df %>% select(-contains("_dy01")); names(df)
df <- df %>% select(-contains("_wk01")); names(df)
df <- df %>% select(-contains("_mo01")); names(df)
df <- df %>% select(-contains("_mo03")); names(df)
df <- df %>% select(-contains("_mo06")); names(df)
df <- df %>% select(-contains("_yr03")); names(df)
df <- df %>% select(-contains("_yr05")); names(df)
df <- df %>% select(-contains("_yr10")); names(df)
df <- df %>% select(-contains("_allb4")); names(df)
df <- df %>% select(-contains("_yr10yr01")); names(df)
df <- df %>% select(-contains("_15_min")); names(df)
df <- df %>% select(-contains("_30_min")); names(df)
df <- df %>% select(-contains("1_deg.")); names(df)
df <- df %>% select(-contains("7daymax")); names(df)
df <- df %>% select(-contains("jplmur")); names(df)

names(df) <- gsub(".tif", "", names(df)); names(df)
names(df) <- gsub(".nc", "", names(df)); names(df)

names(df) <- gsub("diste.from.port.v20201104f", "distance.from.port", names(df)); names(df)
names(df) <- gsub("gpw_v4_population_density_rev11_2pt5_min", "population_density", names(df)); names(df)
names(df) <- gsub("sed_export", "sedimentation", names(df)); names(df)

names(df) <- gsub("_daily", "", names(df)); names(df)
names(df) <- gsub("_monthly", "", names(df)); names(df)

names(df) <- gsub("_crw_", "_", names(df)); names(df)
# names(df) <- gsub("_viirs_", "_", names(df)); names(df)
names(df) <- gsub("_noaa_", "_", names(df)); names(df)
names(df) <- gsub("_yr01", "", names(df)); names(df)

# df <- df %>% mutate(bathymetry = ifelse(bathymetry <= -30, NA, bathymetry))
df <- df %>% filter(!is.na(sedimentation))
df <- df %>% filter(!is.na(bathymetry))
df <- df %>% filter(!is.na(population_density))

# merge data sources for chla, kd490, par, and kd

# Define the statistical summaries
summaries <- c("mean", "sd", "q05", "q95", "mean_annual_range")

# Define variables and their corresponding data sources
variable_sources <- list(
  "chlorophyll_a" = c("esa_oc_cci_v6.0", "npp_viirs"),
  "kd490" = c("esa_oc_cci", "viirs"),
  "par" = c("aqua_modis", "nasa_viirs")
  # "kdpar" = c("viirs")
)

# Function to merge multiple columns
merge_columns <- function(...) {
  cols <- list(...)
  non_na_counts <- rowSums(!is.na(do.call(cbind, cols)))
  merged_col <- rowMeans(do.call(cbind, cols), na.rm = TRUE)
  merged_col[non_na_counts == 0] <- NA
  return(merged_col)
}

# Loop over variables and summaries to create merged columns
for (variable in names(variable_sources)) {
  
  sources <- variable_sources[[variable]]
  
  for (stat in summaries) {
    
    # Initialize a list to hold columns from different sources
    cols_list <- list()
    
    for (source in sources) {
      
      # Construct the column name
      col_name <- paste0(stat, "_", variable, "_", source)
      if (col_name %in% names(df)) {
        # Add the column to the list
        cols_list[[source]] <- df[[col_name]]
      } else {
        # If the column doesn't exist, fill with NA
        cols_list[[source]] <- rep(NA, nrow(df))
      }
    }
    
    # Merge the columns using the merging function
    merged_col <- merge_columns(cols_list[[1]], cols_list[[2]])
    
    # Name of the merged column
    merged_col_name <- paste0(stat, "_", variable)
    
    # Add the merged column to the data frame
    df[[merged_col_name]] <- merged_col
    
  }
}

# Optional: Remove the original columns
original_cols <- c()
for (variable in names(variable_sources)) {
  sources <- variable_sources[[variable]]
  for (stat in summaries) {
    for (source in sources) {
      col_name <- paste0(stat, "_", variable, "_", source)
      original_cols <- c(original_cols, col_name)
    }
  }
}

df <- df %>% select(-all_of(original_cols))

# List of columns with more than 10% NAs
columns_with_na <- df %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  select(where(~ . > 0.01)) %>%
  names()

# Loop through each column and plot
cols_with_na_in_summary <- c()

# Loop through each column and plot
for (col in columns_with_na) {
  
  df_summary <- df %>%
    group_by(lon, lat) %>%
    summarise(y = mean(.data[[col]], na.rm = TRUE))
  
  if (sum(is.na(df_summary$y)) > 0) {
    print(col)
    # Append the column name to the vector
    cols_with_na_in_summary <- c(cols_with_na_in_summary, col)
  }
  
  p <- ggplot(df_summary, aes(lon, lat, fill = y)) +
    geom_raster(show.legend = FALSE) +
    scale_fill_viridis_c() + 
    labs(title = col, fill = col)
  
  print(p)
}

df <- df %>% select(-all_of(cols_with_na_in_summary))

# df <- df %>% select(where(~ mean(is.na(.)) <= 0.1)); names(df)
# df <- df %>% select(where(~ all(!is.na(.)))); names(df)
df <- df %>% select(where(~ n_distinct(.) > 1)); names(df)

# Define the columns to keep in a specific order
key_columns <- c("lon", "lat", "year", "month", "day", "bathymetry", 
                 "distance.from.port", "population_density", "sedimentation")

# Get the rest of the columns, sorted alphabetically
other_columns <- setdiff(names(df), key_columns) %>% sort()

# Reorder the columns by combining the key columns with the alphabetically sorted ones
df <- df %>% select(all_of(key_columns), all_of(other_columns))

# View the reordered column names
names(df)

visdat::vis_miss(df, warn_large_data = F)

# rastertize for maxent 

eds <- rast()

for (v in 6:ncol(df)) {
  
  var_name <- colnames(df)[v]
  
  p <- df %>%
    dplyr::select(lon, lat, !!sym(var_name)) %>% 
    group_by(lon, lat) %>% 
    summarise(mean_value = mean(!!sym(var_name), na.rm = TRUE)) %>% 
    ungroup()
  
  colnames(p)[3] = var_name
  
  p = p %>% 
    rast(type = "xyz")
  
  eds <- c(eds, p)
  
}

plot(eds)

# Subset environmental data based on VIF results
load("output/vif.RData")
eds <- raster::subset(eds, v@results$Variables)
names(eds)

eds = raster::stack(eds)
eds = raster::readAll(eds)
save(eds, file = "data/eds.rdata")

#interpolate to 5m
load("data/tutuila_hybrid_5m_bathymetry.rdata")
load("data/tutuila_hybrid_10m_bathymetry.rdata")

# Convert RasterStack to SpatRaster
# eds_terra <- rast(eds)
eds_terra <- eds
df_terra <- rast(df)

rm(eds, df)

# If eds_terra does not have a CRS, assign it (assuming same as df_terra):
crs(eds_terra) <- crs(df_terra)

# Initialize an empty SpatRaster to store the processed layers
eds_clipped_stack <- rast()

# Loop through each layer in eds_terra
for (i in 2:nlyr(eds_terra)) {
  
  # i = 10
  
  # Extract single layer
  eds_layer <- eds_terra[[i]]
  # plot(eds_layer)
  
  # If needed, reproject (and resample) the layer to match df_terra
  # If CRS is the same, this step will just resample to match resolution and extent
  eds_layer_projected <- project(eds_layer, df_terra, method = "bilinear")
  
  # Mask the projected layer using df_terra, retaining only where df_terra is not NA
  eds_layer_clipped <- mask(eds_layer_projected, df_terra)
  # plot(eds_layer_clipped)
  
  # Add the processed layer to our output stack
  eds_clipped_stack <- c(eds_clipped_stack, eds_layer_clipped)
  print(i)
  gc()
  
}

# eds_clipped_stack now contains all layers from eds, 
# reprojected/resampled and clipped by df’s NA pattern.
plot(eds_clipped_stack)

names(df_terra) = "bathymetry"

eds = c(df_terra, eds_clipped_stack)

# Run VIF step for variable selection
# note, this VIF selection was done at EDS outputs resolution : 0.005, 0.005
# v <- usdm::vifstep(terra::rast(eds), th = 3)
v <- usdm::vifstep(eds, th = 3)
save(v, file = "output/vif.RData")
load("output/vif.RData")

# Subset environmental data based on VIF results
eds <- raster::subset(eds, v@results$Variables)
names(eds)

# Plot the environmental data
plot(eds, col = colorRamps::matlab.like(100))

eds = raster::stack(eds)
eds = raster::readAll(eds)

save(eds, file = "data/eds.rdata")

# for Kira's environmental layer maps 0-18m
eds[["bathymetry"]][eds[["bathymetry"]] <= -18] <- NA
eds <- mask(eds, eds[["bathymetry"]])
plot(eds)
eds = raster::readAll(eds)
save(eds, file = "/Users/Kisei.Tanaka/Desktop/eds.rdata")

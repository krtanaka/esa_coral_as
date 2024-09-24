# Load required libraries
library(raster)      # For working with raster data
library(dplyr)       # For data manipulation
library(readr)       # For reading CSV files
library(colorRamps)  # For color ramp functions
library(ggrepel)     # For repelling overlapping text labels in ggplot2
library(sf)
library(sp)

# Clear the workspace
rm(list = ls())

# Set the target spatial resolution in meters
spatial_resolution = 1000

dat = st_read("N:/GIS/Projects/LBSP/Sediment/Output/AS2a10_PourPt_Output.shp")
# Load necessary libraries
library(terra)
library(sf)

# Assuming 'dat' is your sf object
# If it's not already loaded, you can read it using:
# dat <- st_read("path_to_your_shapefile.shp")

# Reproject the data to WGS84 latitude and longitude
dat_latlon <- st_transform(dat, crs = 4326)  # EPSG code 4326 corresponds to WGS84 lat/lon

# Convert the sf object to a terra SpatVector
dat_spatvect <- vect(dat_latlon)

# Set the desired raster resolution in degrees
# For example, 0.0001 degrees (~11 meters at the equator)
resolution <- 0.01  # Adjust as needed

# Create a raster template with the extent and CRS of your data
raster_template <- rast(
  extent = ext(dat_spatvect),
  resolution = resolution,
  crs = crs(dat_spatvect)
)

# Rasterize the 'sed_export' field onto the raster template
sed_export_raster <- rasterize(
  x = dat_spatvect,
  y = raster_template,
  field = "sed_export",
  fun = mean,  # Use mean if multiple points fall into the same cell
  background = NA  # Set background value for cells with no data
)

# Optionally, plot the raster to visualize
plot(sed_export_raster, main = "Rasterized sed_export in Lat/Lon")

# Save the raster to a file if needed
writeRaster(sed_export_raster, "data/sed_export.tif", overwrite = TRUE)

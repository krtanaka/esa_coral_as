library(dplyr)
library(terra)
library(readr)

rm(list = ls())

file_paths <- c(
  "N:/GIS/Projects/CommonMaps/Bathymetry/tut_dball.asc",
  "N:/GIS/Projects/CommonMaps/Bathymetry/tut_dbmb.asc",
  "N:/GIS/Projects/CommonMaps/Bathymetry/Tutuila_5m.asc",
  "N:/GIS/Projects/CommonMaps/Bathymetry/cudem_ninth_amsam/cudem_ninth_amsam_J1078896.tif",
  "N:/GIS/Projects/CommonMaps/Bathymetry/2022_ngs_topobathy_dem_tutuila/2022_ngs_topobathy_dem_tutuila_J1078893_000_000.tif",
  "N:/GIS/Projects/CommonMaps/Bathymetry/2022_ngs_topobathy_dem_tutuila/2022_ngs_topobathy_dem_tutuila_J1078893_001_000.tif"
)

rasters <- lapply(file_paths, rast)
# lapply(rasters, function(r) { plot(r); res(r) })

base <- rasters[[3]]

rasters_resampled <- lapply(rasters, function(r) resample(r, base, method = "near"))

fine_topo <- mean(do.call(c, rasters_resampled), na.rm = TRUE)
fine_topo[fine_topo >= 0] <- NA
fine_topo[fine_topo <= -100] <- NA
plot(fine_topo)

# Assign the UTM Zone 2S CRS to fine_topo
crs(fine_topo) <- "EPSG:32702"

# Reproject to WGS84 (lat/lon)
fine_topo_latlon <- project(fine_topo, "EPSG:4326")

plot(fine_topo_latlon)

output_path <- "G:/Environmental_Data_Summary/ESA/Bathymetry_cudem_ngs_pacioos_merged/Bathymetry_Cudem_NGS_PacIoos_merged_5m_all_units.nc"
writeCDF(fine_topo_latlon, output_path, overwrite = TRUE)

# Load required libraries
library(tidyverse)
library(rgbif)
library(sf)
library(raster)
library(geodata)
library(ENMeval)
library(ecospat)
library(rJava)
library(dismo)
library(countrycode)
library(readr)
library(colorRamps)
library(ggmap)
library(ggspatial)
library(doParallel)

# Clear workspace
rm(list = ls())

select = dplyr::select

# Register Google Maps API key
ggmap::register_google("AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

# Load custom functions (if any)
source("script/functions.R")

# Define species list and select the species
species <- c("Acropora globiceps", "Isopora crateriformis")[2]

# Use paste() to construct file paths dynamically based on the selected species
ncrmp <- read_csv(paste0("data/occurances_", species, "_ncrmp.csv")) %>%
  select(Longitude, Latitude, Scientific.Name, Source)

gbif <-  read_csv(paste0("data/occurances_", species, "_gbif_obis.csv")) %>%
  select(Longitude, Latitude, Scientific.Name, Source)

nps <-  read_csv(paste0("data/occurances_", species, "_nps.csv")) %>%
  select(Longitude, Latitude, Scientific.Name, Source)

cb <-  read_csv(paste0("data/occurances_", species, "_coralbelt.csv")) %>%
  select(Longitude, Latitude, Scientific.Name, Source)

# Combine datasets
occ_df <- bind_rows(ncrmp, gbif, nps, cb)
# occ_df <- bind_rows(ncrmp, cb)

# map = ggmap::get_map(location = c(-170.7231, -14.30677),
#                      maptype = "satellite",
#                      zoom = 11,
#                      # color = "bw",
#                      force = T)
# 
# ggmap(map) +
#   geom_spatial_point(data = occ_df, aes(Longitude, Latitude, fill = Source, color = Source),
#                      size = 3,
#                      shape = 21, alpha = 0.8, crs = 4326) +
#   annotate("text", x = -170.8375, y = -14.23714,
#            label = species,
#            hjust = 0, vjust = 1, size = 5, color = "white", fontface = "bold") +
#   scale_y_continuous(limits = c(-14.36541, -14.23714), "") +
#   scale_x_continuous(limits = c(-170.8375, -170.5498), "") +
#   theme_minimal() +
#   theme(legend.position = c(0.92, 0.22),
#         legend.background = element_blank(),             # Transparent background
#         legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
#         legend.text = element_text(color = "white", face = "bold"),  # White and bold text
#         legend.title = element_text(color = "white", face = "bold"))
#
# ggsave(last_plot(), filename =  file.path(paste0("output/combined_occurances_", species, ".png")), width = 9)

# Check how many occurrences subset for each spp.
table(occ_df$Scientific.Name)

occ_df = occ_df %>% 
  filter(Longitude >= min(ncrmp$Longitude), Longitude <= max(ncrmp$Longitude),
         Latitude >= min(ncrmp$Latitude), Latitude <= max(ncrmp$Latitude)) %>% 
  select(Longitude, Latitude, Scientific.Name) %>% 
  distinct()

load("data/eds.rdata")

v = usdm::vifstep(terra::rast(eds), th = 5)
eds = raster::subset(eds, v@results$Variables); names(eds)
# eds <- dropLayer(eds, "sed_export"); names(eds)

plot(eds, col = matlab.like(100))

maxent_results = run_maxent(occ_df, eds)
beepr::beep(2)

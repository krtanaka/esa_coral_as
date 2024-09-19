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
species_list <- c(
  "Acropora globiceps",
  "Isopora crateriformis"
)[2]

# Select species-specific data and process accordingly
if (species_list == "Acropora globiceps") {
  
  # Load NCRMP data
  ncrmp <- read_csv("data/A_globiceps_AS.csv") %>%
    filter(ISLAND == "Tutuila", AdColDen > 0) %>%
    dplyr::select(LONGITUDE, LATITUDE) %>%
    mutate(Scientific.Name = species_list, Source = "NCRMP") %>%
    rename(Longitude = LONGITUDE, Latitude = LATITUDE)
  
  # Load GBIF data
  gbif <- read_csv("data/gbif_occurances_Acropora_globiceps.csv") %>%
    mutate(Scientific.Name = species_list, Source = "GBIF") %>%
    select(Longitude, Latitude, Scientific.Name, Source)
  
  # Load NPS data
  nps <- read_csv("data/nps_occurances_Acropora globiceps.csv") %>%
    mutate(Scientific.Name = species_list, Source = "NPS") %>%
    select(Longitude, Latitude, Scientific.Name, Source)
  
  # Combine both datasets
  occ_df <- bind_rows(ncrmp, gbif, nps)
  # occ_df <- bind_rows(ncrmp)
  
} else {
  
  # Load NCRMP data for Isopora crateriformis
  ncrmp <- read_csv("data/I_craterformis_AS.csv") %>%
    filter(ISLAND == "Tutuila", AdColDen > 0) %>%
    dplyr::select(LONGITUDE, LATITUDE) %>%
    mutate(Scientific.Name = species_list, Source = "NCRMP") %>%
    rename(Longitude = LONGITUDE, Latitude = LATITUDE)
  
  # Load GBIF data for Isopora crateriformis
  gbif <- read_csv("data/gbif_occurances_Isopora_crateriformis.csv") %>%
    mutate(Scientific.Name = species_list, Source = "GBIF") %>%
    select(Longitude, Latitude, Scientific.Name, Source)
  
  # Load NPS data
  nps <- read_csv("data/nps_occurances_Isopora crateriformis.csv") %>%
    mutate(Scientific.Name = species_list, Source = "NPS") %>%
    select(Longitude, Latitude, Scientific.Name, Source)
  
  # Combine both datasets
  occ_df <- bind_rows(ncrmp, gbif, nps)
  # occ_df <- bind_rows(ncrmp)
  
}

map = ggmap::get_map(location = c(mean(ncrmp$Longitude), mean(ncrmp$Latitude)),
                     maptype = "satellite",
                     zoom = 11,
                     # color = "bw",
                     force = T)

ggmap(map, darken = c(0.5, "black")) +
  geom_spatial_point(data = occ_df, aes(Longitude, Latitude, fill = Source, color = Source),
                     size = 3,
                     shape = 21, alpha = 0.8, crs = 4326) +
  annotate("text", x = -170.8375, y =  -14.23039, 
           label = species_list,
           hjust = 0, vjust = 1, size = 6, color = "white", fontface = "bold") +
  scale_y_continuous(limits = c(-14.37032, -14.23039), "") +
  scale_x_continuous(limits = c(-170.8375, -170.508), "") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.2),
        legend.background = element_blank(),             # Transparent background
        legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
        legend.text = element_text(color = "white", face = "bold"),  # White and bold text
        legend.title = element_text(color = "white", face = "bold"))

ggsave(last_plot(), filename =  file.path(paste0("output/ncrmp_gbif_nps", species_list, ".png")), height = 4, width = 8.5)

# Check how many occurrences subset for each spp.
table(occ_df$Scientific.Name)

load("data/eds.rdata")

v = usdm::vifstep(terra::rast(eds), th = 5)
eds = raster::subset(eds, v@results$Variables); names(eds)

plot(eds, col = matlab.like(100))

maxent_results = run_maxent(occ_df, eds)
beepr::beep(2)

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

# Load necessary dplyr function
select <- dplyr::select

# Register Google Maps API key
ggmap::register_google(key = "AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

# Get map for the given coordinates
map = ggmap::get_map(location = c(-170.6927, -14.29396),
                     maptype = "satellite",
                     zoom = 11,
                     color = "bw",
                     force = TRUE)

# Load custom functions (if any)
source("script/functions.R")

# Load environmental dataset
load("data/eds.rdata")
  
# Define species list and select the species
species <- c("Acropora globiceps", "Isopora crateriformis")[1]
survey <- c("ncrmp", "combined", "no_nps")[2]

# Load NCRMP occurrences
ncrmp <- read_csv(paste0("data/occurances_", species, "_ncrmp_exp.csv"))

# Set file paths based on the survey type
file_paths <- switch(survey,
                     
                     ncrmp = list(ncrmp = paste0("data/occurances_", species, "_ncrmp_exp.csv")),
                     
                     combined = list(
                       ncrmp = paste0("data/occurances_", species, "_ncrmp_exp.csv"),
                       gbif = paste0("data/occurances_", species, "_gbif_obis.csv"),
                       nps = paste0("data/occurances_", species, "_nps.csv"),
                       crag = paste0("data/occurances_", species, "_crag.csv")),
                     
                     no_nps = list(
                       ncrmp = paste0("data/occurances_", species, "_ncrmp_exp.csv"),
                       gbif = paste0("data/occurances_", species, "_gbif_obis.csv"),
                       crag = paste0("data/occurances_", species, "_crag.csv"))
)

# Initialize an empty list to store the data frames
data_list <- list()

# Read the files if they exist
for (dataset in names(file_paths)) {
  if (file.exists(file_paths[[dataset]])) {
    data <- read_csv(file_paths[[dataset]]) %>%
      select(Longitude, Latitude, Scientific.Name, Source)
    data_list[[dataset]] <- data
  }
}

# Combine all available datasets
occ_df <- bind_rows(data_list) %>%
  filter(Latitude >= -14.38, Latitude <= -14.22,
         Longitude >= -170.85, Longitude <= -170.53)

# plot(occ_df$Longitude, occ_df$Latitude)
# maps::map(add = T)

# Check the number of occurrences for each species
sum_p = nrow(occ_df)

ggmap(map) +
  geom_spatial_point(data = occ_df, aes(Longitude, Latitude, 
                                        fill = Source   ,  
                                        color = Source),
                     alpha = 0.8,  # Map size to y
                     shape = 21, size = 5, crs = 4326, show.legend = T) +
  coord_sf(crs = 4326) +    # Use coord_sf to address the warning
  scale_y_continuous(limits = c(-14.38, -14.22), "") +
  scale_x_continuous(limits = c(-170.85, -170.53), "") +
  scale_fill_discrete("") + 
  scale_color_discrete("") + 
  annotate("text", x = -170.85, y = -14.22, 
           label = paste0(species, "\nn = ", sum_p), 
           hjust = 0, vjust = 1, size = 6, color = "white", fontface = "bold") + 
  theme(
    legend.position = c(0.8, 0.22),
    # legend.position = c(0.1, 0.75),
    legend.background = element_blank(),             # Transparent background
    legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
    legend.text = element_text(color = "white", face = "bold", size = 18),  # White and bold text
    legend.title = element_text(color = "white", face = "bold")
  )

ggsave(last_plot(), file = paste0("data/occurances_", species, "_", survey, ".png"),  width = 8, height = 4.2)

# Filter occurrences based on the NCRMP data
occ_df <- occ_df %>%
  filter(Longitude >= min(ncrmp$Longitude), Longitude <= max(ncrmp$Longitude),
         Latitude >= min(ncrmp$Latitude), Latitude <= max(ncrmp$Latitude)) %>%
  select(Longitude, Latitude, Scientific.Name) %>%
  distinct()

# Run Maxent model
maxent_results <- run_maxent(occ_df, eds, survey)
beepr::beep(2)

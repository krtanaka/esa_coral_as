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
ggmap::register_google("AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

# Load custom functions (if any)
source("script/functions.R")

# Define species list and select the species
species <- c("Acropora globiceps", "Isopora crateriformis", "Genus Tridacna")[1]
survey = c("ncrmp", "combined")[1]

# Load NCRMP occurrences
ncrmp <- read_csv(paste0("data/occurances_", species, "_ncrmp_exp.csv"))

# Set file paths based on the model
file_paths <- if (survey == "ncrmp") {
  
  list(ncrmp = paste0("data/occurances_", species, "_ncrmp_exp.csv"))
  
} else {
  
  list(
    ncrmp = paste0("data/occurances_", species, "_ncrmp_exp.csv"),
    gbif = paste0("data/occurances_", species, "_gbif_obis.csv"),
    nps = paste0("data/occurances_", species, "_nps.csv"),
    crag = paste0("data/occurances_", species, "_crag.csv")
  )
}

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

# Generate the map
map <- ggmap::get_map(location = c(-170.7231, -14.30677),
                      maptype = "satellite",
                      zoom = 11,
                      force = TRUE)
n = dim(occ_df)[1]

ggmap(map) +
  geom_spatial_point(data = occ_df, aes(Longitude, Latitude, fill = Source, color = Source),
                     size = 3, shape = 21, alpha = 0.8, crs = 4326) +
  annotate("text", x = -170.85, y = -14.22,
           label = paste0(species, "\nn = ", n), hjust = 0, vjust = 1, size = 6, color = "white", fontface = "bold") +
  scale_fill_discrete("") + 
  scale_color_discrete("") + 
  scale_y_continuous(limits = c(-14.38, -14.22), "") +
  scale_x_continuous(limits = c(-170.85, -170.53), "") +
  ggdark::dark_mode() + 
  theme(legend.position = c(0.85, 0.22),
        legend.background = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        legend.text = element_text(color = "white", size = 12, face = "bold"), 
        legend.title = element_text(color = "white", size = 12, face = "bold"))

# # Save the plot
# ggsave(last_plot(), filename = file.path(paste0("output/occurances_", species, "_", survey, ".png")), width = 9)

# Check the number of occurrences for each species
table(occ_df$Scientific.Name)

# Filter occurrences based on the NCRMP data
occ_df <- occ_df %>%
  filter(Longitude >= min(ncrmp$Longitude), Longitude <= max(ncrmp$Longitude),
         Latitude >= min(ncrmp$Latitude), Latitude <= max(ncrmp$Latitude)) %>%
  select(Longitude, Latitude, Scientific.Name) %>%
  distinct()

# Load environmental dataset
load("data/eds.rdata")

# Run VIF step for variable selection
v <- usdm::vifstep(terra::rast(eds), th = 3)

# Subset environmental data based on VIF results
eds <- raster::subset(eds, v@results$Variables)
names(eds)

# # Plot the environmental data
# plot(eds, col = matlab.like(100))

# Run Maxent model
maxent_results <- run_maxent(occ_df, eds, survey)
beepr::beep(2)

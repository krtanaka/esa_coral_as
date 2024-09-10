require(tidyverse)
require(rgbif)
require(sf)
require(raster)
require(geodata)
require(ENMeval)
require(ecospat)
require(rJava)
require(dismo)
require(countrycode)
library(readr)
library(colorRamps)
library(ggmap)
library(ggspatial)
library(doParallel)

rm(list = ls())

select = dplyr::select

source("script/functions.R")

species_list <- c(
  "Acropora globiceps",
  "Isopora crateriformis"
)[2]

if (species_list == "Acropora globiceps") {
  
  ncrmp = read_csv("data/A_globiceps_AS.csv") %>%
    filter(ISLAND == "Tutuila", AdColDen > 0) %>%
    dplyr::select(LONGITUDE, LATITUDE) %>%
    mutate(Scientific.Name = species_list,
           Source = "NCRMP") %>%
    rename(Longitude = LONGITUDE, Latitude = LATITUDE)
  
  gbif = read_csv("data/gbif_occurances_Acropora_globiceps.csv") %>% 
    mutate(Scientific.Name = species_list,
           Source = "GBIF") %>% 
    select(Longitude, Latitude, Scientific.Name, Source)
  
  plot(ncrmp$Longitude, ncrmp$Latitude, col = 2, pch = 20); maps::map(add = T)
  points(gbif$Longitude, gbif$Latitude, col = 4, pch = 20)
  
  occ_df = rbind(ncrmp, gbif)

} else {
  
  ncrmp = read_csv("data/I_craterformis_AS.csv") %>%
    filter(ISLAND == "Tutuila", AdColDen > 0) %>%
    dplyr::select(LONGITUDE, LATITUDE) %>%
    mutate(Scientific.Name = species_list,
           Source = "NCRMP") %>%
    rename(Longitude = LONGITUDE, Latitude = LATITUDE)
  
  gbif = read_csv("data/gbif_occurances_Isopora_crateriformis.csv") %>% 
    mutate(Scientific.Name = species_list,
           Source = "GBIF") %>% 
    select(Longitude, Latitude, Scientific.Name, Source)
  
  occ_df = rbind(ncrmp, gbif)
  
  # map = ggmap::get_map(location = c(mean(ncrmp$Longitude), mean(ncrmp$Latitude)),
  #                      maptype = "satellite",
  #                      zoom = 11,
  #                      color = "bw",
  #                      force = T)
  # 
  # ggmap(map, darken = c(0.5, "black")) +
  #   geom_spatial_point(data = occ_df, aes(Longitude, Latitude, fill = Source, color = Source),
  #                      size = 4,
  #                      shape = 21, alpha = 0.7, crs = 4326) +
  #   annotate("text", x = min(ncrmp$Longitude), y = max(ncrmp$Latitude), label = "NCRMP & GBIF",
  #            hjust = 0, vjust = 1, size = 6, color = "white", fontface = "bold") +
  #   scale_y_continuous(limits = c(min(ncrmp$Latitude), max(ncrmp$Latitude)), "") +
  #   scale_x_continuous(limits = c(min(ncrmp$Longitude), max(ncrmp$Longitude)), "") +
  #   theme_minimal() +
  #   theme(legend.position = c(0.9, 0.2),
  #         legend.background = element_blank(),             # Transparent background
  #         legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
  #         legend.text = element_text(color = "white", face = "bold"),  # White and bold text
  #         legend.title = element_text(color = "white", face = "bold"))
  # 
  # ggsave(last_plot(), filename =  file.path("output/ncrmp_gbif.png"), height = 5, width = 10)
  
  
}

# Check how many occurrences subset for each spp.
table(occ_df$Scientific.Name)

load("data/eds.rdata")

v = usdm::vifstep(terra::rast(eds), th = 5)
eds = raster::subset(eds, v@results$Variables); names(eds)

plot(eds, col = matlab.like(100))

maxent_results = run_maxent(occ_df, eds)
beepr::beep(2)

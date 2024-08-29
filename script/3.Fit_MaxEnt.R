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
)[1]

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
  
  plot(ncrmp$Longitude, ncrmp$Latitude, col = 2, pch = 20); maps::map(add = T)
  points(gbif$Longitude, gbif$Latitude, col = 4, pch = 20)
  
  occ_df = rbind(ncrmp, gbif)
  
}

# Check how many occurrences subset for each spp.
table(occ_df$Scientific.Name)

load("data/eds.rdata")

v = usdm::vifstep(terra::rast(eds), th = 5)
eds = raster::subset(eds, v@results$Variables); names(eds)

plot(eds, col = matlab.like(100))

maxent_results = run_maxent(occ_df, eds)
beepr::beep(2)

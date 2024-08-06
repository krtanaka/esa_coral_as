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

source("functions.R")

species_list <- c(
  "Acropora globiceps",
  "Isopora crateriformis"
)[1]

if (species_list == "Acropora globiceps") df = read_csv("A_globiceps_AS.csv") 
if (species_list == "Isopora crateriformis") df = read_csv("I_craterformis_AS.csv") 

occ_df = df %>% 
  filter(ISLAND == "Tutuila") %>% 
  filter(AdColDen > 0) %>% 
  dplyr::select(LONGITUDE, LATITUDE) %>% 
  mutate(Scientific.Name = species_list) %>% 
  as.data.frame()

colnames(occ_df)[1:2] = c("Longitude", "Latitude")

# Check how many occurrences subset for each spp.
table(occ_df$Scientific.Name)

load("eds.rdata")

v = usdm::vifstep(terra::rast(eds), th = 10)
eds = raster::subset(eds, v@results$Variables); names(eds)

plot(eds, col = matlab.like(100))

maxent_results = run_maxent(occ_df, eds)
beepr::beep(2)

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

# Load environmental dataset
load("data/eds.rdata")

# Define species list and select the species
species <- c("Acropora globiceps", "Isopora crateriformis", "Genus Tridacna")[1]
survey <- c("ncrmp", "combined", "no_nps")[3]

# Load NCRMP occurrences
ncrmp <- read_csv(paste0("data/occurances_", species, "_ncrmp_exp.csv"))

# Set file paths based on the survey type
file_paths <- switch(survey,
                     
                     ncrmp = list(ncrmp = paste0("data/occurances_", species, "_ncrmp_exp.csv")),
                     
                     combined = list(
                       ncrmp = paste0("data/occurances_", species, "_ncrmp_exp.csv"),
                       gbif = paste0("data/occurances_", species, "_gbif_obis.csv"),
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

plot(occ_df$Longitude, occ_df$Latitude)
maps::map(add = T)

# Check the number of occurrences for each species
table(occ_df$Scientific.Name)

# Filter occurrences based on the NCRMP data
occ_df <- occ_df %>%
  filter(Longitude >= min(ncrmp$Longitude), Longitude <= max(ncrmp$Longitude),
         Latitude >= min(ncrmp$Latitude), Latitude <= max(ncrmp$Latitude)) %>%
  select(Longitude, Latitude, Scientific.Name) %>%
  distinct()

# Run Maxent model
maxent_results <- run_maxent(occ_df, eds, survey)
beepr::beep(2)

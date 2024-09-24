# https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
# https://data-blog.gbif.org/post/gbif-filtering-guide/

# install.packages("usethis")
# usethis::edit_r_environ()

library(rgbif)
library(dplyr)
library(CoordinateCleaner)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(robis)
library(patchwork)

rm(list = ls())

ggmap::register_google(key = "AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

map = ggmap::get_map(location = c(-170.7325, -14.3258),
                     maptype = "satellite",
                     zoom = 11,
                     # color = "bw",
                     force = TRUE)

species_list <- c("Acropora globiceps",
                  "Isopora crateriformis")

world_sf <- ne_countries(scale = "medium", returnclass = "sf")

pb <- txtProgressBar(min = 0, max = length(species_list), style = 3)

for (s in 1:length(species_list)) {
  
  # s = 1
  
  max_uncertainty <- 5000  # Adjust this value based on your needs
  
  species = species_list[s]
  species = gsub(" ", "_", species)
  
  setTxtProgressBar(pb, s)
  
  taxonkey <- name_backbone(species_list[s])$usageKey
  
  # set up gbif credentials first
  # https://docs.ropensci.org/rgbif/articles/gbif_credentials.html
  
  gbif_download <- occ_download(
    pred("taxonKey", taxonkey),
    pred("hasCoordinate", TRUE), 
    pred("hasGeospatialIssue", FALSE), # remove GBIF default geospatial issues
    format = "SIMPLE_CSV") 
  
  occ_download_wait(gbif_download) 
  
  # filtering GBIF  
  gbif <- gbif_download %>%
    occ_download_get() %>%  # Download the occurrence data from GBIF
    occ_download_import() %>%  # Import the downloaded data
    setNames(tolower(names(.))) %>%  # Set all column names to lowercase for consistency with CoordinateCleaner
    filter(occurrencestatus == "PRESENT") %>%  # Keep only records where species occurrence is marked as "PRESENT"
    filter(!basisofrecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")) %>%  # Exclude fossil and living specimens
    filter(year >= 1980) %>%  # Filter records from the year 2000 or later
    filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>%  # Keep records with high coordinate precision or where precision is not available
    filter(coordinateuncertaintyinmeters < 10000 | is.na(coordinateuncertaintyinmeters)) %>%  # Keep records with uncertainty below 10,000 meters or where it's not available
    filter(!coordinateuncertaintyinmeters %in% c(301, 3036, 999, 9999)) %>%  # Remove records with specific known issues in coordinate uncertainty
    filter(!decimallatitude == 0 | !decimallongitude == 0) %>%  # Exclude records where both latitude and longitude are 0
    mutate(decimalLongitude = decimallongitude,  # Standardize column names for longitude
           decimalLatitude = decimallatitude) %>%  # Standardize column names for latitude
    cc_cen(buffer = 2000) %>%  # Remove points near country centroids within 2km 
    cc_cap(buffer = 2000) %>%  # Remove points near capital city centroids within 2km
    cc_inst(buffer = 2000) %>%  # Remove points near known institutions (zoos, herbaria, etc.) within 2km
    cc_zero() %>%  # Remove records with zero coordinates (invalid locations)
    cc_equ() %>%  # Identify and remove records with identical latitude/longitude pairs
    cc_val() %>%  # Remove points with invalid coordinates (outside realistic bounds)
    # cc_outl() %>%  # Remove spatial outliers that deviate significantly from other points
    distinct(decimallongitude, decimallatitude, .keep_all = TRUE)  # Keep only unique records by longitude, latitude, species, dataset, and country
  
  # filtering OBIS data
  obis <- occurrence(
    scientificname = species_list[s],
    startdate = "1980-01-01") 
  
  if (nrow(obis) == 0) {
    
    obis = gbif[1,]
    
  }else{
    
    obis = obis %>% 
      filter(!is.na(decimalLongitude) & !is.na(decimalLatitude))
    
    if ("coordinateUncertaintyInMeters" %in% names(obis)) {
      
      obis = obis %>% 
        mutate(coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters)) %>% 
        filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%  # Remove records with missing coordinates
        filter(coordinateUncertaintyInMeters <= max_uncertainty | is.na(coordinateUncertaintyInMeters)) # Filter by uncertainty
      
    }
    
    obis = obis %>% 
      cc_cen(buffer = 2000) %>%  # Remove points near country centroids within 2km 
      cc_cap(buffer = 2000) %>%  # Remove points near capital city centroids within 2km
      cc_inst(buffer = 2000) %>%  # Remove points near known institutions (zoos, herbaria, etc.) within 2km
      cc_zero() %>%  # Remove records with zero coordinates (invalid locations)
      cc_equ() %>%  # Identify and remove records with identical latitude/longitude pairs
      cc_val() %>%  # Remove points with invalid coordinates (outside realistic bounds)
      # cc_outl() %>%  # Remove spatial outliers that deviate significantly from other points
      distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
    
  }
  
  # flag non-ocean points 
  obis$non_ocean = cc_sea(obis, value = "flagged")
  gbif$non_ocean = cc_sea(gbif, value = "flagged")
  
  obis <- obis %>% 
    dplyr::mutate(Scientific.Name = species,
                  Longitude  = decimalLongitude, 
                  Latitude  = decimalLatitude) %>%
    dplyr::select(Scientific.Name, Longitude, Latitude) %>%
    na.omit() %>% 
    distinct() %>% 
    mutate(Source = "OBIS")
  
  gbif <- gbif %>% 
    dplyr::rename(Scientific.Name = species,
                  Longitude  = decimallongitude, 
                  Latitude  = decimallatitude) %>%
    dplyr::select(Scientific.Name, Longitude, Latitude) %>%
    na.omit() %>% 
    distinct() %>% 
    mutate(Source = "GBIF")
  
  df = rbind(gbif, obis) %>% 
    distinct()
  
  dfi = df %>% 
    dplyr::select(Longitude, Latitude, Source) %>% 
    distinct() %>% 
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
  
  source_counts <- table(dfi$Source)
  
  labels <- c(
    'GBIF' = paste0('GBIF (n=', source_counts['GBIF'], ')'),
    'OBIS' = paste0('OBIS (n=', source_counts['OBIS'], ')')
  )
  
  facet_labels <- c(
    "GBIF" = paste0("GBIF (n=", source_counts['GBIF'], ")"), 
    "OBIS" = paste0("OBIS (n=", source_counts['OBIS'], ")")
  )
  
  p1 = ggplot() +
    geom_sf(data = world_sf, fill = "grey80", color = "grey60") +  # World map
    geom_sf(data = dfi, fill = "red", size = 3, shape = 21, alpha = 0.5) +  # Your data points
    coord_sf(crs = "+proj=robin") +  # Robinson projection
    theme_minimal() +
    facet_wrap(~Source, ncol = 2, labeller = labeller(Source = facet_labels)) +  # Custom facet titles
    theme(legend.position = "bottom")
  
  p2 = ggmap(map) +
    geom_spatial_point(data = df, 
                       aes(Longitude, Latitude, 
                           fill = Source,  
                           color = Source),
                       shape = 21, crs = 4326,
                       size = 5) +
    coord_sf(crs = 4326) +    # Use coord_sf to address the warning
    annotate("text", x = -170.8375, y = -14.23855, 
             label = paste0(species, "\n1980-2024\nSource: GBIF & OBIS"), 
             hjust = 0, vjust = 1, size = 4, color = "white", fontface = "bold") + 
    scale_y_continuous(limits = c(-14.36940, -14.23855), "") +
    scale_x_continuous(limits = c(-170.8375, -170.5450), "") +
    theme_minimal() + 
    theme(legend.position = c(0.92, 0.22),
          legend.background = element_blank(),             # Transparent background
          legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
          legend.text = element_text(color = "white", face = "bold"),  # White and bold text
          legend.title = element_text(color = "white", face = "bold"))
  
  p1
  
  ggsave(last_plot(), file = paste0("data/occurances_", species, "_gbif_obis_global.png"), height = 5)
  
  p2
  
  ggsave(last_plot(), file = paste0("data/occurances_", species, "_gbif_obis_local.png"), height = 5)
  
  readr::write_csv(df, file = paste0("data/occurances_", species, "_gbif_obis.csv"))
  
}

close(pb)

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)
library(ggmap)
library(ggspatial)

# Clear the workspace
rm(list = ls())

# Species list
species_list <- "Genus Tridacna"

for (s in 1:length(species_list)) {
  
  # s = 1
  
  species = species_list[s]
  
  if (s == 1) df = read_csv("data/original_data/AMSAM_giant_clams.csv")

  dfi = df %>% 
    filter(ISLAND == "Tutuila") %>% 
    mutate(y = 1,
           Longitude = CENTROIDLON,
           Latitude = CENTROIDLAT,
           Year = OBS_YEAR) %>% 
    group_by(Longitude, Latitude, Year) %>%
    summarise(y = mean(y, na.rm = T),
              y = ifelse(y > 0, "present", "absent")) %>% 
    mutate(Scientific.Name = species, 
           Source = "RAMP/NCRMP") %>% 
    dplyr::select(Longitude, Latitude, Year, Scientific.Name, Source) %>% 
    na.omit()
  
  readr::write_csv(dfi, file = paste0("data/occurances_", species, ".csv"))
  
  ggmap::register_google(key = "AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")
  
  # Get map for the given coordinates
  map = ggmap::get_map(location = c(mean(dfi$Longitude), mean(dfi$Latitude)),
                       maptype = "satellite",
                       zoom = 11,
                       # color = "bw",
                       force = TRUE)
  
  sum_p = dim(dfi)[1] %>% as.numeric()
  
  dfi$Year = as.character(dfi$Year)
  
  # Create plot
  ggmap(map) +
    geom_spatial_point(data = dfi, aes(Longitude, Latitude, fill = Year, color = Year),  # Map size to y
                       shape = 21, crs = 4326, alpha = 0.8, size = 4) +
    coord_sf(crs = 4326) +
    scale_y_continuous(limits = range(dfi$Latitude), "") +
    scale_x_continuous(limits = range(dfi$Longitude), "") +
    annotate("text", x = min(dfi$Longitude), y = max(dfi$Latitude), 
             label = paste0(species, "\n2002-2015\nSource: RAMP/NCRMP\nn = ", sum_p), 
             hjust = 0, vjust = 1, size = 4, color = "white", fontface = "bold") + 
    theme_minimal() + 
    theme(legend.position = c(0.92, 0.22),
          legend.background = element_blank(), # Transparent background
          legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
          legend.text = element_text(color = "white", face = "bold"),  # White and bold text
          legend.title = element_text(color = "white", face = "bold"))
  
  ggsave(last_plot(), file = paste0("data/occurances_", species, "_ncrmp.png"), width = 8)
  
}

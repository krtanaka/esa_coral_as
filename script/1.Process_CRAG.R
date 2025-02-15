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
species_list <- c("Acropora globiceps", "Isopora crateriformis")[2]

for (s in 1:length(species_list)) {
  
  # s = 1
  
  species = species_list[s]
  
  if (s == 1) df = read_csv("data/original_data/crag_isopora_occurences.csv")
  # if (s == 2) df = read_csv("data/original_data/I_craterformis_AS.csv")
  
  dfi = df %>% 
    filter(island == "Tutuila") %>% 
    mutate(y = 1,
           Longitude = longitude,
           Latitude = latitude) %>% 
    group_by(Longitude, Latitude) %>%
    summarise(y = mean(y, na.rm = T),
              y = ifelse(y > 0, "present", "absent")) %>% 
    mutate(Scientific.Name = species, 
           Source = "CRAG") %>% 
    dplyr::select(Longitude, Latitude, Scientific.Name, Source)
  
  readr::write_csv(dfi, file = paste0("data/occurances_", species, "_crag.csv"))
  
  ggmap::register_google(key = "AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")
  
  map = ggmap::get_map(location = c(mean(dfi$Longitude), mean(dfi$Latitude)),
                       maptype = "satellite",
                       zoom = 11,
                       # color = "bw",
                       force = TRUE)
  
  sum_p = dim(dfi)[1] %>% as.numeric()
  
  ggmap(map) +
    geom_spatial_point(data = dfi, aes(Longitude, Latitude),
                       shape = 21, crs = 4326, show.legend = F, color = "red", fill = "red", alpha = 0.5, size = 4) +
    coord_sf(crs = 4326) +
    scale_y_continuous(limits = c(-14.38, -14.22), "") +
    scale_x_continuous(limits = c(-170.85, -170.53), "") +
    annotate("text", x = -170.85, y = -14.22, 
             label = paste0(species, "\n2016-2020\nSource: CRAG\nn = ", sum_p), 
             hjust = 0, vjust = 1, size = 4, color = "white", fontface = "bold") + 
    theme_minimal() +
    theme(legend.position = c(0.92, 0.22),
          legend.background = element_blank(), 
          legend.key = element_rect(colour = NA, fill = NA), 
          legend.text = element_text(color = "white", face = "bold"),
          legend.title = element_text(color = "white", face = "bold"))
  
  ggsave(last_plot(), file = paste0("data/occurances_", species, "_crag.png"), width = 9, height = 4.8)
  
}

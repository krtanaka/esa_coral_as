# https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
# https://data-blog.gbif.org/post/gbif-filtering-guide/

# install.packages("usethis")
# usethis::edit_r_environ()

library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)
library(ggmap)
library(ggspatial)

rm(list = ls())

species_list <- c("Acropora globiceps",
                  "Isopora crateriformis")

for (s in 1:length(species_list)) {
  
  # s = 1
  
  species = species_list[s]
  
  load("data/npsa_benthic_data.rdata")
  
  dfi = df %>% 
    mutate(y = ifelse(Taxon_Name == species_list[s], 1, 0),
           island = ifelse(Longitude < -170.5, "Tutuila", "Ofu")) %>% 
    filter(island == "Tutuila") %>% 
    group_by(Longitude, Latitude) %>%
    summarise(y = mean(y, na.rm = T),
              y = ifelse(y > 0, "present", "absent"))
  
  ggmap::register_google(key = "AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")
  
  # Get map for the given coordinates
  map = ggmap::get_map(location = c(mean(dfi$Longitude), mean(dfi$Latitude)),
                       maptype = "satellite",
                       zoom = 13,
                       # color = "bw",
                       force = TRUE)
  
  sum_p = table(dfi$y)[2] %>% as.numeric()
  
  # Create plot
  ggmap(map) +
    geom_spatial_point(data = dfi, aes(Longitude, Latitude, 
                                       fill = y,  
                                       color = y,
                                       alpha = y,
                                       size = y),  # Map size to y
                       shape = 21, crs = 4326, show.legend = F) +
    scale_fill_manual(values = c("absent" = "white", "present" = "red")) +  # Set colors
    scale_color_manual(values = c("absent" = "white", "present" = "red")) + 
    scale_alpha_manual(values = c("absent" = 0.4, "present" = 0.9)) +  # Dim absent, full opacity for present
    scale_size_manual(values = c("absent" = 2, "present" = 4)) +  # Smaller for absent, larger for present
    coord_sf(crs = 4326) +    # Use coord_sf to address the warning
    scale_y_continuous(limits = c(-14.28128 , -14.22946), "") +
    scale_x_continuous(limits = c(-170.7243, -170.6528), "") +
    annotate("text", x = -170.7243, y = -14.22946, label = paste0(species, "\n2007-2019\nSource = NPS\nn = ", sum_p), 
             hjust = 0, vjust = 1, size = 5, color = "white", fontface = "bold") + 
    theme_minimal() + 
    theme(legend.position = c(0.92, 0.22),
          legend.background = element_blank(),             # Transparent background
          legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
          legend.text = element_text(color = "white", face = "bold"),  # White and bold text
          legend.title = element_text(color = "white", face = "bold"))
  
  ggsave(last_plot(), file = paste0("data/nps_occurances_", species, ".png"), height = 6, width = 8)
  readr::write_csv(df, file = paste0("data/nps_occurances_", species, ".csv"))
  
}

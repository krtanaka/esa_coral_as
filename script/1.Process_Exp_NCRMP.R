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
species_list <- c("Acropora globiceps", "Isopora crateriformis")

# Annotation:
# Expanded occurrence data from CoralBelt_Adults_raw_CLEANED_2023.csv
# Includes AGLO and conspecifics (n = 132 occurrences) for the years 2015, 2018, 2020, and 2023.
# Mean depths range from 2.1m to 23.7m.
# ICRA data (n = 1475 occurrences) shows consistent scoring in 2023, 2020, and 2018 (~98%-100%), 
# but a significant drop to 22% in 2015, suggesting possible misclassification at the genus level in 2015.
# The 2015 data has been filtered to include ISSP colonies with ICRA-compatible morphology (e.g., "Plating", "Encrusting", "Laminar"), 
# based on observations of knobby or branching columnar morphologies for IPAL in 2023.

for (s in 1:length(species_list)) {
  
  # s = 2
  
  species = species_list[s]
  
  if (s == 1) df = read_csv("data/original_data/TUT_AGLO & conspecifics.csv")
  if (s == 2) df = read_csv("data/original_data/TUT_ICRA expanded 2015 filter.csv")
  
  dfi = df %>% 
    mutate(y = ifelse(TAXONNAME == species_list[s], 1, 0),
           Longitude = LONGITUDE,
           Latitude = LATITUDE) %>% 
    group_by(Longitude, Latitude) %>%
    summarise(y = mean(y, na.rm = T),
              y = ifelse(y > 0, "present", "absent"))
  
  df = dfi %>%
    filter(y == "present") %>% 
    mutate(Scientific.Name = species, 
           Source = "Exp.NCRMP") %>% 
    dplyr::select(Longitude, Latitude, Scientific.Name, Source)
  
  readr::write_csv(df, file = paste0("data/occurances_", species, "_ncrmp_exp.csv"))
  
  ggmap::register_google(key = "AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")
  
  # Get map for the given coordinates
  map = ggmap::get_map(location = c(mean(dfi$Longitude), mean(dfi$Latitude)),
                       maptype = "satellite",
                       zoom = 11,
                       # color = "bw",
                       force = TRUE)
  
  sum_p = table(dfi$y)[2] %>% as.numeric()
  
  ggmap(map) +
    geom_spatial_point(data = dfi, aes(Longitude, Latitude, 
                                       fill = y,  
                                       color = y,
                                       alpha = y,
                                       size = y),  
                       shape = 21, crs = 4326, show.legend = F) +
    scale_fill_manual(values = c("absent" = "yellow", "present" = "red")) +  
    scale_color_manual(values = c("absent" = "yellow", "present" = "red")) + 
    scale_alpha_manual(values = c("absent" = 0.5, "present" = 0.8)) + 
    scale_size_manual(values = c("absent" = 2, "present" = 3)) + 
    coord_sf(crs = 4326) + 
    scale_y_continuous(limits = c(-14.38, -14.22), "") +
    scale_x_continuous(limits = c(-170.85, -170.53), "") +
    annotate("text", x = -170.85, y = -14.22, 
             label = paste0(species, "\n2015-2023\nSource: NCRMP\nn = ", sum_p), 
             hjust = 0, vjust = 1, size = 4, color = "white", fontface = "bold") + 
    theme_minimal() + 
    theme(legend.position = c(0.92, 0.22),
          legend.background = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.text = element_text(color = "white", face = "bold"), 
          legend.title = element_text(color = "white", face = "bold"))
  
  ggsave(last_plot(), file = paste0("data/occurances_", species, "_ncrmp_exp.png"), width = 9, height = 4.8)
  
}

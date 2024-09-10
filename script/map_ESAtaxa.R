# Load necessary libraries
library(dplyr)
library(ggplot2)
library(colorRamps)
library(patchwork)
library(sf)
library(ggmap)
library(tidyr)

# Clear the environment
rm(list = ls())

# Define species list
sp_list = c("I_craterformis", "A_globiceps")

# Initialize empty list to store plots
p = list()

# Loop through species list
for (s in sp_list) {
  
  # Read in data for each species
  df = read.csv(paste0("data/", s, "_AS.csv"))
  
  # Create histogram of AdColDen
  hist(df$AdColDen)
  
  # Filter dataset based on specific conditions
  df = df %>% filter(ISLAND %in% c("Tutuila", "Ofu & Olosega"))
  df = df %>% filter(ISLAND %in% c("Tutuila"))
  df = df %>% filter(OBS_YEAR != 2020)
  df$lon = df$LONGITUDE
  df$lat = df$LATITUDE
  
  # Summarize data by lon and lat, calculating the mean of AdColCount
  df = df %>% 
    group_by(lon, lat) %>% 
    summarise(AdColCount = mean(AdColCount, na.rm = TRUE))
  
  # Register Google API key (replace with your own API key)
  ggmap::register_google(key = "AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")
  
  # Get map for the given coordinates
  map = ggmap::get_map(location = c(mean(df$lon), mean(df$lat)),
                       maptype = "satellite",
                       zoom = 11,
                       color = "bw",
                       force = TRUE)
  
  # Create plot
  pp = ggmap(map) +
    geom_spatial_point(data = df, aes(lon, lat, 
                                      fill = AdColCount,  # Map fill to AdColCount
                                      color = AdColCount, # Map color to AdColCount
                                      size = AdColCount), 
                       shape = 21, alpha = 0.8, crs = 4326) + 
    annotate("text", x = min(df$lon), y = max(df$lat), label = s, 
             hjust = 0, vjust = 1, size = 6, color = "white", fontface = "bold") + 
    scale_y_continuous(limits = c(min(df$lat), max(df$lat)), "") +
    scale_x_continuous(limits = c(min(df$lon), max(df$lon)), "") +
    scale_color_gradient(low = "yellow", high = "red", guide = "legend") +  # Define color gradient
    scale_fill_gradient(low = "yellow", high = "red", guide = "legend") +   # Define fill gradient
    theme_minimal() + 
    theme(legend.position = c(0.9, 0.2),
          legend.background = element_blank(),             # Transparent background
          legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
          legend.text = element_text(color = "white", face = "bold"),  # White and bold text
          legend.title = element_text(color = "white", face = "bold"))
  
  print(pp)
  
  ggsave(last_plot(), filename =  file.path(paste0("output/ncrmp_", s, ".png")), height = 5, width = 11)
  
  p[[s]] = pp
}

# Display the list of plots
p[[1]] / p[[2]]

ggsave(last_plot(), filename =  file.path("output/ncrmp.png"), height = 10, width = 11)

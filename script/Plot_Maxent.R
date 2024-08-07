library(readr)
library(raster)
library(terra)
library(dplyr)
library(colorRamps)
library(ggmap)
library(ggspatial)

# run "2.Prep_Prediction_Layers.R" first

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

load("maxent_result_Isopora crateriformis.rda")
load("maxent_result_Acropora globiceps.rda")

plot(maxent_result$model)

r <- predict(maxent_result$model, eds)
plot(r, col = matlab.like(100))
r = rasterToPoints(r) %>% as.data.frame()

# use ggmap
ggmap::register_google("AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

# Get the coordinates of the cell centers
coords <- coordinates(eds %>% stack())

# Calculate the mean latitude and longitude
mean_lat <- mean(coords[, 2], na.rm = TRUE)
mean_lon <- mean(coords[, 1], na.rm = TRUE)

map = ggmap::get_map(location = c(mean_lon, mean_lat),
                     maptype = "satellite",
                     zoom = 11,
                     force = T)
ggmap(map) +
  geom_spatial_point(data = r, aes(x, y, fill = layer, color = layer), 
                     size = 8,
                     shape = 22, alpha = 0.8, crs = 4326) + 
  scale_fill_gradientn(colors = matlab.like(100), "Predicted \nOccupancy \n(0-1)") + 
  scale_color_gradientn(colors = matlab.like(100), "Predicted \nOccupancy \n(0-1)") + 
  # ggtitle("Spatial distribution of U. stolonifera predicted habitat suitability") +
  theme(legend.position = c(0.92, 0.81),
        legend.background = element_blank(), # Makes the legend background transparent
        legend.box.background = element_blank(), # Makes the legend box background transparent
        legend.text = element_text(color = "white"), # Makes the legend text white
        legend.title = element_text(color = "white") # Makes the legend title white
  )

ggsave(last_plot(), filename =  file.path("output/SDM_output.png"), height = 5.5, width = 5.5)

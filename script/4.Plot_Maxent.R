library(readr)
library(raster)
library(terra)
library(dplyr)
library(colorRamps)
library(ggmap)
library(ggspatial)

# run "2.Prep_Prediction_Layers.R" first

source("script/functions.R")

species_list <- c(
  "Acropora globiceps",
  "Isopora crateriformis"
)[1]

load(paste0("output/maxent_result_", species_list, ".rda"))

maxent_result$model@results %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("rowname") %>% 
  filter(grepl("contribution", rowname)) %>% 
  mutate(rowname = gsub(".contribution", "", rowname)) %>% 
  tibble::column_to_rownames("rowname") %>% 
  rownames_to_column("rowname") %>% 
  arrange(V1) %>%
  mutate(rowname = factor(rowname, levels = rowname)) %>% 
  ggplot(aes(V1, rowname, fill = V1)) + 
  labs(x = "%", y = "", title = paste0("Variable Contribution for ", species_list)) +
  geom_point(shape = 21, size = 5, show.legend = F) + 
  scale_fill_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt")

load("data/eds.rdata")

r <- predict(maxent_result$model, eds); plot(r, col = matlab.like(100))
r <- rasterToPoints(r) %>% as.data.frame()

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
                     size = 4,
                     shape = 22, alpha = 0.7, crs = 4326) + 
  scale_fill_gradientn(colors = matlab.like(100), "", limits = c(0,1)) + 
  scale_color_gradientn(colors = matlab.like(100), "", limits = c(0,1)) + 
  ggtitle(paste0("Predicted probability of presence for ", species_list)) +
  theme(legend.position = c(0.1, 0.85),
        legend.background = element_blank(), # Makes the legend background transparent
        legend.box.background = element_blank(), # Makes the legend box background transparent
        legend.text = element_text(color = "white"), # Makes the legend text white
        legend.title = element_text(color = "white") # Makes the legend title white
  )

ggsave(last_plot(), filename =  file.path(paste0("output/maxent_map_", species_list, ".png")), height = 5.5, width = 5.5)

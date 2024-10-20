library(readr)
library(raster)
library(terra)
library(dplyr)
library(colorRamps)
library(ggmap)
library(ggspatial)
library(tidyverse)

# run "2.Prep_Prediction_Layers.R" first

source("script/functions.R")

species_list <- c("Acropora globiceps", "Isopora crateriformis", "Genus Tridacna")[1]

load(paste0("output/maxent_result_", species_list, ".rda"))

var_list = maxent_result$model@results %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("rowname") %>% 
  filter(grepl("contribution", rowname)) %>% 
  mutate(rowname = gsub(".contribution", "", rowname)) %>% 
  tibble::column_to_rownames("rowname") %>% 
  rownames_to_column("rowname") %>% 
  arrange(V1) %>%
  mutate(rowname = factor(rowname, levels = rowname)) 

var_list %>% 
  filter(V1 > 0) %>% 
  ggplot(aes(V1, rowname, fill = V1, color = V1)) + 
  labs(x = "Variable Contribution %", y = "", title = species_list) +
  geom_point(shape = 21, size = 5, show.legend = F) + 
  scale_x_log10("log10 Variable Contribution (%)") +
  ggdark::dark_theme_classic(base_size = 15) +
  scale_fill_gradientn(colors = colorRamps::matlab.like(100), trans = "log10") + 
  scale_color_gradientn(colors = colorRamps::matlab.like(100), trans = "log10")

# var_list %>% 
#   filter(V1 > 0) %>% 
#   ggplot(aes(V1, rowname, fill = V1, color = V1)) + 
#   labs(x = "Variable Contribution %", y = "", title = species_list) +
#   geom_point(shape = 21, size = 5, show.legend = F) + 
#   ggdark::dark_theme_classic(base_size = 12) +
#   scale_fill_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") + 
#   scale_color_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt")

ggsave(last_plot(),
       filename =  file.path(paste0("output/maxent_var_", species_list, ".png")), 
       height = 6, width = 5.5)

env_vars <- names(maxent_result$model@presence)

response_list <- list()

for (var in env_vars) {
  
  # var = "population_density"
  
  response_data <- response(maxent_result$model, var = var, expand = TRUE, plot = FALSE)
  
  response_df <- data.frame(x = response_data[, 1], y = response_data[, 2])
  
  response_df = response_df %>% 
    mutate(x = round(x, 2)) %>% 
    group_by(x) %>% 
    summarise(y = mean(y)) %>% 
    mutate(variable = var)
  
  # Append the data frame to the list
  response_list[[var]] <- response_df
  
  # p = response_df %>% 
  #   ggplot(aes(x = x, y = y, fill = y)) +
  #   geom_point(shape = 21, size = 5, show.legend = F) +
  #   theme_classic() +
  #   scale_fill_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") + 
  #   labs(x = var, y = "Predicted Suitability",
  #        title = species_list)
  # 
  # response_list[[var]] <- p
  
}

response_combined <- bind_rows(response_list)

var_list = var_list %>% filter(V1 > 0) %>% select(rowname)

response_combined = response_combined %>% filter(variable %in% var_list$rowname)

reversed_var_order <- rev(var_list$rowname)

response_combined$variable <- factor(response_combined$variable, 
                                     levels = reversed_var_order)

response_combined %>% 
  ggplot(aes(x = x, y = y, fill = y, color = y)) +
  geom_point(shape = 21, size = 1, show.legend = F) +
  facet_wrap(~variable, scales = "free") +
  scale_fill_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") +
  scale_color_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") +
  ggdark::dark_theme_minimal() + 
  labs(x = "Environmental Variable", 
       y = "Predicted Suitability",
       title = paste(species_list, "maxent response curves"))

ggsave(last_plot(), 
       filename =  file.path(paste0("output/maxent_response_", species_list, ".png")),
       height = 8, width = 13)

response_combined %>% 
  filter(variable == "bathymetry") %>% 
  ggplot(aes(x = x, y = y, fill = y, color = y)) +
  geom_point(shape = 21, size = 1, show.legend = F) +
  scale_fill_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") +
  scale_color_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") +
  ggdark::dark_theme_minimal() + 
  labs(x = "Bathymetry (m)", 
       y = "Predicted Suitability",
       title = species_list)

ggsave(last_plot(), 
       filename =  file.path(paste0("output/maxent_response_bathymetry_", species_list, ".png")), 
       height = 3, width = 4)

auc = maxent_result$model@results %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("rowname") %>% 
  filter(grepl("Training.AUC", rowname)) %>% 
  select(V1) %>% 
  as.numeric() %>% 
  round(2)

load("data/eds.rdata")

r <- predict(maxent_result$model, eds); plot(r, col = matlab.like(100))
r <- readAll(r); save(r, file = paste0(paste0("output/maxent_raster_", species_list, ".rdata")))
r <- rasterToPoints(r) %>% as.data.frame()

# use ggmap
ggmap::register_google("AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

# Get the coordinates of the cell centers
coords <- coordinates(eds %>% stack())

# Calculate the mean latitude and longitude
mean_lat <- mean(coords[, 2], na.rm = TRUE)
mean_lon <- mean(coords[, 1], na.rm = TRUE)
min_lat <- min(coords[, 2])
max_lat <- max(coords[, 2])
min_lon <- min(coords[, 1])
max_lon <- max(coords[, 1])

map = ggmap::get_map(location = c(mean_lon, mean_lat),
                     maptype = "satellite",
                     zoom = 11,
                     color = "bw",
                     force = T)

ggmap(map, darken = c(0.5, "black")) +
  geom_spatial_point(data = r, aes(x, y, fill = layer, color = layer), 
                     size = 0.5,
                     shape = 22, alpha = 0.8, crs = 4326) + 
  annotate("text", x = min(r$x), y = max(r$y),
           label = paste0(species_list, "\nPred. Occ. Prob.\nAUC = ",auc),
           hjust = 0, vjust = 1, size = 4, color = "white", fontface = "bold") +
  scale_fill_gradientn(colors = matlab.like(100), "", limits = c(0,1)) + 
  scale_color_gradientn(colors = matlab.like(100), "", limits = c(0,1)) + 
  scale_y_continuous(limits = range(r$y), "") +
  scale_x_continuous(limits = range(r$x), "") +
  ggdark::dark_mode() + 
  theme(legend.position = c(0.95, 0.25),
        legend.background = element_blank(), 
        legend.box.background = element_blank(), 
        legend.text = element_text(color = "white"), 
        legend.title = element_text(color = "white"))

ggsave(last_plot(), filename =  file.path(paste0("output/maxent_map_", species_list, ".png")), width = 10)

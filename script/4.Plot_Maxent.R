library(readr)
library(raster)
library(terra)
library(dplyr)
library(colorRamps)
library(ggmap)
library(ggspatial)
library(tidyverse)
library(dismo)
library(patchwork)
library(cowplot)
library(ggpubr)

rm(list = ls())

select = dplyr::select

load("data/eds.rdata")

# run "2.Prep_Prediction_Layers.R" first

source("script/functions.R")

species_list <- c("Acropora globiceps", "Isopora crateriformis")[1]
survey_list <- c("ncrmp", "combined", "no_nps")[2]

ggmap::register_google("AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

map1 = ggmap::get_map(location = c(-170.705, -14.294),
                      maptype = "satellite",
                      zoom = 11,
                      force = T)

map2 = ggmap::get_map(location = c(-170.68846, -14.25177),
                      maptype = "satellite",
                      zoom = 13,
                      force = T)

map3 = ggmap::get_map(location = c(-170.705, -14.294),
                      maptype = "satellite",
                      zoom = 11,
                      color = "bw",
                      force = T)

map4 = ggmap::get_map(location = c(-170.68846, -14.25177),
                      maptype = "satellite",
                      zoom = 13,
                      color = "bw",
                      force = T)

for (species in species_list) {
  for (survey in survey_list) {
    
    # species = "Isopora crateriformis"
    # species = "Acropora globiceps"
    # survey = "combined"
    
    file_path <- paste0("output/maxent_result_", species, "_", survey, ".rda")
    
    if (file.exists(file_path)) {
      
      load(file_path)
      
      var_list <- maxent_result$model@results %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("rowname") %>% 
        filter(grepl("contribution", rowname)) %>% 
        mutate(rowname = gsub(".contribution", "", rowname)) %>% 
        tibble::column_to_rownames("rowname") %>% 
        rownames_to_column("rowname") %>% 
        mutate(rowname = factor(rowname, levels = rowname)) %>%
        arrange(-V1)
      
      var_list %>%
        mutate(rowname = gsub("^dhw\\.", "", rowname),
               rowname = gsub("_degree_heating_weeks", "_dhw", rowname),
               rowname = gsub("_chlorophyll_a", "_chla", rowname),
               rowname = gsub("meandur_", "mean_dur_", rowname),
               rowname = gsub("_viirs", "", rowname),
               rowname = gsub("distance.from.port", "distance_from_port", rowname)) %>% 
        filter(V1 > 0) %>%
        ggplot(aes(V1, reorder(rowname, V1), fill = V1)) +  # Reorder rowname by V1
        labs(x = "log10 Variable Contribution (%)", y = "") +
        geom_point(shape = 21, size = 5, show.legend = FALSE) + 
        scale_x_log10(labels = scales::label_number(), limits = c(0.001, 100)) +
        theme_classic(base_size = 20) + 
        scale_fill_gradientn(colors = colorRamps::matlab.like(100), trans = "log10") +
        scale_color_gradientn(colors = colorRamps::matlab.like(100), trans = "log10") +
        annotate("text", x = 0.001, y = Inf, 
                 label = species, 
                 hjust = 0, vjust = 1,
                 size = 6, 
                 # color = "white", 
                 fontface = "bold") + 
        # annotate("text", x = Inf, y = -Inf, 
        #          label = paste0("survey data = ", survey), 
        #          hjust = 1, vjust = -0.5, size = 6, 
        #          # color = "white", 
        #          fontface = "bold") + 
        theme(#axis.text.y = element_text(color = "white"),
          plot.margin = margin(5, 30, 5, 5))
      
      ggsave(last_plot(),
             filename = file.path(paste0("output/maxent_vart_", species, "_", survey, ".png")),
             height = 5, width = 8)
      
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
        
        response_list[[var]] <- response_df
        
      }
      
      response_combined <- bind_rows(response_list)
      
      var_order = var_list %>%
        filter(V1 > 0) %>%
        select(rowname) %>% 
        mutate(rowname = gsub("^dhw\\.", "", rowname),
               rowname = gsub("_degree_heating_weeks", "_dhw", rowname),
               rowname = gsub("_chlorophyll_a", "_chla", rowname),
               rowname = gsub("meandur_", "mean_dur_", rowname),
               rowname = gsub("_viirs", "", rowname),
               rowname = gsub("distance.from.port", "distance_from_port", rowname))
      
      response_combined = response_combined %>% 
        mutate(variable = gsub("^dhw\\.", "", variable),
               variable = gsub("_degree_heating_weeks", "_dhw", variable),
               variable = gsub("_chlorophyll_a", "_chla", variable),
               variable = gsub("meandur_", "mean_dur_", variable),
               variable = gsub("_viirs", "", variable),
               variable = gsub("distance.from.port", "distance_from_port", variable))
      
      response_combined = response_combined %>% 
        filter(variable %in% var_order$rowname) 
      
      var_order <- var_order$rowname
      
      response_combined$variable <- factor(response_combined$variable, levels = var_order)
      
      response_combined  %>%
        ggplot(aes(x = x, y = y, fill = y, color = y)) +
        geom_point(shape = 21, size = 3, show.legend = F, alpha = 0.5) +
        facet_wrap(~variable, scales = "free", ncol = 3) +
        scale_fill_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") +
        scale_color_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") +
        # ggdark::dark_mode() + 
        theme_pubr() + 
        labs(x = "", 
             y = "Predicted Suitability",
             title = species
             # title = paste(species, "MaxEnt response curves")
             )
      
      ggsave(last_plot(), 
             filename =  file.path(paste0("output/maxent_response_", species, "_", survey, ".png")),
             height = 9, width = 8)
      
      response_combined %>% 
        filter(variable == "bathymetry") %>% 
        ggplot(aes(x = x, y = y, fill = y, color = y)) +
        geom_point(shape = 21, size = 3, show.legend = F, alpha = 0.5) +
        facet_wrap(~variable, scales = "free") +
        scale_fill_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") +
        scale_color_gradientn(colors = colorRamps::matlab.like(100), trans = "sqrt") +
        # ggdark::dark_mode() + 
        theme_pubr() + 
        labs(x = "Bathymetry (m)", 
             y = "Predicted Suitability",
             title = species_list)
      
      ggsave(last_plot(), 
             filename =  file.path(paste0("output/maxent_response_bathymetry_", species, "_", survey, ".png")),
             height = 3, width = 4)
      
      auc = maxent_result$model@results %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("rowname") %>% 
        filter(grepl("Training.AUC", rowname)) %>% 
        select(V1) %>% 
        as.numeric() %>% 
        round(2)
      
      r <- predict(maxent_result$model, eds)#; plot(r, col = matlab.like(100))
      r <- rast(r) %>% terra::as.data.frame(xy = T)
      # r <- readAll(r)
      save(r, file = paste0(paste0("output/maxent_raster_", species, "_", survey, ".rdata")))
      
      # Load NCRMP occurrences
      ncrmp <- read_csv(paste0("data/occurances_", species, "_ncrmp_exp.csv"))
      
      # Set file paths based on the model
      file_paths <- if (survey == "ncrmp") {
        
        list(ncrmp = paste0("data/occurances_", species, "_ncrmp_exp.csv"))
        
      } else {
        
        list(
          ncrmp = paste0("data/occurances_", species, "_ncrmp_exp.csv"),
          gbif = paste0("data/occurances_", species, "_gbif_obis.csv"),
          nps = paste0("data/occurances_", species, "_nps.csv"),
          crag = paste0("data/occurances_", species, "_crag.csv")
        )
      }
      
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
      
      p1 = ggmap(map1) +
        geom_raster(data = r, aes(x = x, y = y, fill = layer)) +
        annotate("text", x = -170.85, y = -14.22,
                 label = paste0(species, "\nAUC = ", auc, "\nsurvey = ", survey),
                 hjust = 0, vjust = 1, size = 4, color = "white", fontface = "bold") +
        scale_fill_gradientn(colors = matlab.like(100), "Predicted Habitat Suitability", limits = c(0,1),
                             breaks = c(0, 0.5, 1), guide = guide_colorbar(direction = "horizontal",
                                                                           title.position = "top",
                                                                           barwidth = 9.5, barheight = 1)) +
        scale_color_gradientn(colors = matlab.like(100), "Predicted Habitat Suitability", limits = c(0,1),
                              breaks = c(0, 0.5, 1), guide = guide_colorbar(direction = "horizontal",
                                                                            title.position = "top",
                                                                            barwidth = 9.5, barheight = 1)) +
        scale_y_continuous(limits = c(-14.38, -14.22), "") +
        scale_x_continuous(limits = c(-170.85, -170.53), "") +
        theme_minimal(base_size = 10) + 
        theme(legend.position = c(0.8, 0.12),
              legend.background = element_blank(), 
              legend.box.background = element_blank(), 
              legend.text = element_text(color = "white", size = 10, face = "bold"), 
              legend.title = element_text(color = "white", face = "bold")) + 
        coord_sf(crs = 4326)
      
      p2 = ggmap(map2) +
        geom_raster(data = r, aes(x = x, y = y, fill = layer)) +
        scale_fill_gradientn(colors = matlab.like(100), limits = c(0,1), 
                             breaks = c(0, 0.5, 1), guide = "none") + 
        scale_color_gradientn(colors = matlab.like(100), limits = c(0,1), 
                              breaks = c(0, 0.5, 1), guide = "none") + 
        scale_y_continuous(limits = c(-14.28128, -14.22946), "") +
        scale_x_continuous(limits = c(-170.7243, -170.6528), "") +
        theme_minimal(base_size = 10) + 
        coord_sf(crs = 4326)
      
      combined_plot <- p1 + p2
      
      ggsave(plot = combined_plot,
             filename =  file.path(paste0("output/maxent_map_", species, "_", survey, ".png")), 
             width = 15.5, 
             height = 4.9,
             limitsize = FALSE,
             bg = "transparent")
      
      p1 = ggmap(map3) +
        geom_raster(data = r, aes(x = x, y = y, fill = layer), alpha = 0.8) +
        geom_point(data = occ_df, aes(Longitude, Latitude), fill = "green", color = "green", shape = 21, size = 2, alpha = 0.8) +
        annotate("text", x = -170.85, y = -14.22,
                 label = paste0(species, "\nAUC = ", auc, "\nsurvey = ", survey),
                 hjust = 0, vjust = 1, size = 4, color = "white", fontface = "bold") +
        scale_fill_viridis_c("Predicted Habitat Suitability", limits = c(0,1), option = "magma",
                             breaks = c(0, 0.5, 1), guide = guide_colorbar(direction = "horizontal",
                                                                           title.position = "top",
                                                                           barwidth = 9.5, barheight = 1)) +
        scale_color_viridis_c("Predicted Habitat Suitability", limits = c(0,1),
                              breaks = c(0, 0.5, 1), guide = guide_colorbar(direction = "horizontal",
                                                                            title.position = "top",
                                                                            barwidth = 9.5, barheight = 1)) +
        scale_y_continuous(limits = c(-14.38, -14.22), "") +
        scale_x_continuous(limits = c(-170.85, -170.53), "") +
        theme_minimal(base_size = 10) + 
        theme(legend.position = c(0.8, 0.12),
              panel.background = element_rect(fill = "gray80"),
              panel.grid.major = element_blank(),
              legend.background = element_blank(), 
              legend.box.background = element_blank(), 
              legend.text = element_text(color = "white", size = 10, face = "bold"), 
              legend.title = element_text(color = "white", face = "bold")) + 
        coord_sf(crs = 4326)
      
      p2 = ggmap(map4) +
        geom_raster(data = r, aes(x = x, y = y, fill = layer), alpha = 0.8) +
        geom_point(data = occ_df, aes(Longitude, Latitude), fill = "green", color = "green", shape = 21, size = 3, alpha = 0.8) +
        scale_fill_viridis_c(limits = c(0,1), option = "magma",
                             breaks = c(0, 0.5, 1), guide = "none") + 
        scale_color_viridis_c(limits = c(0,1), option = "magma",
                              breaks = c(0, 0.5, 1), guide = "none") + 
        scale_y_continuous(limits = c(-14.28128, -14.22946), "") +
        scale_x_continuous(limits = c(-170.7243, -170.6528), "") +
        theme_minimal(base_size = 10) +
        coord_sf(crs = 4326)
      
      combined_plot <- p1 + p2
      
      ggsave(plot = combined_plot,
             filename =  file.path(paste0("output/maxent_map_", species, "_", survey, "_", "surveypoints.png")), 
             width = 15.5, 
             height = 4.9,
             limitsize = FALSE,
             bg = "transparent")
      
    } else {
      
      cat("File not found for species:", species, ", survey:", survey, "\n")
      
    }
  }
}

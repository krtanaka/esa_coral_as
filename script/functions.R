# ---- Run MaxEnt Batch ----

run_maxent = function(occ_sf, env, survey, partition) {
  
  set.seed(2024)
  
  occ_sf = occ_df
  env = eds
  
  # survey = "ncrmp"
  
  # partition = "randomkfold"
  # partition = "block"
  
  species = unique(occ_sf$Scientific.Name)
  
  cat("Running MaxEnt calibration for species", 
      species, "with", 
      survey, "survey data, using", 
      partition, "partitioning...\n")
  
  # Loop through each species
  for(sp in species) {
    
    # sp = unique(occ_sf$Scientific.Name)
    
    # Subset occurrence data for the species
    occ_sp = occ_sf[occ_sf$Scientific.Name == sp, c("Longitude", "Latitude")]
    
    # crs(env) = "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
    crs(env) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
    
    # plot(env[["bathymetry"]]); points(occ_sp, col='red', pch = 20) # Plot the first environmental layer
    
    tryCatch({
      
      if (partition == "randomkfold") {
        
        enmeval_results = ENMevaluate(occ_sp, env,
                                      bg = NULL,
                                      tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
                                      partitions = "randomkfold",
                                      partition.settings = list(kfolds = 5),
                                      algorithm = "maxnet",
                                      parallel = T,
                                      numCores = detectCores() / 2,
                                      updateProgress = T,
                                      taxon.name = sp)
        
      } else if (partition == "block") {
        
        enmeval_results = ENMevaluate(occ_sp, env,
                                      bg = NULL,
                                      tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
                                      partitions = "block",
                                      algorithm = "maxnet",
                                      parallel = T,
                                      numCores = detectCores() / 2,
                                      updateProgress = TRUE,
                                      taxon.name = sp)
        
      }
      
    }, error = function(e) {
      
      closeAllConnections()
      cat("Error occurred during MaxEnt calibration: ", e$message, "\nConnections closed.\n")
      stop(e)
      
    })
    
    enmeval_df = enmeval_results@results
    
    readr::write_csv(enmeval_results@results %>%
                       mutate(across(where(is.numeric), ~ round(.x, 3))), 
                     file = paste0("output/maxent_result_enmeval_df_", sp, "_", survey, "_", partition, ".csv"))
    
    # Subset the ENMeval results to get the best model
    enmeval_bestm = subset(enmeval_df, delta.AICc == 0)
    
    # Decode the features
    maxent_feats = as.character(enmeval_bestm$fc)[1]
    maxent_rm = as.character(enmeval_bestm$rm)[1]
    
    # Print out the results
    cat(paste("Best features for", sp, ":", maxent_feats, "\n"))
    cat(paste("Best regularization multiplier for", sp, ":", maxent_rm, "\n"))
    
    # Run MaXent SDM
    sp_maxent_model = dismo::maxent(env, as.matrix(occ_sp), features = maxent_feats, betamultiplier = maxent_rm)
    
    # Store maxent model
    maxent_result = list(enm = enmeval_df, model = sp_maxent_model)
    save(maxent_result, file = paste0("output/maxent_result_", sp, "_", survey, "_", partition, ".rda"))
    
  }
  
  # return(list(enm = enm_results, models = maxent_models))
}

# ---- Clip Prediction Raster by Species -+---

spp_clip_raster <- function(spp, env_rs, domain, depth_cutoff) {
  
  # domain <- "Tutuila"
  
  box <- readr::read_csv("data/Bounding_Boxes.csv")
  
  box_i <- box %>% 
    group_by(unit) %>% 
    summarise(ymin = min(ymin), 
              ymax = max(ymax),
              xmin = min(xmin),
              xmax = max(xmax)) 
  
  box_r <- box %>% 
    group_by(region) %>% 
    summarise(ymin = min(ymin), 
              ymax = max(ymax),
              xmin = min(xmin),
              xmax = max(xmax)) %>% 
    rename(unit = region)
  
  box <- rbind(box_i, box_r) %>% filter(unit == domain)
  
  # Define the extent from island_boxes
  extent_box <- extent(box$xmin, box$xmax, box$ymin, box$ymax)
  
  # Clip env_rs using the defined extent
  clipped_raster <- crop(env_rs, extent_box)
  
  clipped_raster[["Bathymetry.Min"]][ clipped_raster[["Bathymetry.Min"]] <= -depth_cutoff] <- NA
  clipped_raster = mask(rast(clipped_raster), rast(clipped_raster[["Bathymetry.Min"]]))
  
  return(clipped_raster)
}
# ---- Run MaxEnt Batch ----

run_maxent = function(occ_sf, env, survey) {
  
  set.seed(2024)
# 
# occ_sf = occ_df
# env = eds

  # survey = "ncrmp"
  
  species = unique(occ_sf$Scientific.Name)
  
  cat("Running MaxEnt calibration for species", 
      species, "with", 
      survey, "survey data...\n")
  
  # Loop through each species
  for(sp in species) {
    
    # sp = unique(occ_sf$Scientific.Name)
    
    # Subset occurrence data for the species
    occ_sp = occ_sf[occ_sf$Scientific.Name == sp, c("Longitude", "Latitude")]
    
    # crs(env) = "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
    crs(env) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
    
    # plot(env[["bathymetry"]]); points(occ_sp, col='red', pch = 20) # Plot the first environmental layer
    
    enmeval_results <- tryCatch({
      
      # Try block partitioning first
      ENMevaluate(occ_sp, env,
                  bg = NULL,
                  tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
                  partitions = "block",
                  algorithm = "maxnet",
                  parallel = TRUE,
                  numCores = detectCores() / 2,
                  updateProgress = TRUE,
                  taxon.name = sp)
      
    }, error = function(e) {
      
      cat("Error in block partitioning: ", e$message, "\nTrying random k-fold partitioning...\n")
      closeAllConnections()
      
      # If block partitioning fails, try random k-fold partitioning
      ENMevaluate(occ_sp, env,
                  bg = NULL,
                  tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
                  partitions = "randomkfold",
                  partition.settings = list(kfolds = 5),
                  algorithm = "maxnet",
                  parallel = TRUE,
                  numCores = detectCores() / 2,
                  updateProgress = TRUE,
                  taxon.name = sp)
    })
    
    enmeval_df = enmeval_results@results
    
    readr::write_csv(enmeval_results@results %>%
                       mutate(across(where(is.numeric), ~ round(.x, 3))), 
                     file = paste0("output/maxent_result_enmeval_df_", sp, "_", survey, ".csv"))
    
    # Subset the ENMeval results to get the best model
    enmeval_bestm = subset(enmeval_df, delta.AICc == 0)
    enmeval_bestm = enmeval_bestm[1,]
    
    # Decode the features
    maxent_feats = as.character(enmeval_bestm$fc)
    maxent_rm = as.character(enmeval_bestm$rm)
    
    # Print out the results
    cat(paste("Best features for", sp, ":", maxent_feats, "\n"))
    cat(paste("Best regularization multiplier for", sp, ":", maxent_rm, "\n"))
    
    # Function to map feature classes to MaxEnt arguments
    get_maxent_feature_args <- function(fc) {
      # Initialize all features to 'false'
      features <- list(
        linear = "linear=false",
        quadratic = "quadratic=false",
        product = "product=false",
        threshold = "threshold=false",
        hinge = "hinge=false"
      )
      
      # Map abbreviations to features
      if (grepl("L", fc)) features$linear <- "linear=true"
      if (grepl("Q", fc)) features$quadratic <- "quadratic=true"
      if (grepl("P", fc)) features$product <- "product=true"
      if (grepl("T", fc)) features$threshold <- "threshold=true"
      if (grepl("H", fc)) features$hinge <- "hinge=true"
      
      # Return the feature arguments as a vector
      return(unlist(features))
    }
    
    # Get the MaxEnt feature arguments
    feature_args <- get_maxent_feature_args(maxent_feats)
    
    # Combine all MaxEnt arguments
    maxent_args <- c(
      paste0("betamultiplier=", maxent_rm),
      feature_args,
      "responsecurves=true",
      "jackknife=true"
    )
    
    # Run MaxEnt SDM
    sp_maxent_model <- dismo::maxent(
      x = env,
      p = as.matrix(occ_sp),
      args = maxent_args,
      removeDuplicates = T
    )
    
    # Run MaXent SDM
    sp_maxent_model = dismo::maxent(env,
                                    as.matrix(occ_sp),
                                    features = maxent_feats,
                                    betamultiplier = maxent_rm,
                                    removeDuplicates = T)
    
    # Store maxent model
    maxent_result = list(enm = enmeval_df, model = sp_maxent_model)
    save(maxent_result, file = paste0("output/maxent_result_", sp, "_", survey, ".rda"))
    
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
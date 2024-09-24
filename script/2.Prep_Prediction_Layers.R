library(terra)
library(dplyr)

rm(list = ls())
select = dplyr::select

df <- readRDS("data/eds_grid_for_prediction.rds") %>% filter(unit == "Tutuila")

names(df) <- gsub("Daily", "daily", names(df)); names(df)
names(df) <- gsub("Weekly", "weekly", names(df)); names(df)
names(df) <- gsub("_8Day_", "_8days_", names(df)); names(df)
names(df) <- gsub("Monthly", "monthly", names(df)); names(df)

names(df) <- gsub("DHW.", "dhw.", names(df)); names(df)
names(df) <- gsub("Major_", "major_", names(df)); names(df)
names(df) <- gsub("Np10y_", "np10y_", names(df)); names(df)
names(df) <- gsub("MeanMax_", "mean_max_", names(df)); names(df)
names(df) <- gsub("CI95Max_", "ci95_max_", names(df)); names(df)
names(df) <- gsub("MeanDur_", "mean_durnal_", names(df)); names(df)
names(df) <- gsub("MaxMax_", "max_max_", names(df)); names(df)

names(df) <- gsub("_DY01", "_01dy", names(df)); names(df)
names(df) <- gsub("_WK01", "_01wk", names(df)); names(df)
names(df) <- gsub("_MO01", "_01mo", names(df)); names(df)
names(df) <- gsub("_MO03", "_03mo", names(df)); names(df)
names(df) <- gsub("_MO06", "_06mo", names(df)); names(df)
names(df) <- gsub("_YR10YR01", "_10yr_01yr", names(df)); names(df)
names(df) <- gsub("_YR01", "_01yr", names(df)); names(df)
names(df) <- gsub("_YR03", "_03yr", names(df)); names(df)
names(df) <- gsub("_YR05", "_05yr", names(df)); names(df)
names(df) <- gsub("_YR10", "_10yr", names(df)); names(df)
names(df) <- gsub("_ALLB4", "_all_before", names(df)); names(df)

df <- df %>% select(-contains("kd490")); names(df)
df <- df %>% select(-contains("hi_otp")); names(df)
df <- df %>% select(-contains("mhi")); names(df)
df <- df %>% select(-contains("biweekly_range")); names(df)
df <- df %>% select(-contains("monthly_range")); names(df)
df <- df %>% select(-contains("_dy01")); names(df)
df <- df %>% select(-contains("_wk01")); names(df)
df <- df %>% select(-contains("_mo01")); names(df)
df <- df %>% select(-contains("_mo03")); names(df)
df <- df %>% select(-contains("_mo06")); names(df)
df <- df %>% select(-contains("_yr03")); names(df)
df <- df %>% select(-contains("_yr05")); names(df)
df <- df %>% select(-contains("_yr10")); names(df)
df <- df %>% select(-contains("_allb4")); names(df)
df <- df %>% select(-contains("_yr10yr01")); names(df)
df <- df %>% select(-contains("_15_min")); names(df)
df <- df %>% select(-contains("_30_min")); names(df)
df <- df %>% select(-contains("1_deg.")); names(df)
df <- df %>% select(-contains("7daymax")); names(df)
df <- df %>% select(-contains("jplmur")); names(df)

df <- df %>% select(where(~ all(!is.na(.)))); names(df)
df <- df %>% select(where(~ n_distinct(.) > 1)); names(df)

names(df) <- gsub("_daily", "", names(df)); names(df)
names(df) <- gsub("_monthly", "", names(df)); names(df)

names(df) <- gsub("_crw_", "_", names(df)); names(df)
names(df) <- gsub("_viirs_", "_", names(df)); names(df)
names(df) <- gsub("_noaa_", "_", names(df)); names(df)
names(df) <- gsub("_yr01", "", names(df)); names(df)

visdat::vis_miss(df, warn_large_data = F)

eds <- rast()

for (v in 6:36) {
  
  var_name <- colnames(df)[v]

  p <- df %>%
    dplyr::select(lon, lat, !!sym(var_name)) %>% 
    group_by(lon, lat) %>% 
    summarise(mean_value = mean(!!sym(var_name), na.rm = TRUE)) %>% 
    ungroup()
  
  colnames(p)[3] = var_name
  
  p = p %>% 
    rast(type = "xyz")
  
  eds <- c(eds, p)

}

eds = raster::stack(eds)
plot(eds[[1:12]])
plot(eds[[13:24]])
plot(eds[[25:31]])

save(eds, file = "data/eds.rdata")


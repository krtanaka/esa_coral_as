library(sdmTMB)
library(sdmTMBextra)
library(dplyr)
library(ggplot2)
library(rgdal)
library(colorRamps)
library(patchwork)
library(raster)
library(sf)
library(readr)
library(ggpubr)
library(tidyr)

df = read.csv("I_craterformis_AS.csv")
df = read.csv("A_globiceps_AS.csv")

df = df %>% filter(ISLAND %in% c("Tutuila", "Ofu & Olosega"))
df = df %>% filter(ISLAND %in% c("Tutuila"))
df = df %>% filter(OBS_YEAR != 2020)

zone <- (floor((df$LONGITUDE[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("LONGITUDE", "LATITUDE")]), paste0("+proj=utm +units=km +zone=", zone))))
colnames(xy_utm) = c("X", "Y")
df = cbind(df, xy_utm)
plot(xy_utm, pch = ".", bty = 'n')
rm(xy_utm)

load('/Users/kisei.tanaka/pifsc_efh/data/MHI_islands_shp.RData')
crs(ISL_bounds) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
ISL_this = ISL_bounds[which(ISL_bounds$ISLAND %in% toupper(df$ISLAND)),]
ISL_this_utm = spTransform(ISL_this,CRS(paste0("+proj=utm +units=km +zone=",zone)))
ISL_this_sf = st_transform(st_as_sf(ISL_this), crs = paste0("+proj=utm +units=km +zone=",zone))

rea_spde <- make_mesh(df, c("X", "Y"), cutoff = 0.1) # search
rea_spde_coast = add_barrier_mesh(rea_spde , ISL_this_sf)

plot(rea_spde_coast$mesh, asp = 1, main = ""); axis(1); axis(2)
plot(ISL_this_utm, add = TRUE)
points(rea_spde_coast$loc_xy,col = "green", pch = ".", cex = 5)
bar_i = rea_spde_coast$barrier_triangles
norm_i = rea_spde_coast$normal_triangles
points(rea_spde_coast$spde$mesh$loc[,1], rea_spde_coast$spde$mesh$loc[,2], pch = ".", col = "black")
points(rea_spde_coast$mesh_sf$V1[bar_i], rea_spde_coast$mesh_sf$V2[bar_i], col = "red", pch = 20, cex = 0.5)
points(rea_spde_coast$mesh_sf$V1[norm_i], rea_spde_coast$mesh_sf$V2[norm_i], col = "blue", pch = 20, cex = 0.5)

df %>% 
  ggplot(aes(LONGITUDE, LATITUDE, fill = AdColCount, size = AdColCount)) + 
  geom_point(shape = 21,  alpha = 0.5) + 
  facet_grid(~ANALYSIS_YEAR)

# df = add_utm_columns(df, ll_names = c("LONGITUDE", "LATITUDE"))
df$depth = (df$MAX_DEPTH_M + df$MIN_DEPTH_M)/2

fit <- sdmTMB(
  AdColCount ~  0 + as.factor(OBS_YEAR) + s(depth, k = 5),
  # reml = T,
  # time = "OBS_YEAR",
  data = df, 
  mesh = rea_spde_coast,
  silent = F,
  control = sdmTMBcontrol(newton_loops = 1),
  family = tweedie(link = "log")
)

fit
sanity(fit)
tidy(fit, conf.int = TRUE)
tidy(fit, effects = "ran_par", conf.int = TRUE)

grid = raster("M:/Environmental_Data_Summary/Data_Download/Bathymetry_NGDC_AmericanSamoa_90m/Block_Level_Data/Tutuila_Bathymetry_NGDC_AmericanSamoa_90m.nc")

grid = raster("M:/Environmental_Data_Summary/Data_Download/Bathymetry_NGDC_AmericanSamoa_90m/Block_Level_Data/Ofu_&_Olosega_Bathymetry_NGDC_AmericanSamoa_90m.nc")

grid = raster("M:/Environmental_Data_Summary/Data_Download/Bathymetry_NGDC_AmericanSamoa_90m/Bathymetry_NGDC_AmericanSamoa_90m_all_units.nc")

grid[grid <= -30] <- NA

grid = grid %>% rasterToPoints() %>% as.data.frame()
colnames(grid) = c("LONGITUDE", "LATITUDE", "depth")
# grid = add_utm_columns(grid, ll_names = c("LONGITUDE", "LATITUDE"))
xy_utm = as.data.frame(cbind(utm = project(as.matrix(grid[, c("LONGITUDE", "LATITUDE")]), paste0("+proj=utm +units=km +zone=", zone))))
colnames(xy_utm) = c("X", "Y")
grid = cbind(grid, xy_utm)
plot(xy_utm, pch = ".", bty = 'n')
rm(xy_utm)
grid$depth = grid$depth * -1

grid_year = NULL

for (y in 1:length(unique(df$OBS_YEAR))) {
  
  grid$OBS_YEAR = unique(df$OBS_YEAR)[y]
  
  grid_year = rbind(grid_year, grid)
  
}

# Predict on the fitted data; see ?predict.sdmTMB
p <- predict(fit)

p %>% 
  ggplot(aes(AdColCount, exp(est), fill = AdColCount)) + 
  geom_point(shape = 21, show.legend = F, size = 3, alpha = 0.8) + 
  labs(x = "Observation", y = "Prediction") + 
  scale_fill_gradientn(colors = matlab.like(100)) + 
  coord_fixed(ratio = 1)

# Predict on new data:
p <- predict(fit, newdata = grid_year)
p$est = exp(p$est)
p = p %>% subset(est <= max(df$AdColCount))
head(p)

p %>% 
  ggplot(aes(LONGITUDE, LATITUDE, fill = est)) + 
  geom_raster() + 
  scale_fill_gradientn(colors = matlab.like(100), trans = "sqrt", "Est. AdColCount") + 
  facet_grid(~OBS_YEAR) + 
  ggtitle("I_craterformis_AS") + 
  # coord_fixed() + 
  theme(panel.background = element_rect(fill = "gray10"),
        panel.grid = element_line(color = "gray15"))



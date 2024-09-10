library(sdmTMB)
library(sdmTMBextra)
library(dplyr)
library(ggplot2)
# library(rgdal)
library(colorRamps)
library(patchwork)
library(raster)
library(sf)
library(readr)
library(ggpubr)
library(tidyr)

rm(list = ls())

sp = c("I_craterformis", "A_globiceps")[1]

# Use paste to construct the file path dynamically
df = read.csv(paste0("data/", sp, "_AS.csv"))

df = df %>% filter(ISLAND %in% c("Tutuila", "Ofu & Olosega"))
df = df %>% filter(ISLAND %in% c("Tutuila"))
df = df %>% filter(OBS_YEAR != 2020)
df$lon = df$LONGITUDE
df$lat = df$LATITUDE

zone <- (floor((df$lon[1] + 180)/6) %% 60) + 1
coords <- df[,c("lon", "lat")]
coordinates(coords) <- ~lon + lat
sf_coords <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
sf_coords_utm <- st_transform(sf_coords, crs = sprintf("+proj=utm +zone=%d +datum=WGS84 units=km", zone))
df$X <- st_coordinates(sf_coords_utm)[,1]
df$Y <- st_coordinates(sf_coords_utm)[,2]

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
coords <- df[,c("LONGITUDE", "LATITUDE")]
coordinates(coords) <- ~LONGITUDE + LATITUDE
sf_coords <- st_as_sf(grid, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
sf_coords_utm <- st_transform(sf_coords, crs = sprintf("+proj=utm +zone=%d +datum=WGS84 units=km", zone))
grid$X <- st_coordinates(sf_coords_utm)[,1]
grid$Y <- st_coordinates(sf_coords_utm)[,2]
# grid = add_utm_columns(grid, ll_names = c("LONGITUDE", "LATITUDE"))
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

dfi = p %>% 
  # mutate(lon = round(LONGITUDE, 2),
  # lat = round(LATITUDE, 2)) %>%
  mutate(lon = round(LONGITUDE, 3),
         lat = round(LATITUDE, 3),
         year = OBS_YEAR) %>%
  group_by(lon, lat, year) %>% 
  summarise(est = mean(est, na.rm = T)
            # AdColCount = log(AdColCount + 1)
  )

p1 = ggmap(map) +
  geom_spatial_point(data = dfi %>% filter(year == "2015"), aes(lon, lat, 
                                                  fill = est,
                                                  color = est), 
                     shape = 21, 
                     size = 0.5,
                     alpha = 0.8, crs = 4326) + 
  annotate("text", x = min(df$lon), y = max(df$lat), label = paste0(sp, " 2015"), 
           hjust = 0, vjust = 1, size = 6, color = "white", fontface = "bold") + 
  scale_y_continuous(limits = c(min(df$lat), max(df$lat)), "") +
  scale_x_continuous(limits = c(min(df$lon), max(df$lon)), "") +
  scale_color_gradientn(colours = matlab.like(100), trans = "sqrt", "Est. AdColCount", limits = c(0, 150))  +
  scale_fill_gradientn(colors = matlab.like(100), trans = "sqrt", "Est. AdColCount", limits = c(0, 150)) + 
  theme_minimal() + 
  theme(legend.position = c(0.9, 0.2),
        legend.background = element_blank(),             # Transparent background
        legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
        legend.text = element_text(color = "white", face = "bold"),  # White and bold text
        legend.title = element_text(color = "white", face = "bold"))

p2 = ggmap(map) +
  geom_spatial_point(data = dfi %>% filter(year == "2023"), aes(lon, lat, 
                                                                fill = est,
                                                                color = est), 
                     shape = 21, 
                     size = 0.5,
                     alpha = 0.8, crs = 4326) + 
  annotate("text", x = min(df$lon), y = max(df$lat), label = paste0(sp, " 2023"),  
           hjust = 0, vjust = 1, size = 6, color = "white", fontface = "bold") + 
  scale_y_continuous(limits = c(min(df$lat), max(df$lat)), "") +
  scale_x_continuous(limits = c(min(df$lon), max(df$lon)), "") +
  scale_color_gradientn(colours = matlab.like(100), trans = "sqrt", "Est. AdColCount", limits = c(0, 150))  +
  scale_fill_gradientn(colors = matlab.like(100), trans = "sqrt", "Est. AdColCount", limits = c(0, 150)) + 
  theme_minimal() + 
  theme(legend.position = c(0.9, 0.2),
        legend.background = element_blank(),             # Transparent background
        legend.key = element_rect(colour = NA, fill = NA), # Transparent key background
        legend.text = element_text(color = "white", face = "bold"),  # White and bold text
        legend.title = element_text(color = "white", face = "bold"))

p1 / p2

ggsave(last_plot(), filename =  file.path("output/sdmtmb.png"), height = 10, width = 11)

# p %>% 
#   ggplot(aes(LONGITUDE, LATITUDE, fill = est)) + 
#   geom_raster() + 
#   scale_fill_gradientn(colors = matlab.like(100), trans = "sqrt", "Est. AdColCount") + 
#   facet_grid(~OBS_YEAR) + 
#   ggtitle("I_craterformis_AS") + 
#   # coord_fixed() + 
#   theme(panel.background = element_rect(fill = "gray10"),
#         panel.grid = element_line(color = "gray15"))



library(dismo)
library(raster)
library(ecospat)
library(maps)
library(readr)
library(dplyr)
library(terra)
library(colorRamps)
library(ggplot2)
library(ggmap)

rm(list = ls())

load("data/eds.rdata")

presence <- read.csv("data/original_data/crag_isopora_occurences.csv")

coordinates_df <- presence %>%
  dplyr::select(x = longitude, y = latitude) %>%
  mutate(x = as.numeric(x), y = as.numeric(y))

predictors <- eds

load(paste0("output/maxent_result_Isopora crateriformis_combined.rda"))

prediction <- predict(maxent_result$model, predictors)
prediction <- rast(prediction)

pred_df <- as.data.frame(prediction, xy = TRUE)

ggmap::register_google(key = "AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

map = ggmap::get_map(location = c(mean(pred_df$x), mean(pred_df$y)),
                     maptype = "satellite",
                     zoom = 11,
                     color = "bw",
                     force = TRUE)

p1 = ggmap(map) +
  geom_raster(data = pred_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c() +
  geom_point(data = presence, aes(x = longitude, y = latitude), color = "red", size = 5, shape = 18) +
  coord_sf(xlim = c(-170.85, -170.75), ylim = c(-14.375, -14.3), expand = FALSE) +
  labs(
    fill = "Predicted Habitat Suitability",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    legend.position = c(0.2, 0.1),  # Adjust legend position inside the plot
    legend.direction = "horizontal",  # Horizontal legend
    legend.background = element_blank(),  # Transparent background
    legend.text = element_text(color = "white"),  # White legend text
    legend.title = element_text(color = "white", face = "bold")  # White bold title
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",  # Move title to the top of the legend
    title.hjust = 0.5       # Center the title horizontally
  ))

p2 = ggmap(map) +
  geom_raster(data = pred_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c() +  
  geom_point(data = presence, aes(x = longitude, y = latitude), color = "red", size = 5, shape = 18) +
  coord_sf(xlim = c(-170.66, -170.5536), ylim = c(-14.33, -14.24124), expand = FALSE) +
  labs(
    fill = "Predicted Habitat Suitability",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    legend.position = c(0.8, 0.1),  # Adjust legend position inside the plot
    legend.direction = "horizontal",  # Horizontal legend
    legend.background = element_blank(),  # Transparent background
    legend.text = element_text(color = "white"),  # White legend text
    legend.title = element_text(color = "white", face = "bold")  # White bold title
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",  # Move title to the top of the legend
    title.hjust = 0.5       # Center the title horizontally
  ))

p1 + p2

# Extract predicted values at presence locations
# Ensure 'presence' is a 'SpatVector' with the correct CRS
presence_vect <- vect(presence, geom = c("longitude", "latitude"), crs = crs(prediction))

# Extract predicted values using terra's extract function
pres_values <- terra::extract(prediction, presence_vect)$layer

# Remove NA values
pres_values <- pres_values[!is.na(pres_values)]

# Get all predicted values
all_values <- values(prediction, mat = FALSE)
all_values <- all_values[!is.na(all_values)]

# Summary of predicted values at presence locations
summary(pres_values)

# Summary of all predicted values
summary(all_values)

# Calculate the Boyce Index
breakpoints <- quantile(all_values, probs = seq(0, 1, length.out = 10))

boyce_result <- ecospat.boyce(fit = all_values,
                              obs = pres_values,
                              PEplot = F,
                              # , res = 50
                              # , nclass = 10
                              # , window.w = breakpoints
)

cat("Boyce Index:", boyce_result$cor, "\n")

# Initialize an empty data frame
boyce_results_df <- data.frame(nclass = integer(), Boyce_Index = numeric())

# Define the nclass values
nclass_values <- seq(5, 100, 10)

# Loop through nclass values and calculate the Boyce Index
for (n in nclass_values) {
  boyce_result <- ecospat.boyce(
    fit = all_values,
    obs = pres_values,
    nclass = n,
    PEplot = TRUE
  )
  
  # Store the result in the data frame
  boyce_results_df <- rbind(boyce_results_df, data.frame(nclass = n, Boyce_Index = boyce_result$cor))
}

# View the data frame
print(boyce_results_df)

# The Boyce Index ranges from -1 to +1:
# Values close to +1 indicate predictions consistent with the distribution of presences.
# Values around 0 indicate no difference from a random model.
# Values close to -1 indicate an incorrect model.

# Predicted-to-Expected Ratio ($F.ratio)
# $F.ratio: A vector of ratios representing the frequency of observed presences in each habitat suitability class relative to the expected frequency if presences were randomly distributed.

# Increasing Trend: The $F.ratio generally increases with higher habitat suitability classes ($HS), indicating that higher predicted suitability corresponds to higher observed presences.
# Values > 1: Ratios greater than 1 indicate that the observed frequency of presences is higher than expected under a random distribution in those habitat suitability classes.

# Plot the P/E curve

plot(boyce_result$HS, boyce_result$F.ratio, type = "b", pch = 19, col = "blue",
     xlab = "Habitat Suitability Class Midpoints",
     ylab = "Predicted-to-Expected Ratio",
     main = "Boyce Index P/E Curve")

# Add a horizontal line at y = 1
abline(h = 1, lty = 2, col = "red")

data.frame(HS = boyce_result$HS, 
           Ratio = boyce_result$F.ratio) %>% 
  ggplot(aes(HS, Ratio, fill = Ratio)) + 
  geom_line() + 
  geom_point(shape = 21, size = 5, show.legend = F) + 
  scale_fill_viridis_c("") + 
  labs(x = "Habitat Suitability Class Midpoints",
       y = "Predicted-to-Expected Ratio") +
  ggtitle("Boyce Index P/E Curve") + 
  theme_classic()

# P/E Ratio > 1: Indicates that presences are more frequent than expected in those habitat suitability classes.
# Trend: An increasing P/E ratio with higher suitability classes confirms that higher suitability predictions correspond to higher observed presences.

# Output the Boyce Index
cat("Boyce Index for Hawaii:", boyce_result$cor, "\n")

all_df <- data.frame(
  Suitability = all_values,
  Data = "All Predicted Values"
)

pres_df <- data.frame(
  Suitability = pres_values,
  Data = "Presence Locations"
)

p1 = ggplot() +
  geom_histogram(
    data  = all_df,
    aes(x = Suitability, fill = Data),
    bins  = 50,
    alpha = 0.6,
    color = "black"
  ) +
  geom_dotplot(
    data       = pres_df,
    aes(x      = Suitability, fill = Data),
    # binwidth   = 0.02,    # adjust if needed
    stackgroups= TRUE,
    stackratio = 1.2,
    # dotsize    = 2,
    alpha = 0.8,
    shape      = 21        # shape 8 is a star
  ) +
  scale_fill_manual(
    name   = "",
    values = c(
      "All Predicted Values" = "#1f77b4",
      "Presence Locations"   = "#ff7f0e"
    )
  ) +
  labs(
    x = "Normalized Suitability",
    y = "Frequency"
  ) +
theme_pubr() + 
  theme(legend.position = c(0.7, 0.8))

boyce_results_df <- boyce_results_df %>%
  mutate(mean_boyce = mean(Boyce_Index, na.rm = TRUE),
         species = "I.crateriformis")

p2 = ggplot(boyce_results_df, aes(x = species, y = Boyce_Index)) +
  geom_boxplot(aes(fill = mean_boyce), show.legend = F) +
  geom_jitter(aes(fill = Boyce_Index), shape = 21, size = 5, alpha = 0.5, width = 0.05, height = 0, show.legend = F) +  
  labs(
    x = "",
    y = "Boyce Index",
    fill = "Mean Boyce Index"
  ) +
  # scale_y_continuous(limits = c(-1, 1)) +
  theme_pubr() + 
  theme(
    axis.text.x = element_blank(),
  )

p1 + p2
ggsave(last_plot(), file = "output/Boyce_Index.png", height = 5, width = 8)

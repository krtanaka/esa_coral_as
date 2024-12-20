library(dismo)  # For working with MaxEnt models
library(pROC)   # For ROC curves and AUC calculations
library(PresenceAbsence)  # For confusion matrix and related metrics
library(ggplot2)

rm(list = ls())

species <- c("Acropora globiceps", "Isopora crateriformis", "Genus Tridacna")[1]
load(paste0("output/maxent_raster_", species, "_ncrmp.rdata"))
# load(paste0("output/maxent_raster_", species, "_no_nps.rdata"))
# load(paste0("output/maxent_raster_", species, "_combined.rdata"))

predicted_suitability = rasterFromXYZ(r)

species_csv <- list(
  "Acropora globiceps" = "A_globiceps_AS.csv",
  "Isopora crateriformis" = "I_craterformis_AS.csv"
)

# use NPS data because CRAG and gbif are occurrence only
load("data/npsa_benthic_data.rdata")

validation_data = df %>% 
  filter(#Loc_Type == "Fixed",
    Latitude >= -14.38, Latitude <= -14.22,
    Longitude >= -170.85, Longitude <= -170.53) %>%
  mutate(y = as.integer(Taxon_Name == species)) %>%
  group_by(Longitude, Latitude) %>%
  summarize(presence = if_else(any(y > 0, na.rm = TRUE), 1, 0), .groups = "drop") %>%
  distinct() %>%
  na.omit() %>%
  rename_with(tolower)

table(validation_data$presence)

# Extract predicted suitability values at the locations of validation points
coordinates <- cbind(validation_data$longitude, validation_data$latitude)
validation_data$predicted_suitability <- raster::extract(predicted_suitability, coordinates)

ggmap::register_google(key = "AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

# Get map for the given coordinates
map = ggmap::get_map(location = c(mean(validation_data$longitude), mean(validation_data$latitude)),
                     maptype = "satellite",
                     zoom = 13,
                     color = "bw",
                     force = TRUE)

p1 = ggmap(map)+ 
  geom_point(data = validation_data, aes(longitude, latitude, fill = predicted_suitability, color = predicted_suitability), shape = 21, size = 3, alpha = 0.5) + 
  geom_point(data = validation_data %>% filter(presence == 1), aes(longitude, latitude), fill = "red", shape = 21, size = 3) + 
  scale_y_continuous(limits =  range(validation_data$latitude), name = NULL) +
  scale_x_continuous(limits =  range(validation_data$longitude), name = NULL) +
  scale_fill_viridis_c("Predicted Habitat Suitability", limits = c(0,1),
                       breaks = c(0, 0.5, 1), guide = guide_colorbar(direction = "horizontal",
                                                                     title.position = "top",
                                                                     barwidth = 10, barheight = 0.8)) +
  scale_color_viridis_c("Predicted Habitat Suitability", limits = c(0,1),
                        breaks = c(0, 0.5, 1), guide = guide_colorbar(direction = "horizontal",
                                                                      title.position = "top",
                                                                      barwidth = 10, barheight = 0.8)) +
  labs(title = species) +
  theme(legend.position = c(0.22, 0.9),
        legend.background = element_blank(), 
        legend.box.background = element_blank(), 
        legend.text = element_text(color = "white", size = 10, face = "bold"), 
        legend.title = element_text(color = "white", face = "bold")) + 
  coord_sf(crs = 4326)

ggsave(plot = p1,
       filename =  file.path(paste0("output/maxent_map_", species, "_suitability_extracted.png")), 
       width = 6, 
       height = 4.5,
       limitsize = FALSE,
       bg = "transparent")

# Calculate ROC curve and AUC
roc_obj <- roc(validation_data$presence, validation_data$predicted_suitability)

# Print the ROC object to check its contents
print(roc_obj)

# Check levels and direction of ROC object
roc_obj$levels
roc_obj$direction

# Explicitly set levels
roc_obj <- roc(validation_data$presence, validation_data$predicted_suitability, levels = c(0, 1), direction = "<")

# Manually calculate AUC as a check (though roc() should work)
auc_manual <- sum((roc_obj$sensitivities[-1] + roc_obj$sensitivities[-length(roc_obj$sensitivities)]) * 
                    diff(roc_obj$specificities) / 2)

print(auc_manual)

auc_value <- auc_manual

# Generate ROC object
roc_obj <- roc(validation_data$presence, validation_data$predicted_suitability)

# Convert ROC object to a data frame for ggplot
roc_data <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)

# Plot ROC curve with ggplot

ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line() +  # ROC curve line
  geom_abline(linetype = "dashed") +  # Diagonal line for random performance
  labs(
    title = paste("ROC Curve (AUC =", round(auc_manual, 2), ")"),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) + coord_fixed()

ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(size = 1) +  # ROC curve line
  geom_abline(linetype = "dashed") +  # Diagonal line for random performance
  labs(
    title = paste(species, "\nROC Curve (AUC =", round(auc_manual, 2), ")"),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  # ggdark::dark_theme_minimal(base_size = 15) + 
  theme_classic(base_size = 15) + 
  coord_fixed()

# Save the plot with a transparent background
ggsave(last_plot(), file = paste0("output/roc_", species, ".png"), height = 5, bg = "transparent")

# Define a threshold (e.g., 0.5) to classify presence/absence
threshold <- 0.5

# Classify presence/absence based on predicted suitability
validation_data$predicted_presence <- ifelse(validation_data$predicted_suitability >= threshold, 1, 0)

# Convert both columns to factors with the same levels (0 and 1)
validation_data$predicted_presence <- factor(validation_data$predicted_presence, levels = c(0, 1))
validation_data$presence <- factor(validation_data$presence, levels = c(0, 1))

# Load caret package (if not already loaded)
library(caret)

# Now run the confusionMatrix function
confusion_matrix <- confusionMatrix(validation_data$predicted_presence, validation_data$presence)

# Display the confusion matrix and associated metrics (accuracy, sensitivity, specificity)
print(confusion_matrix)

# Create a confusion matrix as a data frame
confusion_matrix_df = confusion_matrix$table %>% as.data.frame()

confusion_matrix_df <- confusion_matrix_df %>%
  mutate(Freq = Freq / sum(Freq),
         Freq = round(Freq, 2))

# Plot the confusion matrix using ggplot
ggplot(confusion_matrix_df, aes(x = as.factor(Prediction), y = as.factor(Reference), fill = Freq)) +
  geom_tile(color = "white", show.legend = FALSE) +  # Create the heatmap tiles
  geom_text(
    aes(label = Freq), 
    size = 10, 
    color = "white", 
    stroke = 0.5, 
    fontface = "bold", 
    fill = "black"
  ) +  # Add text with white outline and black fill
  labs(
    title = species,
    x = "Prediction",
    y = "Reference"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 15) + 
  theme(panel.grid.major = element_blank())

ggsave(last_plot(), file = paste0("output/confusion_", species, ".png"), height = 5, bg = "transparent")

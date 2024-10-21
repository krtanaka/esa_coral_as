library(dismo)  # For working with MaxEnt models
library(pROC)   # For ROC curves and AUC calculations
library(PresenceAbsence)  # For confusion matrix and related metrics
library(ggplot2)

rm(list = ls())

species <- c("Acropora globiceps", "Isopora crateriformis", "Genus Tridacna")[2]
load(paste0("output/maxent_raster_", species, ".rdata"))
predicted_suitability <- r

species_csv <- list(
  "Acropora globiceps" = "A_globiceps_AS.csv",
  "Isopora crateriformis" = "I_craterformis_AS.csv"
)

# Dynamically load the appropriate CSV based on the selected species
validation_data <- read_csv(paste0("data/original_data/", species_csv[[species]])) %>%
  filter(ISLAND == "Tutuila") %>% 
  mutate(longitude = LONGITUDE, latitude = LATITUDE) %>%
  group_by(longitude, latitude) %>%
  summarise(presence = ifelse(AdColCount > 0, 1, 0)) %>% 
  na.omit()

# Extract predicted suitability values at the locations of validation points
coordinates <- cbind(validation_data$longitude, validation_data$latitude)
validation_data$predicted_suitability <- raster::extract(predicted_suitability, coordinates)

ggplot()+ 
  geom_point(data = validation_data, aes(longitude, latitude, fill = predicted_suitability), shape = 21) + 
  geom_point(data = validation_data %>% filter(presence == 1), aes(longitude, latitude, fill = presence), shape = 22)

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
  ggdark::dark_theme_minimal(base_size = 15) + 
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

# Plot the confusion matrix using ggplot
ggplot(confusion_matrix_df, aes(x = as.factor(Prediction), y = as.factor(Reference), fill = Freq)) +
  geom_tile(color = "white", show.legend = F) +  # Create the heatmap tiles
  geom_text(aes(label = Freq), size = 10) +  # Add the counts as text on the tiles
  scale_fill_gradientn(colors = matlab.like(100)) +  # Color gradient for the counts
  labs(
    title = species,
    x = "Prediction",
    y = "Reference"
  ) +
  coord_fixed() + 
  ggdark::dark_theme_minimal(base_size = 15)

ggsave(last_plot(), file = paste0("output/confusion_", species, ".png"), height = 5, bg = "transparent")

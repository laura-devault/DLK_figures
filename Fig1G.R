library(here)
library(tidyverse)
library(dplyr)
library(ggplot2)

# Read summary_0915.csv
data_0915 <- read.csv("summary_0915.csv")

# Read summary_0927.csv
data_0927 <- read.csv("summary_0927.csv")

# Read summary_0929.csv
data_0929 <- read.csv("summary_0929.csv")

# Extract the last four characters from the file name and convert to a date
data_0915$date <- as.Date(substr("summary_0915.csv", nchar("summary_0915.csv") - 3, nchar("summary_0915.csv")), format = "%m%d")
data_0927$date <- as.Date(substr("summary_0927.csv", nchar("summary_0927.csv") - 3, nchar("summary_0927.csv")), format = "%m%d")
data_0929$date <- as.Date(substr("summary_0929.csv", nchar("summary_0929.csv") - 3, nchar("summary_0929.csv")), format = "%m%d")
library(dplyr)

# For data_0915
normfactor_0915 <- data_0915 %>%
  filter(str_detect(Slice, "A")) %>%
  summarise(normfactor = mean(X.Area, na.rm = TRUE)) %>%
  pull(normfactor)

data_0915 <- data_0915 %>%
  mutate(normalized = X.Area / normfactor_0915)

# For data_0927
normfactor_0927 <- data_0927 %>%
  filter(str_detect(Slice, "A")) %>%
  summarise(normfactor = mean(X.Area, na.rm = TRUE)) %>%
  pull(normfactor)

data_0927 <- data_0927 %>%
  mutate(normalized = X.Area / normfactor_0927)

# For data_0929
normfactor_0929 <- data_0929 %>%
  filter(str_detect(Slice, "A")) %>%
  summarise(normfactor = mean(X.Area, na.rm = TRUE)) %>%
  pull(normfactor)

data_0929 <- data_0929 %>%
  mutate(normalized = X.Area / normfactor_0929)

# Combine data frames into a single table
combined_data <- bind_rows(
  data_0915 %>% mutate(dataset = "data_0915"),
  data_0927 %>% mutate(dataset = "data_0927"),
  data_0929 %>% mutate(dataset = "data_0929")
)

# Reset row names
rownames(combined_data) <- NULL

# Replace values in the "Slice" column with A, B, C, or D
combined_data$Slice <- gsub("[^ABCD]", "", combined_data$Slice)

combined_data <- combined_data %>%
  mutate(pretreatment = ifelse(Slice %in% c("A", "B"), "dmso", "dlki"))

combined_data <- combined_data %>%
  mutate(treatment = ifelse(Slice %in% c("A", "C"), "dmso", "noc"))

combined_data <- combined_data %>%
  filter(combined_data$Total.Area >=75)

# Convert treatment labels to uppercase
combined_data$treatment <- toupper(combined_data$treatment)

# Reorder the levels of the "pretreatment" column
combined_data$pretreatment <- factor(combined_data$pretreatment, levels = c("dmso", "dlki"))

# Create the plot with updated legend label
plot <- ggplot(combined_data, aes(x = pretreatment, y = normalized, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_point(shape = 16, position = position_dodge(width = 0.8), alpha = 0.6) +
  labs(x = "Pretreatment", y = "Normalized Fragmentation") +
  scale_fill_discrete(name = "Treatment") +  
  scale_x_discrete(labels = c("dmso" = "DMSO", "dlki" = "DLKi")) +  # Update X-axis labels# Change the legend label
  theme_bw() +
  ggtitle("Microtubule Integrity") +
  theme(legend.position = "bottom", text = element_text(size = 24))
# Display the updated plot
print(plot)

# Fit the two-way ANOVA model
model <- aov(normalized ~ pretreatment * treatment, data = combined_data)

# Perform the ANOVA analysis
anova_result <- summary(model)
# Print the ANOVA table
print(anova_result)

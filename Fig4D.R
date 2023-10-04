library(here)
library(tidyverse)
library(readxl)

GC_area <- read_excel(here("GC_area.xlsx"))

GC_area$Condition <- toupper(GC_area$Condition)


# Mutate the data frame to create the "Pretreatment" column
GC_area <- GC_area %>%
  mutate(Pretreatment = ifelse(Condition %in% c('DMSO', 'NOC'), "DMSO", ifelse(Condition %in% c('DLKI', 'DLKINOC'), "DLKi", NA)))

# Mutate the data frame to create the "Treatment" column
GC_area <- GC_area %>%
  mutate(Treatment = ifelse(Condition %in% c('DMSO', 'DLKI'), "DMSO", ifelse(Condition %in% c('NOC', 'DLKINOC'), "NOC", NA)))


# Reorder Pretreatment column as a factor with custom levels
GC_area$Pretreatment <- factor(GC_area$Pretreatment, levels = c("DMSO", "DLKi"))

# Create a box plot with separate groups based on Pretreatment and Treatment
plot <- ggplot(GC_area, aes(x = Pretreatment, y = `%Area`, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_point(shape = 16, position = position_dodge(width = 0.8), alpha = 0.6) +
  labs(x = "Pretreatment", y = "Normalized Area") +
  theme_bw() +
  ggtitle("Growth Cone Area") +
  theme(legend.position = "bottom", text = element_text(size = 24))


print(plot)
# Perform two-way ANOVA on the data
anova_result <- aov(`%Area` ~ Treatment * Pretreatment, data = GC_area)

# Print the ANOVA table
summary(anova_result)

# Extract average %Area with confidence intervals for each item in the Condition column
means <- aggregate(GC_area$`%Area`, by = list(Condition = GC_area$Condition), FUN = mean)
ci <- tapply(GC_area$`%Area`, GC_area$Condition, FUN = function(x) t.test(x)$conf.int)
average_area <- cbind(means, ci)

library(here)
library(tidyverse)

mydata <- read_csv(here("quant_western/MKK4.csv"))


ggplot(mydata, aes(x = factor(Pretreatment, levels = c("DMSO", "DLKi")), y = mydata$`pMKK4/MKK4`, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_point(shape = 16, position = position_dodge(width = 0.8), alpha = 0.6) +
  labs(x = "Pretreatment", y = "pMKK4/MKK4", fill = "Treatment") +
  theme_bw()+
  theme(legend.position = "bottom", text = element_text(size = 24))

# Perform two-way ANOVA
model <- aov(`pMKK4/MKK4` ~ Pretreatment * Treatment, data = mydata)

# Print the ANOVA table
summary(model)

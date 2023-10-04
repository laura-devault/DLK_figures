library(readxl)
library(ggplot2)


filopodia <- read_excel("Filopodia.xlsx")


# Create a scatter plot using ggplot()
ggplot(filopodia, aes(x = factor(Pretreatment, levels = c("DMSO", "DLKi")), y = `filopodia/length`, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_point(shape = 16, position = position_dodge(width = 0.8), alpha = 0.6) +
  labs(x = "Pretreatment", y = "Filopodia per 1000um length") +
  ggtitle("Filopodia length by Treatment and Pretreatment") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 24))  # Increase label size to 24
#ANOVA

# Perform two-way ANOVA
model <- aov(`filopodia/length` ~ Pretreatment * Treatment, data = filopodia)

# Create the 'Condition' column
filopodia$Condition <- paste(filopodia$Pretreatment, filopodia$Treatment, sep = " + ")


# Print the ANOVA table
summary(model)

# Extract average %Area with confidence intervals for each item in the Condition column
means <- aggregate(`filopodia/length` ~ Condition, data = filopodia, FUN = mean)

ci <- tapply(filopodia$`filopodia/length`, filopodia$Condition, 
             FUN = function(x) t.test(x)$conf.int)

average_area <- cbind(means, ci)

# Print the result
print(average_area)



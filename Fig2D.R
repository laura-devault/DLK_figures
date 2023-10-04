library(here)
library(tidyverse)

MF <- read.csv("webgestalt/cell_weighted.csv", header = TRUE)

MF$NES <- -MF$NES

# Create a new column "color" based on the sign of NES
MF$color <- ifelse(MF$NES >= 0, "Positive", "Negative")

# Order the data frame by NES in descending order
MF <- MF %>% arrange(desc(NES))

# Create the horizontal bar chart
ggplot(data = MF, aes(x = NES, y = fct_reorder(Description, NES), fill = color)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Downregulated", "Upregulated")) +
  labs(x = "Normalized Enrichment Score", y = "Gene Set") +
  theme_classic() +
  theme(legend.position = "bottom",
        text = element_text(size = 24))  # Set the font size here (e.g., size = 12)

title ="GO_enrichment"

#produce expression density plots for DLK-dependent transcripional program

library(tidyverse)
library(Mus.musculus)
library(ggrepel)
library(here)
library(dplyr)
library(ggplot2)

source(here('R/map_go_to_gene.R'))

noc <- read.csv(here("data/dlki_noc_diantonio_results/effect_of_dlki_noc.csv"))

# Flip the sign of log2FoldChange column
noc$log2FoldChange <- -noc$log2FoldChange

go_terms =  c("GO:0007010", "GO:0030036")

term_list = go_to_gene(go_terms)

#add column with GO term

noc <- noc %>%
  mutate(cytoskeleton_organization = gene %in% term_list$`cytoskeleton organization`,
         actin_cytoskeleton_organization = gene %in% term_list$`actin cytoskeleton organization`) %>%
  mutate(pvalue_adjusted = pmin(-log10(pvalue), 30),
         gene_expression_change = pmax(-5, log2FoldChange))

# Create a copy of the dataframe to avoid modifying the original data
noc_enriched <- noc

# Update the 'GO_enrichment' column based on the conditions
noc_enriched$GO_enrichment <- ifelse(noc_enriched$actin_cytoskeleton_organization == TRUE, "actin",
                                     ifelse(noc_enriched$cytoskeleton_organization == TRUE, "cytoskeleton", FALSE))
# Update the 'sig' column based on the conditions
noc_enriched$sig <- ifelse(abs(noc_enriched$gene_expression_change) > 0.5, noc_enriched$GO_enrichment, FALSE)

# Filter the data and select the top 10 genes with the greatest gene_expression_change values
top_genes <- noc_enriched %>%
  filter(sig != "FALSE") %>%
  top_n(12, abs(gene_expression_change))


noc_plot <- ggplot() +
  geom_point(data = subset(noc_enriched, sig == "FALSE"), aes(x = gene_expression_change, y = pvalue_adjusted), color = "grey", size = 2) +
  geom_point(data = subset(noc_enriched, sig %in% c("cytoskeleton", "actin")), aes(x = gene_expression_change, y = pvalue_adjusted, color = sig), size = 2) +
  geom_label_repel(data = top_genes, aes(x = gene_expression_change, y = pvalue_adjusted, label = gene),
                   size = 9, color = "black", fill = "white",
                   box.padding = 0.4, box.colour = "black",
                   segment.color = "black", max.overlaps = Inf, nudge_y = -0.3) +
  scale_color_manual(values = c("cytoskeleton" = "blue", "actin" = "red"), labels = c("cytoskeleton" = "cytoskeleton organization", "actin" = "actin cytoskeleton organization")) +
  labs(title = "GO Enrichment in DLKi in Noc", x = "log2(Fold Change)", y = "-log10(p-value)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16), color = c("red", "blue")), title = "Terms", labels = c("Cytoskeleton", "Actin")))

noc_plot

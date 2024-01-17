library(readxl)
library(tidyverse)
library(ggrepel)
library(here)
library(Mus.musculus)

# Set the working directory to the "noc_treatment" folder
setwd(here("noc_treatment"))

source(here('R/map_go_to_gene.R'))


# Read in the data from the xlsx files
effect_of_noc <- read_csv("data/dlki_noc_diantonio_results/effect_of_noc.csv")
effect_of_dlki <- read_csv("data/dlki_noc_diantonio_results/effect_of_dlki.csv")
effect_of_dlki_noc <- read_csv("data/dlki_noc_diantonio_results/effect_of_dlki_noc.csv")

# Flip the sign of log2FoldChange column
effect_of_dlki_noc$log2FoldChange <- -effect_of_dlki_noc$log2FoldChange


MAPK =  go_to_gene("GO:0000165")

MAPK_genes <- effect_of_dlki_noc$gene[abs(effect_of_dlki_noc$log2FoldChange) > 1 & effect_of_dlki_noc$gene %in% MAPK$`MAPK cascade`]# Print the resulting list of selected genes


# Define the genes to label
genes_to_label <- c("Dusp1", "Dusp10", "Dusp14", "Dusp16", "Dusp19", "Dusp3", "Dusp4", "Dusp5", "Dusp6", "Dusp8")

# Define the TFs to label
tfs_to_label <- c("Jun", "Atf3", "Egr1")

gene_list <- c("Atf3", "Jun", "Egr1", "Dusp1")

# Assuming 'effect_of_dlki_noc' is your dataset
# Create a volcano plot using ggplot2 with point coloring for TFs and genes_to_label
volcano_plot <- ggplot(effect_of_dlki_noc, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = filter(effect_of_dlki_noc, gene %in% MAPK_genes),
             aes(color = factor("MAPK_genes", levels = c("Others", "MAPK_genes", "genes_to_label", "TFs"))), size = 4, alpha = 1) +
  geom_point(data = filter(effect_of_dlki_noc, !(gene %in% c(genes_to_label, tfs_to_label, MAPK_genes))),
             aes(color = factor("Others", levels = c("Others", "MAPK_genes", "genes_to_label", "TFs"))), size = 2, alpha = 0.5, color = "grey50") +
  geom_point(data = filter(effect_of_dlki_noc, gene %in% tfs_to_label),
             aes(color = factor("TFs", levels = c("Others", "MAPK_genes", "genes_to_label", "TFs"))), size = 4, alpha = 1) +
  geom_point(data = filter(effect_of_dlki_noc, gene %in% genes_to_label),
             aes(color = factor("genes_to_label", levels = c("Others", "MAPK_genes", "genes_to_label", "TFs"))), size = 4, alpha = 1) +
  labs(x = "log2 Fold Change", y = "-log10(padj)",
       title = "DLK-dependent effects of nocodazole") +
  theme_bw() +
  scale_x_continuous(limits = c(-5, 5), oob = scales::squish) +
  scale_y_continuous(limits = c(0, 50), oob = scales::squish) +
  scale_color_manual(values = c("Others" = "grey50", "MAPK_genes" = "orange", "genes_to_label" = "blue", "TFs" = "red"),
                     labels = c("genes_to_label" = "DUSPs", "TFs" = "top Transcription Factors", "MAPK_genes" = "MAPK cascade", "Others" = "Others"),
                     guide = guide_legend(title = "Legend")) +
  theme(legend.text = element_text(size = 12))

# Add labels for gene_list using geom_label_repel
volcano_plot +
  guides(color = guide_legend(override.aes = list(size = 5)))+
  geom_label_repel(data = filter(effect_of_dlki_noc, gene %in% gene_list),
                   aes(x = log2FoldChange, y = -log10(padj), label = gene),
                   size = 8, color = "black", fill = "white",
                   box.padding = 0.4,
                   segment.color = "black", max.overlaps = Inf, nudge_y = -0.3)


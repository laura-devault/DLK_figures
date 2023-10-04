library(tidyverse)
library(Mus.musculus)
library(ggrepel)
library(here)
library(ggplot2)
library(patchwork)
library(dplyr)

source(here("R/volcano_plot_source_facet.R"))

# Load data frames
df_list <- list(
  diantonio_dlki_noc = read_csv(here("data/dlki_noc_diantonio_results/effect_of_dlki_noc.csv")),
  watkins_dlki = read_csv(here("data/watkins/dlki_vs_dmso_ngf_minus_unfltr_watkins_20220814.csv"))
)

# Iterate through each dataframe in df_list and flip the sign of log2FoldChange column
df_list <- lapply(df_list, function(df) {
  df$log2FoldChange <- -df$log2FoldChange
  return(df)
})


# Define GO terms of interest
go_terms <- c("GO:0030029")

# Retrieve gene symbols for GO terms
term_list <- mapIds(Mus.musculus,
                    keys = go_terms,
                    keytype = "GOALL",
                    column = "SYMBOL",
                    multiVals = "list"
) %>%
  map(unique)

# Map GOIDs to names
term_list_to_name <- mapIds(Mus.musculus,
                            keys = go_terms,
                            keytype = "GOID",
                            column = "TERM"
)

# Assign names to term_list
names(term_list) <- term_list_to_name

# Combine data frames and assign source labels
data_df <- df_list %>%
  bind_rows(.id = "source") %>%
  mutate(source = factor(source,
                         levels = c(
                           "diantonio_dlki_noc",
                           "watkins_dlki"
                         ),
                         labels = c(
                           "Noc",
                           "NGF"
                         )
  ))

# Iterate over each item in term_list and add corresponding columns to data_df
for (term_name in names(term_list)) {
  data_df <- data_df %>%
    mutate(!!term_name := gene %in% term_list[[term_name]])
}

# Change the column names
colnames(data_df)[colnames(data_df) == "actin filament-based process"] <- "AF"


# Set limits to collapse the -log10(padj) and log2Foldchange
y_limit <- 30
x_limit <- 4

# Filter data_df for source "Noc"
data_noc <- subset(data_df, source == "Noc")
data_ngf <- subset(data_df, source == "NGF")

# Adjust points beyond y-axis limit to 30
data_noc <- within(data_noc, pvalue_adjusted <- pmin(-log10(pvalue), 30))
data_ngf <- within(data_ngf, pvalue_adjusted <- pmin(-log10(pvalue), 30))

##label actin genes that overlap between the data sets.
# Filter data_df based on log2FoldChange and AF
filtered_df <- subset(data_df, abs(log2FoldChange) >= 1 & AF == TRUE) %>%
  group_by(gene) %>%
  filter(n() > 1)

# Create a logical function to check if an item appears in filtered_df
is_both_sig <- function(source, gene) {
  any(filtered_df$source == source & filtered_df$gene == gene)
}

# Add the "both_sig" column to data_ngf
data_ngf <- data_ngf %>%
  rowwise() %>%
  mutate(both_sig = is_both_sig(source, gene)) %>%
  ungroup()

# Add the "both_sig" column to data_noc
data_noc <- data_noc %>%
  rowwise() %>%
  mutate(both_sig = is_both_sig(source, gene)) %>%
  ungroup()

#
# Create the volcano plot for Source "Noc" with red points for genes with TRUE values in column "AF" and blue points (50% transparent) for genes with TRUE values in column "actin_org"
plot_noc <- data_noc %>%
  mutate(highlight = ifelse(AF, 'Actin filament-based process', 'Background')) %>%
  mutate(highlight = factor(highlight, levels = c('Background', 'Actin filament-based process'))) %>%
  ggplot(aes(x = log2FoldChange,
             y = pvalue_adjusted,
             alpha = highlight,
             size = highlight,
             shape = highlight,
             color = highlight)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') +
  geom_vline(xintercept = 0.5, linetype='dashed') +
  geom_vline(xintercept = -0.5, linetype='dashed') +
  scale_color_manual(values = c('Background' = 'grey',
                                'Actin filament-based process' = 'red')) +
  scale_alpha_manual(values = c('Background' = 0.4,
                                'Actin filament-based process' = 1)) +
  scale_shape_manual(values = c('Background' = 1,
                                'Actin filament-based process' = 16)) +
  scale_size_manual(values = c('Background' = 1,
                               'Actin filament-based process' = 2)) +
  geom_label_repel(data = ~filter(., both_sig == TRUE),
                   aes(label = gene),
                   size = 5,
                   color = "black",
                   fill = "white",
                   box.padding = 0.4,
                   box.color = "black",
                   segment.color = "black") +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(-3,3), ylim = c(0,30))+
  labs(title = "Low Dose Noc", x = "log2(Fold Change)", y = "-log10(p-value)")+
  guides(alpha = "none",color = guide_legend(override.aes = list(alpha = 1)))


plot_ngf <- data_ngf %>%
  mutate(highlight = ifelse(AF, 'Actin filament-based process', 'Background')) %>%
  mutate(highlight = factor(highlight, levels = c('Background', 'Actin filament-based process'))) %>%
  ggplot(aes(x = log2FoldChange,
             y = pvalue_adjusted,
             alpha = highlight,
             size = highlight,
             shape = highlight,
             color = highlight)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') +
  geom_vline(xintercept = 0.5, linetype='dashed') +
  geom_vline(xintercept = -0.5, linetype='dashed') +
  scale_color_manual(values = c('Actin filament-based process' = 'red')) +  # Remove 'Background' from color scale
  scale_alpha_manual(values = c('Background' = 0.4, 'Actin filament-based process' = 1), guide = "none") +
  scale_shape_manual(values = c('Background' = 1, 'Actin filament-based process' = 16), guide = "none") +
  scale_size_manual(values = c('Background' = 1, 'Actin filament-based process' = 2), guide = "none") +
  geom_label_repel(data = . %>% filter(both_sig == TRUE),
                   aes(label = gene),
                   size=5,
                   color = "black",
                   fill = "white",
                   box.padding = 0.4,
                   box.color = "black",
                   segment.color = "black") +
  theme_bw() +
  coord_cartesian(xlim = c(-3,3), ylim = c(0,30)) +
  labs(title = "NGF Deprivation", x = "log2(Fold Change)", y = "-log10(p-value)") +
  guides(color = guide_legend(title = "Legend",
                              override.aes = list(alpha = 1),
                              keywidth = 0.7,
                              keyheight = 0.7,
                              title.theme = element_text(
                                size = 12,
                                face='bold')),
         size = FALSE,  # Remove size legend
         shape = FALSE)  # Remove shape legend

# Combine the plots side by side
combined_plot <- plot_noc + plot_ngf +
  plot_layout(ncol = 2,
              guides = 'collect')

# Display the combined plots with labeled genes in white boxes and modified legend
combined_plot

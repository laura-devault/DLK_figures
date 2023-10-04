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
                         ))) %>%
  group_by(gene) %>%
  mutate(both_sig_lesser_thres = ifelse(all(abs(log2FoldChange) > 0.5)
                                        & all(padj < 0.05)
                                        & n() > 1,
                                        TRUE, FALSE),
         both_sig_greater_thres = ifelse(all(abs(log2FoldChange) > 2)
                                         & all(padj < 0.05)
                                         & n() > 1,
                                         TRUE, FALSE),
         pvalue_adjusted = pmin(-log10(pvalue), 30),
         gene_expression_change = pmax(-5, log2FoldChange))


# Zscan2 is a gene where one of two sources is significant, and the other is
# not. In that case, both sould be False
stopifnot(!all(data_df %>%
                 filter(gene == 'Zscan2') %>%
                 pull(both_sig_lesser_thres)))
# we expect 13 gene pairs to be sig at the higher thres
stopifnot(sum(data_df$both_sig_greater_thres) == 26)

data_df_split = data_df %>%
  ungroup() %>%
  group_by(source) %>%
  group_split() %>%
  map(droplevels)

names(data_df_split) = levels(unlist(map(data_df_split,
                                        ~unique(pull(.,source)))))

# NOTE!!! THIS IS HARD CODED -- MUST CHECK THE ORDER OF THE NAMES
# BEFORE RUNNING THIS
names(data_df_split) = c('Low Dose Noc', 'NGF Deprivation')

# Set limits to collapse the -log10(padj) and log2Foldchange
y_limit <- 30
x_limit <- 4

plot_split = function(df, plot_title, remove_legend = FALSE){
  stopifnot(is.data.frame(df))
  # Create the volcano plot for Source "Noc" with red points for genes with TRUE values in column "AF" and blue points (50% transparent) for genes with TRUE values in column "actin_org"
  df %>%
    mutate(highlight = ifelse(both_sig_lesser_thres == TRUE, 'shared significant gene', 'Background')) %>%
    mutate(highlight = factor(highlight, levels = c('Background', 'shared significant gene'))) %>%
    # create the plt object
    ggplot(aes(x     = gene_expression_change,
               y     = pvalue_adjusted,
               alpha = highlight,
               size  = highlight,
               shape = highlight,
               color = highlight)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05),
               linetype='dashed') +
    geom_vline(xintercept = 0.5, linetype='dashed') +
    geom_vline(xintercept = -0.5, linetype='dashed') +
    scale_color_manual(values  = c('Background' = 'black',
                                  'shared significant gene' = "blue")) +
    scale_alpha_manual(values  = c('Background' = 0.4,
                                  'shared significant gene' = 1)) +
    scale_shape_manual(values  = c('Background' = 1,
                                  'shared significant gene' = 16)) +
    scale_size_manual(values   = c('Background' = 1,
                                 'shared significant gene' = 2)) +
    geom_label_repel(data = ~filter(., both_sig_greater_thres),
                     aes(label = gene,
                         x = gene_expression_change,
                         y = pvalue_adjusted),
                     size = 5,
                     color = "black",
                     fill = "white",
                     box.padding = 0.4,
                     box.colour = "black",
                     segment.color = "black",
                     max.overlaps = Inf) +
    theme_bw() +
    coord_cartesian(xlim = c(-5,5), ylim = c(0,30)) +
    labs(title = plot_title, x = "log2(Fold Change)", y = "-log10(p-value)") +
    guides(alpha = "none",
           color = guide_legend(override.aes = list(alpha = 1)))
    # + theme(legend.position = 'none')
  # this is how you remove the legend
  #
}

plot_list = map2(data_df_split, names(data_df_split), plot_split)

plot_list$`Low Dose Noc` = plot_list$`Low Dose Noc` +
  theme(legend.position = 'none')

#Combine the plots side by side
combined_plot <- plot_list$`Low Dose Noc` + plot_list$`NGF Deprivation` +
  plot_layout(ncol = 2,
              guides = 'collect')+
  plot_annotation(tag_levels = "A")


# Display the combined plots with labeled genes in white boxes and legend
combined_plot

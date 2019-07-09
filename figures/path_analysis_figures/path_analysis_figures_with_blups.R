#script to merge and plot direct and indirect effects
library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/7_path_analysis/output/blups_analysis/")

#looking at fri and flc individually
path_data_FRI_FLC_ind <- rbindlist(lapply(list.files(pattern = "path_eff_blups.csv"), 
                                          function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE)))
path_data_FRI_FLC_ind$pop <- revalue(path_data_FRI_FLC_ind$pop, c("sha" = "Sha", "oy" = "Oy-0", "jea" = "Jea", "ita" = "Ita-0", "cvi" = "Cvi-0", "bur" = "Bur-0", "blh" = "Blh-1"))
path_data_FRI_FLC_ind$pop_facet <- factor(path_data_FRI_FLC_ind$pop, levels = c("sha" = "Sha", "oy" = "Oy-0", "jea" = "Jea", "ita" = "Ita-0", "cvi" = "Cvi-0", "bur" = "Bur-0", "blh" = "Blh-1"))

path_data_FRI_FLC_ind$trait_facet <- factor(path_data_FRI_FLC_ind$trait, levels = c("bd", "rdry", "h3h1", "idry"))
path_data_FRI_FLC_ind$effect_type <- revalue(path_data_FRI_FLC_ind$effect_type, c("dir" = "Direct", "ind" = "Indirect"))

geom_point_size <- 30

path_data_diff_B4 <- subset(path_data_FRI_FLC_ind, env == "diff" & allele_comb == "B4")
path_data_diff_B4[path_data_diff_B4 == 0] <- NA
path_data_diff_B4_ggplot <- ggplot(data = path_data_diff_B4, aes(x = est, y = pop_facet, color = effect_type)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.25, size = 2) +
    geom_vline(xintercept = 0, 
               color = "black", size = 0.25) + 
    facet_wrap(~ trait_facet, nrow = 1) + 
    xlab("QTL effect") + 
    ylab("Population") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none") + 
    scale_x_continuous(
        labels = scales::number_format(accuracy = 0.1)) +
    ggtitle("m_4_41028")

path_data_diff_B42 <- subset(path_data_FRI_FLC_ind, env == "diff" & allele_comb == "B42")
path_data_diff_B42[path_data_diff_B42 == 0] <- NA
path_data_diff_B42_ggplot <- ggplot(data = path_data_diff_B42, aes(x = est, y = pop_facet, color = effect_type)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.25, size = 2) +
    geom_vline(xintercept = 0, 
               color = "black", size = 0.25) + 
    facet_wrap(~ trait_facet, nrow = 1) + 
    xlab("QTL effect") + 
    ylab("Population") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.title = element_text(size = 20),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20)) + 
    ggtitle("m_4_9222034")

path_data_diff_B5 <- subset(path_data_FRI_FLC_ind, env == "diff" & allele_comb == "B5")
path_data_diff_B5[path_data_diff_B5 == 0] <- NA
path_data_diff_B5_ggplot <- ggplot(data = path_data_diff_B5, aes(x = est, y = pop_facet, color = effect_type)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.25, size = 2) +
    geom_vline(xintercept = 0, 
               color = "black", size = 0.25) + 
    facet_wrap(~ trait_facet, nrow = 1) + 
    xlab("QTL effect") + 
    ylab("Population") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.title = element_text(size = 20),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20)) + 
    scale_x_continuous(breaks = round(seq(min(path_data_diff_B5$est, na.rm = TRUE), max(path_data_diff_B5$est, na.rm = TRUE), by = 0.3), 1)) +
    ggtitle("m_5_3799350")

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/path_analysis_figures/")

ggarrange(path_data_diff_B4_ggplot, path_data_diff_B5_ggplot,
          vjust = 1.1,
          hjust = c(-4.5, -0.5),
          labels = c("A", "B"),
          font.label = list(size = 30))
ggsave("SAR4_SAR5_effect_comparison.png", device = "png", width = 15, height = 6)

#no legend and yaxis text
path_data_diff_B42_ggplot <- ggplot(data = path_data_diff_B42, aes(x = est, y = pop_facet, color = effect_type)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.25, size = 2) +
    geom_vline(xintercept = 0, 
               color = "black", size = 0.25) + 
    facet_wrap(~ trait_facet, nrow = 1) + 
    xlab("QTL effect") + 
    ylab("Population") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.title = element_text(size = 20),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20))  +
    ggtitle("m_4_9222034")

#yaxis text and legend
path_data_diff_B5_ggplot <- ggplot(data = path_data_diff_B5, aes(x = est, y = pop_facet, color = effect_type)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.25, size = 2) +
    geom_vline(xintercept = 0, 
               color = "black", size = 0.25) + 
    facet_wrap(~ trait_facet, nrow = 1) + 
    xlab("QTL effect") + 
    ylab("Population") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.title = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none") + 
    #scale_x_continuous(breaks = round(seq(min(path_data_diff_B5$est, na.rm = TRUE), max(path_data_diff_B5$est, na.rm = TRUE), by = 0.3), 1)) +
    ggtitle("m_5_3799350")

ggarrange(path_data_diff_B4_ggplot, path_data_diff_B42_ggplot, path_data_diff_B5_ggplot, 
          vjust = 1.1,
          #hjust = c(-4.5, , -4),
          labels = c("A", "B", "C"),
          font.label = list(size = 30))
ggsave("SAR4_SAR42_SAR5_effect_blups_comparison.png", device = "png", width = 18, height = 12)

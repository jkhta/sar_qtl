#script to merge and plot direct and indirect effects
library(data.table)
library(ggplot2)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/7_path_analysis/data/")
path_data <- rbindlist(lapply(list.files(pattern = "path_eff_blups.csv"), 
                              function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE)))
path_data$pop_facet <- factor(path_data$pop, levels = c("sha", "oy", "jea", "ita", "cvi", "bur", "blh"))
path_data$trait_facet <- factor(path_data$trait, levels = c("bd", "rdry", "h3h1", "idry"))

for (i in unique(path_data$env)) {
    for (j in unique(path_data$allele_comb)) {
        path_data_allele_subset <- subset(path_data, env == i & allele_comb == j)
        ggplot(data = path_data_allele_subset, aes(x = est, y = pop_facet, color = effect_type)) +
            geom_point() +
            geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.1) +
            geom_vline(xintercept = 0, 
                       color = "black", size = 0.25) + 
            facet_wrap(~ trait_facet, nrow = 1) + 
            xlab("QTL effect") + 
            ylab("population") +
            theme(legend.text = element_text(size = 20),
                  legend.title = element_blank(), 
                  axis.title = element_text(size = 20),
                  axis.text.y = element_text(size = 20)) + 
            ggtitle(i)
        ggsave(paste(i, j, "FRI_FLC_ind_blups_path_effects.png", sep = "_"), width = 6, height = 4, units = "in")
    }
}


#script to merge and plot direct and indirect effects
library(data.table)
library(ggplot2)

rm(list = ls())

setwd("/Users/jkhta/Desktop/nam_cam_fixing/34 - pop_path_analysis/output/")
path_data <- rbindlist(lapply(list.files(pattern = "path_eff.csv"), 
                              function(x) fread(x, sep = ",", header = TRUE, stringsAsFactors = FALSE)))
path_data$pop_facet <- factor(path_data$pop, levels = c("sha", "oy", "jea", "ita", "cvi", "bur", "blh"))
path_data$trait_facet <- factor(path_data$trait, levels = c("bd", "rdry", "h3h1", "idry"))

for (i in unique(path_data$allele_comb)) {
  path_data_allele_subset <- subset(path_data, allele_comb == i & env == "diff")
  ggplot(data = path_data_allele_subset, aes(x = est, y = pop_facet, color = effect_type)) +
    geom_point() +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.2) +
    geom_vline(xintercept = 0, 
               color = "black", size = 0.25) + 
    facet_wrap(~ trait_facet, nrow = 2) + 
    xlab("QTL effect") + 
    ylab("population") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.title = element_text(size = 20)) +
    ggtitle(i)
  ggsave(paste(i, "pop_path_effects.png", sep = "_"), width = 6, height = 4, units = "in")
}

for (i in unique(path_data$allele_comb)) {
  path_data_allele_subset <- subset(path_data, allele_comb == i & env == "diff")
  path_data_allele_subset$trait_facet <- factor(path_data_allele_subset$trait_facet, levels = c("idry", "h3h1", "rdry", "bd"))
  path_data_allele_subset$pop_facet <- factor(path_data_allele_subset$pop_facet, levels = c("blh", "bur", "cvi", "ita", "jea", "oy", "sha"))
  ggplot(data = path_data_allele_subset, aes(x = est, y = trait_facet, color = effect_type)) +
    geom_point() +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.2) +
    geom_vline(xintercept = 0, 
               color = "black", size = 0.25) + 
    facet_wrap(~ pop_facet, nrow = 2) + 
    xlab("QTL effect") + 
    ylab("trait") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.title = element_text(size = 20)) +
    ggtitle(i)
  ggsave(paste(i, "pop_path_effects_by_pop.png", sep = "_"), width = 8, height = 4, units = "in")
}

for (i in unique(path_data$allele_comb)) {
  path_data_allele_subset <- subset(path_data, allele_comb == i & env == "sun")
  ggplot(data = path_data_allele_subset, aes(x = est, y = pop_facet, color = effect_type)) +
    geom_point() +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.2) +
    geom_vline(xintercept = 0, 
               color = "black", size = 0.25) + 
    facet_wrap(~ trait_facet, nrow = 2) + 
    xlab("QTL effect") + 
    ylab("population") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.title = element_text(size = 20)) +
    ggtitle(i)
  ggsave(paste(i, "sun_pop_path_effects.png", sep = "_"), width = 6, height = 4, units = "in")
}

for (i in unique(path_data$allele_comb)) {
  path_data_allele_subset <- subset(path_data, allele_comb == i & env == "shade")
  path_data_allele_subset$trait_facet <- factor(path_data_allele_subset$trait_facet, levels = c("idry", "h3h1", "rdry", "bd"))
  path_data_allele_subset$pop_facet <- factor(path_data_allele_subset$pop_facet, levels = c("blh", "bur", "cvi", "ita", "jea", "oy", "sha"))
  ggplot(data = path_data_allele_subset, aes(x = est, y = trait_facet, color = effect_type)) +
    geom_point() +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.2) +
    geom_vline(xintercept = 0, 
               color = "black", size = 0.25) + 
    facet_wrap(~ pop_facet, nrow = 2) + 
    xlab("QTL effect") + 
    ylab("trait") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_blank(), 
          axis.title = element_text(size = 20)) +
    ggtitle(i)
  ggsave(paste(i, "shade_pop_path_effects_by_pop.png", sep = "_"), width = 8, height = 4, units = "in")
}

test <- subset(path_data, pop == "sha" & env == "shade" & allele_comb == "BA")
test$trait <- factor(test$trait, levels = c("bd", "rdry", "h3h1", "idry"))
test[label == "BA_bd_ind_shade"]$est <- NA
test[label == "BA_h3h1_dir_shade"]$est <- NA
test[label == "BA_idry_dir_shade"]$est <- NA

ggplot(test, aes(x = allele_comb, y = est, group = effect_type, fill = effect_type)) + 
  geom_bar(aes(y = est), position = position_dodge(width = 1), size = 3, stat = "identity") + 
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.4, position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0) +
  facet_grid(~ trait) +
  theme(legend.text = element_text(size = 25),
        legend.title = element_blank(), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  #axis.title.x = element_blank(), 
  #strip.text = element_text(size = 20)) +
  ylab("QTL effects") +
  xlab("Trait")

ggplot(test, aes(x = allele_comb, y = est, fill = effect_type)) +
  geom_bar(aes(x = allele_comb, y = est, fill = effect_type), position = position_dodge(width = 0), size = 5, stat = "identity") + 
  geom_errorbar(aes(ymin = (est - se), ymax = (est + se)), width = 0.4, position = position_dodge(width = 0)) +
  geom_hline(yintercept = 0) +
  facet_grid(~ trait, space = "free") + 
  #facet_grid(~ effect_type) +
  #facet_grid(~ path, margins = c("rdry", "h3h1", "idry")) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  #facet_grid(~ trait) + 
  #axis.title.x = element_blank(), 
  #strip.text = element_text(size = 20)) +
  ylab("Effect estimates") +
  xlab(NULL) +
  labs(color = "Treatment")


path_data_allele_subset <- subset(path_data, allele_comb == "BA" & env == "diff")
ggplot(data = path_data_allele_subset, aes(x = est, y = pop_facet, color = effect_type)) +
  geom_point() +
    geom_errorbarh(aes(xmax = est + se, xmin = est - se), height = 0.2) +
  geom_vline(xintercept = 0, 
             color = "black", size = 0.25) + 
  facet_wrap(~ trait_facet, nrow = 2) + 
  xlab("QTL effect") + 
  ylab("population") +
  theme(legend.text = element_text(size = 30),
        legend.title = element_blank(), 
        axis.title = element_text(size = 30)) +
  ggtitle(i)


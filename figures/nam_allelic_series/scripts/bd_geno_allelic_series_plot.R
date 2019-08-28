#this script will plot the allelic series for each QTL and trait across populations
library(data.table)
library(ggplot2)
library(ggpubr)

rm(list = ls())

#reading in the allelic series data
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/data/")

nam_trait_lme4qtl_fixef <- fread("trait_lme4qtl_fixef_no_cov.csv",
                                 sep = ",",
                                 header = TRUE, 
                                 stringsAsFactors = FALSE)

#this is so the allele represents substituting from Col-0 allele to alternative
#parent allele
nam_trait_lme4qtl_fixef$fixef <- nam_trait_lme4qtl_fixef$fixef * -1

bd_geno_fixef <- subset(nam_trait_lme4qtl_fixef, trait == "bd_geno")

#changing QTL names 
bd_geno_fixef$qtl <- mapvalues(bd_geno_fixef$qtl, from = c("m_4_407208", "m_5_3799350"), to = c("BD4_1", "BD5_2"))

bd_geno_fixef_ggplot_m_4_407208 <- ggplot(subset(bd_geno_fixef, qtl == "BD4_1"), aes(x = pop, y = fixef, fill = pop)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
    facet_grid(~ qtl, scale = "free") +
    xlab("Population") +
    ylab("QTL effect") +
    theme(text = element_text(size = 25),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")

bd_geno_fixef_ggplot_m_5_3799350 <- ggplot(subset(bd_geno_fixef, qtl == "BD5_2"), aes(x = pop, y = fixef, fill = pop)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
    facet_grid(~ qtl, scale = "free") +
    xlab("Population") +
    ylab("QTL effects") +
    theme(text = element_text(size = 25),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.title.y = element_blank())

ggarrange(bd_geno_fixef_ggplot_m_4_407208, 
          bd_geno_fixef_ggplot_m_5_3799350,
          nrow = 1, ncol = 2,
          labels = c("A", "B"),
          font.label = list(size = 25))

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/img/")
ggsave("bd_geno_allelic_series.png", device = "png", width = 14, height = 8)

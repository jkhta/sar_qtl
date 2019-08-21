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

#only working with gxe traits
nam_trait_lme4qtl_fixef_gxe <- subset(nam_trait_lme4qtl_fixef, grepl("_gxe", trait))

for (i in unique(nam_trait_lme4qtl_fixef_gxe$trait)) {
  nam_trait_lme4qtl_fixef_gxe_subset <- subset(nam_trait_lme4qtl_fixef_gxe, trait == i)
  nam_trait_lme4qtl_fixef_gxe_subset <- nam_trait_lme4qtl_fixef_gxe_subset[order(nam_trait_lme4qtl_fixef_gxe_subset$qtl), ]
  
  for (j in unique(nam_trait_lme4qtl_fixef_gxe_subset$qtl)) {
    nam_trait_lme4qtl_fixef_gxe_subset_qtl <- subset(nam_trait_lme4qtl_fixef_gxe_subset, qtl == j)
    plot_name <- paste(i, "fixef_ggplot", j, sep = "_")
    assign(plot_name, ggplot(nam_trait_lme4qtl_fixef_gxe_subset_qtl, aes(x = pop, y = fixef_norm, fill = pop)) +
             geom_bar(stat = "identity") +
             geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
             facet_grid(~ qtl, scale = "free") +
             xlab("Population") +
             ylab("QTL effects") +
             ggtitle(i) + 
             theme(text = element_text(size = 20),
                   axis.text = element_text(size = 20),
                   plot.title = element_text(hjust = 0.5),
                   strip.text.x = element_text(size = 20),
                   legend.position = "none"))
  }
}


ggarrange(bd_gxe_fixef_ggplot_m_4_407208, 
          bd_gxe_fixef_ggplot_m_4_9222034,
          bd_gxe_fixef_ggplot_m_5_3799350,
          r_dry_gxe_fixef_ggplot_m_4_407208,
          r_dry_gxe_fixef_ggplot_m_4_9222034,
          r_dry_gxe_fixef_ggplot_m_5_26130021,
          r_dry_gxe_fixef_ggplot_m_5_3799350,
          h3_h1_gxe_fixef_ggplot_m_2_11952697,
          h3_h1_gxe_fixef_ggplot_m_5_4991121,
          i_dry_gxe_fixef_ggplot_m_4_407208,
          i_dry_gxe_fixef_ggplot_m_5_4216890,
          nrow = 4, 
          ncol = 3, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
          font.label = list(size = 25))

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/nam_allelic_series/img/")
ggsave("gxe_allelic_series_larger.png", device = "png", width = 25, height = 10)

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

#reading in the file to change qtl names
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/qtl_table/data/")
trait_qtl_name_match <- fread("sar_gxe_matched_qtl_name.csv",
                              sep = ",", 
                              header = TRUE,
                              stringsAsFactors = FALSE)
trait_qtl_name_match$trait_qtl <- with(trait_qtl_name_match, paste(Trait, `QTL Marker`, sep = "_"))

#this is so the allele represents substituting from Col-0 allele to alternative
#parent allele
nam_trait_lme4qtl_fixef$fixef <- nam_trait_lme4qtl_fixef$fixef * -1

#only working with gxe traits
nam_trait_lme4qtl_fixef_gxe <- subset(nam_trait_lme4qtl_fixef, grepl("_gxe", trait))
nam_trait_lme4qtl_fixef_gxe$fixef_norm <- with(nam_trait_lme4qtl_fixef_gxe, fixef/abs(trt_eff))
nam_trait_lme4qtl_fixef_gxe$fixef_std_norm <- with(nam_trait_lme4qtl_fixef_gxe, fixef_std/abs(trt_eff))

#changing trait names
new_trait_name_df <- data.frame(old_name = c("bd_gxe", "h3_h1_gxe", "i_dry_gxe", "r_dry_gxe"), 
                                new_name = c("BD", "IG", "IB", "RB"),
                                stringsAsFactors = FALSE)
nam_trait_lme4qtl_fixef_gxe$trait <- mapvalues(nam_trait_lme4qtl_fixef_gxe$trait, 
                                               from = new_trait_name_df$old_name,
                                               to = new_trait_name_df$new_name)
nam_trait_lme4qtl_fixef_gxe$trait_qtl <- with(nam_trait_lme4qtl_fixef_gxe, paste(trait, qtl, sep = "_"))
nam_trait_lme4qtl_fixef_gxe$QTL <- mapvalues(nam_trait_lme4qtl_fixef_gxe$trait_qtl,
                                             from = trait_qtl_name_match$trait_qtl,
                                             to = trait_qtl_name_match$QTL)

for (i in unique(nam_trait_lme4qtl_fixef_gxe$trait)) {
  nam_trait_lme4qtl_fixef_gxe_subset <- subset(nam_trait_lme4qtl_fixef_gxe, trait == i)
  nam_trait_lme4qtl_fixef_gxe_subset <- nam_trait_lme4qtl_fixef_gxe_subset[order(nam_trait_lme4qtl_fixef_gxe_subset$qtl), ]
  
  for (j in unique(nam_trait_lme4qtl_fixef_gxe_subset$QTL)) {
    nam_trait_lme4qtl_fixef_gxe_subset_qtl <- subset(nam_trait_lme4qtl_fixef_gxe_subset, QTL == j)
    plot_name <- paste(i, "fixef_ggplot", j, sep = "_")
    assign(plot_name, ggplot(nam_trait_lme4qtl_fixef_gxe_subset_qtl, aes(x = pop, y = fixef_norm, fill = pop)) +
             geom_bar(stat = "identity") +
             geom_errorbar(aes(ymin = fixef_norm - fixef_std_norm, ymax = fixef_norm + fixef_std_norm), width = 0.2) +
             facet_grid(~ QTL, scale = "free") +
             xlab("Population") +
             ylab("Plasticity change") +
             theme(text = element_text(size = 20),
                   axis.text = element_text(size = 20),
                   plot.title = element_text(hjust = 0.5),
                   strip.text.x = element_text(size = 20),
                   legend.position = "none"))
  }
}

#plotting all plots
ggarrange(BD_fixef_ggplot_BD_1,
          BD_fixef_ggplot_BD_2,
          BD_fixef_ggplot_BD_3,
          RB_fixef_ggplot_RB_1,
          RB_fixef_ggplot_RB_2,
          RB_fixef_ggplot_RB_3,
          RB_fixef_ggplot_RB_4,
          IG_fixef_ggplot_IG_1,
          IG_fixef_ggplot_IG_2,
          IB_fixef_ggplot_IB_1,
          IB_fixef_ggplot_IB_2,
          nrow = 4, 
          ncol = 3, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
          font.label = list(size = 25))

#choosing the plots that look like they have allelic series
#BD_1, RB_2, RB_3, IG_1, IG_2, IB_2

#now modifying the plots
RB_fixef_ggplot_RB_2 <- RB_fixef_ggplot_RB_2 + theme(axis.title.y = element_blank())
RB_fixef_ggplot_RB_3 <- RB_fixef_ggplot_RB_3 + theme(axis.title.y = element_blank())
IG_fixef_ggplot_IG_2 <- IG_fixef_ggplot_IG_2 + theme(axis.title.y = element_blank())
IB_fixef_ggplot_IB_2 <- IB_fixef_ggplot_IB_2 + theme(axis.title.y = element_blank())

ggarrange(BD_fixef_ggplot_BD_1,
          RB_fixef_ggplot_RB_2,
          RB_fixef_ggplot_RB_3,
          IG_fixef_ggplot_IG_1,
          IG_fixef_ggplot_IG_2,
          IB_fixef_ggplot_IB_2,
          nrow = 2, 
          ncol = 3, 
          labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 25))

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/img/")
ggsave("gxe_allelic_series_normalized_subset.png", device = "png", width = 18, height = 10)

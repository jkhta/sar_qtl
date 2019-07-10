#this script will plot the allelic series for each QTL and trait across populations
library(data.table)
library(ggplot2)
library(ggpubr)

rm(list = ls())

#reading in the allelic series data
setwd("/Users/James/Documents/GitHub/sar_qtl/figures/nam_allelic_series/output/")

nam_trait_lme4qtl_fixef <- fread("trait_lme4qtl_fixef.csv",
                                 sep = ",",
                                 header = TRUE, 
                                 stringsAsFactors = FALSE)

#this is so the allele represents substituting from Col-0 allele to alternative
#parent allele
nam_trait_lme4qtl_fixef$fixef <- nam_trait_lme4qtl_fixef$fixef * -1

bd_gxe_fixef <- subset(nam_trait_lme4qtl_fixef, trait == "bd_gxe")

bd_gxe_fixef_ggplot_m_4_407208 <- ggplot(subset(bd_gxe_fixef, qtl == "m_4_407208"), aes(x = pop, y = fixef, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
  facet_grid(~ qtl, scale = "free") +
  xlab("Population") +
  ylab("QTL effects") +
  ggtitle("Bolting time") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 20),
        legend.position = "none")

bd_gxe_fixef_ggplot_m_4_9222034 <- ggplot(subset(bd_gxe_fixef, qtl == "m_4_9222034"), aes(x = pop, y = fixef, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
  facet_grid(~ qtl, scale = "free") +
  xlab("Population") +
  ylab("QTL effects") +
  ggtitle("Bolting time") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 20),
        legend.position = "none")

bd_gxe_fixef_ggplot_m_5_3799350 <- ggplot(subset(bd_gxe_fixef, qtl == "m_5_3799350"), aes(x = pop, y = fixef, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
  facet_grid(~ qtl, scale = "free") +
  xlab("Population") +
  ylab("QTL effects") +
  ggtitle("Bolting time") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 20),
        legend.position = "none")

h3h1_gxe_fixef <- subset(nam_trait_lme4qtl_fixef, trait == "h3_h1_gxe")

h3h1_gxe_fixef_ggplot_m_2_11952697 <- ggplot(subset(h3h1_gxe_fixef, qtl == "m_2_11952697"), aes(x = pop, y = fixef, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
  facet_grid(~ qtl, scale = "free") +
  xlab("Population") +
  ylab("QTL effects") +
  ggtitle("Inflorescence growth") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 20),
        legend.position = "none")

h3h1_gxe_fixef_ggplot_m_5_4991121 <- ggplot(subset(h3h1_gxe_fixef, qtl == "m_5_4991121"), aes(x = pop, y = fixef, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
  facet_grid(~ qtl, scale = "free") +
  xlab("Population") +
  ylab("QTL effects") +
  ggtitle("Inflorescence growth") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 20),
        legend.position = "none")

r_dry_gxe_fixef <- subset(nam_trait_lme4qtl_fixef, trait == "r_dry_gxe")
r_dry_gxe_fixef_ggplot <- ggplot(r_dry_gxe_fixef, aes(x = pop, y = fixef, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
  facet_grid(~ qtl, scale = "free") +
  xlab("Population") +
  ylab("QTL effects") +
  ggtitle("Rosette biomass") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 20),
        legend.position = "none")
i_dry_gxe_fixef <- subset(nam_trait_lme4qtl_fixef, trait == "i_dry_gxe")
i_dry_gxe_fixef_ggplot_m_5_20260494 <- ggplot(subset(i_dry_gxe_fixef, qtl == "m_5_20260494"), aes(x = pop, y = fixef, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
  facet_grid(~ qtl, scale = "free") +
  xlab("Population") +
  ylab("QTL effects") +
  ggtitle("Inflorescence biomass") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 20),
        legend.position = "none")

i_dry_gxe_fixef_ggplot_m_5_4216890 <- ggplot(subset(i_dry_gxe_fixef, qtl == "m_5_4216890"), aes(x = pop, y = fixef, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
  facet_grid(~ qtl, scale = "free") +
  xlab("Population") +
  ylab("QTL effects") +
  ggtitle("Inflorescence biomass") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 20),
        legend.position = "none")

ggarrange(bd_gxe_fixef_ggplot_m_4_407208, 
          bd_gxe_fixef_ggplot_m_4_9222034,
          bd_gxe_fixef_ggplot_m_5_3799350,
          r_dry_gxe_fixef_ggplot, 
          h3h1_gxe_fixef_ggplot_m_2_11952697,
          h3h1_gxe_fixef_ggplot_m_5_4991121,
          i_dry_gxe_fixef_ggplot_m_5_20260494,
          i_dry_gxe_fixef_ggplot_m_5_4216890,
          nrow = 4, 
          ncol = 2, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          font.label = list(size = 25))

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/nam_allelic_series/output/")
ggsave("gxe_allelic_series_larger.png", device = "png", width = 25, height = 10)

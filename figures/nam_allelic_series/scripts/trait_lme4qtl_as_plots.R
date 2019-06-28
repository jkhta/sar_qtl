#this script will plot the allelic series for each QTL and trait across populations
library(data.table)
library(ggplot2)

rm(list = ls())

#reading in the allelic series data
setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/output/")

nam_trait_lme4qtl_fixef <- fread("trait_lme4qtl_fixef.csv",
                                 sep = ",",
                                 header = TRUE, 
                                 stringsAsFactors = FALSE)

#this is so the allele represents substituting from Col-0 allele to alternative
#parent allele
nam_trait_lme4qtl_fixef$fixef <- nam_trait_lme4qtl_fixef$fixef * -1

#for each trait, plot out the allelic series
for (i in unique(nam_trait_lme4qtl_fixef$trait)) {
  trait_fixef <- subset(nam_trait_lme4qtl_fixef, trait == i)
  
  ggplot(trait_fixef, aes(x = pop, y = fixef, fill = pop)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = fixef - fixef_std, ymax = fixef + fixef_std), width = 0.2) +
    facet_grid(~ qtl, scale = "free") +
    xlab("Population") +
    ylab("QTL effects") +
    theme(text = element_text(size = 20),
          legend.position = "none")
    
  ggfilename <- paste(i, "as_with_std.png")
  
  ggsave(ggfilename, device = "png", width = 20, height = 10)
}

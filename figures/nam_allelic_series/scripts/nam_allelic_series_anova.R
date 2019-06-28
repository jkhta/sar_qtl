library(data.table)
library(HH)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/output/")

#reading in the fixed effect data
trait_qtl_fixef <- fread("trait_lme4qtl_fixef.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

#grabbing only the FRI fixed effects
FRI_subset <- subset(trait_qtl_fixef, qtl == "m_4_407208" & trait == "bd_geno")
FRI_subset$pop <- as.factor(FRI_subset$pop)

#running the anova with sufficient summary statistics
FRI_aov <- aovSufficient(fixef ~ pop, 
                         data = FRI_subset, 
                         weights = FRI_subset$count, 
                         sd = FRI_subset$fixef_std*sqrt(FRI_subset$count)/2)
summary(FRI_aov)

#multiple comparisons between groups
FRI_mmc <- mmc(FRI_aov,
               linfct = mcp(pop = "Tukey"),
               df = FRI_aov$df.residual,
               vcov. = vcovSufficient)

as.data.frame(anova(FRI_aov))

#glht for cld for groupings 
FRI_glht <- glht(FRI_aov, 
                 linfct = mcp(pop = "Tukey"), 
                 df = FRI_aov$df.residual,
                 vcov. = vcovSufficient)

summary(FRI_glht)
cld(FRI_glht, level = 0.05/11)
plot(cld(FRI_glht, level = 0.05/11))

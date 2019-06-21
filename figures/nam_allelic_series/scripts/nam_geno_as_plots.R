library(data.table)
library(tidyr)
library(ggplot2)
library(plyr)
library(ggpubr)

rm(list = ls())

setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/input/")

pop_order <- names(readRDS("nam_all_traits_ind_pop_pheno_geno_proximal_james.RDS"))
pop_order <- revalue(pop_order, c("21RV_21RV" = "Blh-1", 
                                  "20RV_20RV" = "Bur-0", 
                                  "8RV_8RV" = "Cvi-0", 
                                  "29RV_29RV" = "Ita-0", 
                                  "28RV_28RV" = "Jea", 
                                  "27RV_27RV" = "Oy-0", 
                                  "13RV_13RV" = "Sha"))
#reading in only the gxe data
trait_models <- lapply(list.files()[grepl("geno", list.files()) & grepl("GridLMM_stepwise_model", list.files())], function(x) readRDS(x))
trait_names <- sapply(strsplit(list.files()[grepl("geno", list.files()) & grepl("GridLMM_stepwise_model", list.files())], split = "_james"), function(x) x[1])

snp_beta_extractor <- function(trait_model) {
  #grabbing sig qlt
  trait_sig_qtl <- trait_model$found_qtl
  
  #grabbing the individual pop betas
  pop_betas <- trait_model$ind_pop_qtl_betas
  
  #grabbing only the snp betas, separate from the intercept and trait covariates
  beta_ind <- tail(2:ncol(pop_betas), length(trait_sig_qtl))
  pop_snp_betas <- pop_betas[,..beta_ind]
  
  #relabeling snps to their marker names
  colnames(pop_snp_betas) <- trait_sig_qtl
  
  #adding the population column
  pop_snp_betas <- cbind(data.table(pop = pop_order), pop_snp_betas)
  return(pop_snp_betas)
}

#convertin data to long format
trait_snp_betas <- lapply(trait_models, function(x) snp_beta_extractor(x))
names(trait_snp_betas) <- trait_names

bd_geno_long <- gather(trait_snp_betas[["bd_geno"]], qtl, qtl_effect, 2:ncol(trait_snp_betas[["bd_geno"]]))
bd_geno_long$qtl_effect <- bd_geno_long$qtl_effect * -1
bd_geno_as <- ggplot(bd_geno_long, aes(x = pop, y = qtl_effect, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, 
             color = "black", size = 1) + 
  facet_wrap(~ qtl, nrow = 1) +
  xlab("Population") + 
  ylab("QTL Effect") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 15))

r_dry_geno_long <- gather(trait_snp_betas[["r_dry_geno"]], qtl, qtl_effect, 2:ncol(trait_snp_betas[["r_dry_geno"]]))
r_dry_geno_long$qtl_effect <- r_dry_geno_long$qtl_effect * -1
r_dry_geno_as <- ggplot(r_dry_geno_long, aes(x = pop, y = qtl_effect, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, 
             color = "black", size = 1) + 
  facet_wrap(~ qtl, nrow = 1) +
  xlab("Population") + 
  ylab("QTL Effect") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        legend.position = "none")

h3_h1_geno_long <- gather(trait_snp_betas[["h3_h1_geno"]], qtl, qtl_effect, 2:ncol(trait_snp_betas[["h3_h1_geno"]]))
h3_h1_geno_long$qtl_effect <- h3_h1_geno_long$qtl_effect * -1
h3_h1_geno_as <- ggplot(h3_h1_geno_long, aes(x = pop, y = qtl_effect, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, 
             color = "black", size = 1) + 
  facet_wrap(~ qtl, nrow = 1) +
  xlab("Population") + 
  ylab("QTL Effect") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 15))

i_dry_geno_long <- gather(trait_snp_betas[["i_dry_geno"]], qtl, qtl_effect, 2:ncol(trait_snp_betas[["i_dry_geno"]]))
i_dry_geno_long$qtl_effect <- i_dry_geno_long$qtl_effect * -1
i_dry_geno_as <- ggplot(i_dry_geno_long, aes(x = pop, y = qtl_effect, fill = pop)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, 
             color = "black", size = 1) + 
  facet_wrap(~ qtl, nrow = 1) +
  xlab("Population") + 
  ylab("QTL Effect") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        legend.position = "none")

ggarrange(bd_geno_as, 
          h3_h1_geno_as, 
          r_dry_geno_as, 
          i_dry_geno_as, 
          nrow = 4, 
          ncol = 1, 
          labels = c("A", "B", "C", "D"))


setwd("/Users/jkhta/Documents/GitHub/sar_qtl/figures/nam_allelic_series/output/")
ggsave("geno_as_plots.png", device = "png", width = 18, height = 20)



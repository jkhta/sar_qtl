#this script will generate a plot for the linkage map comparing the covariate vs. no covariate methods
library(ggpubr)
library(png)

rm(list = ls())

setwd("/Users/James/Documents/GitHub/sar_qtl/figures/path_sum_model/img/")

sum_path_sun <- readPNG("sig_paths_sum_sun_model.png")
sum_path_shade <- readPNG("sig_paths_sum_shade_model.png")

im_A <- ggplot() + 
    background_image(sum_path_sun) +
    # This ensures that the image leaves some space at the edges
    theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent"))

im_B <- ggplot() + background_image(sum_path_shade) + 
    theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent"))

ggarrange(im_A, im_B,
          labels = c("A", "B"), 
          nrow = 2,
          font.label = list(size = 55))

ggsave("sun_vs_shade_sum_path_model.png", device = "png", height = 20, width = 12, bg = "transparent")

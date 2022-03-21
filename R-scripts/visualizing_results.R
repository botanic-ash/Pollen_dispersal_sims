#load in results from main_loop.R
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")
load("prop_alleles_capt.Rdata")

#libraries 
library(dplyr)
library(tidyr)
library(ggplot2)

#Preparing data for plotting
combined_data = rbind(prop_capt_all_same[,,1], prop_capt_all_eligible[,,1], prop_capt_skewed[,,1])
combined_data = as.data.frame(combined_data)
filtered = combined_data %>% filter(maternal_trees == 10)

ggplot(data=filtered, aes(x=total_seeds, y=prop_capt, color=donor_type)) +
  geom_line() +
  geom_point()

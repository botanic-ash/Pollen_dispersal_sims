#Code written by Kaylee Rosenberger

#libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(wesanderson)

#load in results from main_loop.R
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")
###MAKE SURE TO LOAD IN THE CORRECT DATA
load("prop_alleles_capt_new.Rdata")
sample_type = "ideal"
#load("prop_alleles_capt_skewed_new.Rdata")
#sample_type = realistic

#defining arrays to store the results in for easy plotting 
same = array(dim=c(50,2)) #50 rows for 50 replicates, 2 columns for prop_capt and donor type
eligible = array(dim=c(50,2))
skewed = array(dim=c(50,2))

#defining which scenario (row) you want to pull out to plot 
mat_trees=100
num_seeds=500

#looping over simulation replicates to pull equivalent scenarios to comapre among donor modes 
for(i in 1:50){
  temp = as.data.frame(prop_capt_all_same[,,i])
  temp_mod = filter(temp, total_seeds==num_seeds, maternal_trees==mat_trees)
  same[i,1] = temp_mod[,1]#getting row 1, the prop. alleles conserved
  same[i,2] = temp_mod[,5]#getting row 5, the donor type 
  
  temp = as.data.frame(prop_capt_all_eligible[,,i])
  temp_mod = filter(temp, total_seeds==num_seeds, maternal_trees==mat_trees)
  eligible[i,1] = temp_mod[,1] 
  eligible[i,2] = temp_mod[,5]
  
  temp = as.data.frame(prop_capt_skewed[,,i])
  temp_mod = filter(temp, total_seeds==num_seeds, maternal_trees==mat_trees)
  skewed[i,1] = temp_mod[,1]
  skewed[i,2] = temp_mod[,5]
}

#data processing for plotting
equal_comparison = rbind(same, eligible, skewed)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle(paste(num_seeds, "seeds from", mat_trees, "trees")) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(.3,1) +
  theme(legend.position = "none")
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/Figures")
ggsave(paste(num_seeds, "_", mat_trees, "_", sample_type, ".png", sep=""), height=5, width=5, units="in")
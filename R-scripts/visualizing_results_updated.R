#Code written by Kaylee Rosenberger

#libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(wesanderson)

#load in results from main_loop.R
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")#
###MAKE SURE TO LOAD IN THE CORRECT DATA
#load("prop_alleles_capt_new_727.Rdata")
#sample_type = "ideal"
load("prop_alleles_capt_skewed_new_727.Rdata")
sample_type = "realistic"

##########################################################################################
#First type of plot:
#defining arrays to store the results in for easy plotting 
same = array(dim=c(50,2)) #50 rows for 50 replicates, 2 columns for prop_capt and donor type
eligible = array(dim=c(50,2))
skewed = array(dim=c(50,2))

#defining which scenario (row) you want to pull out to plot 
#NOTE: can make these a vector or list to loop over scenarios
mat_trees=25
num_seeds=400

#looping over simulation replicates to pull equivalent scenarios to comapre among donor modes 
for(i in 1:50){
  same_df = as.data.frame(prop_capt_all_same[,,i])
  filter_same_df = filter(same_df, total_seeds==num_seeds, maternal_trees==mat_trees)
  same[i,1] = filter_same_df[,1]#getting row 1, the prop. alleles conserved
  same[i,2] = filter_same_df[,5]#getting row 5, the donor type 
  
  eligible_df = as.data.frame(prop_capt_all_eligible[,,i])
  filter_eligible_df = filter(eligible_df, total_seeds==num_seeds, maternal_trees==mat_trees)
  eligible[i,1] = filter_eligible_df[,1] 
  eligible[i,2] = filter_eligible_df[,5]
  
  skewed_df = as.data.frame(prop_capt_skewed[,,i])
  filter_skewed_df = filter(skewed_df, total_seeds==num_seeds, maternal_trees==mat_trees)
  skewed[i,1] = filter_skewed_df[,1]
  skewed[i,2] = filter_skewed_df[,5]
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
  ylim(0.25,1) +
  theme(legend.position = "none")
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/Figures")
ggsave(paste(num_seeds, "_", mat_trees, "_", sample_type, ".png", sep=""), height=5, width=5, units="in")

#############################################################
#Now we are making plots directly comparing the scenarios where an equal number of seeds are sampled per tree
#vs. scenarios where an unequal number of seeds are sampled
#Here, we plot the two types of scenarios on the same plots 

#load data
load("prop_alleles_capt_new.Rdata")
#renaming data to differentiate 
prop_capt_all_same_equal = prop_capt_all_same
prop_capt_all_eligible_equal = prop_capt_all_eligible
prop_capt_skewed_equal = prop_capt_skewed

#load data
load("prop_alleles_capt_skewed_new.Rdata")
prop_capt_all_same_skewed = prop_capt_all_same
prop_capt_all_eligible_skewed = prop_capt_all_eligible
prop_capt_skewed_skewed = prop_capt_skewed

#defining arrays to store the results in for easy plotting - EQUAL
same_equal = array(dim=c(50,3)) #50 rows for 50 replicates, 2 columns for prop_capt and donor type
eligible_equal = array(dim=c(50,3))
skewed_equal = array(dim=c(50,3))

#defining arrays to store the results in for easy plotting- SKEWED
same_skewed = array(dim=c(50,3)) #50 rows for 50 replicates, 2 columns for prop_capt and donor type
eligible_skewed = array(dim=c(50,3))
skewed_skewed = array(dim=c(50,3))

####################
#50 maternal trees 500 seeds
for(i in 1:50) {
  same_equal[i,1] = prop_capt_all_same_equal[10,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_equal[i,2] = prop_capt_all_same_equal[10,5,i]
  same_equal[i,3] = "equal"
  
  same_skewed[i,1] = prop_capt_all_same_skewed[5,1,i]
  same_skewed[i,2] = prop_capt_all_same_skewed[5,5,i]
  same_skewed[i,3] = "skewed"
  
  eligible_equal[i,1] = prop_capt_all_eligible_equal[10,1,i]
  eligible_equal[i,2] = prop_capt_all_eligible_equal[10,5,i]
  eligible_equal[i,3] = "equal"
  
  eligible_skewed[i,1] = prop_capt_all_eligible_skewed[5,1,i]
  eligible_skewed[i,2] = prop_capt_all_eligible_skewed[5,5,i]
  eligible_skewed[i,3] = "skewed"
  
  skewed_equal[i,1] = prop_capt_skewed_equal[10,1,i]
  skewed_equal[i,2] = prop_capt_skewed_equal[10,5,i]
  skewed_equal[i,3] = "equal"
  
  skewed_skewed[i,1] = prop_capt_skewed_skewed[5,1,i]
  skewed_skewed[i,2] = prop_capt_skewed_skewed[5,5,i]
  skewed_skewed[i,3] = "skewed"
}

#data processing for plotting
equal_comparison = rbind(same_equal, eligible_equal, skewed_equal, same_skewed, eligible_skewed, skewed_skewed)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type", "scenario_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=donor_type, y=as.numeric(prop_capt), fill=scenario_type)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  labs(color="donor_type", fill="scenario_type") +
  ggtitle("50 maternal trees, 500 seeds total") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.85,1)
ggsave("comparison_50_500.png", height=5, width=5, units="in")
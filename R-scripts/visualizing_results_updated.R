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

#############################################################################################
#Now we are making plots directly comparing the scenarios where an equal number of seeds are sampled per tree
#vs. scenarios where an unequal number of seeds are sampled
#Here, we plot the two types of scenarios on the same plots 

#directory of Rdata
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")

#renaming data to differentiate between sample types
#load in ideal scenario data
load("prop_alleles_capt_new_727.Rdata")
prop_capt_all_same_ideal = prop_capt_all_same
prop_capt_all_eligible_ideal = prop_capt_all_eligible
prop_capt_skewed_ideal = prop_capt_skewed

#load in realistic scenario data
load("prop_alleles_capt_skewed_new_727.Rdata")
prop_capt_all_same_realistic = prop_capt_all_same
prop_capt_all_eligible_realistic = prop_capt_all_eligible
prop_capt_skewed_realistic = prop_capt_skewed

#defining arrays to store the results in for easy plotting - EQUAL
same_ideal = array(dim=c(50,3)) #50 rows for 50 replicates, 2 columns for prop_capt and donor type
eligible_ideal = array(dim=c(50,3))
skewed_ideal = array(dim=c(50,3))

#defining arrays to store the results in for easy plotting- SKEWED
same_realistic = array(dim=c(50,3)) #50 rows for 50 replicates, 2 columns for prop_capt and donor type
eligible_realistic = array(dim=c(50,3))
skewed_realistic = array(dim=c(50,3))

#defining which scenario to pull
mat_trees = 100
num_seeds = 200

for(i in 1:50) {
  #SAME IDEAL
  same_df_ideal = as.data.frame(prop_capt_all_same_ideal[,,i])
  filter_same_df_ideal = filter(same_df_ideal, total_seeds==num_seeds, maternal_trees==mat_trees)
  same_ideal[i,1] = filter_same_df_ideal[,1] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_ideal[i,2] = filter_same_df_ideal[,5]
  same_ideal[i,3] = "ideal"
  
  #SAME REALISTIC
  same_df_realistic = as.data.frame(prop_capt_all_same_realistic[,,i])
  filter_same_df_realistic = filter(same_df_realistic, total_seeds==num_seeds, maternal_trees==mat_trees)
  same_realistic[i,1] = filter_same_df_realistic[,1]
  same_realistic[i,2] = filter_same_df_realistic[,5]
  same_realistic[i,3] = "realistic"
  
  #ELIGIBLE IDEAL
  eligible_df_ideal = as.data.frame(prop_capt_all_eligible_ideal[,,i])
  filter_eligible_df_ideal = filter(eligible_df_ideal, total_seeds==num_seeds, maternal_trees==mat_trees)
  eligible_ideal[i,1] = filter_eligible_df_ideal[,1]
  eligible_ideal[i,2] = filter_eligible_df_ideal[,5]
  eligible_ideal[i,3] = "ideal"
  
  #ELIGIBLE REALISTIC
  eligible_df_realistic = as.data.frame(prop_capt_all_eligible_realistic[,,i])
  filter_eligible_df_realistic = filter(eligible_df_realistic, total_seeds==num_seeds, maternal_trees==mat_trees)
  eligible_realistic[i,1] = filter_eligible_df_realistic[,1]
  eligible_realistic[i,2] = filter_eligible_df_realistic[,5]
  eligible_realistic[i,3] = "realistic"
  
  #SKEWED IDEAL
  skewed_df_ideal = as.data.frame(prop_capt_skewed_ideal[,,i])
  filter_skewed_df_ideal = filter(skewed_df_ideal, total_seeds==num_seeds, maternal_trees==mat_trees)
  skewed_ideal[i,1] = filter_skewed_df_ideal[,1]
  skewed_ideal[i,2] = filter_skewed_df_ideal[,5]
  skewed_ideal[i,3] = "ideal"
  
  #SKEWED REALISTIC
  skewed_df_realistic = as.data.frame(prop_capt_skewed_realistic[,,i])
  filter_skewed_df_realistic = filter(skewed_df_realistic, total_seeds==num_seeds, maternal_trees==mat_trees)
  skewed_realistic[i,1] = filter_skewed_df_realistic[,1]
  skewed_realistic[i,2] = filter_skewed_df_realistic[,5]
  skewed_realistic[i,3] = "realistic"
}

#data processing for plotting
equal_comparison = rbind(same_ideal, eligible_ideal, skewed_ideal, same_realistic, eligible_realistic, skewed_realistic)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type", "scenario_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=donor_type, y=as.numeric(prop_capt), fill=scenario_type)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", hide.ns = TRUE) +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  labs(color="donor_type", fill="scenario_type") +
  ggtitle(paste(num_seeds, "seeds from", mat_trees, "trees")) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.25,1)
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/Figures")
ggsave(paste(num_seeds, "_", mat_trees, "_comparison", ".png", sep=""), height=5, width=5, units="in")

##########################
#comparing proportion of alleles captured between ideal and realistic scnearios 
#pick a scenario, and get proportion of alleles captured for all replicates
#then, average prop. alleles captured across all replicates 

#taking the average of all scenarios 
mean(as.numeric(same_ideal[,1]))
mean(as.numeric(same_realistic[,1]))
  
mean(as.numeric(eligible_ideal[,1]))
mean(as.numeric(eligible_realistic[,1]))
  
mean(as.numeric(skewed_ideal[,1]))
mean(as.numeric(skewed_realistic[,1]))

#load in results from main_loop.R
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")
###MAKE SURE TO LOAD IN THE CORRECT DATA
load("prop_alleles_capt.Rdata")
##load("prop_alleles_capt_skewed.Rdata")

#libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)

########################################################################################################
#Data preparation
#Here we average the results across replicates and make new data frames to store this averaged data

#arrays to save averaged results
prop_capt_all_same_avg = array(dim=c(465,5))#465 scenarios and 5 columns needed to store data
prop_capt_all_eligible_avg = array(dim=c(465,5))
prop_capt_skewed_avg = array(dim=c(465,5))
#naming columns for data frame
colnames(prop_capt_all_same_avg) = c("prop_capt", "total_seeds", "maternal_trees", "num_donors", "donor_type")
colnames(prop_capt_all_eligible_avg) = c("prop_capt", "total_seeds", "maternal_trees", "num_donors", "donor_type")
colnames(prop_capt_skewed_avg) = c("prop_capt", "total_seeds", "maternal_trees", "num_donors", "donor_type")

#looping over each scenario (row = 465) to get the average prop. alleles captured in each replicate (slice ~640)
#and for each donor type (why the code is repeated 3x)
for(i in 1:465) {
  #ALL SAME
  avg_prop_all_capt = mean(as.numeric(prop_capt_all_same[i,1,]), na.rm = TRUE) #calculates the means of all replicates for one scenario
  #saving results
  prop_capt_all_same_avg[i,] = prop_capt_all_same[i,,1] #copying over data from previous array
  prop_capt_all_same_avg[i,1] = avg_prop_all_capt
  
  #ALL ELIGIBLE
  avg_prop_all_capt = mean(as.numeric(prop_capt_all_eligible[i,1,]), na.rm = TRUE) #calculates the means of all replicates for one scenario
  #saving results
  prop_capt_all_eligible_avg[i,] = prop_capt_all_eligible[i,,1] #copying over data from previous array
  prop_capt_all_eligible_avg[i,1] = avg_prop_all_capt
  
  #SKEWED
  avg_prop_all_capt = mean(as.numeric(prop_capt_skewed[i,1,]), na.rm = TRUE) #calculates the means of all replicates for one scenario
  #saving results
  prop_capt_skewed_avg[i,] = prop_capt_skewed[i,,1] #copying over data from previous array
  prop_capt_skewed_avg[i,1] = avg_prop_all_capt
}

#Preparing data for plotting
#combining each of the separate arrays into one large array for plotting
combined_data = rbind(prop_capt_all_same_avg, prop_capt_all_eligible_avg, prop_capt_skewed_avg)
combined_data = as.data.frame(combined_data) #converting array to dataframe to use in ggplot

##################################################################################################
#The first type of plot shows averaged results across replicates
#X axis is number of seeds sampled and y axis is the (averaged) proportion of alleles captured
#First, looking into total seeds vs. prop alleles capt for each donor type
#
#the plots shown below vary by number of maternal trees sampled. A common question is how many unique
#trees should be sampled, so we can compare that here
#Additionally, when we try to plot all numbers of maternal trees on one plot, it gets too busy
#so we need to have some constants

#filtering the data--looking at 10 maternal trees only
filtered = combined_data %>% filter(maternal_trees == 10)
ggplot(data=filtered, aes(x=as.numeric(total_seeds), y=as.numeric(prop_capt), group=donor_type, color=donor_type)) +
  geom_line() +
  xlab("Total seeds sampled") +
  ylab("Proportion of alleles captured") +
  ggtitle("Diversity captured from sampling 10 maternal trees")

#looking at 50 maternal trees
filtered = combined_data %>% filter(maternal_trees == 50)
ggplot(data=filtered, aes(x=as.numeric(total_seeds), y=as.numeric(prop_capt), group=donor_type, color=donor_type)) +
  geom_line() +
  xlab("Total seeds sampled") +
  ylab("Proportion of alleles captured") +
  ggtitle("Diversity captured from sampling 50 maternal trees")

#1 maternal tree
filtered = combined_data %>% filter(maternal_trees == 1)
ggplot(data=filtered, aes(x=as.numeric(total_seeds), y=as.numeric(prop_capt), group=donor_type, color=donor_type)) +
  geom_line() +
  xlab("Total seeds sampled") +
  ylab("Proportion of alleles captured") +
  ggtitle("Diversity captured from sampling 1 maternal trees")

##################################################################################################
#Next, looking into diversity captured on y axis vs. number of maternal trees sampled on x axis
#Again, we need to have a constant--total seeds sampled will be constant
#This is a follow up to the last group of plots to more directly compare how sampling more
#maternal trees impacts the diversity captured

#first, looking at a total size of 100 seeds
filtered = combined_data %>% filter(total_seeds == 100)
ggplot(data=filtered, aes(x=as.numeric(maternal_trees), y=as.numeric(prop_capt), group=donor_type, color=donor_type)) +
  geom_line() +
  xlab("Number of maternal trees sampled") +
  ylab("Proportion of alleles captured") +
  ggtitle("Diversity captured from sampling 100 total seeds from varying maternal trees")

#next looking at 200 seeds
filtered = combined_data %>% filter(total_seeds == 200)
ggplot(data=filtered, aes(x=as.numeric(maternal_trees), y=as.numeric(prop_capt), group=donor_type, color=donor_type)) +
  geom_line() +
  xlab("Number of maternal trees sampled") +
  ylab("Proportion of alleles captured") +
  ggtitle("Diversity captured from sampling 200 total seeds from varying maternal trees")

#looking at 50 seeds
filtered = combined_data %>% filter(total_seeds == 50)
ggplot(data=filtered, aes(x=as.numeric(maternal_trees), y=as.numeric(prop_capt), group=donor_type, color=donor_type)) +
  geom_line() +
  xlab("Number of maternal trees sampled") +
  ylab("Proportion of alleles captured") +
  ggtitle("Diversity captured from sampling 50 total seeds from varying maternal trees")

###############################################################################################
#Now, we are comparing "equivalent scenarios"
#We define these equivalent scenarios as the same nimber of maternal trees sampled and the same
#total number of seeds sampled from each tree. 
#With these plots, we can directly compare how the pollen donation mode impacts the diversity captured
#by the same sampling strategy.
#Here, we pull directly from the results array without averaging, so we can note the variation across 
#simulation replicates
#When all plots below are taken together, we can compare across different sample sizes, and across
#different number of trees sampled. 
#NOTE: This section is for the original set of data (where an equal number of seeds is taken from each tree)

#defining arrays to store the results in for easy plotting 
same_25_10 = array(dim=c(500,2)) #500 rows for 500 replicates, 2 columns for prop_capt and donor type
eligible_25_10 = array(dim=c(500,2))
skewed_25_10 = array(dim=c(500,2))

#################
#comparing sampling 25 maternal trees 10 seeds from each tree
#looping over simulation replicates to pull equivalent scenarios to comapre among donor modes 
#We pull directly from the results array by hard coding (which should be avoided, but is in progress)
for(i in 1:500) {
  same_25_10[i,1] = prop_capt_all_same[265,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[265,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[265,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[265,5,i]
    
  skewed_25_10[i,1] = prop_capt_skewed[265,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[265,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("50 seeds from 1 maternal tree") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0,1) +
  theme(legend.position = "none")

###########################

#looping over simulation replicates to pull equivalent scenarios to comapre among donor modes 
#5 seeds 10 maternal tree
for(i in 1:500) {
  same_25_10[i,1] = prop_capt_all_same[20,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[20,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[20,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[20,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[20,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[20,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("5 seeds from 10 maternal trees") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0,1) +
  theme(legend.position = "none")

###########################

#2 seeds 25 maternal tree
for(i in 1:500) {
  same_25_10[i,1] = prop_capt_all_same[7,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[7,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[7,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[7,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[7,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[7,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("2 seeds from 25 maternal trees") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0,1) +
  theme(legend.position = "none")

#######################

#1 seeds 50 maternal tree
for(i in 1:500) {
  same_25_10[i,1] = prop_capt_all_same[1,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[1,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[1,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[1,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[1,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[1,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("1 seed from 50 maternal trees") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw()+
  theme(axis.title = element_blank()) +
  ylim(0,1) +
  theme(legend.position = "none")

#########################

#50 maternal trees 5 seeds per
for(i in 1:500) {
  same_25_10[i,1] = prop_capt_all_same[5,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[5,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[5,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[5,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[5,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[5,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("5 seeds from 50 maternal trees") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0,1) +
  theme(legend.position = "none")

######################

#25 maternal trees 10 seeds per
for(i in 1:500) {
  same_25_10[i,1] = prop_capt_all_same[15,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[15,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[15,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[15,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[15,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[15,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("10 seeds from 25 maternal trees") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0,1) +
  theme(legend.position = "none")

###########################

#10 maternal trees 25 seeds per
for(i in 1:500) {
  same_25_10[i,1] = prop_capt_all_same[40,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[40,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[40,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[40,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[40,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[40,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("25 seeds from 10 maternal trees") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw()  +
  theme(axis.title = element_blank()) +
  ylim(0,1) +
  theme(legend.position = "none")

######################

#1 maternal trees 250 seeds per
for(i in 1:500) {
  same_25_10[i,1] = prop_capt_all_same[465,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[465,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[465,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[465,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[465,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[465,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("250 seeds from 1 maternal trees") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0,1) +
  theme(legend.position = "none")

#####################################################################################################
#SKEWED
#comparing equivalent scenarios, as we did in the previous section, but these plots are for the 
#scenarios in which an unequal number of seeds is taken from each tree 
#Again, these plots do not average results, so we can note the variation across simulation replicates
#We can directly compare how the proportion of alleles captured varies when pollen donation mode varies
#for a given sampling strategy (number of maternal trees and total sample size)

#defining the results arrays to store filtered data
same_25_10 = array(dim=c(50,2)) #50 rows for 50 replicates, 2 columns for prop_capt and donor type
eligible_25_10 = array(dim=c(50,2))
skewed_25_10 = array(dim=c(50,2))

#################### 

#looping over simulation replicates to pull equivalent scenarios to comapre among donor modes 
#50 maternal trees, 100 total seeds
for(i in 1:50) {
  same_25_10[i,1] = prop_capt_all_same[1,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[1,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[1,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[1,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[1,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[1,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("50 maternal trees, 100 seeds total") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.2,1) +
  theme(legend.position = "none")

#####################

#50 maternal trees, 200 total seeds
for(i in 1:50) {
  same_25_10[i,1] = prop_capt_all_same[2,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[2,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[2,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[2,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[2,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[2,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("50 maternal trees, 200 seeds total") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.2,1) +
  theme(legend.position = "none")

#########################

#25 maternal trees, 100 total seeds
for(i in 1:50) {
  same_25_10[i,1] = prop_capt_all_same[4,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[4,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[4,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[4,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[4,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[4,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("25 maternal trees, 100 seeds total") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.2,1) +
  theme(legend.position = "none")

######################

#25 maternal trees 200 seeds
for(i in 1:50) {
  same_25_10[i,1] = prop_capt_all_same[6,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[6,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[6,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[6,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[6,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[6,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("25 maternal trees, 200 seeds total") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.2,1) +
  theme(legend.position = "none")

#######################

#10 maternal trees 100 seeds
for(i in 1:50) {
  same_25_10[i,1] = prop_capt_all_same[11,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[11,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[11,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[11,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[11,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[11,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("10 maternal trees, 100 seeds total") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.2,1) +
  theme(legend.position = "none")

##########################

#10 maternal trees 200 seeds
for(i in 1:50) {
  same_25_10[i,1] = prop_capt_all_same[16,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[16,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[16,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[16,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[16,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[16,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("10 maternal trees, 200 seeds total") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.2,1) +
  theme(legend.position = "none")

#########################

#2 maternal trees 100 seeds
for(i in 1:50) {
  same_25_10[i,1] = prop_capt_all_same[66,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[66,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[66,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[66,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[66,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[66,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("2 maternal trees, 100 seeds total") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.2,1) +
  theme(legend.position = "none")

#########################

#2 maternal trees 200 seeds
for(i in 1:50) {
  same_25_10[i,1] = prop_capt_all_same[86,1,i] #this is hard coded to pull the scenario--having issues filtering data from a 3D array
  same_25_10[i,2] = prop_capt_all_same[86,5,i]
  
  eligible_25_10[i,1] = prop_capt_all_eligible[86,1,i]
  eligible_25_10[i,2] = prop_capt_all_eligible[86,5,i]
  
  skewed_25_10[i,1] = prop_capt_skewed[86,1,i]
  skewed_25_10[i,2] = prop_capt_skewed[86,5,i]
}

#data processing for plotting
equal_comparison = rbind(same_25_10, eligible_25_10, skewed_25_10)
equal_comparison = as.data.frame(equal_comparison)
colnames(equal_comparison) = c("prop_capt", "donor_type")
#plotting using a boxplot, includes all simulation replicates to note variation 
ggplot(data=equal_comparison, aes(x=(donor_type), y=as.numeric(prop_capt), group=donor_type, fill=donor_type)) +
  geom_boxplot() +
  ylab("Proportion of alleles captured") +
  xlab("Pollen donation type") +
  ggtitle("2 maternal trees, 200 seeds total") +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  ylim(0.2,1) +
  theme(legend.position = "none")

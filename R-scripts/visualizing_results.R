#load in results from main_loop.R
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")
load("prop_alleles_capt.Rdata")

#libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)

#arrays to save averaged results
prop_capt_all_same_avg = array(dim=c(465,5))
prop_capt_all_eligible_avg = array(dim=c(465,5))
prop_capt_skewed_avg = array(dim=c(465,5))
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

###############################################################
#First, looking into total seeds vs. prop alleles capt for each donor type

#filtering the data--looking at 10 maternal trees only
#when we try to plot all numbers of maternal trees on one plot, it gets too busy
#we need to have some constants
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

filtered = combined_data %>% filter(maternal_trees == 1)
ggplot(data=filtered, aes(x=as.numeric(total_seeds), y=as.numeric(prop_capt), group=donor_type, color=donor_type)) +
  geom_line() +
  xlab("Total seeds sampled") +
  ylab("Proportion of alleles captured") +
  ggtitle("Diversity captured from sampling 1 maternal trees")

################################################################
#Next looking into diversity captured vs. number of maternal trees sampled
#need to have a constant--total seeds sampled will be constant

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

#next looking at 50 seeds
filtered = combined_data %>% filter(total_seeds == 50)
ggplot(data=filtered, aes(x=as.numeric(maternal_trees), y=as.numeric(prop_capt), group=donor_type, color=donor_type)) +
  geom_line() +
  xlab("Number of maternal trees sampled") +
  ylab("Proportion of alleles captured") +
  ggtitle("Diversity captured from sampling 50 total seeds from varying maternal trees")

###########################################################
#comparing equivalent scenarios 
#comparing sampling 25 maternal trees 10 seeds from each tree (row 15)
same_25_10 = array(dim=c(500,2)) #500 rows for 500 replicates, 2 columns for prop_capt and donor type
eligible_25_10 = array(dim=c(500,2))
skewed_25_10 = array(dim=c(500,2))

#looping over simulation replicates to pull equivalent scenarios to comapre among donor modes 
#50 seeds 1 maternal tree
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

###################################################################################
#SKEWED
#comparing equivalent scenarios 
#comparing sampling 25 maternal trees 10 seeds from each tree (row 15)
same_25_10 = array(dim=c(50,2)) #50 rows for 50 replicates, 2 columns for prop_capt and donor type
eligible_25_10 = array(dim=c(50,2))
skewed_25_10 = array(dim=c(50,2))

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

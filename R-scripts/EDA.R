#Exploratory data analysis for one population simulation of ideal sampling scenarios 

#libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())

#load in data 
load("R-scripts/Rdata/alleles_capt_ideal_onepop.Rdata")

############################################################################################
#GET DATA IN PROPER FORMAT FOR PLOTTING AN ANALYSIS
#We need all the data which is currently in 3 separate 3D matrices to be combined into one large matrix in tidy format

#1. Converting 3D results matrices to 2D for plotting and analysis purposes
same_long=NULL
eligible_long=NULL
skewed_long=NULL
j=1
for(i in 1:50) {
    temp = rbind(prop_capt_all_same[,,j], prop_capt_all_same[,,(j+1)])
    same_long = rbind(same_long, temp)
    j = j+2
}
j=1
for(i in 1:50) {
    temp = rbind(prop_capt_all_eligible[,,j], prop_capt_all_eligible[,,(j+1)])
    eligible_long = rbind(eligible_long, temp)
    j = j+2
}
j=1
for(i in 1:50) {
    temp = rbind(prop_capt_skewed[,,j], prop_capt_skewed[,,(j+1)])
    skewed_long = rbind(skewed_long, temp)
    j = j+2
}
rm(temp)

#2. Binding all data frames together into one huge data frame 
#This is considered tidy format? 
tidy_df = rbind(same_long, eligible_long, skewed_long)
tidy_df = as.data.frame(tidy_df)

############################################################################################
#EXPLORATORY DATA ANALYSIS 

#Plot of total alleles present in simulations--to see the distribution of variation between simulation replicates 
ggplot(tidy_df, aes(x=total_alleles)) +
    geom_bar()

#Plotting curves 
tidy_df %>% 
    ggplot(aes(x=total_seeds, y=prop_capt, color=donor_type)) +
    geom_point(alpha=0.25) +
    geom_jitter() +
    theme(axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank()) +
    facet_wrap(vars(maternal_trees))

tidy_df %>% 
    filter(total_seeds==100) %>% 
    ggplot(aes(x=total_seeds, y=prop_capt, color=donor_type)) +
    geom_point(alpha=0.25) +
    geom_jitter() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    facet_wrap(vars(maternal_trees))
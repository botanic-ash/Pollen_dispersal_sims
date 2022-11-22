#Exploratory data analysis for one population simulation of ideal sampling scenarios 

#libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())

#load in data 
load("alleles_capt_ideal_onepop.Rdata")

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
#Just adding a column indicating the number of failures (number of alleles not captured) for binomial data
tidy_df$num_fails = c(as.numeric(tidy_df$total_alleles) - as.numeric(tidy_df$num_capt))

###########################################################################################
#EXPLORATORY DATA ANALYSIS 
#1. Examine the data
View(tidy_df)

#2. Visualizing the data 
#Plot of total alleles present in simulation replicates to see the distribution of variation between replicates 
#Each bar represents a different count of alleles
#The 3 increments on the plot show the allele count being present in one, two, or three different simulations by chance 
#However, most replicates model a different number of total alleles due to the stochastic simulations 
ggplot(tidy_df, aes(x=as.numeric(total_alleles))) +
    geom_bar() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

#Showing a table of all instances of total alleles present in simulations 
#Occurence values of 2805 = present in 1 independent simulation (935 scenarios x pollen donors = 1 simulation)
#Occurence values of 5610 = present in 2 independent simulations 
#Occurence values of 8415 = present in 3 independent simulations 
#Also, the total alleles simulated in simulations ranges from 235 to 288
#--this results in variation during sampling as well, for the same sample size of 50 seeds, 
#you may capture more diversity in the simulation with 235 alleles than 288 total alleles
table(tidy_df$total_alleles)

#Plotting all the data! 
#Plot proportion of alleles captured vs total number of seeds sampled for each number of maternal trees
#Here we see that in scenarios with fewer maternal trees sampled (facet 1, 2), there are greater differences in the proportion of alleles captured between pollen donor scenarios (large difference between the curves)
#In scenarios with more maternal trees sampled (facet 50, 100), similar proportion of alleles are captured across donor types (though there is less information here with fewer scenarios). 
#Additionally, comparing a given pollen donor type across varying number of trees sampled (compare a color across facets), we see that when more maternal trees are sampled, 
#much more diversity is captured (curve gets higher--greater proportion of alleles captured ).
#Lastly, when more seeds are sampled per tree (going along the x-axis), we see a slight increase in the diversity captured (see curves upward)
tidy_df %>% 
    ggplot(aes(x=as.numeric(total_seeds), y=as.numeric(prop_capt), color=donor_type)) +
    geom_point(alpha=0.25) +
    facet_wrap(vars(maternal_trees)) +
    ylim(0,1) +
    theme(axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank()) 

#Inspecting scenarios more closely--100 seeds total sampled, faceted by number of maternal trees sampled again and proportion of alleles captured on the y-axis (plotted with jitter to spread points across x-axis)
#Again, we see the trend of pollen donor type appearing to have more influence on the proportion of alleles captured in scenarios with fewer maternal trees sampled, since most of the diversity would be coming from pollen donors in these cases. 
tidy_df %>% 
    filter(total_seeds==100) %>% 
    ggplot(aes(x=as.numeric(total_seeds), y=as.numeric(prop_capt), color=donor_type)) +
    geom_point(alpha=0.25) +
    facet_wrap(vars(maternal_trees)) +
    ylim(0,1) +
    geom_jitter() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

#Boxplots of specific scenarios for each # of maternal tree sampled, 100 total seeds sampled
#Again, this shows the variation in alleles captured for each pollen donor type in scenarios with fewer maternal trees
#When many maternal trees are sampled, there is no difference in diversity captured 
tidy_df %>% 
    filter(total_seeds==100) %>% 
    ggplot(aes(x=as.numeric(total_seeds), y=as.numeric(prop_capt), color=donor_type)) +
    geom_boxplot(alpha=0.25) +
    facet_wrap(vars(maternal_trees)) +
    ylim(0,1) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 

#Plotting a scenario of 1 maternal tree sampled to better view the curve of the data 
tidy_df %>% 
    filter(maternal_trees==1) %>% 
    ggplot(aes(x=as.numeric(total_seeds), y=as.numeric(prop_capt), color=donor_type)) +
    geom_point(alpha=0.25) +
    facet_wrap(vars(maternal_trees)) +
    ylim(0,1) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
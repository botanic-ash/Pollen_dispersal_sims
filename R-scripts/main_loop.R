#Main loop
#This script has the main processing loop, which runs the functions defined in the previous script
#assinging mating, pollen dispersal, and seed collecting in the population
#Then, the results of sampling are saved

##################################################################################
#Library functions
library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(hierfstat)
library(poppr)
library(dplyr)

#including R-script containing functions used for import, conversions, and sampling
source("C:\\Users\\kayle\\Documents\\Pollen_dispersal_sims\\R-scripts\\import_seed_functions.R")

#defining the working directory containing simulation files
mydir = "C:\\Users\\kayle\\Documents\\Pollen_dispersal_sims\\Simulations\\two_pop_2500"
setwd(mydir)

#importing and converting arlequin files to genepop files
import_arp2gen_files(mydir,".arp$")

#importing and converting genepop files to genalex
import_gen2genalex_files(mydir, ".gen$")

##############################################################################################################
#DEFINING VARIABLES 

num_loci = 20 #number of loci simulated, needed to make a dataframe to save the data

total_seeds = 250 #total seeds to be sampled 

load("combined_list_params.Rdata") #loading in function parameters defined in defining_function_parameters.R script
#this Rdata file contains the three list for all_same, all_eligible, and skewed scenarios 

#defining array to store seeds that collectors have 'sampled'
#first we need to create column names depending on how many loci are present in simulations
#then, define the matrix, convert to dataframe, and rename the columns to label the data
#this dataframe keeps track of the alleles that are captured during sampling
loci_names = c()
for(i in 1:num_loci){
  loci_names = c(loci_names, paste("locus", i, "a", sep=""))
  loci_names = c(loci_names, paste("locus", i, "b", sep=""))
}

#creating a container to store the results of sampling and other important data
# four columns to save the proportion of alleles captured, number of seeds sampled,
#number of trees sampled, and number of pollen donors (4 cols)
#each row indicates the scenario (465)
#the third dimension is the simulation replicate (10)
prop_capt_all_same = array(dim=c(465,4,10))
prop_capt_all_eligible = array(dim=c(465,4,10))
prop_capt_skewed = array(dim=c(465,4,10))


#######################################################################################################
#Main processing loop
#first, import a single simulation replicate as genalex and convert to a dataframe 

#list of genalex files for all simulation replicates--genalex files end in .csv
genalex_list = list.files(mydir, ".csv$")

#for every simulation replicate, process data to be usable for the function, then call function
#finally, save results (prop. alleles capt, number seeds sampled, number trees sampled, and number pollen donors)
for(i in 1:length(genalex_list)) {
  #first import and process the data
  #import genalex files as dataframe
  data = read.csv(genalex_list[[i]], header=FALSE)
  #cut off first 2 rows in dataframe -- the population data is not required for this
  data = data[-2,]
  data = data[-1,]
  #giving the dataframe columns new names
  names(data) = c("Ind", "Pop", loci_names)
  data = data[-1,] #removing the first row -- repeat of now column headers
  
  #call sampling function here--save result in 3D matrix? (third dim. is for replicates?)
  #for each element in scenario list--for 'all same' sampling (defined in parameters.R file)
  for(x in 1:length(all_same_params)) {
    #call the function using that scenario and save data
    temp = sample_seed(data, all_same_params[[x]][[1]], all_same_params[[x]][[2]], all_same_params[[x]][[3]], all_same_params[[x]][[4]])
    #save these results--save proportion of alleles captured in temp
    #n_distinct counts all distinct values in a column -- we want to sum this across multiple columns, so that's where sum and sapply come in
    captured = sum(sapply(temp[1:40], n_distinct))
    total = sum(sapply(data[3:42], n_distinct))
    prop_capt_all_same[x,1,i] = (captured/total)#proportion of alleles captured
    prop_capt_all_same[x,2,i] = ((all_same_params[[x]][[1]])*(all_same_params[[x]][[2]][[1]]))#total seeds sampled --if possible, change all_same_params[[x]][[2]][[1]] hard coding 
    prop_capt_all_same[x,3,i] = (all_same_params[[x]][[1]]) #number of trees sampled 
    prop_capt_all_same[x,4,i] = (all_same_params[[x]][[3]]) #number of pollen donors
  }
  
  #'all eligible' sampling
  for(x in 1:length(all_eligible_params)) {
    temp = sample_seed(data, all_eligible_params[[x]][[1]], all_eligible_params[[x]][[2]], all_eligible_params[[x]][[3]], all_eligible_params[[x]][[4]])
    #save these results--save proportion of alleles captured in temp
    #n_distinct counts all distinct values in a column -- we want to sum this across multiple columns, so that's where sum and sapply come in
    captured = sum(sapply(temp[1:40], n_distinct)) #hard coded --try to remove this
    total = sum(sapply(data[3:42], n_distinct)) #hard coded -- try to remove this
    prop_capt_all_eligible[x,1,i] = (captured/total) #proportion of alleles captured
    prop_capt_all_eligible[x,2,i] = ((all_eligible_params[[x]][[1]])*(all_eligible_params[[x]][[2]][[1]]))#total seeds sampled
    prop_capt_all_eligible[x,3,i] = (all_eligible_params[[x]][[1]]) #number of trees sampled 
    prop_capt_all_eligible[x,4,i] = (all_eligible_params[[x]][[3]]) #number of pollen donors
  }
  
  for(x in 1:length(skewed_params)) {
    temp = sample_seed(data, skewed_params[[x]][[1]], skewed_params[[x]][[2]], skewed_params[[x]][[3]], skewed_params[[x]][[4]])
    #save these results--save proportion of alleles captured in temp
    #n_distinct counts all distinct values in a column -- we want to sum this across multiple columns, so that's where sum and sapply come in
    captured = sum(sapply(temp[1:40], n_distinct)) #hard coded --try to remove this
    total = sum(sapply(data[3:42], n_distinct)) #hard coded -- try to remove this
    prop_capt_skewed[x,1,i] = (captured/total)#proportion of alleles captured
    prop_capt_skewed[x,2,i] = ((skewed_params[[x]][[1]])*(skewed_params[[x]][[2]][[1]]))#total seeds sampled
    prop_capt_skewed[x,3,i] = (skewed_params[[x]][[1]]) #number of trees sampled 
    prop_capt_skewed[x,4,i] = (skewed_params[[x]][[3]]) #number of pollen donors
  }
  
}

#saving results to Rdata file
save(prop_capt_all_same, prop_capt_all_eligible, prop_capt_skewed, file="prop_alleles_capt.Rdata")
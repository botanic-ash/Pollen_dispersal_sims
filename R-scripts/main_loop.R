#Code written by Kaylee Rosenberger
#Main loop
#This script has the main processing loop, which runs the functions defined in the previous script
#assinging mating, pollen dispersal, and seed collecting in the population
#Then, the results of sampling are saved

##################################################################################
#Library functions
library(adegenet)
library(diveRsity)
library(poppr)
library(tidyr)
library(dplyr)

setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")
###MAKE SURE TO LOAD IN THE CORRECT DATA
load("combined_list_params.Rdata") #loading in function parameters defined in defining_function_parameters.R script
#this Rdata file contains the three list for all_same, all_eligible, and skewed scenarios 

#including R-script containing functions used for import, conversions, and sampling
source("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts/import_seed_functions.R")
#including edited arp2gen function 
source("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts/arp2gen_edit.R")

#defining the working directory containing simulation files
mydir = "C:/Users/kayle/Documents/Pollen_dispersal_sims/Simulations/one_pop_2500"
setwd(mydir)

#importing and converting arlequin files to genepop files
import_arp2gen_files(mydir,".arp$")

#importing and converting genepop files to genalex
import_gen2genalex_files(mydir, ".gen$")

##############################################################################################################
#DEFINING VARIABLES 

num_loci = 20 #number of loci simulated, needed to make a dataframe to save the data

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
# five columns to save the proportion of alleles captured, number of seeds sampled,
#number of trees sampled, and number of pollen donors, and the donor type (5 cols)
#each row indicates the scenario (465)
#the third dimension is the simulation replicate (50)
#saving each of the pollen donor scenarios in different arrays--we can combine them later if we need to! 
  #this is just easier because it's less filtering and more organized
#Idealized scenarios
prop_capt_all_same = array(dim=c(215,5,50))
prop_capt_all_eligible = array(dim=c(215,5,50))
prop_capt_skewed = array(dim=c(215,5,50))

#realistic scenarios--we have different numbers of scenarios for each
#prop_capt_all_same = array(dim=c(215,5,50))
#prop_capt_all_eligible = array(dim=c(215,5,50))
#prop_capt_skewed = array(dim=c(215,5,50))


#######################################################################################################
#Main processing loop
#goes through each simulation replicate, calling the function with the parameters defined in the 
#defining_function_parameters.R script.
#results are saved in a 3D array 

#list of genalex files for all simulation replicates--genalex files end in .csv
#these have the simulated genetic data
genalex_list = list.files(mydir, ".csv$")
setwd(mydir)

#Main loop overview:
  #for every simulation replicate, process data to be usable for the function
  #then, we have three separate loops, that loop over the parameters created in defining_function_parameters.R
  #there are three separate loops for the three pollen donor scnearios (skewed, all eligble, and all same)
  #calculate proportion of alleles captured by 
  #finally, save results (prop. alleles capt, number seeds sampled, number trees sampled, and number pollen donors)
for(i in 1:length(genalex_list)) {
  print(paste("replicate number ", i))
  #first import and process the data
  #import genalex files as dataframe
  data = read.csv(genalex_list[[i]], header=FALSE)
  #cut off first 2 rows in dataframe -- this is the population data, which is not required for our purposes
  data = data[-2,]
  data = data[-1,]
  #giving the dataframe columns new names
  names(data) = c("Ind", "Pop", loci_names)
  data = data[-1,] #removing the first row -- repeat of now column headers
  
  #All same pollen scenarios
  #for each element in scenario list--for 'all same' sampling (defined in defining_function_parameters.R file)
  for(x in 1:length(all_same_params)) {
    #call the function using that scenario and save the function return in temp
    temp = sample_seed(data, all_same_params[[x]][[1]], all_same_params[[x]][[2]], all_same_params[[x]][[3]], all_same_params[[x]][[4]])
    #calculating proportion of alleles captured 
    #Counting the total number of alleles present in the parental dataset (total)
    total = 0 #sum to keep track of total alleles 
    k=3
    for(i in 1:num_loci) {
      allele_list = unique(data[,k:k+1]) #getting unique values in column k for locus i (locusiA)
      total = total + n_distinct(allele_list) #getting the number of distinct values for locus 1 to count alleles 
      k = k+2 #increment k for next loop iteration
    }
    #counting the number of alleles in the seeds sampled dataset (captured)
    captured = 0
    k=1
    for(i in 1:num_loci) {
      allele_list = unique(temp[,k:k+1]) #getting unique values in column k for locus i (locusiA)
      captured = captured + n_distinct(allele_list) #getting the number of distinct values for locus 1 to count alleles 
      k = k+2 #increment k for next loop iteration
    }
    prop_capt_all_same[x,1,i] = (captured/total)#proportion of alleles captured = captured/total, save these results
    prop_capt_all_same[x,2,i] = sum(all_same_params[[x]][[2]])#total seeds sampled --if possible, change all_same_params[[x]][[2]][[1]] hard coding 
    prop_capt_all_same[x,3,i] = (all_same_params[[x]][[1]]) #number of trees sampled 
    prop_capt_all_same[x,4,i] = (all_same_params[[x]][[3]]) #number of pollen donors
    prop_capt_all_same[x,5,i] = "all_same"
  }
  
  #'all eligible' sampling
  for(x in 1:length(all_eligible_params)) {
    #call the function using that scenario and save the function return in temp
    temp = sample_seed(data, all_eligible_params[[x]][[1]], all_eligible_params[[x]][[2]], all_eligible_params[[x]][[3]], all_eligible_params[[x]][[4]])
    #calculating proportion of alleles captured 
    #calculating proportion of alleles captured 
    #Counting the total number of alleles present in the parental dataset (total)
    total = 0 #sum to keep track of total alleles 
    k=3 #starting from 3rd column
    for(i in 1:num_loci) {
      allele_list = unique(data[,k:k+1]) #getting unique values in column k for locus i (locusiA)
      total = total + n_distinct(allele_list) #getting the number of distinct values for locus 1 to count alleles 
      k = k+2 #increment k for next loop iteration
    }
    #counting the number of alleles in the seeds sampled dataset (captured)
    captured = 0
    k=1
    for(i in 1:num_loci) {
      allele_list = unique(temp[,k:k+1]) #getting unique values in column k for locus i (locusiA)
      captured = captured + n_distinct(allele_list) #getting the number of distinct values for locus 1 to count alleles 
      k = k+2 #increment k for next loop iteration
    }
    prop_capt_all_eligible[x,1,i] = (captured/total) #proportion of alleles captured= captured/total
    prop_capt_all_eligible[x,2,i] = sum(all_eligible_params[[x]][[2]])#total seeds sampled
    prop_capt_all_eligible[x,3,i] = (all_eligible_params[[x]][[1]]) #number of trees sampled 
    prop_capt_all_eligible[x,4,i] = (all_eligible_params[[x]][[3]]) #number of pollen donors
    prop_capt_all_eligible[x,5,i] = "all_eligible"
  }
  
  #skewed sampling
  for(x in 1:length(skewed_params)) {
    temp = sample_seed(data, skewed_params[[x]][[1]], skewed_params[[x]][[2]], skewed_params[[x]][[3]], skewed_params[[x]][[4]])
    #calculating proportion of alleles captured 
    #calculating proportion of alleles captured 
    #Counting the total number of alleles present in the parental dataset (total)
    total = 0 #sum to keep track of total alleles 
    k=3 #starting from 3rd column
    for(i in 1:num_loci) {
      allele_list = unique(data[,k:k+1]) #getting unique values in column k for locus i (locusiA)
      total = total + n_distinct(allele_list) #getting the number of distinct values for locus 1 to count alleles 
      k = k+2 #increment k for next loop iteration
    }
    #counting the number of alleles in the seeds sampled dataset (captured)
    captured = 0
    k=1
    for(i in 1:num_loci) {
      allele_list = unique(temp[,k:k+1]) #getting unique values in column k for locus i (locusiA)
      captured = captured + n_distinct(allele_list) #getting the number of distinct values for locus 1 to count alleles 
      k = k+2 #increment k for next loop iteration
    }
    prop_capt_skewed[x,1,i] = (captured/total)#proportion of alleles captured = captured/ total
    prop_capt_skewed[x,2,i] = sum(skewed_params[[x]][[2]])#total seeds sampled
    prop_capt_skewed[x,3,i] = (skewed_params[[x]][[1]]) #number of trees sampled 
    prop_capt_skewed[x,4,i] = (skewed_params[[x]][[3]]) #number of pollen donors
    prop_capt_skewed[x,5,i] = "skewed"
  }
  
}

colnames(prop_capt_all_same) = c("prop_capt", "total_seeds", "maternal_trees", "num_donors", "donor_type")
colnames(prop_capt_all_eligible) = c("prop_capt", "total_seeds", "maternal_trees", "num_donors", "donor_type")
colnames(prop_capt_skewed) = c("prop_capt", "total_seeds", "maternal_trees", "num_donors", "donor_type")

#saving EQUAL results to Rdata file
#setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")
#save(prop_capt_all_same, prop_capt_all_eligible, prop_capt_skewed, file="prop_alleles_capt_new.Rdata")

#saving SKEWED results to Rdata file
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")
save(prop_capt_all_same, prop_capt_all_eligible, prop_capt_skewed, file="prop_alleles_capt_skewed_new.Rdata")

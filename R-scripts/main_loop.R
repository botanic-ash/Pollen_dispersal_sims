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

setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims")

###MAKE SURE TO LOAD IN THE CORRECT DATA
#load("R-scripts/combined_list_params_ideal.Rdata") #loading in function parameters defined in defining_function_parameters.R script
load("R-scripts/combined_list_params_realistic.Rdata")
#this Rdata file contains the three list for all_same, all_eligible, and skewed scenarios 

#including R-script containing functions used for import, conversions, and sampling
source("R-scripts/import_seed_functions.R")
#including edited arp2gen function 
source("R-scripts/arp2gen_edit.R")

#Defining this because arp2gen didn't want to use relative file paths 
mydir = "C:/Users/kayle/Documents/Pollen_dispersal_sims/Simulations/one_pop_2500/"

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
# five columns to save the proportion of alleles captured, 
#number of seeds sampled, number of trees sampled, pollen donor type (6 cols),
#number of alleles captured, total and number of alleles present, and simulation replicate 
#each row indicates the scenario (465)
#the third dimension is the simulation replicate (50)
#saving each of the pollen donor scenarios in different arrays--we can combine them later if we need to! 
  #this is just easier because it's less filtering and more organized
#Idealized scenarios
#prop_capt_all_same = array(dim=c(935,7,50))
#prop_capt_all_eligible = array(dim=c(935,7,50))
#prop_capt_skewed = array(dim=c(935,7,50))

#realistic scenarios--we have different numbers of scenarios for each
prop_capt_all_same = array(dim=c(217,7,50))
prop_capt_all_eligible = array(dim=c(217,7,50))
prop_capt_skewed = array(dim=c(217,7,50))


#######################################################################################################
#Main processing loop
#goes through each simulation replicate, calling the function with the parameters defined in the 
#defining_function_parameters.R script.
#results are saved in a 3D array 

#list of genalex files for all simulation replicates--genalex files end in .csv
#these have the simulated genetic data
genalex_list = list.files(mydir, ".csv$")

#Main loop overview:
  #for every simulation replicate, process data to get in the format required for the function sample_seed()
  #then, we have three separate loops that loop over the parameters created in defining_function_parameters.R for each pollen donor type
  #there are three separate loops for the three pollen donor scnearios (skewed, all eligble, and all same)
  #calculate proportion of alleles captured by dividing number of alleles sampled / total alleles present in the simulation
  #finally, save results (prop. alleles capt, number seeds sampled, number trees sampled, and number pollen donors)
for(i in 1:length(genalex_list)) {
  print(paste("replicate number ", i))
  #first import and process the data
  #import genalex files as dataframe
  genetic_data = read.csv(paste(mydir, genalex_list[[i]], sep=""), header=FALSE)
  #cut off first 2 rows in dataframe -- this is the population data, which is not required for our purposes
  genetic_data = genetic_data[-2,]
  genetic_data = genetic_data[-1,]
  #giving the dataframe columns new names
  names(genetic_data) = c("Ind", "Pop", loci_names)
  genetic_data = genetic_data[-1,] #removing the first row -- repeat of now column headers
  
  #calculating the number of alleles present in the parental dataset--we only need to do this once for each replicate
  total = 0  #sum to keep track of total alleles
  k=3 #counter variable for column (locus) of parental dataset
  for(z in 1:num_loci){
    parental_allele_list = table(c(as.matrix(genetic_data[,k:(k+1)]))) #getting alleles and their frequencies for locus i in parental dataframe
    #parental_allele_list = parental_allele_list[parental_allele_list>3] #subsetting parental data to only include alleles with frequency greater than 3
    total_names = names(parental_allele_list)
    total = total + n_distinct(total_names) #getting the number of distinct values for locus 1 to count alleles 
    
    k = k+2 #increment k for next loop iteration
  }
  
  #All same pollen scenarios
  #for each element in scenario list--for 'all same' sampling (defined in defining_function_parameters.R file)
  for(x in 1:length(all_same_params)) {
    #call the function using that scenario and save the function return in temp
    temp = sample_seed(genetic_data, all_same_params[[x]][[1]], all_same_params[[x]][[2]], all_same_params[[x]][[3]], all_same_params[[x]][[4]])
    #calculating proportion of alleles captured 
    captured = 0 #sum to keep track of alleles captured by sampling
    j=1 #counter variable to column (locus) of seed dataset 
    for(z in 1:num_loci) {
      seed_allele_list = table(c(as.matrix(temp[,j:(j+1)]))) #getting unique values for locus i in seed dataframe 
      captured_names = names(seed_allele_list)#getting the names of the alleles captured from sampling
      captured = captured + n_distinct(captured_names)#making sure none of the super rare alleles excluded from parental dataset are included here
      
      j = j+2 #increment j for next loop iteration 
    }
    prop_capt_all_same[x,1,i] = (captured/total)#proportion of alleles captured = captured/total, save these results
    prop_capt_all_same[x,2,i] = sum(all_same_params[[x]][[2]])#total seeds sampled --if possible, change all_same_params[[x]][[2]][[1]] hard coding 
    prop_capt_all_same[x,3,i] = (all_same_params[[x]][[1]]) #number of trees sampled 
    prop_capt_all_same[x,5,i] = captured
    #clearing out containers
    rm(captured_names)
    rm(seed_allele_list)
    rm(temp)
  }
  #defining these columns with vectors to speed up code 
  prop_capt_all_same[,4,i] = "all_same"
  prop_capt_all_same[,6,i] = total
  prop_capt_all_same[,7,i] = i
  
  #'all eligible' sampling
  for(x in 1:length(all_eligible_params)) {
    #call the function using that scenario and save the function return in temp
    temp = sample_seed(genetic_data, all_eligible_params[[x]][[1]], all_eligible_params[[x]][[2]], all_eligible_params[[x]][[3]], all_eligible_params[[x]][[4]])
    #calculating proportion of alleles captured 
    #calculating proportion of alleles captured 
    captured = 0 #sum to keep track of alleles captured by sampling
    j=1 #counter variable to column (locus) of seed dataset 
    for(z in 1:num_loci) {
      seed_allele_list = table(c(as.matrix(temp[,j:(j+1)]))) #getting unique values for locus i in seed dataframe 
      captured_names = names(seed_allele_list)#getting the names of the alleles captured from sampling
      captured = captured + n_distinct(captured_names)#making sure none of the super rare alleles excluded from parental dataset are included here
      
      j = j+2 #increment j for next loop iteration 
    }
    prop_capt_all_eligible[x,1,i] = (captured/total) #proportion of alleles captured= captured/total
    prop_capt_all_eligible[x,2,i] = sum(all_eligible_params[[x]][[2]])#total seeds sampled
    prop_capt_all_eligible[x,3,i] = (all_eligible_params[[x]][[1]]) #number of trees sampled 
    prop_capt_all_eligible[x,5,i] = captured
    #clearing out containers
    rm(captured_names)
    rm(seed_allele_list)
    rm(temp)
  }
  prop_capt_all_eligible[,4,i] = "all_eligible"
  prop_capt_all_eligible[,6,i] = total
  prop_capt_all_eligible[,7,i] = i
  
  
  #skewed sampling
  for(x in 1:length(skewed_params)) {
    temp = sample_seed(genetic_data, skewed_params[[x]][[1]], skewed_params[[x]][[2]], skewed_params[[x]][[3]], skewed_params[[x]][[4]])
    #calculating proportion of alleles captured 
    captured = 0 #sum to keep track of alleles captured by sampling
    j=1 #counter variable to column (locus) of seed dataset 
    for(z in 1:num_loci) {
      seed_allele_list = table(c(as.matrix(temp[,j:(j+1)]))) #getting unique values for locus i in seed dataframe 
      captured_names = names(seed_allele_list)#getting the names of the alleles captured from sampling
      captured = captured + n_distinct(captured_names)#making sure none of the super rare alleles excluded from parental dataset are included here
      
      j = j+2 #increment j for next loop iteration 
    }
    prop_capt_skewed[x,1,i] = (captured/total)#proportion of alleles captured = captured/ total
    prop_capt_skewed[x,2,i] = sum(skewed_params[[x]][[2]])#total seeds sampled
    prop_capt_skewed[x,3,i] = (skewed_params[[x]][[1]]) #number of trees sampled 
    prop_capt_skewed[x,5,i] = captured
    #clearing out containers
    rm(captured_names)
    rm(seed_allele_list)
    rm(temp)
  }
  prop_capt_skewed[,4,i] = "skewed"
  prop_capt_skewed[,6,i] = total
  prop_capt_skewed[,7,i] = i
  
  
  rm(total_names)
  rm(parental_allele_list)
}

colnames(prop_capt_all_same) = c("prop_capt", "total_seeds", "maternal_trees", "donor_type", "num_capt", "total_alleles", "replicate")
colnames(prop_capt_all_eligible) = c("prop_capt", "total_seeds", "maternal_trees",  "donor_type", "num_capt", "total_alleles", "replicate")
colnames(prop_capt_skewed) = c("prop_capt", "total_seeds", "maternal_trees", "donor_type", "captured", "total_alleles", "replicate")

#saving ideal results to Rdata file
#save(prop_capt_all_same, prop_capt_all_eligible, prop_capt_skewed, file="../../R-scripts/Rdata/alleles_capt_ideal_onepop.Rdata")

#saving realistic results to Rdata file
save(prop_capt_all_same, prop_capt_all_eligible, prop_capt_skewed, file="../../R-scripts/Rdata/alleles_capt_realistic_onepop.Rdata")

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

load("combined_list_params.Rdata")

#defining array to store seeds that collectors have 'sampled'
#first we need to create column names depending on how many loci are present in simulations
#then, define the matrix, convert to dataframe, and rename the columns to label the data
#this dataframe keeps track of the alleles that are captured during sampling
loci_names = c()
for(i in 1:num_loci){
  loci_names = c(loci_names, paste("locus", i, "a", sep=""))
  loci_names = c(loci_names, paste("locus", i, "b", sep=""))
}

seeds_sampled = matrix(nrow = total_seeds, ncol = ((2*num_loci)+1))
seeds_sampled = as.data.frame(seeds_sampled)
names(seeds_sampled) = c(loci_names, "Pop")


#######################################################################################################
#Main processing loop
#first, import a single simulation replicate as genalex and convert to a dataframe 

#list of genalex files for all simulation replicates--genalex files end in .csv
genalex_list = list.files(mydir, ".csv$")

#for every simulation replicate, process data, call function
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
  #for each element in scenario list (defined in parameters.R file)
    #call the function using that scenario and save data
}

#call sample seed function
#args: data, num trees, num seeds, num donors, probability
#x = sample_seed(data, trees_to_sample, seeds_to_sample, num_pollen_donors, pollen_probability)
temp = sample_seed(data, combined_list_params[[1]][[1]], combined_list_params[[1]][[2]], combined_list_params[[1]][[3]], combined_list_params[[1]][[4]])
#third parameter is currently causing issues--assuming the next two will as well 



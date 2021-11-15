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
source("C:\\Users\\kayle\\Documents\\Pollen_dispersal_sims\\R-scripts\\defining_functions.R")

#defining the working directory containing simulation files
mydir = "C:\\Users\\kayle\\Documents\\Pollen_dispersal_sims\\Simulations\\example_population"
setwd(mydir)

#importing and converting arlequin files to genepop files
import_arp2gen_files(mydir,".arp$")

#importing and converting genepop files to genalex
import_gen2genalex_files(mydir, ".gen$")

##############################################################################################################
#DEFINING VARIABLES 

#number of loci modeled in simulations
num_loci = 1

#number of trees in population
total_trees = 100
#number of trees in population that collectors will sample from
trees_to_sample = 10
#number of seeds collectors will sample from each tree
seeds_to_sample = 10

#number of pollen donors per maternal tree
num_pollen_donors = 1 #defined as 1 for simplicity
pollen_probability = c(1) #defined as 1 since there is 1 father per mother
#can change this vector as number of fathers increases, will be proportions, eg., 0.2, 0.2, 0.6 (but should add up to 1!)

#defining array to store seeds that collectors have 'sampled'
#first we need to create column names depending on how many loci are present in simulations
#then, define the matrix, convert to dataframe, and rename the columns to label the data
#this dataframe keeps track of the alleles that are captured during sampling
loci_names = c()
for(i in 1:num_loci){
  loci_names = c(loci_names, paste("locus", i, "a", sep=""))
  loci_names = c(loci_names, paste("locus", i, "b", sep=""))
}

seeds_sampled = matrix(nrow = total_trees, ncol = (2*num_loci))
seeds_sampled = as.data.frame(seeds_sampled)
names(seeds_sampled) = c(loci_names)


#######################################################################################################
#Main processing loop
#first, import a single simulation replicate as genalex and convert to a dataframe 

#list of genalex files for all simulation replicates
genalex_list = list.files(mydir, ".csv$")
#for number of simulation replicates
for(i in 1:length(genalex_list)) {
  #first import and process the data
  #import genalex files as dataframe
  data = read.csv(genalex_list[[i]], header=FALSE)
  #cut off first 2 rows in dataframe -- the population data is not required for this
  data = data[-2,]
  data = data[-1,]
  data = data[,-5] #getting rid of the empty column
  #giving the dataframe columns new names
  names(data) = c("Ind", "Pop", loci_names)
  data = data[-1,] #removing the first row -- repeat of now column headers
  
  #call sampling function here--save result in 3D matrix? (third dim. is for replicates?)
}
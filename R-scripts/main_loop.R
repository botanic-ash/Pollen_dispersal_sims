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

#working directory containing genalex files
mydir = "C:\\Users\\kayle\\Documents\\Pollen_dispersal_sims\\Simulations\\example_population"
setwd(mydir)

#number of loci simulated
num_loci = 1

#list of genalex files for all simulation replicates
genalex_list = list.files(mydir, ".csv$")

#######################################################################################################
#Main processing loop
#first, import a single simulation replicate as genalex and convert to a dataframe 

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
}
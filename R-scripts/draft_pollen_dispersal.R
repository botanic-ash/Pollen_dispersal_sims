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

mydir = "C:\\Users\\kayle\\Documents\\Pollen_dispersal_sims\\Simulations\\example_population"
setwd(mydir)
##################################################################################
#FILE IMPORT AND CONVERSIONS

#Defining an import function
#converts all arlequin simulation files in a folder to genepop files
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(z in 1:length(temp_list_1)){temp_list_2[[z]]=arp2gen(temp_list_1[z])}
  temp_list_2
} 

#importing and converting files
import_arp2gen_files(mydir,".arp$")

##############################################################################################################
#make a list of genind objects here for each simulation replicate

#just using one simulation replicate for now.. 
genind_obj = read.genepop("example_population_0.gen", ncode=3)
#since we will mostly be accessing the genind object in tab, defining this for simplicity 
gen_data = genind_obj@tab
##############################################################################################################
#DEFINING VARIABLES 

#number of trees in population
total_trees = 100
#number of seeds carried by each tree
#a simplification, maybe a place holder? 
seeds_per_tree = 100
#array to store the total seed population
#array dimensions are total number of mothers x total number of seeds per tree
seed_population = array(dim=c(total_trees, seeds_per_tree))

#number of trees that collectors will sample from
trees_to_sample = 10
#number of seeds collectors will sample on each tree
seeds_to_sample = 10
#array to store seeds that collectors have 'sampled'
#dimensions are number of trees sampled from x number of seeds per tree
seeds_sampled = array(dim=c(trees_to_sample, seeds_to_sample))

#number of pollen donors per maternal tree
num_pollen_donors = 1 #defined as 1 for simplicity
pollen_probability = c(1) #defined as 1 since there is 1 father per mother -- can change this as number of fathers increases, will be a proportion

#############################################################################################################
#Defniing the function to create sets of seeds and sample them

# move function to another script once complete with testing
make_seed = function() {
  for(i in 1:total_trees) {
    #choose mother -- a single row representing one maternal tree
    #choose father(s) randomly -- a vector of rows representing potential pollen donors  
    for(j in 1:seeds_per_tree) {
      #choose maternal allele randomly -- need to sample() non-zero values only 
      #choose father based on probability vector 
        #then choose allele from that father randomly -- need to sample() non-zero values only 
      #save maternal allele + paternal allele as a seed
      #EX.) seed_array[i,j] = list(maternal_allele, paternal_allele)
    }
  }
}


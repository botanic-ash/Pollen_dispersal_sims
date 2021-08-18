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

#################################################################################
#make a list of genind objects here for each simulation replicate

#just using one simulation replicate for now.. 
genind_obj = read.genepop("example_population_0.gen", ncode=3)
gen_data = genind_obj@tab
#################################################################################
#DEFINING VARIABLES 

#how many trees in entire population - required to build all seed set (potential trees to sample from)
total_trees = 100

#defining an array to store seeds on mother trees 
all_seeds = array(c(100,100)) #size = size of mothers x number of seeds they carry (simplified)

#number of trees collectors will sample seed from 
n = 10
#number of seeds collectors will sample from each tree 
m = 10
#array to store seeds that have been collected 
collected_seeds = array(c(n,m)) #size = number of trees sampled x number of seeds sampled per tree 

#number of pollen donors
#changes depending on the scenario 
num_fathers = 1 

#defining a pollen probability vector 
#changes depending on the scenario 
pollenation_probability = vector(length = num_fathers)
pollenation_probability[1] = 1
#could find some way to automate this, so we don't have to edit everytime

#rows to sample
#defining which rows or 'individuals' to sample from in the genind object 
#(rows correspond to individuals)
rows_to_sample_mothers = n
rows_to_sample_fathers = num_fathers

################################################################################
#Defniing the function to create sets of seeds and sample them
# move function to another script once done 

build_seed = function(num_mothers, num_fathers, pollen_probability_vector, num_seeds) {
      for(i in 1:n) {
        rows_to_sample_mothers = sample(total_trees, 1) #sample from the entire dataset, 1 mother (repeat n times for n number mothers)
        for(j in 1:num_fathers) { #choosing which individuals to assign as fathers, for a given mother (one mother tree can have multiple pollen donors)
          rows_to_sample_fathers[j] = sample(total_trees, num_fathers)
        }
        m_alleles = sample(gen_data[rows_to_sample_mothers,], 1)#error here, it can select an allele that the individual doesn't have
        #need to filter out data = 0 somehow 
        p_alleles = sample(gen_data[rows_to_sample_fathers,], 1)
        seed_genotype = list(m_allele, p_allele)
      }
}

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
pollen_probability = c(1) #defined as 1 since there is 1 father per mother
#can change this vector as number of fathers increases, will be proportions, eg., 0.2, 0.2, 0.6 (but should add up to 1!)

#############################################################################################################
#Defniing the function to create sets of seeds and sample them

# move functions to another script once complete with testing

#
#make_seed()
#args:
#total_trees: total trees in population that have seeds
#seeds_per_tree: number of seeds present on each tree (a simplification)
#num_pollen_donors: the number of trees contributing pollen to a single mother tree

#Function creates sets of seeds for all trees in the population
#A seed is created by choosing a mother tree and a father tree
#(or sets of father trees), and selecting alleles randomly 
#from each parent for the seed to 'inherit'
#If there are multiple father trees (pollen donors), then
#The fathers are selected, and each one donates pollen (alleles)
#to seeds based on a probability vector
#Note: all trees are hermaphroditic and have the potential to self,
#since trees are selected randomly with sample() 
make_seed = function(total_trees, seeds_per_tree, num_pollen_donors) {
  for(i in 1:total_trees) {
    #choose mother -- a single row representing one maternal tree
    #Ex.) x <- sample(total_trees, 1)
    #choose father(s) randomly -- a vector of rows representing potential pollen donors  
    #Ex.) vec <- sample(total_trees, num_pollen_donors)
    for(j in 1:seeds_per_tree) {
      #choose maternal allele randomly -- need to sample() non-zero values only 
      #Ex.) sample(gen_data[x,], 1) -- figure out how to filter data that is == 0
      #choose father based on probability vector 
      #Ex.) y <- sample(vec, 1, prob = pollen_probability)
        #then choose allele from that father randomly -- need to sample() non-zero values only 
        #Ex.) sample(gen_data[y,], 1)
      #save the genotype of the seed (maternal allele + paternal allele)
      #EX.) seed_array[i,j] = list(maternal_allele, paternal_allele)
    }
  }
  #return seed_array
}
#Note: currently we are using just 1 locus for simplicity (modeled in Simcoal)
#it will be more difficult to implement multiple alleles 



#
#sample_seed()
#args: 
#trees_to_sample: number of trees collectors will sample seed from
#seeds_to_sample: number of seeds collectors will sample from each tree

#Function samples seeds from the array created in the function above
#This represents collectors sampling seeds from trees within the population
sample_seed = function(trees_to_sample, seeds_to_sample) {
  for(i in 1:trees_to_sample) { #for each tree...
    #sample a row from seed_array -- corresponding to a mother
    #Ex.) sample(total_trees, 1) -> row 25
    for(j in 1:seeds_to_sample) { #sample this many seeds
      #sample a 'seed' from the row
      #Ex.) sample(seeds_per_tree, 1) -> column 60
      #save the seed's genotype to seeds_sampled array
      #Ex.) seeds_sampled[25,60] = seed_array[25,60]
    }
  }
  #return seeds_sampled
}


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

mydir = "C:\\Users\\kayle\\Documents\\Pollen_dispersal_sims\\Simulations\\example_population"
setwd(mydir)
##################################################################################
#FILE IMPORT AND CONVERSIONS

#Defining an import function
#converts all arlequin simulation files in a folder to genepop files
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath) #directory containing files
  temp_list_1 = list.files(mypath, mypattern) #list of files with .arp extension
  temp_list_2 = list(length = length(temp_list_1)) #empty list with the same length as the first list
  for(z in 1:length(temp_list_1)){temp_list_2[[z]]=arp2gen(temp_list_1[z])} # for every item in list 1, convert it to .gen file and save to list 2
  temp_list_2
} 

#importing and converting files
import_arp2gen_files(mydir,".arp$")

#creating a function to convert genepop files to genalex files
#there 2 main steps--first you need to convert genepop files to genind objects
#then you convert genind objects to genalex files
import_gen2genalex_files = function(mypath, mypattern) {
  setwd(mypath) #directory containing files
  temp_list_1 = list.files(mypath, mypattern) #list of all files with .gen extension
  temp_list_2 = list(length = length(temp_list_1)) #empty list with same length as the first
  temp_list_3 = list(length = length(temp_list_1))#another empty list with same length as the first
  for(z in 1:length(temp_list_1)){
    temp_list_2[[z]]=read.genepop(temp_list_1[z], ncode = 3)#converting each genepop file to a genind object
    temp_list_3[[z]]=genind2genalex(temp_list_2[[z]], filename = (paste("example_population_",z,".csv", sep="")), overwrite = TRUE) #converting each genind object to a genalex file
  }
  temp_list_3 #retun
}


#importing and converting files to genalex
import_gen2genalex_files(mydir, ".gen$")

#NOTE: Will need to import genalex files as a dataframe (likely inside the main processing loop)
#also, will need to do some processing on the dataframe (cut out first 2 rows -- population data)

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

#array to store seeds that collectors have 'sampled'
#dimensions are number of trees sampled from x number of seeds per tree
seeds_sampled = array(dim=c(trees_to_sample, seeds_to_sample))

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
sample_seed_0 = function(trees_to_sample, seeds_to_sample) {
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

##########################################################################################
#New merged function--combining make_seed and sample_seed() from above
#Function creates seeds for collectors to sample from trees. First, a maternal tree is defined. Then father tree(s) is/are 
#defined. From there, alleles are selected from the mother and father randomly to create seeds.
#arugments:
  #num_trees_to_sample: number of trees in the population collectors will sample seed from 
  # num_seeds_to_sample: number of seeds collectors will take from each tree
  #num_pollen_donors: the number of trees that donate pollen to a single mother tree
  #pollen_probability: the probability of pollen donors to pollinate a mother tree (it can be equal for all donors, or skewed)
  #num_loci: number of loci simulated 
sample_seed = function(num_trees_to_sample, num_seeds_to_sample, num_pollen_donors, pollen_probability, num_loci) {
  #loop over x number of trees collectors are sampling from
  #Ex.) for(x in 1:trees_to_sample)
    #choose mother
    #Ex.) x <- sample(total_trees, 1) #selecting a row from the dataframe 
    #choose father(s) randomly -- a vector of rows representing potential pollen donors  
    #Ex.) vec <- sample(total_trees, num_pollen_donors)
    
    #create x number of seeds that collectors are sampling per tree
    #Ex.) for(x in 1:seeds_to_sample)
      #choose father based on probability vector 
      #Ex.) y <- sample(vec, 1, prob = pollen_probability)
  
      #Loop over number of loci simulated, in order to save the data
      #Ex.) for(x in 1:num_loci)
      #then select allele from father randomly, save
      #select allele randomly for mother, save
}


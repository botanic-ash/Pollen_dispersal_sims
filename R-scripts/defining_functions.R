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

##########################################################################################
#sample_seed() function
#Function creates seeds for collectors to sample from trees. First, a maternal tree is defined. Then father tree(s) is/are 
#defined. From there, alleles are selected from the mother and father randomly to create seeds.
#arugments:
  #data: the data imported from genalex files
  #num_trees_to_sample: number of trees in the population collectors will sample seed from 
  #num_seeds_to_sample: number of seeds collectors will take from each tree *CHANGE THIS TO BE A VECTOR OF VARYING SAMPLE SIZES PER TREE
  #num_pollen_donors: the number of trees that donate pollen to a single mother tree
  #pollen_probability: the probability of pollen donors to pollinate a mother tree (it can be equal for all donors, or skewed)
  #num_loci: number of loci simulated *TAKE THIS VARIABLE OUT, calculate by number of columns-2/2

sample_seed = function(data, num_trees_to_sample, num_seeds_to_sample, num_pollen_donors, pollen_probability) {
  #loop over x number of trees collectors are sampling from
  i=1 #simple counter variable to keep track of the current row (individual) 
  #defining the number of loci simulated
  num_loci = ((ncol(data)-2)/2)
  #defining the number of populations simulated
  dat_uniq = unique(data$Pop)
  num_pops = length(dat_uniq)
  for(w in 1:num_pops) {
    #subset data to for the current population group
    temp_data = data[which(data$Pop==paste("pop", w, sep="")),]
    total_trees = nrow(temp_data)
    for(x in 1:num_trees_to_sample) {
      #choose mother
      mother <- sample(total_trees, 1) #selecting a row from the dataframe 
      #choose father(s) randomly -- a vector of rows representing potential pollen donors  
      fathers <- sample(total_trees, num_pollen_donors)
      
      #create y number of seeds that collectors are sampling per tree
      #NOTE: change this loop to not sample an equal number from each tree 
      #SELECT A SIZE FROM THE VECTOR THEN LOOP THROUGH
      for(y in 1:num_seeds_to_sample[x]) {
        #choose father based on probability vector 
        if(length(fathers)>1){
          donor <- sample(fathers, 1, prob = pollen_probability)
        }else if (length(fathers==1)){
          donor = fathers
        }
    
        #Loop over number of loci simulated, in order to save the data
        j=1 #counter variable to keep track of the current column
        for(z in 1:num_loci) {
          #then select allele from father randomly, save
          p_alleles = c(temp_data[donor,paste("locus", z, "a", sep="")], data[donor,paste("locus", z, "b", sep="")])
          seeds_sampled[i,j] = sample(p_alleles,1)
          j=j+1
          #select allele randomly for mother, save
          m_alleles = c(temp_data[mother,paste("locus", z, "a", sep="")], data[mother,paste("locus", z, "b", sep="")])
          seeds_sampled[i,j] = sample(m_alleles,1)
          j=j+1
        }
        seeds_sampled[i,j] = paste("pop", w, sep="")
        i=i+1
      }
    }
  }
  seeds_sampled
}



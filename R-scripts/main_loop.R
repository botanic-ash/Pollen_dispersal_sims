#MSAT sampling code written by Kaylee Rosenberger and edited by Ash Hamilton
#Gibbs sampling + MH hastng code written by Ash Hamilton
#Main loop

##################################################################################
#Library functions
library(adegenet)
library(diveRsity)
library(poppr)
library(tidyr)
library(dplyr)

setwd("/Users/Ashley/Desktop/grad_school/2023_Winter/fund_comp_bio/Pollen_dispersal_sims")

#including edited arp2gen function 
source("R-scripts/arp2gen_edit.R")

#Defining this because arp2gen didn't want to use relative filepaths, using 1 parental pop with 2500 inds
mydir = "/Users/Ashley/Desktop/grad_school/2023_Winter/fund_comp_bio/Pollen_dispersal_sims/Simulations/one_pop_2500/"

####creating important functions####
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

#the function that samples the parental population and generates children WITH NO MUTATION by randomly choosing one of each allele from each parent for all loci
#output is a list of matrices, the row of each matrix corresponds to information about that child, the 1st matrix is the offspring genotype data, the second is the maternal gentoype, the last is the paternal genotype
sample_seed = function(data, num_trees_to_sample, num_seeds_to_sample, num_pollen_donors, pollen_probability) {
    #loop over x number of trees collectors are sampling from
    i=1 #simple counter variable to keep track of the current row (individual) 
    #defining the number of loci simulated
    num_loci = ((ncol(data)-2)/2)
    #defining the number of populations simulated
    dat_uniq = unique(data$Pop)
    num_pops = length(dat_uniq)
    #defining seeds_sampled to store results
    total_seeds = sum(num_seeds_to_sample)
    seeds_sampled = array(dim=c((num_pops*total_seeds), ((num_loci*2)))) #defining the array to store the results of the function in
    seeds_sampled = as.data.frame(seeds_sampled) #making it a data frame because it's easier to work with
    colnames(seeds_sampled) = colnames(data[,-c(1,2)])
    #number rows = total number of seeds sampled (number of populations*number of trees per pop*number of seeds per tree)
    #number columns = genotype of each seed sampled number of loci*2 (for two alleles) + 1 (to track the population), got rid of the plus 1 bc I don't need to track the population it came from
    w = 1
    seed_father_df <- data[FALSE,]
    seed_mother_df <- data[FALSE,]
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
            for(y in 1:(num_seeds_to_sample[x])) {
                #####################
                # #testing
                #  print(paste("Current tree: ", x))
                #  print(paste("Number seed: ", y))
                #####################
                #choose father based on probability vector 
                if(length(fathers)>1){
                    donor <- sample(fathers, 1, prob = pollen_probability)
                }else if (length(fathers==1)){
                    donor = fathers
                }
                #Loop over number of loci simulated, in order to save the data
                donor_data <- temp_data[donor,]
                rownames(donor_data) = NULL
                mother_data <- temp_data[mother,]
                rownames(mother_data) = NULL
                seed_father_df <- rbind(seed_father_df, donor_data)  
                seed_mother_df <- rbind(seed_mother_df, mother_data)
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
                #seeds_sampled[i,j] = paste("pop", w, sep="")
                
                i=i+1
            }
        }
    }
    return(list(seeds_sampled, seed_mother_df[,-2], seed_father_df[,-2])) #rows correspond to child number, row 1 = child one, child 1's mom and child 1's dad
}


#importing and converting arlequin files to genepop files
import_arp2gen_files(mydir,".arp$")

#importing and converting genepop files to genalex
import_gen2genalex_files(mydir, ".gen$")

##############################################################################################################
#DEFINING VARIABLES 

num_loci = 20 #number of loci simulated, needed to make a dataframe to save the data, starts with 20, maybe cut down later?

#defining array to store seeds that collectors have 'sampled'
#first we need to create column names depending on how many loci are present in simulations
#then, define the matrix, convert to dataframe, and rename the columns to label the data
#this dataframe keeps track of the alleles that are captured during sampling
loci_names = c()
for(i in 1:num_loci){
  loci_names = c(loci_names, paste("locus", i, "a", sep=""))
  loci_names = c(loci_names, paste("locus", i, "b", sep=""))
}


#define scenario params
num_trees_to_sample = 2
num_seeds_to_sample = 5
num_pollen_donors = 2
pollen_probability = c(.5, .5)


#######################################################################################################
#Main processing loop
#goes through each simulation replicate, calling the function with the parameters defined in the 
#defining_function_parameters.R script.
#results are saved in a 3D array 

#list of genalex files for all simulation replicates--genalex files end in .csv
#these have the simulated genetic data
genalex_list = list.files(mydir, ".csv$")


####things I think I need####

#list of genalex files for all simulation replicates--genalex files end in .csv
#these have the simulated genetic data
genalex_list = list.files(mydir, ".csv$")[42] #.... what are the diffs between the example pops?, does it matter which one I choose?, picking a random example populations to run everything with

#cut off first 2 rows in dataframe -- this is the population data, which is not required for our purposes
genetic_data = read.csv(paste(mydir, genalex_list, sep=""), header=FALSE)[-c(1,2),]

#giving the dataframe columns new names
names(genetic_data) = c("Ind", "Pop", loci_names)

genetic_data = genetic_data[-1,] #removing the first row -- repeat of now column headers, #is this maternal or seedling genetic data, i think it might be total parental genetic data

#thinning down number of loci for first runs
genetic_data_10_loci = genetic_data


#first thing is num_trees_to_sample
#second thing is num_seeds_to_sample
#third thing is num_pollen_donors
#fourth thing is pollen_probability

for(x in 1:length(skewed_params_AMH)) {
    temp = sample_seed(genetic_data, skewed_params_AMH[[x]][[1]], skewed_params_AMH[[x]][[2]], skewed_params_AMH[[x]][[3]], skewed_params_AMH[[x]][[4]]) #skewed params coming from the loaded r data that is saved from the defining functions script, seems fine, I will just pic my favorite?
    #output dfs are child, mother, father where row corresponds to child
    #don't really think I need this in a loop bc I only really need it over a single param set.... which I can make
}


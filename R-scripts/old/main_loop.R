#MSAT sampling code written by Kaylee Rosenberger and edited by Ash Hamilton
#Gibbs sampling + MH hastng code written by Ash Hamilton

##################################################################################
#Librarys
library(adegenet)
#library(diveRsity) #do I need this?
library(poppr)
library(tidyr)
library(dplyr)

#setwd("/Users/Ashley/Desktop/grad_school/2023_Winter/fund_comp_bio/final_project/Pollen_dispersal_sims")

#Defining this because arp2gen didn't want to use relative filepaths, using 1 parental pop with 2500 inds
mydir = "/Users/Ashley/Desktop/grad_school/2023_Winter/fund_comp_bio/final_project/Pollen_dispersal_sims/Simulations/one_pop_2500/"

mydir = "./one_pop_2500/"

####creating important functions####

#create a function that is the opposite of %in%
`%notin%` <- Negate(`%in%`)

#not sure what this does but I think I need it 
arp2gen<- function (infile)
{
    flForm <- strsplit(infile, split = "\\.")[[1]]
    if (substr(infile, 1, 2) == "./") {
        flForm <- flForm[-1]
    }
    else if (substr(infile, 1, 3) == "../") {
        flForm <- flForm[-(1:2)]
    }
    if (length(flForm) > 3) {
        stop("There were multiple '.' characters in your file name!")
    }
    tstfile <- paste(flForm[1], ".gen", sep = "")
    if (!file.exists(tstfile)) {
        fastScan <- function(fname) {
            s <- file.info(fname)$size
            buf <- readChar(fname, s, useBytes = TRUE)
            if (length(grep("\r", buf)) != 0L) {
                buf <- gsub("\r", "\n", buf)
                buf <- gsub("\n\n", "\n", buf)
            }
            return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
        }
        dat <- fastScan(infile)
        dat <- gsub("^\\s+|\\s+$", "", dat)
        dataType <- grep("*datatype=*", tolower(dat))
        if (strsplit(dat[dataType], "=")[[1]][2] != "MICROSAT") {
            stop("Data are not in 'MICROSAT' format!")
        }
        missDataLine <- grep("*missingdata=*", tolower(dat))
        missData <- noquote(substr(dat[missDataLine], nchar(dat[missDataLine]) -
                                       1, nchar(dat[missDataLine]) - 1))
        sampSizeLine <- grep("*samplesize=*", tolower(dat))
        # if (length(sampSizeLine) > 1) {
        sampNpos <- sapply(sampSizeLine, function(i) {
            return(regexpr("=", dat[i])[1])
        })
        # }
        popSizes <- as.numeric(substr(dat[sampSizeLine], start = sampNpos +
                                          1, stop = nchar(dat[sampSizeLine])))
        npops <- length(popSizes)
        sampStrt <- grep("*sampledata=*", tolower(dat))
        strts <- sapply(sampStrt, function(x) {
            if (dat[(x + 1)] == "") {
                return(x + 2)
            }
            else {
                return(x + 1)
            }
        })
        ends <- strts + ((popSizes * 2) - 1)
        nloci <- length(strsplit(dat[strts[1]], split = "\\s+")[[1]]) -
            2
        popGeno <- lapply(seq_along(strts), function(i) {
            return(dat[strts[i]:ends[i]])
        })
        popSzcheck <- sapply(popGeno, function(x) length(x)/2)
        if (!all(identical(popSzcheck, popSizes))) {
            stop("Failed! Please make sure that your file is formatted correctly.")
        }
        popIdx <- lapply(popGeno, function(x) {
            return(seq(1, length(x), 2))
        })
        popList <- lapply(seq_along(popGeno), function(i) {
            al1 <- matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]],
                                          split = "\\s+")), nrow = popSizes[i], byrow = TRUE)[,
                                                                                              -(1:2)]
            al2 <- matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +
                                                            1)], split = "\\s+")), nrow = popSizes[i], byrow = TRUE)
            tst <- matrix(paste(al1, al2, sep = ""), nrow = popSizes[i])
            tst <- cbind(paste(rep("pop", nrow(tst)), i, " ,",
                               sep = ""), tst)
            rm(al1, al2)
            z <- gc()
            rm(z)
            if (nchar(tst[1, 2]) == 4) {
                tst[tst == paste(missData, missData, sep = "")] <- "0000"
            }
            else {
                tst[tst == paste(missData, missData, sep = "")] <- "000000"
            }
            out <- apply(tst, 1, function(x) {
                return(paste(x, collapse = "\t"))
            })
            out <- c("POP", out)
            rm(tst)
            z <- gc()
            rm(z)
            return(out)
        })
        outfile <- strsplit(infile, "\\.")[[1]]
        if (length(outfile) >= 2) {
            outfile <- paste(outfile[-length(outfile)], collapse = ".")
        }
        else {
            outfile <- outfile[1]
        }
        loci <- paste("locus", 1:nloci, sep = "")
        loci <- c(paste(outfile, "_gen_converted", sep = ""),
                  loci)
        of <- c(loci, unlist(popList))
        out <- file(paste(outfile, ".gen", sep = ""), "w")
        for (i in 1:length(of)) {
            cat(of[i], "\n", file = out, sep = "")
        }
        close(out)
        return(TRUE)
    }
    else {
        return(NULL)
    }
}




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
num_seeds_to_sample = c(5, 5)
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
genetic_data_10_loci = genetic_data[,c(1:22)]

set.seed(2023)

#output dfs are child, mother, father where row corresponds to child
simple_data <- sample_seed(genetic_data_10_loci, num_trees_to_sample = num_trees_to_sample, num_seeds_to_sample = num_seeds_to_sample, num_pollen_donors = num_pollen_donors, pollen_probability = pollen_probability) 

getwd()
save(simple_data, file="../../R-scripts/Rdata/simple_data.Rdata")

#######################################################
# All code below written by Ash Hamilton

# Make sure I have the same dataset every time I reload this script
getwd()
# If loading from /Users/Ashley/Desktop/grad_school/2023_Winter/fund_comp_bio/final_project/Pollen_dispersal_sims
load(file="./R-scripts/Rdata/simple_data.Rdata")
# If loading from my_dir
load(file="../../R-scripts/Rdata/simple_data.Rdata")

children_data <- simple_data[[1]] # Making child data set its own object
mothers <- simple_data[[2]] # Making mother data set its own object
fathers <- simple_data[[3]] # Making father data set its own object

children_data$mom_ID <- mothers[,1] # Adding maternal IDs to child data for parentage assignment
children_data$dad_ID <- fathers[,1] # Adding paternal IDs to child data so I can test parentage assignment
children_data$child_ID <- seq.int(nrow(children_data))

# Creating a data set that contains all of the parents all mixed up, I feel this is likely to rep the data, maybe should set something to make sure no mom is also a dad....
parents <- c(unique(children_data$mom_ID), unique(children_data$dad_ID))
parental_data <- subset(genetic_data_10_loci, Ind %in% parents)[, -2]

# Making the parental data df that isn't dep on the entire script being run top to bottom
parental_data <- rbind(unique(mothers), unique(fathers)) # I think this differs from the old paretal df bc it doesn't have the pop column

# Now I have a children df and a parental df with all the loci! woot!

#####making the model####


### Calculating log likelihood of each possible being the father of each known mother-offspring pair
num_loci<- 10

ll_D_given_OM <- matrix(nrow = nrow(dad_data), ncol = nrow(children_data), 
                        dimnames = list( c(dad_data$Ind), c(paste0(children_data$mom_ID, "-", children_data$child_ID))))


pr_one_D_given_OM_at_each_locus <- matrix(nrow = num_loci, ncol = nrow(children_data), 
                                      dimnames = list(seq(1:num_loci),c(paste0(children_data$mom_ID, "-", children_data$child_ID))))

each_D_pr_one_D_given_OM_at_each_locus <- list(pr_one_D_given_OM_at_each_locus, pr_one_D_given_OM_at_each_locus, pr_one_D_given_OM_at_each_locus)


for(child in 1:nrow(children_data)){
    
all_child_loci <- children_data[child, -c(21:23)] #Pulling a single child and dropping the various IDs so just locus data remains, unsure that I necessarily need to do this??
all_mom_loci<- mothers[child,-1] #Pulling the mother corresponds to the child pulled and dropping the ID so just locus data remains

    for(dad in 1:nrow(dad_data)){
    all_dad_loci <- dad_data[dad,-1] #Pulling a single father and dropping the ID and pop columns so just locus data remains
    
    pr_D_given_OM_at_each_locus <- rep(0, num_loci) # Making a vector to hold the prob of this dad for this mother off pair at each locus
    print(paste0(children_data$mom_ID[child], "-", children_data$child_ID[child]))
    
        for(locus in 1:num_loci){
            
            # Pulling both alleles at a single locus from the mom, child, and dad
            mom_alleles <- all_mom_loci[, c(2*locus -1, 2*locus)]
            child_alleles <- all_child_loci[, c(2*locus -1, 2*locus)]
            dad_alleles <- all_dad_loci[, c(2*locus -1, 2*locus)]
            
            
            # Making a punnet square from the known mother and the potential father
            
            possible_children <- matrix(0, nrow = 4, ncol = 2) # Making my punnet square matrix of possible children for the current mother and father
            
            k <- 1  # Start a ticker so the child genotypes go in the right row of the matrix
            for(i in mom_alleles[1,]){
                for(j in dad_alleles[1,]){
                    possible_children[k,1] <- i # Put allele from mom in col 1, probably could do this seperate from the dads and just overwrite the second column over each dad.. would be faster at least
                    possible_children[k,2] <- j # Put allele from dad in col 2
                    k <- k+1
                }
            }
        
            allele_order_1 <- child_alleles # Getting the allele order from the child
            allele_order_2 <- child_alleles[2:1] # Reversing the allele order from the child
            
            real_children <- possible_children[,1] == allele_order_1[,1] & possible_children[,2] == allele_order_1[,2] # Creating a Boolean vector of 4 True/False's that will give True only if both the real childs alleles match the ones in possible_children 
            real_children_2 <- possible_children[,1] == allele_order_2[,1] & possible_children[,2] == allele_order_2[,2] # Creating a second Boolean vector of 4 True/False's that will give True only if both the real childs alleles in reverse order match the ones in possible_children 
            any_real_child <- real_children + real_children_2 # Returns number of matches per parental allele combo regardless of order in child
            any_real_child[any_real_child == 2] <- 1 # Overwrites any 2s (from homozygotes) with 1s 
            any_real_child_probs <- sum(any_real_child)*.25 # Sum the number of matches and multiply by .25 which is the prob of any of the 4 possible children in the punnet sqaure being made  
            
            pr_D_given_OM_at_each_locus[locus] <- any_real_child_probs # Saving the prob of that father at each locus 
        }
        
        each_D_pr_one_D_given_OM_at_each_locus[[dad]][,child] <- pr_D_given_OM_at_each_locus # Saving each of the probabilities at each locus to see if I can ID where the issue is casuing the returning of 0's even when a child literally has to come from that dad-mom-combo
        
        
        if(sum(pr_D_given_OM_at_each_locus == 0) > 1){ # Technically if any loci return a 0, the dad isn't a valid option so ll should be 0 BUT THIS WILL CHANGE WITH ERROR!!!
            ll_D_given_OM[dad, child] <- 0
        }
        else{
            ll_D_given_OM[dad, child] <- sum(log(pr_D_given_OM_at_each_locus)) # Calculating the log likelihoods of each dad by multiplying probabilities across loci, with enough loci will only leave a single father
        }
    }
}

###Next incorporate probability of paternal ind with allele freqs from the population
#It seems to me IBD is only really used when other parent is unknown, because in this instance I know that 1 of the alleles is from the mother and not randomly from the population, it seems to me that scaling the likelihood of parentage using population level allele frequencies can be done with just the addition of the frequency of the potential paternal allele (if it can be either paternal allele then .5 * freq allele + 5 * freq allele)
#need to ask Matthias about this, 

#there only seems to be options where probabilities assume one parent is known and both are sampled or both are unknown and only one is sampled

#I think where I multiply by 1 or 0 depending on which child is present is what the multinomial does in the Neff 2001 paper




###Incorporate probability of paternal ind with error 

error_2_est = .1
error_1_est = .05

#rates of errors
error_2 = error_2_est
error_1 = error_1_est
e_2 = error_2 / (1 + error_2) #this is wrong, should be k + error_2 where k is the number of codominant loci
e_1 = error_1 / (1 + error_1)


# Function for obtaining allele frequencies at each locus
# Input: Desired locus number for which you will obtain allele frequency data
# Output: A list-   Element 1) vector of allele frequencies at that locus
#                   Element 2) matrix with 2 rows, row 1= allele 
get_alleles_and_freqs <- function(locus){
    
    # Pulling all alleles at a single locus from all parents in the population
    all_alleles <- sort(unique(c(parental_data[, c(2*locus)], parental_data[, 2*locus + 1])))
    
    # Getting the frequency of each allele, I think this should be input somewhere???? currently just going to make everything 1/k 
    allele_freqs <- matrix(c(all_alleles, rep(1/length(all_alleles), length(all_alleles))), ncol = length(all_alleles), nrow = 2, byrow = T)
    return(list(all_alleles,allele_freqs))
}


# Function that takes a vector containing the two alleles of the true genotype and a vector of the two alleles of the  observed genotype and outputs the probability of observing the oberved type given the true type
calc_true_geno_prob <- function(ind_geno, true_geno, allele_freqs){
    #browser()
    # Setting the Kronecker delta-variable of the obs geno
    if(ind_geno[1] == ind_geno[2]){  
        delta_uv <- 1 # If observed genotype = homo, delta_uv = 1
    } else {
        delta_uv <- 0 # If observed genotype = het, delta_uv = 0
    }
    
    # Setting the Kronecker delta-variable of the true geno
    if(true_geno[1] == true_geno[2]){  
        delta_xw <- 1 # If observed genotype = homo, delta_uv = 1
    } else {
        delta_xw <- 0 # If observed genotype = het, delta_uv = 0
    }
    
    # Calculate the frequency of the true genotype, will need to be multiplied by probability of obs geno given true geno
    true_geno_freq <- (2 - delta_xw) * as.numeric(allele_freqs[2, c(allele_freqs[1,] == true_geno[1])]) * as.numeric(allele_freqs[2, c(allele_freqs[1,] == true_geno[2])])
    
    # If true geno is a homozygote
    if(true_geno[1] == true_geno[2]){ 
        
        # If mom geno is the same geno as the true geno
        if(ind_geno[1] == ind_geno[2] & ind_geno[1] == ind_geno[1]){
            true_geno_prob = (1 - error_2)^2
            
            # If mom geno has a single allele that is the same as the true geno
        } else if(ind_geno[1] == true_geno[1] | ind_geno[2] == true_geno[1]){
            true_geno_prob = 2*e_2*(1 - error_2)
            
            # If mom geno has no alleles that are the same as the true geno
        } else { 
            true_geno_prob = (2 - delta_uv) * e_2^2 
        }
        
        # If true geno is a heterozygote    
    } else if(true_geno[1] != true_geno[2]){ 
        
        # If mom geno is the same geno as the true geno
        if(all(ind_geno == true_geno) | all(ind_geno[2:1] == true_geno)){
            true_geno_prob = (1 - error_2)^2 + e_2^2 - (2*e_1*(1 - error_2- e_2)^2)
            
            # If mom geno is homozygous and matches one of the alleles in the true geno    
        } else if(ind_geno[1] == ind_geno[2] & ind_geno[1] %in% true_geno){
            true_geno_prob = (e_2*(1 - error_2)) +  (e_1*(1 - error_2- e_2)^2)
            
            # If neither allele of the mom geno matches either of the alleles in the true geno  
        } else if(ind_geno[1] %notin% true_geno & ind_geno[2] %notin% true_geno){
            true_geno_prob = (2 - delta_uv) * e_2^2 
            
            # Otherwise
        } else {
            true_geno_prob = e_2*(1 - error_2- e_2)
        }
    }
    
    true_geno_prob = true_geno_prob * true_geno_freq # Probability needs to be scaled by the frequency of that genotype
    
    return(true_geno_prob)
}


prob_hidden_genos <- function(locus, ind_ID) {
    #browser()
    if(ind_ID %in% parental_data$Ind){
        all_ind_loci<- parental_data[parental_data$Ind == ind_ID,] #Pulling the parent ind that corresponds to the ind ID  
        # set individuals genotype
        ind_geno <- all_ind_loci[, c(2*locus, 2*locus + 1)]
    } else {
        all_ind_loci<- children_data[children_data$child_ID == ind_ID,] #Pulling the child ind corresponds to the ind ID
        # set individuals genotype
        ind_geno <- all_ind_loci[, c(2*locus - 1, 2*locus)]
        }

    # Pulling all alleles at a single locus from all parents in the population
    all_alleles <- get_alleles_and_freqs(locus = locus)
    
    # Making list of all possible genotypes
    all_hets <- combn(all_alleles[[1]], 2) # All heterozygotous genotypes 
    all_homos <- rbind(all_alleles[[1]], all_alleles[[1]]) # All homozygotous genotypes
    all_genos <- cbind(all_hets, all_homos) # All genotypes
    
    # Getting the frequency of each allele, I think this should be input somewhere....???? currently just going to make everything 1/k 
    allele_freqs <- all_alleles[[2]]
    
    # Make a vector which will hold the probabilities of each mother having the true genotype
    ind_true_geno_probs <- rep(0, ncol(all_genos)) # Correlates to column of all_genos
    
    # Iterate through all possible "true" genotypes
    for(geno in 1:ncol(all_genos)){
        true_geno <- all_genos[,geno] # Set current "true genotype"
        
        ind_true_geno_probs[geno] <- calc_true_geno_prob(ind_geno = ind_geno, true_geno = true_geno, allele_freqs = allele_freqs) #There is some issue with this when it is called on locus 2 after looping.... unsure what it is
            
        } # Close the for loop for all possible genotypes
    
    ind_true_geno_probs <- matrix(c(ind_true_geno_probs), nrow = 1) 
    colnames(ind_true_geno_probs) <- paste0(all_genos[1,], "-", all_genos[2,]) # Make the column name the genotypes
    
    #ind_all_geno_all_loci_probs[[locus]] <- ind_true_geno_probs # Each loci = it's own element of the list
    
    return(ind_true_geno_probs)
}


calc_pr_child <- function(locus, mom_ID, dad_ID, child_ID){
    #browser()
    all_child_loci <- children_data[child_ID, -c(21:23)] # Pulling the single child and dropping the various IDs so just locus data remains
    
    # Pulling both alleles at a single locus from the ind
    child_geno <- all_child_loci[, c(2*locus - 1, 2*locus)] # Ind geno based on ind name
    
    allele_freqs <- get_alleles_and_freqs(locus)[[2]] 
    
    ind_geno <- NA
    prob_mom <- prob_hidden_genos(locus = locus, ind_ID = mom_ID) #THIS IS THE CALL THAT ISN'T WORKING
    
    ind_geno <- NA
    prob_dad <- prob_hidden_genos(locus = locus, ind_ID = dad_ID)
    
    c(colnames(prob_mom), colnames(prob_dad))
    
    # Making list of all true genotype parental pairings, true for all inds at each loci
    all_diff <- combn(colnames(prob_mom), 2) # All diff genotype pairings
    all_same <- rbind(colnames(prob_mom), colnames(prob_dad)) # All same genotype pairings
    all_pairings <- cbind(all_diff, all_same) # All genotype pairings
    
    ll_child_given_any_mom_dad <- vector(length = ncol(all_pairings))
    ll_child_given_any_dad_mom <- vector(length = ncol(all_pairings))
    
    children_of_pairing <- list()
    
    for(pairing in 1:ncol(all_pairings)){
        # Making a punnet square from the each possible mother and the father
        
        possible_children <- matrix(0, nrow = 4, ncol = 2) # Making my punnet square matrix of possible children for the current mother and father
        
        k <- 1  # Start a ticker so the child genotypes go in the right row of the matrix
        
        # Getting the unique alleles from each parent into vectors again
        parent1 <- unlist(strsplit(all_pairings[1,pairing], "-"))
        parent2 <- unlist(strsplit(all_pairings[2,pairing], "-"))
        
        
        for(i in parent1){
            for(j in parent2){
                possible_children[k,1] <- i # Put allele from parent 1 in col 1, 
                possible_children[k,2] <- j # Put allele from parent 2 in col 2, 
                k <- k+1
            }
        }
        
        # Each possible parental pairing has the possible children assigned
        children_of_pairing[[pairing]] <- possible_children
        
        # Making vector to hold all probabilities of observing the child genotype given the two parental genotypes
        obs_geno_prob <- rep(0, 4)
        
        
        for(real_geno in 1:nrow(possible_children)){
            
            obs_geno_prob[real_geno] <- calc_true_geno_prob(ind_geno = child_geno, true_geno = possible_children[real_geno,], allele_freqs = allele_freqs)
            
        }
        
        pr_child_given_real_geno <- sum(obs_geno_prob[real_geno] * (1/4)) # Prob observing the real geno given all possible true genos from the designated parental pairing
        
        # Prob observing this child given obs mom given true mom in spot 1 of the pairing and obs dad given true dad in spot 2 of the pairing
        ll_child_given_mom_dad <- sum(log(prob_mom[colnames(prob_mom) == all_pairings[1,pairing]]), log(prob_dad[colnames(prob_dad) == all_pairings[2,pairing]]), log(pr_child_given_real_geno))
        
        # Prob observing this child given obs mom given true mom in spot 2 of the pairing and obs dad given true dad in spot 1 of the pairing
        ll_child_given_dad_mom <- sum(log(prob_mom[colnames(prob_mom) == all_pairings[2,pairing]]), log(prob_dad[colnames(prob_dad) == all_pairings[1,pairing]]), log(pr_child_given_real_geno))
        
        ll_child_given_any_mom_dad[pairing] <- ll_child_given_mom_dad
        
        ll_child_given_any_mom_dad[pairing] <- ll_child_given_dad_mom
    }
    
    total_ll_child_given_mom_dad <- log(sum(exp(ll_child_given_any_mom_dad)))
    total_ll_child_given_dad_mom <- log(sum(exp(ll_child_given_any_mom_dad)))
    total_ll_child_given_parents <- log(sum(exp(total_ll_child_given_dad_mom), exp(total_ll_child_given_mom_dad)))

    return(total_ll_child_given_parents)
}



num_loci<- 10

mom_ID = "1909"
dad_ID = "2479"
locus = 1
child_ID = "1"


ll_O_given_one_DM <- matrix(nrow = nrow(dad_data), ncol = nrow(children_data), 
                        dimnames = list( c(dad_data$Ind), c(paste0(children_data$mom_ID, "-", children_data$child_ID))))

#looping above functions across unique child - dad pairings to calc prob obs off|known mom, putative dad with error incorporated
for(child in 1:nrow(children_data)){
    child_ID <- children_data$child_ID[child] # Pulling the child ID corresponds to the child
    mom_ID <- children_data$mom_ID[child] # Pulling the mother corresponds to the child
    
    for(dad in 1:nrow(dad_data)){
        dad_ID <- dad_data$Ind[dad] # Pulling a single father and dropping the ID and pop columns so just locus data remains
        ll_O_given_DM_at_each_locus <- rep(0, num_loci) # Making a vector to hold the prob of this dad for specific mother off pair at each locus
        
        for(locus in 1:num_loci){
            ll_O_given_DM_at_each_locus[locus] <- calc_pr_child(locus = locus, mom_ID = mom_ID, dad_ID = dad_ID, child_ID = child_ID) # Saving the prob of that child given the known mother and each potential father at each locus 
        }
        ll_O_given_one_DM[dad, child] <- sum(ll_O_given_DM_at_each_locus) # Since already in log likelihood space can sum across all loci (assuming independence) to get the ll of each dad for each mother/off pair
    }
}

total_ll <- sum(log(colSums(exp(ll_O_given_one_DM))) - log(nrow(ll_O_given_one_DM))) # Overall likelihood of the data (mother and offspring pairs) and the parameters, first sum over columns to sum across all fathers then divide by the probability of each father being the correct father which is the same as number of fathers and therefore the number of rows in the ll_O_given_one_DM df (assumes each father has equal prob to contribute to each mother-off pair)


#corresponds to Neffss t matrix
t(ll_O_given_one_DM)

    
#calculate exclusion probability

#### Re-writing above functions so that they can take different allele frequency values and error values####

# Function for obtaining allele frequencies at each locus
# Input: Desired locus number for which you will obtain allele frequency data
# Output: A list-   Element 1) vector of allele frequencies at that locus
#                   Element 2) matrix with 2 rows, row 1= allele, row 2 = frequency of alleles
get_alleles_and_freqs <- function(locus, parental_data){
    
    # Pulling all alleles at a single locus from all parents in the population
    all_alleles <- sort(unique(c(parental_data[, c(2*locus)], parental_data[, 2*locus + 1])))
    
    # Getting the frequency of each allele, I think this should be input somewhere???? currently just going to make everything 1/k 
    allele_freqs <- matrix(c(all_alleles, rep(1/length(all_alleles), length(all_alleles))), ncol = length(all_alleles), nrow = 2, byrow = T)
    return(list(all_alleles,allele_freqs))
}

calc_true_geno_prob <- function(ind_geno, true_geno, allele_freqs, error_1, error_2){
    #browser()
    k = ncol(allele_freqs)
    e_1 = error_1 / (1 + error_1)
    e_2 = error_2 / (k + error_2)
    
    # Setting the Kronecker delta-variable of the obs geno
    if(ind_geno[1] == ind_geno[2]){  
        delta_uv <- 1 # If observed genotype = homo, delta_uv = 1
    } else {
        delta_uv <- 0 # If observed genotype = het, delta_uv = 0
    }
    
    # Setting the Kronecker delta-variable of the true geno
    if(true_geno[1] == true_geno[2]){  
        delta_xw <- 1 # If observed genotype = homo, delta_uv = 1
    } else {
        delta_xw <- 0 # If observed genotype = het, delta_uv = 0
    }
    
    # Calculate the frequency of the true genotype, will need to be multiplied by probability of obs geno given true geno
    true_geno_freq <- (2 - delta_xw) * as.numeric(allele_freqs[2, c(allele_freqs[1,] == true_geno[1])]) * as.numeric(allele_freqs[2, c(allele_freqs[1,] == true_geno[2])])
    
    # maybe change this to an which %in% with all possible satisfactory conditions listed in the right side of the in, and then have the next thing written based on the in result of the in and have a matrix with corresponding equations in the correct indexed position that can then be evaluated with eval(parse(text = matrix[matched index value from which command])), should run faster
    
    
    # If true geno is a homozygote
    if(true_geno[1] == true_geno[2]){ 
        
        # If mom geno is the same geno as the true geno
        if(ind_geno[1] == ind_geno[2] & ind_geno[1] == ind_geno[1]){
            true_geno_prob = (1 - error_2)^2
            
            # If mom geno has a single allele that is the same as the true geno
        } else if(ind_geno[1] == true_geno[1] | ind_geno[2] == true_geno[1]){
            true_geno_prob = 2*e_2*(1 - error_2)
            
            # If mom geno has no alleles that are the same as the true geno
        } else { 
            true_geno_prob = (2 - delta_uv) * e_2^2 
        }
        
        # If true geno is a heterozygote    
    } else if(true_geno[1] != true_geno[2]){ 
        
        # If mom geno is the same geno as the true geno
        if(all(ind_geno == true_geno) | all(ind_geno[2:1] == true_geno)){
            true_geno_prob = (1 - error_2)^2 + e_2^2 - (2*e_1*(1 - error_2- e_2)^2)
            
            # If mom geno is homozygous and matches one of the alleles in the true geno    
        } else if(ind_geno[1] == ind_geno[2] & ind_geno[1] %in% true_geno){
            true_geno_prob = (e_2*(1 - error_2)) +  (e_1*(1 - error_2- e_2)^2)
            
            # If neither allele of the mom geno matches either of the alleles in the true geno  
        } else if(ind_geno[1] %notin% true_geno & ind_geno[2] %notin% true_geno){
            true_geno_prob = (2 - delta_uv) * e_2^2 
            
            # Otherwise
        } else {
            true_geno_prob = e_2*(1 - error_2- e_2)
        }
    }
    
    true_geno_prob = true_geno_prob * true_geno_freq # Probability needs to be scaled by the frequency of that genotype
    
    return(true_geno_prob)
}

prob_hidden_genos <- function(locus, ind_ID, error_1, error_2) {
    #browser()
    if(ind_ID %in% parental_data$Ind){
        all_ind_loci<- parental_data[parental_data$Ind == ind_ID,] #Pulling the parent ind that corresponds to the ind ID  
        # set individuals genotype
        ind_geno <- all_ind_loci[, c(2*locus, 2*locus + 1)]
    } else {
        all_ind_loci<- children_data[children_data$child_ID == ind_ID,] #Pulling the child ind corresponds to the ind ID
        # set individuals genotype
        ind_geno <- all_ind_loci[, c(2*locus - 1, 2*locus)]
    }
    
    # Pulling all alleles at a single locus from all parents in the population
    all_alleles <- get_alleles_and_freqs(locus = locus, parental_data = parental_data)
    
    # Making list of all possible genotypes
    all_hets <- combn(all_alleles[[1]], 2) # All heterozygotous genotypes 
    all_homos <- rbind(all_alleles[[1]], all_alleles[[1]]) # All homozygotous genotypes
    all_genos <- cbind(all_hets, all_homos) # All genotypes
    
    # Getting the frequency of each allele, I think this should be input somewhere....???? currently just going to make everything 1/k 
    allele_freqs <- all_alleles[[2]]
    
    # Make a vector which will hold the probabilities of each mother having the true genotype
    ind_true_geno_probs <- rep(0, ncol(all_genos)) # Correlates to column of all_genos
    
    # Iterate through all possible "true" genotypes
    for(geno in 1:ncol(all_genos)){
        true_geno <- all_genos[,geno] # Set current "true genotype"
        
        ind_true_geno_probs[geno] <- calc_true_geno_prob(ind_geno = ind_geno, true_geno = true_geno, allele_freqs = allele_freqs, error_1 = error_1, error_2 = error_2)
        
    } # Close the for loop for all possible genotypes
    
    ind_true_geno_probs <- matrix(c(ind_true_geno_probs), nrow = 1) 
    colnames(ind_true_geno_probs) <- paste0(all_genos[1,], "-", all_genos[2,]) # Make the column name the genotypes
    
    #outer function might be faster here
    
    
    #ind_all_geno_all_loci_probs[[locus]] <- ind_true_geno_probs # Each loci = it's own element of the list
    
    return(ind_true_geno_probs)
}

calc_pr_child <- function(locus, mom_ID, dad_ID, child_ID, error_1, error_2){
    #browser()
    all_child_loci <- children_data[child_ID, -c(21:23)] # Pulling the single child and dropping the various IDs so just locus data remains
    
    # Pulling both alleles at a single locus from the ind
    child_geno <- all_child_loci[, c(2*locus - 1, 2*locus)] # Ind geno based on ind name
    
    allele_freqs <- get_alleles_and_freqs(locus, parental_data = parental_data)[[2]] 
    
    ind_geno <- NA
    prob_mom <- prob_hidden_genos(locus = locus, ind_ID = mom_ID, error_1 = error_1, error_2 = error_2) 
    
    ind_geno <- NA
    prob_dad <- prob_hidden_genos(locus = locus, ind_ID = dad_ID, error_1 = error_1, error_2 = error_2)
    
    # Making list of all true genotype parental pairings, true for all inds at each loci
    all_diff <- combn(colnames(prob_mom), 2) # All diff genotype pairings
    all_same <- rbind(colnames(prob_mom), colnames(prob_dad)) # All same genotype pairings
    all_pairings <- cbind(all_diff, all_same) # All genotype pairings
    
    ll_child_given_any_mom_dad <- vector(length = ncol(all_pairings))
    ll_child_given_any_dad_mom <- vector(length = ncol(all_pairings))
    
    children_of_pairing <- list()
    
    for(pairing in 1:ncol(all_pairings)){
        # Making a punnet square from the each possible mother and the father
        
        possible_children <- matrix(0, nrow = 4, ncol = 2) # Making my punnet square matrix of possible children for the current mother and father
        
        k <- 1  # Start a ticker so the child genotypes go in the right row of the matrix
        
        # Getting the unique alleles from each parent into vectors again
        parent1 <- unlist(strsplit(all_pairings[1,pairing], "-"))
        parent2 <- unlist(strsplit(all_pairings[2,pairing], "-"))
        
        
        for(i in parent1){
            for(j in parent2){
                possible_children[k,1] <- i # Put allele from parent 1 in col 1, 
                possible_children[k,2] <- j # Put allele from parent 2 in col 2, 
                k <- k+1
            }
        }
        
        # Each possible parental pairing has the possible children assigned
        children_of_pairing[[pairing]] <- possible_children
        
        # Making vector to hold all probabilities of observing the child genotype given the two parental genotypes
        obs_geno_prob <- rep(0, 4)
        
        
        for(real_geno in 1:nrow(possible_children)){
            
            obs_geno_prob[real_geno] <- calc_true_geno_prob(ind_geno = child_geno, true_geno = possible_children[real_geno,], allele_freqs = allele_freqs, error_1 = error_1, error_2 = error_2)
            
        }
        
        pr_child_given_real_geno <- sum(obs_geno_prob * (1/4)) # Prob observing the real geno given all possible true genos from the designated parental pairing
        
        # Prob observing this child given obs mom given true mom in spot 1 of the pairing and obs dad given true dad in spot 2 of the pairing
        ll_child_given_mom_dad <- sum(log(prob_mom[colnames(prob_mom) == all_pairings[1,pairing]]), log(prob_dad[colnames(prob_dad) == all_pairings[2,pairing]]), log(pr_child_given_real_geno))
        
        # Prob observing this child given obs mom given true mom in spot 2 of the pairing and obs dad given true dad in spot 1 of the pairing
        ll_child_given_dad_mom <- sum(log(prob_mom[colnames(prob_mom) == all_pairings[2,pairing]]), log(prob_dad[colnames(prob_dad) == all_pairings[1,pairing]]), log(pr_child_given_real_geno))
        
        ll_child_given_any_mom_dad[pairing] <- ll_child_given_mom_dad
        
        ll_child_given_any_dad_mom[pairing] <- ll_child_given_dad_mom
    }
    
    total_ll_child_given_mom_dad <- log(sum(exp(ll_child_given_any_mom_dad)))
    total_ll_child_given_dad_mom <- log(sum(exp(ll_child_given_any_dad_mom)))
    total_ll_child_given_parents <- log(sum(exp(total_ll_child_given_dad_mom), exp(total_ll_child_given_mom_dad)))
    
    return(total_ll_child_given_parents)
}

#unclear how to incorporate allele frequencies that differ at each loci
get_ll <- function(children_data, parental_data, error_1, error_2) {
    dad_data <- subset(parental_data, parental_data$Ind %notin% children_data$mom_ID) # Assumes a parent can only ever be a mother or a father
    
    # Make matrix to hold the output of the for loop below
    ll_O_given_one_DM <- matrix(nrow = nrow(dad_data), ncol = nrow(children_data), 
                                dimnames = list( c(dad_data$Ind), c(paste0(children_data$mom_ID, "-", children_data$child_ID))))
    
    for(child in 1:nrow(children_data)){
        child_ID <- children_data$child_ID[child] # Pulling the child ID corresponds to the child
        mom_ID <- children_data$mom_ID[child] # Pulling the mother corresponds to the child
        
        for(dad in 1:nrow(dad_data)){
            dad_ID <- dad_data$Ind[dad] # Pulling a single father and dropping the ID and pop columns so just locus data remains
            ll_O_given_DM_at_each_locus <- rep(0, num_loci) # Making a vector to hold the prob of this dad for specific mother off pair at each locus
            
            for(locus in 1:num_loci){
                ll_O_given_DM_at_each_locus[locus] <- calc_pr_child(locus = locus, mom_ID = mom_ID, dad_ID = dad_ID, child_ID = child_ID, error_1 = error_1, error_2 = error_2) # Saving the prob of that child given the known mother and each potential father at each locus 
            }
            ll_O_given_one_DM[dad, child] <- sum(ll_O_given_DM_at_each_locus) # Since already in log likelihood space can sum across all loci (assuming independence) to get the ll of each dad for each mother/off pair
        }
    }
    
    total_ll <- sum(log(colSums(exp(ll_O_given_one_DM))) - log(nrow(ll_O_given_one_DM))) # Overall likelihood of the data (mother and offspring pairs) and the parameters, first sum over columns to sum across all fathers then divide by the probability of each father being the correct father which is the same as number of fathers and therefore the number of rows in the ll_O_given_one_DM df (assumes each father has equal prob to contribute to each mother-off pair)
    return(list(total_ll, ll_O_given_one_DM))
}


error_2_est = .1
error_1_est = .05

#rates of errors
error_2 = error_2_est
error_1 = error_1_est



# Testing my big compiled function
get_ll(children_data = children_data, parental_data = parental_data, error_1 = error_1, error_2 = error_2)[[1]]

get_ll(children_data = children_data, parental_data = parental_data, error_1 = .05, error_2 = .05)[[1]]

get_ll(children_data = children_data, parental_data = parental_data, error_1 = .005, error_2 = .1)[[1]]
#definitely works but does take a long time to work (a little over a minute).... will only be able to run maybe 1000 iterations.... that will take 16 hours at the current run time

# Uniform prior with a log on it (aka both of my error priors)
log_prior = function(p){
    if((p<0) || (p>1)){  # || here means "or", if p is smaller than 0 or bigger than 1 return -inf (which when exponentiated returns 0), otherwise return 0 (which once exponentiated equals 1) when this function is called, basically just makes sure all evaluated p values between 1 and 0 have the same prior odds of 1 and anything outside of that has a prior odds of 0
        return(-Inf)}
    else{
        return(0)}
}

# Making the MH sampler
error_sampler = function(children_data, parental_data, niter, error_startvals, proposalsd){
    error = matrix(0,ncol = niter, nrow = 2)
    error[1,1] = error_startvals[1]
    error[2,1] = error_startvals[2]
    for(i in 2:niter){
        current_error1 = error[1, i-1]
        current_error2 = error[2, i-1]
        new_error1 = current_error1 + rnorm(1,0,proposalsd) #this isn't good because error can't be negative... to make symmetric call the rnorm and then write an if else that asks if the resulting new_error would negative, if yes then subtract the value from it to get back above 1?
        A = exp(log_prior(new_error1) + get_ll(children_data = children_data, parental_data = parental_data, error_1 = new_error1, error_2 = current_error2)[[1]] - log_prior(current_error1) - get_ll(children_data = children_data, parental_data = parental_data, error_1 = current_error1, error_2 = current_error2)[[1]]) 
        if(runif(1)<A){ #the bigger the A value is, the more likely you get to move to the new p value
            error[1,i] = new_error1       # accept move with probability min(1,A)
        } else {
            error[1,i] = current_error1        # otherwise "reject" move, and stay where we are
        }
        current_error1 = error[1,i] # set the new current error1 that takes into account the last half step  
        new_error2 = current_error2 + rnorm(1,0,proposalsd) #this isn't good because error can't be negative... to make symmetric call the rnorm and then write an if else that asks if the resulting new_error would negative, if yes then subtract the value from it to get back above 1?
        A = exp(log_prior(new_error2) + get_ll(children_data = children_data, parental_data = parental_data, error_1 = current_error1, error_2 = new_error2) - log_prior(current_error2) - get_ll(children_data = children_data, parental_data = parental_data, error_1 = current_error1, error_2 = current_error2)) 
        if(runif(1)<A){ #the bigger the A value is, the more likely you get to move to the new p value
            error[2,i] = new_error2       # accept move with probability min(1,A)
        } else {
            error[2,i] = current_error2        # otherwise "reject" move, and stay where we are
        }
    }
    
    return(error)
}


# Running the Gibbs sampler
first_big_try <- error_sampler(children_data = children_data_error, parental_data = parental_data_error, niter = 500, error_startvals = c(.02,.05), proposalsd = .02) # Attempting 500 iterations, should take about 24 hours to run





hist(first_big_try[1,],prob=T)
hist(first_big_try[2,],prob=T)


# Taking the result of MH updates and using it to obtain most likely dad | data 
final_dad_avg_likelihoods <- get_ll(children_data = children_data, parental_data = parental_data, error_1 = mean(first_try[1,]), error_2 = mean(first_try[2,]))[[2]]

dad_index <- apply(final_dad_avg_likelihoods, 2, function(x) which.max(x)) # Finding the index of the row w/ the highest likelihood

dad_ID <- dimnames(final_dad_avg_likelihoods)[[1]][dad_index] # Getting the names of the estimated most likely dads 

proportions(table(dad_ID == children_data$dad_ID)) # Getting proportion of correct dads to see how my algorithm does



####Simulating error and adding it to the simulated reads####
# Error class 1 = allelic dropout aka failure to amplify one of 2 alleles, worse w/ lower DNA quality, only affects hets by making them into false homos
# Error class 2 = stochastic errors excluding allelic dropout, including false alleles, miscalling, and contaminant DNA

simulating_error <- function(data_for_editing, type1_rate, type2_rate, parental_data = NA){
    if(missing(parental_data)){
        parental_data <- data_for_editing
    }
    
    # Pulling just the loci from the df
    locus_data <- data_for_editing[,grepl( "locus" , colnames(data_for_editing))]
    
    # Making the new matrix with the same shape as the old matrix
    locus_data_new <- data.frame(matrix(ncol = ncol(locus_data), nrow = nrow(locus_data)))
    diffs_in_genos <- data.frame(matrix(ncol = ncol(locus_data)/2, nrow = nrow(locus_data)))
    
    for(locus in 1:num_loci){
        alleles <- get_alleles_and_freqs(locus, parental_data = parental_data)[[1]] # Want to get all possible alleles from the parental data
        k = length(alleles) # k = number of alleles at designated locus
        e_1 = error_1 / (1 + type1_rate)
        e_2 = error_2 / (k + type2_rate)
        genos <- locus_data[, c(2*locus - 1, 2*locus)] # Keep only the genotypes at the designated locus

        for(ind in 1:nrow(genos)){
            
            old_geno <- paste0(genos[ind,1], "-", genos[ind,2])
            current_geno <- old_geno
            #Ask if genotype is het, if yes, then add chance of allelic dropout error
            if(genos[ind,1] != genos[ind,2]){
                het <- paste0(genos[ind,1],"-", genos[ind,2])
                homo1 <- paste0(genos[ind,1],"-", genos[ind,1])
                homo2 <- paste0(genos[ind,2],"-", genos[ind,2])
                possible_genos =  c(het, homo1, homo2)
                new_geno <- sample(possible_genos, size = 1, prob = c(1 - 2*e_1, e_1, e_1))
                current_geno <- new_geno
                }
            # Getting the unique alleles into separate vectors again
            allele1 <- unlist(strsplit(current_geno, "-"))[1]
            allele2 <- unlist(strsplit(current_geno, "-"))[2]
            
            # Adding the error that isn't due to alleleic dropout, this class of error typically occurs at lower rate (effects less inds) but will turn the allele to any other allele at an equal rate --> assume these happen after class 1 errors
            new_allele1 <- sample(c(allele1, alleles[alleles %notin% allele1]), size = 1, prob = c(1-e_2, rep(e_2/(k-1), k - 1)))
            new_allele2 <- sample(c(allele2, alleles[alleles %notin% allele2]), size = 1, prob = c(1-e_2, rep(e_2/(k-1), k - 1)))
            
            # Putting the new genos into a data set
            locus_data_new[ind, 2*locus - 1] <- new_allele1
            locus_data_new[ind, 2*locus] <- new_allele2
            
            # Checking if the new genos are *functionally* diff from the old genos
            new_geno_v1 <- paste0(new_allele1, "-", new_allele2)
            new_geno_v2 <- paste0(new_allele1, "-", new_allele2)
            
            # Noting whether or not the genotype has changed and putting answer in the diffs_in_genos matrix
            geno_changed <- "F"
            if (old_geno %notin% c(new_geno_v1, new_geno_v2)) {
                geno_changed <- "T"
            }
            diffs_in_genos[ind, locus] <- geno_changed
            
            }
    }
    diffs_in_genos<- as.data.frame(lapply(diffs_in_genos, as.logical))
   
    return(list(locus_data_new, diffs_in_genos))     
}
    

# Adding error to the parental data
parental_data_error <- simulating_error(data_for_editing = parental_data, type1_rate = .1, type2_rate = .05, parental_data = parental_data)

parental_data_sim_error <- parental_data_error[[1]] # For input into Gibbs sampler

parental_data_track_errors <- parental_data_error[[2]] # For data summarization purposes
colSums(parental_data_track_errors)


# Adding error to the children data
children_data_error <- simulating_error(data_for_editing = children_data, type1_rate = .1, type2_rate = .05, parental_data = parental_data)

children_data_sim_error <- children_data_error[[1]] # For input into Gibbs sampler

children_data_track_errors <- children_data_error[[2]] # For data summarization purposes
colSums(children_data_track_errors)


#if I can't get my Gibbs sampler likelihood function to work fast enough --> switch out my likelihood function for a simpler likelihood function (maybe known dad) with a known posterior to prove that the Gibbs sampler works
    
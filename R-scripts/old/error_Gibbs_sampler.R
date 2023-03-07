#MSAT sampling code written by Kaylee Rosenberger and edited by Ash Hamilton
#Simulating error + Gibbs sampling code written by Ash Hamilton

#Librarys
library(adegenet)
library(poppr)
library(tidyr)
library(dplyr)
library(data.table)

#Best if this is the absolute path to this script
setwd("/Users/Ashley/Desktop/grad_school/2023_Winter/fund_comp_bio/final_project/Pollen_dispersal_sims")

#NEEDS TO BE ABSOLUTE PATH TO PLACE WHERE THE arlequin FILE IS
mydir = "/Users/Ashley/Desktop/grad_school/2023_Winter/fund_comp_bio/final_project/Pollen_dispersal_sims/Simulations/only_used_pop/"

####Creating functions for importing and sampling initial MSAT data####

# Create a function that is the opposite of %in%
`%notin%` <- Negate(`%in%`)

#A function that samples the parental population and generates children:
#   children generated WITH NO MUTATION by randomly choosing one of each allele from each parent for all loci
#   output is a list of matrices, the row of each matrix corresponds to genetic information related to that child in the offspring matrix:
#       the 1st matrix is the offspring genotype data
#       the 2nd matrix is the maternal genotype
#       the 3rd matrix is the paternal genotype
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

####Loading the already simulated population data####
num_loci = 20 #number of loci simulated, needed to make a data frame to save the data, starts with 20 --> gets cut down to 10 later

# Make a vector of all of the loci names
loci_names = c()
for(i in 1:num_loci){
    loci_names = c(loci_names, paste("locus", i, "a", sep=""))
    loci_names = c(loci_names, paste("locus", i, "b", sep=""))
}

#list of genalex files for all simulation replicates--genalex files end in .csv
genalex_list = list.files(mydir, ".csv$")

#cut off first 2 rows in data frame -- this is the population data, which is not required for our purposes
genetic_data = read.csv(paste(mydir, genalex_list, sep=""), header=FALSE)[-c(1,2),]

#giving the data frame columns new names
names(genetic_data) = c("Ind", "Pop", loci_names)

#this data frame now contains all of the data of 20 MSAT loci from each of 2500 individuals (in a single population) as simulated from SIM2COAL
genetic_data = genetic_data[-1,] #removing the first row -- repeat of now column headers

#thinning down to only 2 loci for now
#genetic_data_loci = genetic_data[,c(1:6)]

#thinning down to only 10 loci for now
genetic_data_loci = genetic_data[,c(1:22)]

####Simulating the offspring####

#Define scenario params
num_trees_to_sample = 2
num_seeds_to_sample = rep(5, num_trees_to_sample)
num_pollen_donors = 2
pollen_probability = c(.5, .5)


set.seed(2023)

#Simulate the offspring and parental dfs, outputs are in a list, each output is a df, each row corresponds to a unique child individual that is located at that row in the child data frame
#   First object in list is df with all simulated children alleles 
#   Second object in list is df with all simulated mother alleles
#   Third object in list is df with all simulated father alleles
simple_data <- sample_seed(genetic_data_loci, num_trees_to_sample = num_trees_to_sample, num_seeds_to_sample = num_seeds_to_sample, num_pollen_donors = num_pollen_donors, pollen_probability = pollen_probability) 

getwd()
save(simple_data, file="./simple_data.Rdata")
#save(simple_data, file="./R-scripts/Rdata/simple_data.Rdata")

####Cleaning the simulated data for easier use####

# Make sure I have the same dataset every time I reload this script
#getwd()
# If loading from /Users/Ashley/Desktop/grad_school/2023_Winter/fund_comp_bio/final_project/Pollen_dispersal_sims
#load(file="./R-scripts/Rdata/simpler_data.Rdata")
# If loading from my_dir
#load(file="../../R-scripts/Rdata/simple_data.Rdata")

children_data <- simple_data[[1]] # Making child data set its own object
mothers <- simple_data[[2]] # Making mother data set its own object
fathers <- simple_data[[3]] # Making father data set its own object

children_data$mom_ID <- mothers[,1] # Adding maternal IDs to child data for parentage assignment
children_data$dad_ID <- fathers[,1] # Adding paternal IDs to child data so I can test parentage assignment
children_data$child_ID <- seq.int(nrow(children_data))

# Creating a data set that contains all of the parents all mixed up, I feel this is likely to rep the data, maybe should set something to make sure no mom is also a dad....
parents <- c(unique(children_data$mom_ID), unique(children_data$dad_ID))
parental_data <- rbind(unique(mothers), unique(fathers)) 

children_data_true <- children_data
parental_data_true <- parental_data


# Now you have a children df and a parental df with all the loci! woot!

####Simulating error and adding it to the simulated reads####

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


# Error class 1 = allelic dropout aka failure to amplify one of 2 alleles, worse w/ lower DNA quality, only affects hets by making them into false homos
# Error class 2 = stochastic errors excluding allelic dropout, including false alleles, miscalling, and contaminant DNA

simulating_error <- function(data_for_editing, type1_rate, type2_rate, parental_data = NA){
    #browser()
    # Set data_for_editing to be the parental data if there is no parental data loaded
    if(missing(parental_data)){
        parental_data <- data_for_editing
    }
    
    # Pulling just the loci from the df
    locus_data <- data_for_editing[,grepl( "locus" , colnames(data_for_editing))]
    
    # Getting loci number from data 
    num_loci <- ncol(locus_data)/2
    
    # Making the new matrix with the same shape as the old matrix
    locus_data_new <- data.frame(matrix(ncol = ncol(locus_data), nrow = nrow(locus_data)))
    colnames(locus_data_new) <- colnames(locus_data)
    diffs_in_genos <- data.frame(matrix(ncol = ncol(locus_data)/2, nrow = nrow(locus_data)))
    
    for(locus in 1:num_loci){
        alleles <- get_alleles_and_freqs(locus, parental_data = parental_data)[[1]] # Want to get all possible alleles from the parental data
        k = length(alleles) # k = number of alleles at designated locus
        e_1 = type1_rate / (1 + type1_rate)
        e_2 = type2_rate / (k + type2_rate)
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
    data_for_editing[names(locus_data_new)] <- locus_data_new
    diffs_in_genos<- as.data.frame(lapply(diffs_in_genos, as.logical)) # Making this output as actual TRUE/FALSE values that R recognizes so I can easily count the number of occurrances
    
    return(list(data_for_editing, diffs_in_genos))     
}


# Adding error to the parental data
parental_data_error <- simulating_error(data_for_editing = parental_data, type1_rate = .2, type2_rate = .1, parental_data = parental_data)

parental_data_sim_error <- parental_data_error[[1]] # For input into Gibbs sampler

parental_data_track_errors <- parental_data_error[[2]] # For data summarization purposes
colSums(parental_data_track_errors)


# Adding error to the children data
children_data_error <- simulating_error(data_for_editing = children_data, type1_rate = .2, type2_rate = .1, parental_data = parental_data)

children_data_sim_error <- children_data_error[[1]] # For input into Gibbs sampler

children_data_track_errors <- children_data_error[[2]] # For data summarization purposes
colSums(children_data_track_errors)


#Saving the error data
error_data <- list(parental_data_error, children_data_error)

save(error_data, file="./error_data.Rdata")


####Making the likelihood function####

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
    
    log_true_geno_prob = log1p(true_geno_prob*true_geno_freq) # Probability needs to be scaled by the frequency of that genotype
    #putting log1p here so that a 0 doesn't throw errors
    
    return(log_true_geno_prob)
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
    ind_true_geno_log_probs <- rep(0, ncol(all_genos)) # Correlates to column of all_genos
    
    # Iterate through all possible "true" genotypes
    for(geno in 1:ncol(all_genos)){
        true_geno <- all_genos[,geno] # Set current "true genotype"
        
        ind_true_geno_log_probs[geno] <- calc_true_geno_prob(ind_geno = ind_geno, true_geno = true_geno, allele_freqs = allele_freqs, error_1 = error_1, error_2 = error_2)
        
    } # Close the for loop for all possible genotypes
    
    ind_true_geno_log_probs <- matrix(c(ind_true_geno_log_probs), nrow = 1) 
    colnames(ind_true_geno_log_probs) <- paste0(all_genos[1,], "-", all_genos[2,]) # Make the column name the genotypes
    
    #outer function might be faster here
    
    
    #ind_all_geno_all_loci_probs[[locus]] <- ind_true_geno_probs # Each loci = it's own element of the list
    
    return(ind_true_geno_log_probs)
}

calc_pr_child <- function(locus, mom_ID, dad_ID, child_ID, error_1, error_2){
    #browser()
    all_child_loci <- children_data[child_ID, -c(21:23)] # Pulling the single child and dropping the various IDs so just locus data remains
    
    # Pulling both alleles at a single locus from the ind
    child_geno <- all_child_loci[, c(2*locus - 1, 2*locus)] # Ind geno based on ind name
    
    allele_freqs <- get_alleles_and_freqs(locus, parental_data = parental_data)[[2]] 
    
    ind_geno <- NA
    log_prob_mom <- prob_hidden_genos(locus = locus, ind_ID = mom_ID, error_1 = error_1, error_2 = error_2) 
    
    ind_geno <- NA
    log_prob_dad <- prob_hidden_genos(locus = locus, ind_ID = dad_ID, error_1 = error_1, error_2 = error_2)
    
    # Making list of all true genotype parental pairings, true for all inds at each loci
    all_diff <- combn(colnames(log_prob_mom), 2) # All diff genotype pairings
    all_same <- rbind(colnames(log_prob_mom), colnames(log_prob_dad)) # All same genotype pairings
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
        
        log_pr_child_given_real_geno <- log1p(sum(expm1(obs_geno_prob) * (1/4))) # Prob observing the real geno given all possible true genos from the designated parental pairing, adding the expm1 and log1p here so 0's don't throw errors
        
        # Prob observing this child given obs mom given true mom in spot 1 of the pairing and obs dad given true dad in spot 2 of the pairing
        log_prob_this_mom <- log_prob_mom[colnames(log_prob_mom) == all_pairings[1,pairing]]
        log_prob_this_dad <- log_prob_dad[colnames(log_prob_dad) == all_pairings[2,pairing]]
        
        wrong_ll_child_given_mom_dad <- sum(log_prob_this_mom, log_prob_this_dad, log_pr_child_given_real_geno)
        
        #since each of the above are log(c + 1), they need to be undone before logging again so:
        ll_child_given_mom_dad <- log1p(expm1(wrong_ll_child_given_mom_dad) - (expm1(log_pr_child_given_real_geno) * expm1(log_prob_this_mom)) - (expm1(log_pr_child_given_real_geno) * expm1(log_prob_this_dad)) - (expm1(log_prob_this_dad) * expm1(log_prob_this_mom)) - expm1(log_pr_child_given_real_geno) - expm1(log_prob_this_mom) - expm1(log_prob_this_dad))
        
        
        # Prob observing this child given obs mom given true mom in spot 2 of the pairing and obs dad given true dad in spot 1 of the pairing
        log_prob_this_mom <- log_prob_mom[colnames(log_prob_mom) == all_pairings[2,pairing]]
        log_prob_this_dad <- log_prob_dad[colnames(log_prob_dad) == all_pairings[1,pairing]]
        
        wrong_ll_child_given_dad_mom <- sum(log_prob_this_mom, log_prob_this_dad, log_pr_child_given_real_geno)
        
        #since each of the above are log(c + 1), they need to be undone before logging again so:
        ll_child_given_dad_mom <- log1p(expm1(wrong_ll_child_given_dad_mom) - (expm1(log_pr_child_given_real_geno) * expm1(log_prob_this_mom)) - (expm1(log_pr_child_given_real_geno) * expm1(log_prob_this_dad)) - (expm1(log_prob_this_dad) * expm1(log_prob_this_mom)) - expm1(log_pr_child_given_real_geno) - expm1(log_prob_this_mom) - expm1(log_prob_this_dad))
        
        ll_child_given_any_mom_dad[pairing] <- ll_child_given_mom_dad
        
        ll_child_given_any_dad_mom[pairing] <- ll_child_given_dad_mom
    }
    
    total_ll_child_given_mom_dad <- log1p(sum(expm1(ll_child_given_any_mom_dad)))
    total_ll_child_given_dad_mom <- log1p(sum(expm1(ll_child_given_any_dad_mom)))
    total_ll_child_given_parents <- log1p(sum(expm1(total_ll_child_given_dad_mom), expm1(total_ll_child_given_mom_dad)))
    
    return(total_ll_child_given_parents)
}

#unclear how to incorporate allele frequencies that differ at each loci
get_ll <- function(children_data, parental_data, error_1, error_2) {
    dad_data <- subset(parental_data, parental_data$Ind %notin% children_data$mom_ID) # Assumes a parent can only ever be a mother or a father
    
    num_loci <- (ncol(parental_data)-1)/2
    
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
            ll_O_given_one_DM[dad, child] <- log1p(prod(expm1(ll_O_given_DM_at_each_locus))) # Since 0's might be around, deciding to multiply (instead of sum....)  across all loci (assuming independence) and then logging again to get the ll of each dad for each mother/off pair
        }
    }
    
    total_ll <- log(sum(colSums(expm1(ll_O_given_one_DM)))) - log(nrow(ll_O_given_one_DM)) # Overall likelihood of the data (mother and offspring pairs) and the parameters, first sum over columns to sum across all fathers then divide by the probability of each father being the correct father which is the same as number of fathers and therefore the number of rows in the ll_O_given_one_DM df (assumes each father has equal prob to contribute to each mother-off pair)
    return(list(total_ll, ll_O_given_one_DM))
}

#Practice run of the likelihood function
#get_ll(children_data = children_data_sim_error, parental_data = parental_data_sim_error, error_1 = .2, error_2 = .1)[[1]]


####Making the Gibbs sampler####
# Uniform prior with a log on it (aka both of my error priors)
log_prior = function(p){
    if((p<0) || (p>1)){  # || here means "or", if p is smaller than 0 or bigger than 1 return -inf (which when exponentiated returns 0), otherwise return 0 (which once exponentiated equals 1) when this function is called, basically just makes sure all evaluated p values between 1 and 0 have the same prior odds of 1 and anything outside of that has a prior odds of 0
        return(-Inf)}
    else{
        return(0)}
}


error_sampler = function(children_data, parental_data, niter, error_startvals, proposalsd, testing){
    #browser()

    error = matrix(0,ncol = niter, nrow = 2)
    error[1,1] = error_startvals[1]
    error[2,1] = error_startvals[2]
    for(i in 2:niter){
        current_error1 = error[1, i-1]
        current_error2 = error[2, i-1]
        new_error1 = current_error1 + rnorm(1,0,proposalsd) #this isn't good because error can't be negative... to make symmetric call the rnorm and then write an if else that asks if the resulting new_error would negative, if yes then subtract the value from it to get back above 1?
        # If not in testing mode, run get_ll function to obtain real likliehoods
        if(testing == F | missing(testing)){
            A = exp(log_prior(new_error1) + get_ll(children_data = children_data, parental_data = parental_data, error_1 = new_error1, error_2 = current_error2)[[1]] - log_prior(current_error1) - get_ll(children_data = children_data, parental_data = parental_data, error_1 = current_error1, error_2 = current_error2)[[1]]) 
        }
        # If in testing mode, use a uniform dist'n with no dep on error to obtain likelihoods
        else(A = exp(log_prior(new_error1) + log(1) - log_prior(current_error1) - log(1)) )
        
        
        if(runif(1)<A){ 
            error[1,i] = new_error1 # Accept move with probability min(1,A)
        } else {
            error[1,i] = current_error1 # Otherwise "reject" move, and stay where we are
        }
        current_error1 = error[1,i] # Set the new current error1 that takes into account the last half step  
        new_error2 = current_error2 + rnorm(1,0,proposalsd) 
        
        # If not in testing mode, run get_ll function to obtain real likliehoods
        if(testing == F | missing(testing)){
            A = exp(log_prior(new_error2) + get_ll(children_data = children_data, parental_data = parental_data, error_1 = current_error1, error_2 = new_error2)[[1]] - log_prior(current_error2) - get_ll(children_data = children_data, parental_data = parental_data, error_1 = current_error1, error_2 = current_error2)[[1]])
        }
        # If in testing mode, use a uniform dist'n with no dep on error to obtain likelihoods
        else(A = exp(log_prior(new_error2) + log(1) - log_prior(current_error2) - log(1)) )
        
        
        if(runif(1)<A){ 
            error[2,i] = new_error2 # Accept move with probability min(1,A)
        } else {
            error[2,i] = current_error2 # Otherwise "reject" move, and stay where we are
        }
    }
    return(error)
}


####Testing the Gibbs sampler####
first_testing_Gibbs_try <- error_sampler(children_data = children_data_sim_error, parental_data = parental_data_sim_error, niter = 100000, error_startvals = c(.35,.05), proposalsd = .05, testing = T)

first_testing_Gibbs_try_noburnins <- first_testing_Gibbs_try[,-c(1:10000)]

# Plotting test of the Gibbs machinery!
plot(density(first_testing_Gibbs_try_noburnins[1,]), xlim= c(0.0001,1), xlab = "Error value", main = "90000 iters sampled from a uniform tagret dist'n")
abline(v = mean(first_testing_Gibbs_try_noburnins[1,]), lty = 2)
lines(density(first_testing_Gibbs_try_noburnins[2,]), add = T, col = "grey")
abline(v = mean(first_testing_Gibbs_try_noburnins[2,]), lty = 2, col = "grey")
curve(dunif(x), from = 0, to = 1, add = T, col = "blue")
abline(v = .5, col = "blue", lty = 2)
legend("bottomright", 
       legend = c("Sampler", "Sampler Mean", "Uniform dist'n", "Uniform mean"), 
       col = c("black", "black", "blue", "blue"), 
       lty = c(1, 2, 1, 2), 
       inset = .05)


####Running the Gibbs sampler####
first_big_Gibbs_try <- error_sampler(children_data = children_data_sim_error, parental_data = parental_data_sim_error, niter = 350, error_startvals = c(.35,.05), proposalsd = .05) # Attempting 350 iterations, would take about 24 hours to run on my laptop
save(first_big_Gibbs_try, file="./first_big_Gibbs_try.Rdata")

#real type1_rate = .2, type2_rate = .1



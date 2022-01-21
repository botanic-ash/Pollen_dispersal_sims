###
#This script fills a list containing sets of function parameters (scenarios for 'simulation')
#The loop goes through each number of maternal trees (50, 25, 10, 5, 2, 1) sampling a different
#number of seeds per tree, from 1 to (250/number of maternal trees to sample) 
#each loop iteration, three 'scenarios' of parameters are made--seed sets made by a single father (all same),
#seed sets made by multiple fathers with equal donation (all unique), and seed sets created by mutliple
#fathers with skewed pollen donation (skewed--one donates 80%, the rest donate evenly)
#
###############################################################################################
#Fixed seeds/tree
#parameters are in the order: dataset (genalex file--not included here, will be added in main loop), 
#trees to sample, seeds to sample, pollen donors, pollen probability


combined_list_params = list() #used for saving scenario parameters--ends up being a list of lists
#each scenario can be accessed using combined_list_params[[x]] where x is the scenario desired (goes in order of the values on the table)
num_maternal_trees = c(50, 25, 10, 5, 2, 1) #different number of maternal trees to be sampled for each scenario

x = 1 #x is the list counter variable that names each of the lists
for(i in 1:length(num_maternal_trees)) { #loops over the vector of maternal trees
  for(j in 1:(250/num_maternal_trees[i])) { #loops ofrom 1 to max number number of seeds to sample per maternal tree
    temp = list(num_maternal_trees[i], j, 1, 1) #All Same scenario--all seeds created from 1 pollen donor
    combined_list_params[[x]] = temp #saving the parameters to a list (list of lists)
    x=x+1 #increment counter
    temp = list(num_maternal_trees[i], j, j, c(rep((1/j), j))) #All Unique* scenario--each pollen donor has equal probability to pollinate seeds
    combined_list_params[[x]] = temp #saving parameters
    x=x+1 #increment
    #Skewed scenarios--multiple pollen donors, but 1 donates the majority of the pollen (80%), the rest have equal probability to donate
    if(j==1){
      temp = list(num_maternal_trees[i], j, 1, 1) #if there's only one seed, it can't be skewed
    } else if (j!=1) {
      temp = list(num_maternal_trees[i], j, j, c(0.8, rep((0.2/(j-1)), (j-1)))) # if more than 1, then skew pollen donation
    } 
    combined_list_params[[x]] = temp
    x=x+1
  }
}

#combined_list_params ends up having 1395 elements! 


#***Note: all unique scenario currently is not that each seed sampled from a tree has a different father
#it is just that each pollen donor (pollen donors = number of seeds to sample) has equal probability to 
#donate pollen on the tree. Thus, there could be replicates since the father is chosen with sample()
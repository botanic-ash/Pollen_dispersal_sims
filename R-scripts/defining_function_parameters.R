#Code written by Kaylee Rosenberger#Code written by Kaylee Rosenberger
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
#num maternal trees to sample, num seeds to sample, pollen donors, pollen probability
setwd("C:/Users/kayle/Documents/Pollen_dispersal_sims/R-scripts")

all_same_params = list() #used for saving scenario parameters--ends up being a list of lists
#each scenario can be accessed using all_same_params[[x]] where x is the scenario desired (goes in order of the values on the table)
all_eligible_params = list()
skewed_params = list()
num_maternal_trees = c(100, 50, 25, 10, 5, 2, 1) #different number of maternal trees to be sampled for each scenario

#All same scenarios
x = 1 #x is the list counter variable that names each of the lists
#this makes it easier to save the data, we can't use i because it only goes up to 5, and j varies each time it iterates
#x keeps track of the position in the list across each of the loops
for(i in 1:length(num_maternal_trees)) { #loops over the vector num_maternal_trees
  for(j in 1:(500/num_maternal_trees[i])) { #loops from 1 to max number number of seeds to sample per maternal tree
    temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), 1, 1) #All Same scenario--all seeds created from 1 pollen donor
    if((round(sum(temp[[4]]), 4))!=1){
      print(paste("pollen probability = ", sum(temp[[4]])))
      break
    }
    all_same_params[[x]] = temp #saving the parameters to a list (list of lists)
    x=x+1 #increment counter
  }
}

#All eligible scenarios
x = 1 #x is the list counter variable that names each of the lists
#this makes it easier to save the data, we can't use i because it only goes up to 5, and j varies each time it iterates
#x keeps track of the position in the list across each of the loops
for(i in 1:length(num_maternal_trees)) { #loops over the vector of maternal trees
  for(j in 1:(500/num_maternal_trees[i])) { #loops ofrom 1 to max number number of seeds to sample per maternal tree
    temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), j, c(rep((1/j), j))) #All Unique* scenario--each pollen donor has equal probability to pollinate seeds
    if((round(sum(temp[[4]]), 4))!=1){
      paste("pollen probability = ", sum(temp[[4]]))
      break
    }
    all_eligible_params[[x]] = temp #saving parameters
    x=x+1 #increment
  }
}

#Skewed scenarios
x = 1#x is the list counter variable that names each of the lists
#this makes it easier to save the data, we can't use i because it only goes up to 5, and j varies each time it iterates
#x keeps track of the position in the list across each of the loops
for(i in 1:length(num_maternal_trees)) { #loops over the vector of maternal trees
  for(j in 1:(500/num_maternal_trees[i])) { #loops ofrom 1 to max number number of seeds to sample per maternal tree
    temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), 10, c(0.6, 0.2, rep(0.05, 3), rep(0.01, 5))) 
    if((round(sum(temp[[4]]), 4))!=1){
      paste("pollen probability = ", sum(temp[[4]]))
      break
    }
    skewed_params[[x]] = temp #saving parameters
    x=x+1
  }
}

#465 scenarios in each parameter list
#so 1395 scenarios total are created
#saving the list in an Rdata file 
save(all_same_params, all_eligible_params, skewed_params, file="combined_list_params_711.Rdata")

#***Note: all unique scenario currently is not that each seed sampled from a tree has a different father
#it is just that each pollen donor (pollen donors = number of seeds to sample) has equal probability to 
#donate pollen on the tree. Thus, there could be replicates since the father is chosen with sample()

##################################################################################################
#Skewed seeds per tree
#params: trees to sample, seeds to sample, pollen donors, pollen probability

all_same_params = list() #used for saving scenario parameters--ends up being a list of lists
#each scenario can be accessed using all_same_params[[x]] where x is the scenario desired (goes in order of the values on the table)
all_eligible_params = list()
skewed_params = list()
num_maternal_trees_skewed = c(100, 50, 25, 10, 5, 5, 2)

#100 maternal trees
x = 1
total_seeds = c(seq(200, 400, by=200)) #manually saving the number of seeds total to sample
for(i in 1:length(total_seeds)) {
  #all same
  temp = list(num_maternal_trees_skewed[1], c((0.50*total_seeds[i]), (0.01*total_seeds[i]), rep((0.005*total_seeds[i]), 98)), 1, 1) #All same
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_same_params[[x]] = temp #saving parameters
  #all eligible
  temp = list(num_maternal_trees_skewed[1], c((0.50*total_seeds[i]), (0.01*total_seeds[i]), rep((0.005*total_seeds[i]), 98)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i]))) #All eligible
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_eligible_params[[x]] = temp
  #skewed
  temp = list(num_maternal_trees_skewed[1], c((0.50*total_seeds[i]), (0.01*total_seeds[i]), rep((0.005*total_seeds[i]), 98)), 10, c(0.6, 0.2, rep(0.05, 3), rep(0.01, 5))) #skewed
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  skewed_params[[x]] = temp
  x=x+1
}

#50 maternal trees
total_seeds = c(seq(100, 500, by=100)) #manually saving the number of seeds total to sample
for(i in 1:length(total_seeds)) {
  #all same
  temp = list(num_maternal_trees_skewed[2], c((0.51*total_seeds[i]), rep((0.01*total_seeds[i]), 49)), 1, 1) #All same
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_same_params[[x]] = temp #saving parameters
  #all eligible
  temp = list(num_maternal_trees_skewed[2], c((0.51*total_seeds[i]), rep((0.01*total_seeds[i]), 49)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i]))) #All eligible
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_eligible_params[[x]] = temp
  #skewed
  temp = list(num_maternal_trees_skewed[2], c((0.51*total_seeds[i]), rep((0.01*total_seeds[i]), 49)), 10, c(0.6, 0.2, rep(0.05, 3), rep(0.01, 5))) #skewed
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  skewed_params[[x]] = temp
  x=x+1
}

#25 maternal trees
total_seeds = c(seq(50, 500, by=50))
for(i in 1:length(total_seeds)) {
  temp = list(num_maternal_trees_skewed[3], c((0.52*total_seeds[i]), rep((0.02*total_seeds[i]), 24)), 1, 1) #all same
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_same_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[3], c((0.52*total_seeds[i]), rep((0.02*total_seeds[i]), 24)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i]))) #all unique
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_eligible_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[3], c((0.52*total_seeds[i]), rep((0.02*total_seeds[i]), 24)), 10, c(0.6, 0.2, rep(0.05, 3), rep(0.01, 5))) #skewed
  if((round(sum(temp[[4]]), 4))!=1){
    paste("pollen probability = ", sum(temp[[4]]))
    break
  }
  skewed_params[[x]] = temp
  x=x+1
}

#10 maternal trees
total_seeds = c(seq(20, 500, by=20))
for(i in 1:length(total_seeds)) {
  temp = list(num_maternal_trees_skewed[4], c((0.5*total_seeds[i]), (0.1*total_seeds[i]), rep((0.05*total_seeds[i]), 8)), 1, 1)
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_same_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[4], c((0.5*total_seeds[i]), (0.1*total_seeds[i]), rep((0.05*total_seeds[i]), 8)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i])))
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_eligible_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[4], c((0.5*total_seeds[i]), (0.1*total_seeds[i]), rep((0.05*total_seeds[i]), 8)), 10, c(0.6, 0.2, rep(0.05, 3), rep(0.01, 5)))
  if((round(sum(temp[[4]]), 4))!=1){
    paste("pollen probability = ", sum(temp[[4]]))
    break
  }
  skewed_params[[x]] = temp
  x=x+1
}

#5 maternal trees--max 80% skew
total_seeds = c(seq(20,500,by=20))
for(i in 1:length(total_seeds)) {
  temp = list(num_maternal_trees_skewed[5], c((0.8*total_seeds[i]), rep((0.05*total_seeds[i]), 4)), 1, 1) #all same
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_same_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[5], c((0.8*total_seeds[i]), rep((0.05*total_seeds[i]), 4)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i])))
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_eligible_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[5], c((0.8*total_seeds[i]), rep((0.05*total_seeds[i]), 4)), 10, c(0.6, 0.2, rep(0.05, 3), rep(0.01, 5)))
  if((round(sum(temp[[4]]), 4))!=1){
    paste("pollen probability = ", sum(temp[[4]]))
    break
  }
  skewed_params[[x]] = temp
  x=x+1
}

#5 maternal trees--max 50% skew
total_seeds = c(seq(10,500,by=10))
for(i in 1:length(total_seeds)) {
  temp = list(num_maternal_trees_skewed[6], c((0.5*total_seeds[i]), (0.2*total_seeds[i]), rep((0.1*total_seeds[i]), 3)), 1, 1)#all same
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_same_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[6], c((0.5*total_seeds[i]), (0.2*total_seeds[i]), rep((0.1*total_seeds[i]), 3)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i])))#all unique
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_eligible_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[6], c((0.5*total_seeds[i]), (0.2*total_seeds[i]), rep((0.1*total_seeds[i]), 3)), 10, c(0.6, 0.2, rep(0.05, 3), rep(0.01, 5)))#skewed
  if((round(sum(temp[[4]]), 4))!=1){
    paste("pollen probability = ", sum(temp[[4]]))
    break
  }
  skewed_params[[x]] = temp
  x=x+1
}

#2 maternal trees
total_seeds = c(seq(5,500,by=5))
for(i in 1:length(total_seeds)) {
  temp = list(num_maternal_trees_skewed[7], c((0.8*total_seeds[i]), (0.2*total_seeds[i])), 1, 1) # all same
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_same_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[7], c((0.8*total_seeds[i]), (0.2*total_seeds[i])), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i]))) # all unique
  if((round(sum(temp[[4]]), 4))!=1){
    print(paste("pollen probability = ", sum(temp[[4]])))
    break
  }
  all_eligible_params[[x]] = temp
  temp = list(num_maternal_trees_skewed[7], c((0.8*total_seeds[i]), (0.2*total_seeds[i])),  10, c(0.6, 0.2, rep(0.05, 3), rep(0.01, 5))) # skewed
  if((round(sum(temp[[4]]), 4))!=1){
    paste("pollen probability = ", sum(temp[[4]]))
    break
  }
  skewed_params[[x]] = temp
  x=x+1
}

#combined_list_params_skewed ends up having 261 elements 
#so 261 scenarios are created 
#saving the list in an Rdata file
save(all_same_params, all_eligible_params, skewed_params, file="combined_list_params_skewed_new.Rdata")

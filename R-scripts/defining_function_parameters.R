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


all_same_params = list() #used for saving scenario parameters--ends up being a list of lists
#each scenario can be accessed using all_same_params[[x]] where x is the scenario desired (goes in order of the values on the table)
all_eligible_params = list()
skewed_params = list()
num_maternal_trees = c(50, 25, 10, 5, 2, 1) #different number of maternal trees to be sampled for each scenario

#All same scenarios
x = 1 #x is the list counter variable that names each of the lists
for(i in 1:length(num_maternal_trees)) { #loops over the vector of maternal trees
  for(j in 1:(250/num_maternal_trees[i])) { #loops ofrom 1 to max number number of seeds to sample per maternal tree
    temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), 1, 1) #All Same scenario--all seeds created from 1 pollen donor
    all_same_params[[x]] = temp #saving the parameters to a list (list of lists)
    x=x+1 #increment counter
  }
}

#All eligible scenarios
x = 1
for(i in 1:length(num_maternal_trees)) { #loops over the vector of maternal trees
  for(j in 1:(250/num_maternal_trees[i])) { #loops ofrom 1 to max number number of seeds to sample per maternal tree
    temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), j, c(rep((1/j), j))) #All Unique* scenario--each pollen donor has equal probability to pollinate seeds
    all_eligible_params[[x]] = temp #saving parameters
    x=x+1 #increment
  }
}

#Skewed scenarios
x = 1
for(i in 1:length(num_maternal_trees)) { #loops over the vector of maternal trees
  for(j in 1:(250/num_maternal_trees[i])) { #loops ofrom 1 to max number number of seeds to sample per maternal tree
    #Skewed scenarios--multiple pollen donors, but 1 donates the majority of the pollen (80%), the rest have equal probability to donate
    if(j==1) {
      temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), 1, 1) #if there's only one seed, it can't be skewed
    } else if (j==2) {
      temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), j, c(0.8, 0.2)) # if only 2 pollen donors, then 1 80% and 1 20% 
    } else if (j==3) {
      temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), j, c(0.4, 0.4, 0.2))
    } else if (j==4) {
      temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), j, c(0.4, 0.4, 0.1, 0.1))
    } else if (j==5) {
      temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), j, c(0.4, 0.4, 0.1, 0.05, 0.05))
    } else if (j>=6&&j<14) {
      temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), j, c(0.4, 0.4, 0.5, 0.5, rep((0.1/(j-4)), (j-4)))) 
    } else if (j>=14) {
      temp = list(num_maternal_trees[i], c(rep(j, (num_maternal_trees[i]))), 14, c(0.4, 0.4, 0.5, 0.5, rep((0.1/10), (10)))) #cap on pollen donors is 14 total (such that the lowest probability to donate pollen is 1%)
    }
    skewed_params[[x]] = temp #saving parameters
    x=x+1
  }
}

#465 scenarios in each parameter list
#so 1395 scenarios total are created
#saving the list in an Rdata file 
save(all_same_params, all_eligible_params, skewed_params, file="combined_list_params.Rdata")

#***Note: all unique scenario currently is not that each seed sampled from a tree has a different father
#it is just that each pollen donor (pollen donors = number of seeds to sample) has equal probability to 
#donate pollen on the tree. Thus, there could be replicates since the father is chosen with sample()

##################################################################################################
#Skewed seeds per tree
#params: trees to sample, seeds to sample, pollen donors, pollen probability

combined_list_params_skewed = list()
num_maternal_trees_skewed = c(50, 25, 10, 5, 5, 2)

#50 maternal trees
x = 1
total_seeds = c(seq(100, 200, by=50)) #manually saving the number of seeds total to sample
for(i in 1:3) {
  temp = list(num_maternal_trees_skewed[1], c((0.51*total_seeds[i]), rep((0.01*total_seeds[i]), 49)), 1, 1) #All same
  combined_list_params_skewed[[x]] = temp #saving parameters
  x=x+1
  temp = list(num_maternal_trees_skewed[1], c((0.51*total_seeds[i]), rep((0.01*total_seeds[i]), 49)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i]))) #All unique
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[1], c((0.51*total_seeds[i]), rep((0.01*total_seeds[i]), 49)), total_seeds[i], c(0.8, rep((0.2/total_seeds[i]), (total_seeds[i]-1))))
  combined_list_params_skewed[[x]] = temp
  x=x+1
}

#25 maternal trees
total_seeds = c(seq(50, 200, by=50))
for(i in 1:4) {
  temp = list(num_maternal_trees_skewed[2], c((0.52*total_seeds[i]), rep((0.02*total_seeds[i]), 24)), 1, 1) #all same
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[2], c((0.52*total_seeds[i]), rep((0.02*total_seeds[i]), 24)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i]))) #all unique
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[2], c((0.52*total_seeds[i]), rep((0.02*total_seeds[i]), 24)), total_seeds[i], c(0.8, rep((0.2/total_seeds[i]), (total_seeds[i]-1)))) #skewed
  combined_list_params_skewed[[x]] = temp
  x=x+1
}

#10 maternal trees
total_seeds = c(seq(20, 200, by=20))
for(i in 1:10) {
  temp = list(num_maternal_trees_skewed[3], c((0.5*total_seeds[i]), (0.1*total_seeds[i]), rep((0.05*total_seeds[i]), 8)), 1, 1)
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[3], c((0.5*total_seeds[i]), (0.1*total_seeds[i]), rep((0.05*total_seeds[i]), 8)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i])))
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[3], c((0.5*total_seeds[i]), (0.1*total_seeds[i]), rep((0.05*total_seeds[i]), 8)), total_seeds[i], c(0.8, rep((0.2/total_seeds[i]), (total_seeds[i]-1))))
  combined_list_params_skewed[[x]] = temp
  x=x+1
}

#5 maternal trees--max 80% skew
total_seeds = c(seq(20,200,by=20))
for(i in 1:10) {
  temp = list(num_maternal_trees_skewed[4], c((0.8*total_seeds[i]), rep((0.05*total_seeds[i]), 4)), 1, 1) #all same
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[5], c((0.8*total_seeds[i]), rep((0.05*total_seeds[i]), 4)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i])))
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[5], c((0.8*total_seeds[i]), rep((0.05*total_seeds[i]), 4)), total_seeds[i], c(0.8, rep((0.2/total_seeds[i]), (total_seeds[i]-1))))
  combined_list_params_skewed[[x]] = temp
  x=x+1
}

#5 maternal trees--max 50% skew
total_seeds = c(seq(10,200,by=10))
for(i in 1:20) {
  temp = list(num_maternal_trees_skewed[5], c((0.5*total_seeds[i]), (0.2*total_seeds[i]), rep((0.1*total_seeds[i]), 3)), 1, 1)#all same
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[5], c((0.5*total_seeds[i]), (0.2*total_seeds[i]), rep((0.1*total_seeds[i]), 3)), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i])))#all unique
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[5], c((0.5*total_seeds[i]), (0.2*total_seeds[i]), rep((0.1*total_seeds[i]), 3)), total_seeds[i], c(0.8, rep((0.2/total_seeds[i]), (total_seeds[i]-1))))#skewed
  combined_list_params_skewed[[x]] = temp
  x=x+1
}

#2 maternal trees
total_seeds = c(seq(5,200,by=5))
for(i in 1:40) {
  temp = list(num_maternal_trees_skewed[6], c((0.8*total_seeds[i]), (0.2*total_seeds[i])), 1, 1) # all same
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[6], c((0.8*total_seeds[i]), (0.2*total_seeds[i])), total_seeds[i], c(rep((1/total_seeds[i]), total_seeds[i]))) # all unique
  combined_list_params_skewed[[x]] = temp
  x=x+1
  temp = list(num_maternal_trees_skewed[6], c((0.8*total_seeds[i]), (0.2*total_seeds[i])),  total_seeds[i], c(0.8, rep((0.2/total_seeds[i]), (total_seeds[i]-1)))) # skewed
  combined_list_params_skewed[[x]] = temp
  x=x+1
}

#combined_list_params_skewed ends up having 261 elements 
#so 261 scenarios are created 
#saving the list in an Rdata file
save(combined_list_params_skewed, file="combined_list_params_skewed.Rdata")

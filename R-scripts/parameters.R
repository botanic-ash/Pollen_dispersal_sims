#This script fills the list containing sets of function parameters (scenarios for 'simulation')
###############################################################################################
#Fixed seeds/tree
#parameters are in the order: dataset (genalex file--not included here, will be added in main loop), 
#trees to sample, seeds to sample, pollen donors, pollen probability

#50 maternal trees
x = 1 #x is the list counter variable that names each of the lists
for(i in 1:5) {
  assign(paste("list", x, sep=""), list(50, i, 1, 1)) #all same
  x=x+1
  assign(paste("list", x, sep=""), list(50, i, i, c(rep((1/i), i)))) #all unique*
  x=x+1
  #skewed
  if(i==1){
    assign(paste("list", x, sep=""), list(50, i, 1, 1)) #if there's only one seed, it can't be skewed
  } else if (i!=1) {
    assign(paste("list", x, sep=""), list(50, i, i, c(0.8, rep((0.2/(i-1)), (i-1)))))
  } 
  x=x+1
}

#25 maternal trees
for(i in 1:10){
  assign(paste("list", x, sep=""), list(25,i,1,1)) #all same
  x=x+1
  assign(paste("list", x, sep=""), list(25,i,i, c(rep((1/i),i)))) #all unique*
  x=x+1
  #skewed
  if(i==1){
    assign(paste("list", x, sep=""), list(25, i, 1, 1)) #if there's only one seed, it can't be skewed
  } else if (i!=1) {
    assign(paste("list", x, sep=""), list(25, i, i, c(0.8, rep((0.2/(i-1)), (i-1)))))
  } 
  x=x+1
}

#10 maternal trees
for(i in 1:25){
  assign(paste("list", x, sep=""), list(10,i,1,1)) #all same
  x=x+1
  assign(paste("list", x, sep=""), list(10,i,i, c(rep((1/i),i)))) #all unique*
  x=x+1
  #skewed
  if(i==1){
    assign(paste("list", x, sep=""), list(10, i, 1, 1)) #if there's only one seed, it can't be skewed
  } else if (i!=1) {
    assign(paste("list", x, sep=""), list(10, i, i, c(0.8, rep((0.2/(i-1)), (i-1)))))
  } 
  x=x+1
}

#5 maternal trees
for(i in 1:50){
  assign(paste("list", x, sep=""), list(5,i,1,1)) #all same
  x=x+1
  assign(paste("list", x, sep=""), list(5,i,i, c(rep((1/i),i)))) #all unique*
  x=x+1
  #skewed
  if(i==1){
    assign(paste("list", x, sep=""), list(5, i, 1, 1)) #if there's only one seed, it can't be skewed
  } else if (i!=1) {
    assign(paste("list", x, sep=""), list(5, i, i, c(0.8, rep((0.2/(i-1)), (i-1)))))
  } 
  x=x+1  
}

#2 maternal trees
for(i in 1:125){
  assign(paste("list", x, sep=""), list(2,i,1,1)) #all same
  x=x+1
  assign(paste("list", x, sep=""), list(2,i,i, c(rep((1/i),i)))) #all unique*
  x=x+1
  #skewed
  if(i==1){
    assign(paste("list", x, sep=""), list(2, i, 1, 1)) #if there's only one seed, it can't be skewed
  } else if (i!=1) {
    assign(paste("list", x, sep=""), list(2, i, i, c(0.8, rep((0.2/(i-1)), (i-1)))))
  } 
  x=x+1 
}

#1 maternal tree
for(i in 1:250){
  assign(paste("list", x, sep=""), list(1,i,1,1)) #all same
  x=x+1
  assign(paste("list", x, sep=""), list(1,i,i, c(rep((1/i),i)))) #all unique*
  x=x+1
  #skewed
  if(i==1){
    assign(paste("list", x, sep=""), list(1, i, 1, 1)) #if there's only one seed, it can't be skewed
  } else if (i!=1) {
    assign(paste("list", x, sep=""), list(1, i, i, c(0.8, rep((0.2/(i-1)), (i-1)))))
  } 
  x=x+1 
}


#fixed_scenarios = list(list1, list2, list3...)
#1395 lists of parameters total

###############################################################################################
#skewed seeds/tree
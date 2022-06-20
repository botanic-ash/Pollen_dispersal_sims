#Counting the total number of alleles present in the parental dataset (total) 
#while ignoring super rare alleles (frequency less than 2 by locus)
#We also calculate the number of alleles captured in the seeds sampled dataset here (captured)
#We make sure that none of the super rare alleles ignored in the parental set are included in the seed set 

num_loci=20#number of loci simulated
total = 0 #sum to keep track of total alleles 
captured = 0 #sum to keep track of alleles captured by sampling
k=3 #counter variable for column (locus) of parental dataset
j=1 #counter variable to column (locus) of seed dataset 

for(i in 1:num_loci) {
  parental_allele_list = table(c(as.matrix(data[,k:k+1]))) #getting alleles and their frequencies for locus i in parental dataframe
  parental_allele_list=parental_allele_list[parental_allele_list>3] #subsetting parental data to only include alleles with frequency greater than 3
  total_names = names(parental_allele_list)
  total = total + n_distinct(names(parental_allele_list)) #getting the number of distinct values for locus 1 to count alleles 
  
  seed_allele_list = table(c(as.matrix(temp[,j:j+1]))) #getting unique values for locus i in seed dataframe 
  captured_names = names(seed_allele_list)#getting the names of the alleles captured from sampling
  captured = captured + n_distinct(intersect(captured_names,total_names))#making sure none of the super rare alleles excluded from parental dataset are included here
  
  k = k+2 #increment k for next loop iteration
  j = j+2 #increment j for next loop iteration 
}

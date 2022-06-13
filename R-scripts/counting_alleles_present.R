#Counting the total number of alleles present in the parental dataset (total)
num_loci=20#number of loci simulated
total = 0 #sum to keep track of total alleles 
k=3 #starting from 3rd column
for(i in 1:num_loci) {
  temp1 = unique(data[,k]) #getting unique values in column k for locus i (locusiA)
  k = k+1 #increment k to move to next column for locus i
  temp2 = unique(data[,k]) #getting unique values in column k for locus i (locusiB)
  temp3 = c(temp1,temp2) #concatenating both vectors into one
  total = total + n_distinct(temp3) #getting the number of distinct values for locus 1 to count alleles 
  k = k+1 #increment k for next loop iteration
}


#counting the number of loci in the seeds sampled dataset (captured)
captured = 0
k=1
for(i in 1:num_loci) {
  temp1 = unique(temp[,k]) #getting unique values in column k for locus i (locusiA)
  k = k+1 #increment k to move to next column for locus i
  temp2 = unique(temp[,k]) #getting unique values in column k for locus i (locusiB)
  temp3 = c(temp1,temp2) #concatenating both vectors into one
  captured = captured + n_distinct(temp3) #getting the number of distinct values for locus 1 to count alleles 
  k = k+1 #increment k for next loop iteration
}

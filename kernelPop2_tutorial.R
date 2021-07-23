

library(rlang)
library(devtools)

#downloading github package
install_github("stranda/kernelPop2")

#using package
library(kernelPop2)

###############################################################################################################
#kernelPop2 tutorial 

###FIRST STEP:
#creating a skeleton landscape object 
#landscape object changes over time, each object represents the landscape at a point in time 
land = landscape.new.empty()
#printing out organization structures of the landscape object 
names(land)

###SECOND STEP:
#Initializing parameter values for the landscape 

#intparam describes integer values such as: number of habitats (h), number of demographic stages 
#present in a habitat (s), number of genetic loci to be simulated (locusnum), etc...
#giving the landscape object some integer parameters-- 4 habitats and 6 demographic stages
land = landscape.new.intparam(land, h=4, s=6)

#floatparams describe float type parameters such as selfing rate (s), and seed/pollen dispersal characteristics
#giving the landscape object the default values for floatparams
land = landscape.new.floatparam(land)

#switchparam are boolean parameters about the landscape 
#giving default values for swtichparams
land = landscape.new.switchparam(land)

#demography is comprised of sub-objects that determine survival and reproduction

#localdem is a description of the sub matrices that describe the demography in each habitat
#it's a list of length equal to the number of habitats
#here we have 6x6 matrices, because we defined 6 stages of life (3 stages for each sex--
#seed/zygote, juvenile/subadult, adult)
#S matrix: rate of survival for all stages of life
#R matrix: reproductive output at all stages of life
#M matrix: probability that male gametes come from a particular stage of life

#defining matrices for the landscape object
S <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.18, 
              0, 0, 0, 0, 0, 0, 0.14, 0, 0.26, 0, 0, 0, 0, 0.7, 0, 0.09, 0, 0,
              0, 0, 0.2, 0, 0.18), byrow = T, nrow = 6)
R <- matrix(c(0, 0, 0, 0, 14, 0, 0, 0, 0, 0, 8.5, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0), byrow = T, nrow = 6)
M <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0.75, 0, 0, 0,
              0, 0, 0), byrow = T, nrow = 6)
#assigning values to the landscape object
land = landscape.new.local.demo(land, S, R, M)
#add a new local demography with different reproduction 
R2 <- matrix(c(0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 5.5, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0), byrow = T, nrow = 6)
land = landscape.new.local.demo(land, S, R2, M)

#epochs
#defining empty matrix
zeromat <- matrix(0, nrow = 4 * 6, ncol = 4 * 6)
#defining yearly extinction rate in each habitat (4 habitats)
extnct <- c(0, 0.1, 0, 0.1)
#defining carrying capacity of each habitat (4 habitats)
k <- c(1000, 600, 600, 1000)
ldem <- c(0.5, 0.5)
#defining seed dispersal kernel matrix
sk <- matrix(0, nrow = 4 * 6, ncol = 6)
sk[, 1] <- rep(3, 4 * 6) #column is the seed dispersal kernel
sk[, 2] <- rep(10, 4 * 6) #this column is the scale component ? for kernel component 1
sk[, 3] <- rep(1.1, 4 * 6) #shape parameter for kernel component 1
sk[, 4] <- rep(100, 4 * 6)#scale component for kernel component 2
sk[, 5] <- rep(50, 4 * 6) #shape param for kernel component 2
sk[, 6] <- rep(0.5, 4 * 6) #mixing parameter--range 0 to 1. if 1, dispersal is determined by kernel 1
#if 2, dispersal is determined by kernel 2

#pollen dispersal matrix:
pk <- matrix(0, nrow = 4 * 6, ncol = 6)
pk[, 1] <- rep(3, 4 * 6) #pollen dispersal kernel
pk[, 2] <- rep(5, 4 * 6)
pk[, 3] <- rep(2, 4 * 6)
pk[, 4] <- rep(100, 4 * 6)
pk[, 5] <- rep(50, 4 * 6)
pk[, 6] <- rep(1, 4 * 6) #mixing parameter

#defining habitat locations
#locations are a set of vectors describing the left, right, top and bottom coordinates (N,S,E,W)
lx = c(0,0,800,800) #left x coordinate
rx = c(600,600,1400,1400) #right x coordinate
bty = c(0,800,0,800) #top y coordinate
ty = c(600,1400,600,1400) #bottom y coordinate

#add epoch elements to the landscape object
land = landscape.new.epoch(land, S = zeromat, R = zeromat, M = zeromat,
                           extinct = extnct, carry = k, localprob = ldem, pollen.kernels = pk,
                           seed.kernels = sk, leftx = lx, rightx = rx, boty = bty, topy = ty)

#Loci--describes locus characteristics including type, ploidy, transmission (biparental/maternal),
#mutation rate, etc...
#initializing 3 new loci below
land = landscape.new.locus(land, type = 0, ploidy = 1, transmission = 1, numalleles = 5) #haploid maternally inherited allele, infinite model, 5 alleles
land = landscape.new.locus(land, type = 1, ploidy = 2, transmission = 0, numalleles = 3) #diploid biparental inherited allele, stepwise mutation model, 3 alleles
land = landscape.new.locus(land, type = 1, ploidy = 2, transmission = 0, numalleles = 3) #diploid biparental allele, stepwise mutation model, 3 alleles
length(land$loci)

#Individuals are represented by a matrix with num rows = num individuals 
#num columns = demographic info and genetic info (1 column for haploid loci, 2 columns for diploid loci)
landscape.ploidy(land)
landscape.democol()

vlen <- ((land$intparam$habitats)*(land$intparam$stages)) #vector length is number habitat*number local stages
vec <- rep(100, vlen)
land <- landscape.new.individuals(land, vec)

#sIMULATION
#Simulating ecology and genetics!!
#commands below create 4 replicate simulations of 10 years
l1 = landscape.simulate(land, 10)
l2 = landscape.simulate(land, 10)
l3 = landscape.simulate(land, 10)
l4 = landscape.simulate(land, 10)
#plotting the landscapes
par(mfrow = c(2,2)) #plots 2x2 paneled plots 
landscape.plot.locations(l1)
landscape.plot.locations(l2)
landscape.plot.locations(l3)
landscape.plot.locations(l4)
#saving landscapes--they can be saved as R binary files (.rda)
save(file="l1.rda", l1)

#Alter landscapes throughout time
names(l1$demography$epochs[[1]])

sk <- matrix(0, nrow = 4 * 6, ncol = 6)
sk[, 1] <- rep(3, 4 * 6)
sk[, 2] <- rep(10, 4 * 6)
sk[, 3] <- rep(1.1, 4 * 6)
sk[, 4] <- rep(400, 4 * 6)
sk[, 5] <- rep(100, 4 * 6)
sk[, 6] <- rep(0.5, 4 * 6)
sk

#changing the seed kernel for land replicate 1
l1$demography$epochs[[1]]$seedkern = sk

#simulate another 10 years
l1 = landscape.simulate(l1, 10)
l2 = landscape.simulate(l2, 10)
l3 = landscape.simulate(l3, 10)
l4 = landscape.simulate(l4, 10)

#plotting locations of individuals in each habitat
par(mfrow = c(2, 2))
landscape.plot.locations(l1)
landscape.plot.locations(l2)
landscape.plot.locations(l3)
landscape.plot.locations(l4)
par(mfrow = c(1, 1))

#Extracting information from landscapes
#exact coordinates of every individual and their parents allows us to determine the dispersal distributions
#for zygotes (seeds) and male gametes (pollen)
source("C:\\Users\\kayle\\Documents\\kernelPop2\\R\\distance-functions.R")
#plotting every dispersal event that gives rise to an individual within a landscape 
par(mfrow = c(2, 2))
hist(seed.dist(l1), breaks = 30, xlab = "seed dispersal distance")
hist(seed.dist(l2), breaks = 30, xlab = "seed dispersal distance")
hist(seed.dist(l3), breaks = 30, xlab = "seed dispersal distance")
hist(seed.dist(l4), breaks = 30, xlab = "seed dispersal distance")
par(mfrow = c(1, 1))
#note the longer axis on l1, where there are more LDD events

#pollen dispersal events--pollen distance distributions are dependent on the spatial structure of the plants 
par(mfrow = c(2, 2))
hist(pollination.dist(l1), breaks = 30, xlab = "pollination dispersal distance")
hist(pollination.dist(l2), breaks = 30, xlab = "pollination dispersal distance")
hist(pollination.dist(l3), breaks = 30, xlab = "pollination dispersal distance")
hist(pollination.dist(l4), breaks = 30, xlab = "pollination dispersal distance")
par(mfrow = c(1, 1))
#note the longer axis on l1, where there are more LDD events for pollen

#POPULATIONS--can determine the number of individuals in each population, and which pop. an individual belongs to
table(landscape.populations(l1))

#GENETIC INFORMATION
#genetic information can be inferred from a landscape
#can get allele indices and states for each genotype at each locus and can get summary stats too (fst, allele freqs., etc...)
#landscape.sample() randomly samples the landscape for individuals to further analyze. 
#code below simulates landscape land for 100 generations, saving the state at every 10 gens in a list
#at gen 50, the mean dispersal distance is increased
gland <- land
landlist <- vector("list", 11)
landlist[[1]] <- gland
for (i in 2:11) {
  print(table(landscape.populations(gland)))
  gland <- landscape.simulate(gland, 10)
  landlist[[i]] <- gland
  if (i == 6) {
    sk <- gland$demography$epochs[[1]]$seedkern
    sk[, 4] <- rep(500, dim(sk)[1])
    sk[, 6] <- rep(0.8, dim(sk)[1])
    gland$demography$epochs[[1]]$seedkern <- sk
  }
}
#code takes list created above and samples 24 individuals from each pop. 
plot.ob <- do.call(rbind, lapply(landlist, function(l) {
  c(l$intparam$currentgen, mean(landscape.amova(l, ns = 24)))
  }))
print(xyplot(plot.ob[, 2] ~ plot.ob[, 1], 
             type = c("b", "smooth"), 
             xlab = "Time in years", 
             ylab = "Mean Phi-ST"))

#WRITING OUT FILES FOR OTHER PROGRAMS
source("C:\\Users\\kayle\\Documents\\kernelPop2\\R\\genepopR.R")
landscape.genepop.output(gland) #outputs data in genepop file format
#landscape.write.foreign() also does a lot of stuff
#like converting to arlequin, etc...

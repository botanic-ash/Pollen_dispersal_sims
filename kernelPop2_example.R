library(kernelPop2)

#create empty landscape
land = landscape.new.empty()

#defining integer parameters-- 4 habitats and 6 demographic stages
land = landscape.new.intparam(land, h=4, s=6)

#defining values for float params as default
land = landscape.new.floatparam(land)

#defining switch params as default
land = landscape.new.switchparam(land)

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

#epochs
#defining empty matrix
zeromat <- matrix(0, nrow = 4 * 6, ncol = 4 * 6)
#defining yearly extinction rate in each habitat (4 habitats)
extnct <- c(0, 0.1, 0, 0.1)
#defining carrying capacity of each habitat (4 habitats)
k <- c(1000,600,600,1000)
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

#set working directory
#using source to get this function
source("C:\\Users\\kayle\\Documents\\kernelPop2\\R\\simcoal.R")
setwd("C:\\Users\\kayle\\Documents\\Pollen_dispersal_sims_test\\Simcoal_files\\example_sp")
#input individuals from Simcoal simulations
land <- landscape.coalinput(land, npp=c(10,10,10,10), arlseq = NULL, arlms = "example_sp_0.arp", msmut = 0.0025)

#sIMULATION
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

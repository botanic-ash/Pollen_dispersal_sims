# Pollen_dispersal_sims
Project repository for complex pollen dispersal/mating using Simcoal 2
In collaboration by Kaylee Rosenberger, Sean Hoban, and Emily Schumacher 

#### Overview
Simulations (such as those done in the Simcoal software) have been used to test seed sampling strategies used to inform botanic gardens and arboreta. These informed sampling strategies are one way of ensuring that a genetically diverse and representative sample has been collected, to be conserved for future use in restoration efforts for example. However, the simulation software Simcoal models a simplified version of mating and pollen dispersal within the population. For example, any tree in the population has equal opportunity to donate pollen to any other tree in the population. In reality, for many species, the trees closest to another will donate the majority of the pollen. Furthermore, previous seed sampling strategies tested by simulations have assumed that one seed be sampled from one maternal tree--but in reality, collectors will sample many seeds from a given tree if they are available. While it has been shown that sampling as many unique maternal trees as possible will provide the most genetically diverse collection, this is often not feasible or realistic. We aim to quanitfy the difference in genetic diversity captured using a combination of different pollen dispersal patterns and sampling techniques. 

#### Directory contents
**Simulations:** contains simulation parameter files representing a hypothetical species 
        one_pop_2500: files represent a hypothetical species with one population of size 2500  
        two_pop_2500: files represent a hypothetical species with two populations each of size 2500  
**R-scripts:** contains R scripts used for data importing and processing, and sampling scripts
    functions.R: this script defines the functions used in the main loop script. There are multiple functions--some import and convert data into more usable file types (genalex). The main function in this loop (sample_seed) imports a genetic data and creates new seed sets. The number of pollen donors, number of seeds sampled per tree, and number of trees to sample from are taken as function parameters, making the inputs highly customizable.   
    parameters.R: this script makes lists that containing sets of function parameters to be passed to the sample_seed function. 
    main_loop.R: this script loops over simulation replicates and calls the functions defined in the functions.R script to run the functions with varying inputs, defined in parameters.R  
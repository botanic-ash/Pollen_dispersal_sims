# Modified Pollen_dispersal_sims
A project repository by Ash Hamilton for estimating error rate and parentage from offspring and parental genetic data

This repo is forked from Kaylee Rosenberger; it uses her sim2coal data and code written by her to generate offspring from that data. 


#### Overview
Simulations (such as those done in the Simcoal software) of MSAT data can be used to test sampling strategies to inform botanic gardens and arboreta how best to sample to preserve genetic diversity. One metric gaining popularity to measure diversity in both simulated and real data is percent capture of rare alleles. However, unlike real data, sim2coal data has no error. While MSATs are prone to a lower error rate per loci than SNPs, they are still subject to two different kinds of error: systematic (class 1) and stochastic (class 2). Class 1 systematic errors can be due to various kinds of “null” alleles including allelic dropout (where the shorter allele gets preferentially amplified), mutations to primer sequences (such that a mutation in the allele makes the allele impossible to amplify), or poor DNA template/primer quality (causing some alleles to consistently amplify worse than others). Class 2 Stochastic errors on the other hand can be due to slippage (where an allele randomly has too many or too few repeats due to an error by the TAQ polymerase during amplification) and human scoring error. Due to the different origins of these two error types, they can occur at very different rates, and these rates can vary between loci. Analyses involving the quantification of rare alleles are likely to be particularly heavily influenced by differences in error rates. As such, while point estimates for both MSAT error types can be obtained from likelihood maximization (see Wang 2004 and Wang & Santure 2009), providing confidence intervals on these error rates through a Bayesian framework will enable a better assessment of confidence in results and statistics that come from analyses of MSAT allele capture. In this project, we attempt to use a Gibbs sampler with Metropolis Hastings updates to estimate the posterior distribution of both class 1 and class 2 error rates and then use the mean estimated error rates to assign paternity via maximum likelihood. 

#### Approach
At a high level, the approach is as follows: (1) use a simulated populations from sim2coal to create a genetic data set representing a plant species, (2) apply a sampling function to this data set, such that we can sample differing numbers of maternal individuals and paternal individuals as well as number of offspring per maternal individual, (3) run that data through a Gibbs sampler algorithm which uses standard uniform priors on each error rate and piece wise Metropolis Hastings updates to each error type (w/ proposals being generated from a random walk w/ sd of .05) to generate the posterior distributions of each error type, (4) the posterior distributions are visualized along with trace plots of the values of error to evaluate convergence, (5) paternal individual identity is assigned via maximum likelihood using the means of the posterior distributions of each error type.

#### Function parameters and assumptions
The 

#### Directory contents
**R-scripts:** contains the R script written for this project
**Figures:** contains figures and plots created in R with baseR
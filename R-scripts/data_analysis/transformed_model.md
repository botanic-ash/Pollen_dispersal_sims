
``` r
#Load libraries required for the whole script
library(rstanarm)
```

    ## Loading required package: Rcpp

    ## This is rstanarm version 2.21.3

    ## - See https://mc-stan.org/rstanarm/articles/priors for changes to default priors!

    ## - Default priors may change, so it's safest to specify priors, even if equivalent to the defaults.

    ## - For execution on a local, multicore CPU with excess RAM we recommend calling

    ##   options(mc.cores = parallel::detectCores())

``` r
options(mc.cores = parallel::detectCores())
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(ggplot2)
theme_set(theme_bw())

source("hpdi.R")

#Load in data 
load("tidy_df.Rdata")
```

``` r
#Converting these to numeric--they already should be, but they must have been converted accidentally when making the matrix a dataframe, etc...
tidy_df$prop_capt = as.numeric(tidy_df$prop_capt)
tidy_df$total_seeds = as.numeric(tidy_df$total_seeds)
tidy_df$maternal_trees = as.numeric(tidy_df$maternal_trees)
#Running the model
transformed_model = glm(prop_capt ~ log(total_seeds) * maternal_trees * donor_type,
                          weights=total_seeds, family = binomial(link='logit'), data = tidy_df)
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
# #Save the model since it takes so long to run
save(transformed_model, file = "transformed_model.Rdata")
#Load the model from previously saved run
load("transformed_model.Rdata")
#Model summary! 
summary(transformed_model, digits = 4)
```

    ## 
    ## Call:
    ## glm(formula = prop_capt ~ log(total_seeds) * maternal_trees * 
    ##     donor_type, family = binomial(link = "logit"), data = tidy_df, 
    ##     weights = total_seeds)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -28.9175   -1.4151   -0.1883    1.0569   11.1638  
    ## 
    ## Coefficients:
    ##                                                      Estimate Std. Error
    ## (Intercept)                                        -1.5952723  0.0097853
    ## log(total_seeds)                                    0.6690339  0.0017414
    ## maternal_trees                                      0.0319706  0.0012266
    ## donor_typeall_same                                  0.5632364  0.0134966
    ## donor_typeskewed                                    0.9318726  0.0129512
    ## log(total_seeds):maternal_trees                    -0.0048363  0.0002135
    ## log(total_seeds):donor_typeall_same                -0.6489173  0.0023768
    ## log(total_seeds):donor_typeskewed                  -0.4647758  0.0022867
    ## maternal_trees:donor_typeall_same                   0.1353492  0.0026443
    ## maternal_trees:donor_typeskewed                     0.0233368  0.0023228
    ## log(total_seeds):maternal_trees:donor_typeall_same -0.0019413  0.0004588
    ## log(total_seeds):maternal_trees:donor_typeskewed    0.0070384  0.0004057
    ##                                                     z value Pr(>|z|)    
    ## (Intercept)                                        -163.027  < 2e-16 ***
    ## log(total_seeds)                                    384.193  < 2e-16 ***
    ## maternal_trees                                       26.065  < 2e-16 ***
    ## donor_typeall_same                                   41.732  < 2e-16 ***
    ## donor_typeskewed                                     71.953  < 2e-16 ***
    ## log(total_seeds):maternal_trees                     -22.653  < 2e-16 ***
    ## log(total_seeds):donor_typeall_same                -273.024  < 2e-16 ***
    ## log(total_seeds):donor_typeskewed                  -203.253  < 2e-16 ***
    ## maternal_trees:donor_typeall_same                    51.184  < 2e-16 ***
    ## maternal_trees:donor_typeskewed                      10.047  < 2e-16 ***
    ## log(total_seeds):maternal_trees:donor_typeall_same   -4.231 2.33e-05 ***
    ## log(total_seeds):maternal_trees:donor_typeskewed     17.348  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 10031191  on 140249  degrees of freedom
    ## Residual deviance:   872445  on 140238  degrees of freedom
    ## AIC: 1612659
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
#Creating a new dataframe of values to base predictions on 
newd = data.frame(maternal_trees=(rep(c(1,2,5,10,25,50,100), each=1500)), total_seeds=rep(seq(1,500,1),21), donor_type=factor(rep((rep(c("all_eligible", "all_same", "skewed"), each=500)), 7)))
#Predictions 
pmu = predict(transformed_model, re.form=NA, transform = TRUE, newdata=newd)

#Creating a dataframe to plot in ggplot 
preds <- cbind(newd, pmu)

#Plotting the data
ggplot(data=preds) +
    geom_point(data = tidy_df, aes(x=as.numeric(total_seeds), y=as.numeric(prop_capt), color=donor_type), alpha=0.25) +
    facet_wrap(vars(maternal_trees)) +
    geom_line(mapping = aes(x=total_seeds, y=pmu, lty=donor_type)) +
    ylim(0,1) +
    ylab("Proportion of alleles captured") +
    xlab("Total seeds sampled")
```

    ## Warning: Removed 5517 row(s) containing missing values (geom_path).

![](transformed_model_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

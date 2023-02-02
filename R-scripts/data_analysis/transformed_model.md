
``` r
#Load libraries required for the whole script
library(rstanarm)
```

    ## Warning: package 'rstanarm' was built under R version 4.0.5

    ## Loading required package: Rcpp

    ## Warning: replacing previous import 'lifecycle::last_warnings' by
    ## 'rlang::last_warnings' when loading 'tibble'

    ## Warning: replacing previous import 'lifecycle::last_warnings' by
    ## 'rlang::last_warnings' when loading 'pillar'

    ## This is rstanarm version 2.21.3

    ## - See https://mc-stan.org/rstanarm/articles/priors for changes to default priors!

    ## - Default priors may change, so it's safest to specify priors, even if equivalent to the defaults.

    ## - For execution on a local, multicore CPU with excess RAM we recommend calling

    ##   options(mc.cores = parallel::detectCores())

``` r
options(mc.cores = parallel::detectCores())
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.0.5

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
```

    ## Warning: package 'tidyr' was built under R version 4.0.5

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.0.5

``` r
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
                          weights=total_seeds, data = tidy_df)
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
    ##     donor_type, data = tidy_df, weights = total_seeds)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -11.6722   -0.7678   -0.0784    0.4258    8.1515  
    ## 
    ## Coefficients:
    ##                                                      Estimate Std. Error
    ## (Intercept)                                         4.587e-01  5.168e-03
    ## log(total_seeds)                                    7.697e-02  8.998e-04
    ## maternal_trees                                      4.689e-03  5.205e-04
    ## donor_typeall_same                                 -1.619e-01  7.309e-03
    ## donor_typeskewed                                   -1.047e-01  7.309e-03
    ## log(total_seeds):maternal_trees                    -7.531e-04  8.954e-05
    ## log(total_seeds):donor_typeall_same                -7.145e-02  1.273e-03
    ## log(total_seeds):donor_typeskewed                  -2.583e-02  1.273e-03
    ## maternal_trees:donor_typeall_same                   1.440e-02  7.361e-04
    ## maternal_trees:donor_typeskewed                     7.523e-03  7.361e-04
    ## log(total_seeds):maternal_trees:donor_typeall_same -5.627e-04  1.266e-04
    ## log(total_seeds):maternal_trees:donor_typeskewed   -3.944e-04  1.266e-04
    ##                                                    t value Pr(>|t|)    
    ## (Intercept)                                         88.750  < 2e-16 ***
    ## log(total_seeds)                                    85.539  < 2e-16 ***
    ## maternal_trees                                       9.008  < 2e-16 ***
    ## donor_typeall_same                                 -22.150  < 2e-16 ***
    ## donor_typeskewed                                   -14.327  < 2e-16 ***
    ## log(total_seeds):maternal_trees                     -8.411  < 2e-16 ***
    ## log(total_seeds):donor_typeall_same                -56.145  < 2e-16 ***
    ## log(total_seeds):donor_typeskewed                  -20.299  < 2e-16 ***
    ## maternal_trees:donor_typeall_same                   19.560  < 2e-16 ***
    ## maternal_trees:donor_typeskewed                     10.220  < 2e-16 ***
    ## log(total_seeds):maternal_trees:donor_typeall_same  -4.443 8.87e-06 ***
    ## log(total_seeds):maternal_trees:donor_typeskewed    -3.115  0.00184 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 2.047454)
    ## 
    ##     Null deviance: 2145942  on 140249  degrees of freedom
    ## Residual deviance:  287131  on 140238  degrees of freedom
    ## AIC: -235827
    ## 
    ## Number of Fisher Scoring iterations: 2

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

    ## Warning: Removed 1000 row(s) containing missing values (geom_path).

![](transformed_model_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

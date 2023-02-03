
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
transformed_model = glm(prop_capt ~ log(total_seeds) * log(maternal_trees) * donor_type, data = tidy_df) #removed weights=total_seeds
# #Save the model since it takes so long to run
save(transformed_model, file = "transformed_model.Rdata")
#Load the model from previously saved run
load("transformed_model.Rdata")
#Model summary! 
summary(transformed_model, digits = 4)
```

    ## 
    ## Call:
    ## glm(formula = prop_capt ~ log(total_seeds) * log(maternal_trees) * 
    ##     donor_type, data = tidy_df)
    ## 
    ## Deviance Residuals: 
    ##       Min         1Q     Median         3Q        Max  
    ## -0.265841  -0.018967  -0.000248   0.021007   0.150021  
    ## 
    ## Coefficients:
    ##                                                           Estimate Std. Error
    ## (Intercept)                                              0.2616000  0.0010836
    ## log(total_seeds)                                         0.1126096  0.0002035
    ## log(maternal_trees)                                      0.0926990  0.0011406
    ## donor_typeall_same                                      -0.0251088  0.0015324
    ## donor_typeskewed                                        -0.0264045  0.0015324
    ## log(total_seeds):log(maternal_trees)                    -0.0159334  0.0002116
    ## log(total_seeds):donor_typeall_same                     -0.1091952  0.0002879
    ## log(total_seeds):donor_typeskewed                       -0.0470822  0.0002879
    ## log(maternal_trees):donor_typeall_same                   0.0905537  0.0016130
    ## log(maternal_trees):donor_typeskewed                     0.0489686  0.0016130
    ## log(total_seeds):log(maternal_trees):donor_typeall_same  0.0168630  0.0002993
    ## log(total_seeds):log(maternal_trees):donor_typeskewed    0.0083232  0.0002993
    ##                                                         t value Pr(>|t|)    
    ## (Intercept)                                              241.42   <2e-16 ***
    ## log(total_seeds)                                         553.24   <2e-16 ***
    ## log(maternal_trees)                                       81.27   <2e-16 ***
    ## donor_typeall_same                                       -16.39   <2e-16 ***
    ## donor_typeskewed                                         -17.23   <2e-16 ***
    ## log(total_seeds):log(maternal_trees)                     -75.29   <2e-16 ***
    ## log(total_seeds):donor_typeall_same                     -379.33   <2e-16 ***
    ## log(total_seeds):donor_typeskewed                       -163.56   <2e-16 ***
    ## log(maternal_trees):donor_typeall_same                    56.14   <2e-16 ***
    ## log(maternal_trees):donor_typeskewed                      30.36   <2e-16 ***
    ## log(total_seeds):log(maternal_trees):donor_typeall_same   56.34   <2e-16 ***
    ## log(total_seeds):log(maternal_trees):donor_typeskewed     27.81   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.001212871)
    ## 
    ##     Null deviance: 8097.39  on 140249  degrees of freedom
    ## Residual deviance:  170.09  on 140238  degrees of freedom
    ## AIC: -543720
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

    ## Warning: Removed 961 row(s) containing missing values (geom_path).

![](transformed_model_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

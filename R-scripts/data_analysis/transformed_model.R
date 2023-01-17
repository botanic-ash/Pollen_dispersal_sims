---
output: github_document
---

```{r}
#Load libraries required for the whole script
library(rstanarm)
options(mc.cores = parallel::detectCores())
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())

source("hpdi.R")

#Load in data 
load("tidy_df.Rdata")
```
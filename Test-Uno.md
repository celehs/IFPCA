Intensity Functional Principal Component Analysis
================

## Overview

This package re-implements the first step of the
{[MASTA](https://celehs.github.io/MASTA/)} package to extract features
from longitudinal encounter records. Compared to {MASTA}, the input data
of {IFPCA} is more compact and memory efficent. Click
[HERE](https://github.com/celehs/IFPCA/tree/master/data-raw) to view
input data structure.

## Installation

Install development version from GitHub.

``` r
## install.packages("remotes")
remotes::install_github("celehs/IFPCA", force=TRUE)
```

Load the package into R.

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:HH':
    ## 
    ##     transpose

``` r
library(IFPCA)
```

## Data Example

### Data Preparation

``` r
url <- "https://raw.githubusercontent.com/celehs/IFPCA/master/data-raw/"
follow_up_train <- fread(paste0(url, "follow_up_train.csv"))[1:600]
follow_up_valid <- fread(paste0(url, "follow_up_valid.csv"))

fu_train <- follow_up_train$fu_time
fu_valid <- follow_up_valid$fu_time
names(fu_train) <- follow_up_train$id
names(fu_valid) <- follow_up_valid$id

#--- the number of longitudinal codes
number_of_codes = 3

D=list()
for (i in 1:number_of_codes){
time_code <- fread(paste0(url, "time_code",i,".csv"))
time <- time_code$month
names(time) <- time_code$id
D[[i]]=time
}
```

### Run the 1st step (feature selection)

``` r
TrainFt = c()
ValidFt = c()
TrainPK = c()
ValidPK = c()

for (i in 1:length(D)){
 time=D[[i]]
 ans <- ifpca(time, fu_train, fu_valid)

#--- create an object for fitting --
Ft_name = colnames(ans$TrainFt) ; Ft_name = paste0(Ft_name,i) ; Ft_name ; colnames(ans$TrainFt) = Ft_name
TrainFt = cbind(TrainFt,ans$TrainFt)
Ft_name = colnames(ans$ValidFt) ; Ft_name = paste0(Ft_name,i) ; Ft_name ; colnames(ans$ValidFt) = Ft_name
ValidFt = cbind(ValidFt,ans$ValidFt)
TrainPK = cbind(TrainPK, ans$TrainPK) 
ValidPK = cbind(ValidPK, ans$ValidPK) 
}

colnames(TrainPK)=colnames(ValidPK) = paste0("pred",1:length(D)) 

Z=list()
Z$TrainFt = TrainFt
Z$ValidFt = ValidFt
Z$TrainPK = TrainPK
Z$ValidPK = ValidPK
```

### Data preparation for the 2nd step (fitting)

``` r
url <- "https://raw.githubusercontent.com/celehs/MASTA/master/data-raw/data_org/"
TrainSurv <- fread(paste0(url, "TrainSurv.csv"))
ValidSurv <- fread(paste0(url, "ValidSurv.csv"))
#load("data_org.rda")
TrainSurv <- data.frame(TrainSurv)
ValidSurv <- data.frame(ValidSurv)
colnames(TrainSurv)=colnames(ValidSurv) = c("case","delta","sx","sc",paste0("base_pred",1:3))

Z$nn = length(fu_train)
Z$codes = paste0("pred",1:3) 
Z$Tend <- 1
Z$TrainSurv <- TrainSurv
Z$ValidSurv <- ValidSurv 
Z$TrainSurv_pred_org <- TrainSurv[,-c(1:4)]
Z$ValidSurv_pred_org <- ValidSurv[,-c(1:4)]
str(Z)
```

    ## List of 11
    ##  $ TrainFt           : num [1:600, 1:15] 49.4 13.9 12.6 1 4 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:600] "1" "2" "3" "4" ...
    ##   .. ..$ : chr [1:15] "1stCode1" "Pk1" "ChP1" "1stScore1" ...
    ##  $ ValidFt           : num [1:500, 1:15] 71.7 70.4 14.9 33.8 98.6 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:500] "90001" "90002" "90003" "90004" ...
    ##   .. ..$ : chr [1:15] "1stCode1" "Pk1" "ChP1" "1stScore1" ...
    ##  $ TrainPK           : int [1:600, 1:3] NA NA NA 34 4 14 NA NA 11 NA ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:600] "1" "2" "3" "4" ...
    ##   .. ..$ : chr [1:3] "pred1" "pred2" "pred3"
    ##  $ ValidPK           : int [1:500, 1:3] 9 NA 5 NA NA 1 NA NA 1 2 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:500] "90001" "90002" "90003" "90004" ...
    ##   .. ..$ : chr [1:3] "pred1" "pred2" "pred3"
    ##  $ nn                : int 600
    ##  $ codes             : chr [1:3] "pred1" "pred2" "pred3"
    ##  $ Tend              : num 1
    ##  $ TrainSurv         :'data.frame':  600 obs. of  7 variables:
    ##   ..$ case      : int [1:600] 1 2 3 4 5 6 7 8 9 10 ...
    ##   ..$ delta     : int [1:600] 1 0 0 0 0 1 0 0 0 0 ...
    ##   ..$ sx        : num [1:600] 9.36 13.93 12.55 14.85 80.66 ...
    ##   ..$ sc        : num [1:600] 49.4 13.9 12.6 14.9 80.7 ...
    ##   ..$ base_pred1: int [1:600] 79 81 55 72 83 47 86 40 85 58 ...
    ##   ..$ base_pred2: int [1:600] 1 0 1 1 1 1 0 0 1 0 ...
    ##   ..$ base_pred3: int [1:600] 0 0 1 0 1 0 1 0 0 1 ...
    ##  $ ValidSurv         :'data.frame':  500 obs. of  7 variables:
    ##   ..$ case      : int [1:500] 90001 90002 90003 90004 90005 90006 90007 90008 90009 90010 ...
    ##   ..$ delta     : int [1:500] 0 0 1 1 0 0 0 0 0 0 ...
    ##   ..$ sx        : num [1:500] 71.72 70.37 9.49 17.77 98.63 ...
    ##   ..$ sc        : num [1:500] 71.7 70.4 14.9 33.8 98.6 ...
    ##   ..$ base_pred1: int [1:500] 39 36 68 59 38 53 88 60 40 33 ...
    ##   ..$ base_pred2: int [1:500] 1 0 1 0 1 0 0 1 0 0 ...
    ##   ..$ base_pred3: int [1:500] 0 1 1 0 0 0 1 0 0 0 ...
    ##  $ TrainSurv_pred_org:'data.frame':  600 obs. of  3 variables:
    ##   ..$ base_pred1: int [1:600] 79 81 55 72 83 47 86 40 85 58 ...
    ##   ..$ base_pred2: int [1:600] 1 0 1 1 1 1 0 0 1 0 ...
    ##   ..$ base_pred3: int [1:600] 0 0 1 0 1 0 1 0 0 1 ...
    ##  $ ValidSurv_pred_org:'data.frame':  500 obs. of  3 variables:
    ##   ..$ base_pred1: int [1:500] 39 36 68 59 38 53 88 60 40 33 ...
    ##   ..$ base_pred2: int [1:500] 1 0 1 0 1 0 0 1 0 0 ...
    ##   ..$ base_pred3: int [1:500] 0 1 1 0 0 0 1 0 0 0 ...

### Run the 2nd step

``` r
library(survival)
library(doParallel)
```

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

``` r
library(foreach)
library(data.table)
library(survC1)
library(rootSolve)
library(splines)
library(gglasso)
```

    ## 
    ## Attaching package: 'gglasso'

    ## The following object is masked from 'package:survival':
    ## 
    ##     colon

``` r
library(glmnet)
```

    ## Loaded glmnet 4.0

``` r
library(rpart)
library(rpart.utils)

url <- "https://raw.githubusercontent.com/celehs/MASTA/master/R/"
source(paste0(url,"masta-fit.R"))
source(paste0(url,"masta-fit-bs-grplasso.R"))
source(paste0(url,"masta-fit-npmle.R"))
source(paste0(url,"masta-fit-cheng.R"))
```

``` r
b=masta.fit(Z, cov_group = NULL, thresh = 0.7, PCAthresh = 0.9, seed = 1234, seed2 = 100) 
b
```

    ## $bgbbest_FromChengInit_BFGS
    ##         bgbm.init                       AIC          BIC     AIC.Orig
    ##  [1,] -9.76795676 -10.22834075 -10.24051174 -10.24051174 -10.24051174
    ##  [2,]  4.22627139   4.44098240   4.47773491   4.47773491   4.47773491
    ##  [3,]  4.47949568   5.22358547   5.27465675   5.27465675   5.27465675
    ##  [4,]  5.29371585   5.77434758   5.81392289   5.81392289   5.81392289
    ##  [5,]  4.79469823   5.49945786   5.52435356   5.52435356   5.52435356
    ##  [6,]  6.01937444   6.82891319   6.84167003   6.84167003   6.84167003
    ##  [7,]  5.67552751   6.67219287   6.68302226   6.68302226   6.68302226
    ##  [8,]  5.42277173   6.77315502   6.78634585   6.78634585   6.78634585
    ##  [9,]  6.00757797   7.70758141   7.73266840   7.73266840   7.73266840
    ## [10,]  5.42753433   7.58483121   7.59330594   7.59330594   7.59330594
    ## [11,]  4.02123441   8.26059042   8.30068456   8.30068456   8.30068456
    ## [12,]  0.14286894   0.07144419   0.07707134   0.07707134   0.07707134
    ## [13,] -0.32068091  -0.46241680  -0.48019900  -0.48019900  -0.48019900
    ## [14,] -0.31976383  -0.52124237  -0.49932617  -0.49932617  -0.49932617
    ## [15,]  0.06620263   0.25202982   0.24114425   0.24114425   0.24114425
    ## [16,] -0.38520792  -0.41814454  -0.40425237  -0.40425237  -0.40425237
    ## [17,] -1.68993194  -0.97618622  -0.97818705  -0.97818705  -0.97818705
    ## [18,] -0.56471194  -0.59541043  -0.61291899  -0.61291899  -0.61291899
    ## [19,] -0.01810521  -0.04193981   0.00000000   0.00000000   0.00000000
    ## [20,]  0.30652211   0.03308233   0.00000000   0.00000000   0.00000000
    ## [21,]  0.11137667   0.22156687   0.00000000   0.00000000   0.00000000
    ## [22,] -0.12898221  -0.01838244   0.00000000   0.00000000   0.00000000
    ## [23,]  1.31048145   1.80983898   1.81451705   1.81451705   1.81451705
    ## [24,]  0.57541565   0.07778873   0.07300103   0.07300103   0.07300103
    ## [25,]  0.01036809   0.10816522   0.11268440   0.11268440   0.11268440
    ##           BIC.Orig
    ##  [1,] -10.24051174
    ##  [2,]   4.47773491
    ##  [3,]   5.27465675
    ##  [4,]   5.81392289
    ##  [5,]   5.52435356
    ##  [6,]   6.84167003
    ##  [7,]   6.68302226
    ##  [8,]   6.78634585
    ##  [9,]   7.73266840
    ## [10,]   7.59330594
    ## [11,]   8.30068456
    ## [12,]   0.07707134
    ## [13,]  -0.48019900
    ## [14,]  -0.49932617
    ## [15,]   0.24114425
    ## [16,]  -0.40425237
    ## [17,]  -0.97818705
    ## [18,]  -0.61291899
    ## [19,]   0.00000000
    ## [20,]   0.00000000
    ## [21,]   0.00000000
    ## [22,]   0.00000000
    ## [23,]   1.81451705
    ## [24,]   0.07300103
    ## [25,]   0.11268440
    ## 
    ## $Cstat_BrierSc_ChengInit_BFGS
    ##                  Init       MLE       AIC       BIC  AIC.Orig  BIC.Orig
    ## Cstat       0.8215967 0.8663195 0.8641285 0.8641285 0.8641285 0.8641285
    ## BrierSc.Adj 0.1789903 0.4677978 0.4674519 0.4674519 0.4674519 0.4674519
    ##                 Cheng     NPMLE    TreeTN   TreeAll         logi
    ## Cstat       0.8215967 0.8630978 0.8613730 0.8196111  0.816663184
    ## BrierSc.Adj 0.1795059 0.4355211 0.2329827 0.2518399 -0.005523029
    ## 
    ## $group
    ##  [1] 1 2 3 4 4 4 4 5 5 5 5 6 6 6

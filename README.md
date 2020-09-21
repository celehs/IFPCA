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
# install.packages("remotes")
remotes::install_github("celehs/IFPCA")
```

Load the package into R.

``` r
library(IFPCA)
```

## Data Example

### Data Preparation

``` r
library(data.table)
```

``` r
i <- 3 # 1, 2, 3 
# time for each code (training + validation)
url <- "https://raw.githubusercontent.com/celehs/IFPCA/master/data-raw/"
time_code <- fread(paste0(url, "time_code", i, ".csv"))
time <- time_code$month
names(time) <- time_code$id
# follow up time for training and validation sets
follow_up_train <- fread(paste0(url, "follow_up_train.csv"))
follow_up_valid <- fread(paste0(url, "follow_up_valid.csv"))
fu_train <- follow_up_train$fu_time
fu_valid <- follow_up_valid$fu_time
names(fu_train) <- follow_up_train$id
names(fu_valid) <- follow_up_valid$id
str(fu_train)
```

    ##  Named num [1:20600] 49.4 13.9 12.6 14.9 80.7 ...
    ##  - attr(*, "names")= chr [1:20600] "1" "2" "3" "4" ...

``` r
str(fu_valid)
```

    ##  Named num [1:500] 71.7 70.4 14.9 33.8 98.6 ...
    ##  - attr(*, "names")= chr [1:500] "90001" "90002" "90003" "90004" ...

### Feature Extraction

``` r
system.time(ans <- ifpca(time, fu_train, fu_valid))
```

    ##    user  system elapsed 
    ##   8.275   0.332   8.614

``` r
data.table(ans$TrainFt) # Extracted Features (Training) 
```

    ##         1stCode        Pk       ChP    1stScore      logN
    ##     1: 49.41273 49.412731 49.412731 -1.77557902 0.0000000
    ##     2: 13.93018 13.930185 13.930185 -1.77557902 0.0000000
    ##     3: 12.55031 12.550308 12.550308 -1.77557902 0.0000000
    ##     4: 14.85010 14.850103 14.850103 -1.77557902 0.0000000
    ##     5: 80.65708 80.657084 80.657084 -1.77557902 0.0000000
    ##    ---                                                   
    ## 20596:  2.00000  3.750801  3.595729 -0.05899512 1.3862944
    ## 20597: 14.00000 14.359918 14.140123  0.71721892 0.6931472
    ## 20598:  8.00000  7.638604  7.434908 -1.32323471 1.0986123
    ## 20599: 13.00000 13.059548 12.711294 -1.36421318 1.0986123
    ## 20600: 23.00000 24.274004 23.841643  3.02440941 1.3862944

``` r
data.table(ans$ValidFt) # Extracted Features (Validation)
```

    ##       1stCode       Pk      ChP  1stScore      logN
    ##   1: 71.72074 71.72074 71.72074 -1.775579 0.0000000
    ##   2: 70.37372 70.37372 70.37372 -1.775579 0.0000000
    ##   3: 14.94867 14.01437 13.64066  3.038376 0.6931472
    ##   4: 33.80698 33.80698 33.80698 -1.775579 0.0000000
    ##   5: 98.62834 98.62834 98.62834 -1.775579 0.0000000
    ##  ---                                               
    ## 496: 35.58111 35.58111 35.58111 -1.775579 0.0000000
    ## 497: 27.36756 27.36756 27.36756 -1.775579 0.0000000
    ## 498: 24.60780 24.60780 24.60780 -1.775579 0.0000000
    ## 499: 27.17043 27.17043 27.17043 -1.775579 0.0000000
    ## 500: 29.30595 29.30595 29.30595 -1.775579 0.0000000

## References

  - Wu, S., Müller, H., Zhang, Z. (2013). **Functional Data Analysis for
    Point Processes with Rare Events**. *Statistica Sinica*, 23:1-23.
    <https://doi.org/10.5705/ss.2010.162>

  - Liang, L., Uno, H., Ma, Y., Cai, T. **Robust Approach to Event Time
    Annotation Using Longitudinal Medical Encounters**. *Working Paper*.

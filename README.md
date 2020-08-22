IFPCA
================

``` r
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.3     ✓ dplyr   1.0.1
    ## ✓ tidyr   1.1.1     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ──────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
url <- "https://raw.githubusercontent.com/celehs/MASTA/master/data-raw/long/"
follow_up_train <- read_csv(paste0(url, "follow_up_train.csv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   id = col_double(),
    ##   fu_time = col_double()
    ## )

``` r
follow_up_valid <- read_csv(paste0(url, "follow_up_valid.csv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   id = col_double(),
    ##   fu_time = col_double()
    ## )

``` r
time_code1 <- read_csv(paste0(url, "time_code1.csv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   id = col_double(),
    ##   month = col_double()
    ## )

``` r
time_code2 <- read_csv(paste0(url, "time_code2.csv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   id = col_double(),
    ##   month = col_double()
    ## )

``` r
time_code3 <- read_csv(paste0(url, "time_code3.csv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   id = col_double(),
    ##   month = col_double()
    ## )

``` r
fu_train <- follow_up_train$fu_time
fu_valid <- follow_up_valid$fu_time
names(fu_train) <- follow_up_train$id
names(fu_valid) <- follow_up_valid$id
```

``` r
time1 <- time_code1$month
time2 <- time_code2$month
time3 <- time_code3$month
names(time1) <- time_code1$id
names(time2) <- time_code2$id
names(time3) <- time_code3$id
time <- list(time1, time2, time3)
```

``` r
ifpca <- function(time, fu_train, fu_valid) {
  input <- list(
    time = time, 
    fu_train = fu_train, 
    fu_valid = fu_valid)
  return(str(input))
}
```

``` r
ifpca(time, fu_train, fu_valid)
```

    ## List of 3
    ##  $ time    :List of 3
    ##   ..$ : Named num [1:168374] 1 1 1 2 4 4 4 5 5 5 ...
    ##   .. ..- attr(*, "names")= chr [1:168374] "4" "4" "4" "4" ...
    ##   ..$ : Named num [1:68242] 11 1 2 4 2 11 22 22 25 29 ...
    ##   .. ..- attr(*, "names")= chr [1:68242] "1" "5" "5" "5" ...
    ##   ..$ : Named num [1:21501] 37 37 38 38 39 39 40 40 41 41 ...
    ##   .. ..- attr(*, "names")= chr [1:21501] "6" "6" "6" "6" ...
    ##  $ fu_train: Named num [1:20600] 49.4 13.9 12.6 14.9 80.7 ...
    ##   ..- attr(*, "names")= chr [1:20600] "1" "2" "3" "4" ...
    ##  $ fu_valid: Named num [1:500] 71.7 70.4 14.9 33.8 98.6 ...
    ##   ..- attr(*, "names")= chr [1:500] "90001" "90002" "90003" "90004" ...

``` r
proc.time()
```

    ##    user  system elapsed 
    ##   1.621   0.128   1.779

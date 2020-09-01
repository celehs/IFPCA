Wide to Long Format
================

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.3     ✓ dplyr   1.0.1
    ## ✓ tidyr   1.1.1     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
Code <- bind_rows(
  Train = read_csv("data_org/TrainCode.csv"), 
  Valid = read_csv("data_org/ValidCode.csv"), .id = "split")
```

    ## Parsed with column specification:
    ## cols(
    ##   id = col_double(),
    ##   fu_time = col_double(),
    ##   month = col_double(),
    ##   pred1 = col_double(),
    ##   pred2 = col_double(),
    ##   pred3 = col_double()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   id = col_double(),
    ##   fu_time = col_double(),
    ##   month = col_double(),
    ##   pred1 = col_double(),
    ##   pred2 = col_double(),
    ##   pred3 = col_double()
    ## )

``` r
fu_time <- Code %>% select(split, id, fu_time) %>% distinct()
follow_up_train <- fu_time %>% filter(split == "Train") %>% select(-split)
follow_up_valid <- fu_time %>% filter(split == "Valid") %>% select(-split)
write_csv(follow_up_train, "long/follow_up_train.csv")
write_csv(follow_up_valid, "long/follow_up_valid.csv")
```

``` r
for (i in 1:3) {
  pred <- paste0("pred", i)
  path <- paste0("long/time_code", i, ".csv")
  wide <- Code[c(Code[, pred] > 0), c("id", "month", pred)]
  freq <- wide[, pred][[1]]
  long <- tibble(
    id = rep(wide$id, freq),
    month = rep(wide$month, freq))
  write_csv(long, path)  
  print(path)  
}
```

    ## [1] "long/time_code1.csv"
    ## [1] "long/time_code2.csv"
    ## [1] "long/time_code3.csv"

``` r
proc.time()
```

    ##    user  system elapsed 
    ##   2.503   0.155   2.645

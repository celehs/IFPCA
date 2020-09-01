Intensity Functional Principal Component Analysis (IFPCA)
================

``` r
library(tidyverse)
library(doParallel)
source("R/GetPK.R")
source("R/masta-fpca-main.R")
source("R/masta-fpca-kern.R")
url <- "data-raw/long/"
follow_up_train <- read_csv(paste0(url, "follow_up_train.csv"))
follow_up_valid <- read_csv(paste0(url, "follow_up_valid.csv"))
time_code1 <- read_csv(paste0(url, "time_code1.csv"))
time_code2 <- read_csv(paste0(url, "time_code2.csv"))
time_code3 <- read_csv(paste0(url, "time_code3.csv"))
```

``` r
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

``` r
time1 <- time_code1$month
time2 <- time_code2$month
time3 <- time_code3$month
names(time1) <- time_code1$id
names(time2) <- time_code2$id
names(time3) <- time_code3$id
str(time1)
```

    ##  Named num [1:168374] 1 1 1 2 4 4 4 5 5 5 ...
    ##  - attr(*, "names")= chr [1:168374] "4" "4" "4" "4" ...

``` r
str(time2)
```

    ##  Named num [1:68242] 11 1 2 4 2 11 22 22 25 29 ...
    ##  - attr(*, "names")= chr [1:68242] "1" "5" "5" "5" ...

``` r
str(time3)
```

    ##  Named num [1:21501] 37 37 38 38 39 39 40 40 41 41 ...
    ##  - attr(*, "names")= chr [1:21501] "6" "6" "6" "6" ...

``` r
#' @title Intensity Functional Principal Component Analysis (IFPCA)
#' @description Performs IFPCA to extract features from longitudinal encounter data.
#' @param data input data. See \code{data(data_org)} for example.
#' @param PPIC_K a logical indicating whether you want to use Pseudo-Poisson Information Criterion to choose 
#' the number of principal components K (K.select="PPIC") \code{TRUE} or another criterion to choose 
#' K (K.select="PropVar") \code{FALSE} in the PP_FPCA_Parallel function (hidden). Default is \code{FALSE}.
#' @param n.grid an integer value for grid points used in estimating covariance function g. Default is \code{401}.
#' @param propvar a proportion of variation used to select number of FPCs. Default is \code{0.85}.
#' @param n_core an integer to specify the number of core using for parallel computing. 
#' @export
ifpca <- function(time, fu_train, fu_valid, 
                  PPIC_K = FALSE, n.grid = 401, propvar = 0.85, n_core = NULL) {
  if (is.null(n_core)) n_core <- parallel::detectCores()
  registerDoParallel(cores = n_core)    
  tseq <- 0:floor(max(c(fu_train, fu_valid)))
  PKTS <- GetPK(id = names(time), t = time, tseq = tseq, nn = length(fu_train))
  return(PKTS[!is.na(PKTS)])
}
```

``` r
pkts <- ifpca(time1, fu_train, fu_valid)
str(pkts)
```

    ##  Named int [1:12523] 1 27 95 1 70 3 1 4 9 31 ...
    ##  - attr(*, "names")= chr [1:12523] "4" "5" "6" "9" ...

``` r
library(MASTA)
```

    ## Loading required package: survival

    ## 
    ## Attaching package: 'MASTA'

    ## The following object is masked _by_ '.GlobalEnv':
    ## 
    ##     FPC.Kern.S

``` r
data <- data_org
TrainSurv <- data.frame(data$TrainSurv)
ValidSurv <- data.frame(data$ValidSurv)
TrainCode <- data.frame(data$TrainCode)
ValidCode <- data.frame(data$ValidCode)
colnames(TrainSurv) <- colnames(ValidSurv) <- c(
  "case", "delta", "sx", "sc", paste0("base_pred", 1:(NCOL(TrainSurv) - 4)))
colnames(TrainCode)[1:2] <- colnames(ValidCode)[1:2] <- c("case", "analysisfu")  
TrainSurv_pred_org <- TrainSurv[-(1:4)]
ValidSurv_pred_org <- ValidSurv[-(1:4)]  
codes <- names(TrainCode[, -(1:3)])
TrainN <- rowsum(TrainCode[, -(1:3)], group = TrainCode[, 1])
ValidN <- rowsum(ValidCode[, -(1:3)], group = ValidCode[, 1])
TrainN <- cbind(case = rownames(TrainN), TrainN)
ValidN <- cbind(case = rownames(ValidN), ValidN)  
colnames(TrainN)[-1] <- paste0(colnames(TrainCode)[-(1:3)], "_total")
colnames(ValidN)[-1] <- paste0(colnames(ValidCode)[-(1:3)], "_total")
TrainPatNum <- unique(TrainCode[, 1])
ValidPatNum <- unique(ValidCode[, 1])
nn <- nrow(TrainSurv) # labeled (training)
nnv <- nrow(ValidSurv) # labeled (validation)
NN <- length(TrainPatNum) - nn # unlabeled  
Tend <- 1
TrainCode$monthstd <- TrainCode$month / TrainCode$analysisfu # standardize follow up time
ValidCode$monthstd <- ValidCode$month / ValidCode$analysisfu # standardize follow up time 
TrainFU <- aggregate(TrainCode$analysisfu, list(TrainCode$case), max)  
ValidFU <- aggregate(ValidCode$analysisfu, list(ValidCode$case), max)
TrainFU <- TrainFU[match(TrainPatNum, TrainFU[, 1]), 2]
ValidFU <- ValidFU[match(ValidPatNum, ValidFU[, 1]), 2]
#--TRAINING---
K <- NULL
ft.e <- ft.e.S <- PKTS <- NULL
FPCA <- vector("list", length(codes))
names(FPCA) <- codes
```

``` r
i <- 1
print(paste("training:", codes[i]))
```

    ## [1] "training: pred1"

``` r
tmp2 <- TrainN[, i + 1] > 0
TrainNP <- TrainN[tmp2, i + 1]
### PKs from Two-step procedure
txt <- paste0("t=with(TrainCode,unlist(sapply(seq_along(month),function(j) rep(month[j],",codes[i],"[j]))))")
eval(parse(text = txt))
id <- (1:(nn+NN))[tmp2]
id <- rep(id, TrainNP)
PKTS <- NULL
PKTS <- cbind(PKTS, GetPK(id = id, t = t, tseq = sort(unique(TrainCode$month)), nn = nrow(TrainN)))    
str(PKTS[!is.na(PKTS)])
```

    ##  int [1:12250] 1 4 19 2 2 1 2 33 43 2 ...

``` r
sort(unique(TrainCode$month))
```

    ##   [1]   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
    ##  [19]  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35
    ##  [37]  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53
    ##  [55]  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71
    ##  [73]  72  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89
    ##  [91]  90  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105 106 107
    ## [109] 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125
    ## [127] 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143
    ## [145] 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161
    ## [163] 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179

``` r
proc.time()
```

    ##    user  system elapsed 
    ##  12.680   0.274  12.906

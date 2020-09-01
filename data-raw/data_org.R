library(data.table)

dir <- "data-raw/data_org/"

TrainSurv <- fread(paste0(dir, "TrainSurv.csv"))
ValidSurv <- fread(paste0(dir, "ValidSurv.csv"))
TrainCode <- fread(paste0(dir, "TrainCode.csv"))
ValidCode <- fread(paste0(dir, "ValidCode.csv"))

TrainLong <- ValidLong <- vector("list", 3)
TrainLong[[1]] <- subset(TrainCode[, c(1:3, 4)], pred1 != 0)
TrainLong[[2]] <- subset(TrainCode[, c(1:3, 5)], pred2 != 0)
TrainLong[[3]] <- subset(TrainCode[, c(1:3, 6)], pred3 != 0)
ValidLong[[1]] <- subset(ValidCode[, c(1:3, 4)], pred1 != 0)
ValidLong[[2]] <- subset(ValidCode[, c(1:3, 5)], pred2 != 0)
ValidLong[[3]] <- subset(ValidCode[, c(1:3, 6)], pred3 != 0)

data_org <- list(
  TrainSurv = TrainSurv,
  ValidSurv = ValidSurv,
  TrainCode = TrainCode, # [rowSums(TrainCode[, -(1:3)]) > 0, ],
  ValidCode = ValidCode) # [rowSums(ValidCode[, -(1:3)]) > 0, ])

data_org$TrainCode[rowSums(data_org$TrainCode[, -(1:3)]) > 0, ]

data_org2 <- list(
  TrainLong = TrainLong,
  ValidLong = ValidLong)

usethis::use_data(data_org, overwrite = TRUE)
usethis::use_data(data_org2, overwrite = TRUE)

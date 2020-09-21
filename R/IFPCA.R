#' @importFrom utils head
#' @importFrom graphics hist
#' @importFrom stats IQR approx bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv density dnorm optimize
NULL

VTM <- function(vc, dm) matrix(vc, ncol = length(vc), nrow = dm, byrow = TRUE)

GetPK <- function(id, t, tseq, nn) {
  PK <- rep(NA, nn)
  e <- diff(tseq)
  aa <- tapply(t, id, FUN = function(t) {
    x <- hist(t, breaks = tseq, plot = FALSE)$counts
    avg_sp <- cumsum(x) / tseq[-1]
    avg_sp2 <- c(0, head(avg_sp, -1))
    which.max((avg_sp - avg_sp2) / (avg_sp2 + e))
  })
  PK[unique(id)] <- tseq[aa + 1]
  PK
}

#' @title Intensity Functional Principal Component Analysis (IFPCA)
#' @description Performs IFPCA to extract features from longitudinal encounter data.
#' @param time longitudinal encounter times
#' @param fu_train follow-up time (training)
#' @param fu_valid follow-up time (validation)
#' @param PPIC_K a logical indicating whether you want to use Pseudo-Poisson Information Criterion to choose 
#' the number of principal components K (K.select="PPIC") \code{TRUE} or another criterion to choose 
#' K (K.select="PropVar") \code{FALSE} in the PP_FPCA_Parallel function (hidden). Default is \code{FALSE}.
#' @param n.grid an integer value for grid points used in estimating covariance function g. Default is \code{401}.
#' @param propvar a proportion of variation used to select number of FPCs. Default is \code{0.85}.
#' @export
ifpca <- function(time, fu_train, fu_valid, 
                  PPIC_K = FALSE, n.grid = 401, propvar = 0.85) {
  Tend <- 1.0
  as_int <- TRUE # FALSE
  names_time <- names(time)
  count <- tapply(time, names_time, length) 
  fu <- c(fu_train, fu_valid)
  time_std <- time / fu[names_time]   
  train <- names(fu_train) 
  valid <- names(fu_valid) 
  time_train <- time[names_time %in% train]
  time_valid <- time[names_time %in% valid]
  names_time_train <- names(time_train)  
  names_time_valid <- names(time_valid)
  if (as_int) {
    names_time <- as.integer(names_time)
    names_time_train <- as.integer(names_time_train)
    names_time_valid <- as.integer(names_time_valid)
  }  
  count_train <- count[names(count) %in% train]
  count_valid <- count[names(count) %in% valid]
  count_train_all <- 0 * fu_train
  count_valid_all <- 0 * fu_valid
  count_train_all[names(count_train)] <- count_train  
  count_valid_all[names(count_valid)] <- count_valid
  # TRAINING
  PKTS <- GetPK(
    id = names_time_train, ### INT/CHAR ### 
    t = time_train, 
    tseq = 0:floor(max(fu_train)), 
    nn = length(fu_train))
  time_std_train <- time_std[names(time_std) %in% train] 
  time1_train <- tapply(time_train, names_time_train, min) 
  h1 <- bw.nrd(time_std_train)
  h2 <- bw.nrd(time_std_train)^(5/6)
  if (PPIC_K) {
    tmp <- PP_FPCA(
      time_std_train, 
      h1 = h1, 
      h2 = h2, 
      count_train, 
      bw = "nrd", 
      ngrid = n.grid, 
      Tend = Tend,
      K.select = "PPIC", 
      derivatives = TRUE)
  } else {
    tmp <- PP_FPCA(
      time_std_train, 
      h1 = h1, 
      h2 = h2, 
      count_train, 
      bw = "nrd", 
      ngrid = n.grid, 
      Tend = Tend, 
      K.select = "PropVar", 
      propvar = propvar, 
      derivatives = TRUE)
  }      
  ft.e <- cbind(
    matrix(fu_train, nrow = length(fu_train), ncol = 3), 
    -tmp$baseline[1], log(1 + count_train_all))
  pos <- count_train_all > 0
  ft.e[pos, 1] <- time1_train  
  locm <- unlist(apply(tmp$densities[, 1:sum(pos) + 1], 2, which.max))
  ft.e[pos, 2] <- ft.e[pos, 2] * tmp$densities[locm, 1]
  ft.e[pos, 3] <- ft.e[pos, 3] * tmp$derivatives[
    sapply(1:sum(pos), function(i) {
      which.max(tmp$derivatives[1:locm[i], i + 1])
    }), 1]
  ft.e[pos, 4] <- tmp$scores[, 2] 
  ft.e.S <- cbind(
    ft.e[, 1], VTM(-tmp$baseline[1:4], length(fu_train)), log(1 + count_train_all))
  ft.e.S[pos, 2:5] <- as.matrix(tmp$scores[, 2:5])
  FPCA <- list(
    K = tmp$K,
    scores = tmp$scores, 
    dens = tmp$densities, 
    deriv = tmp$derivatives,
    mean = tmp$mean, 
    basis = tmp$basis, 
    baseline = tmp$baseline)
  colnames(ft.e) <- c("1stCode", "Pk", "ChP", "1stScore", "logN")
  colnames(ft.e.S) <- c("1stCode", "1stScore", "2ndScore", "3rdScore", "4thScore", "logN")
  rownames(ft.e.S) <- rownames(ft.e) <- names(PKTS) <- train
  # VALIDATION
  PKTS2 <- GetPK(
    id = names_time_valid, ### INT/CHAR ### 
    t = time_valid, 
    tseq = 0:floor(max(fu_valid)), 
    nn = length(fu_valid))  
  time_std_valid <- time_std[names(time_std) %in% valid]   
  time1_valid <- tapply(time_valid, names_time_valid, min)  
  tmp <- PP.FPCA.Pred(time_std_valid, count_valid, FPCA$mean, FPCA$basis, FPCA$K)
  FPCA$ValidPred <- tmp
  ft.e2 <- cbind(
    matrix(fu_valid, nrow = length(fu_valid), ncol = 3), 
    -tmp$baseline[1], log(1 + count_valid_all))
  pos <- count_valid_all > 0  
  locm <- unlist(apply(tmp$densities[, 1:sum(pos) + 1], 2, which.max))  
  ft.e2[pos, 2] <- ft.e2[pos, 2] * tmp$densities[locm, 1]
  ft.e2[pos, 3] <- ft.e2[pos, 3] * tmp$derivatives[
    sapply(1:sum(pos), function(i) {
      which.max(tmp$derivatives[1:locm[i], i + 1])
    }), 1]
  ft.e2[pos, 4] <- tmp$scores[, 2]   
  ft.e.S2 <- cbind(
    ft.e2[, 1], VTM(-tmp$baseline[1:4], length(fu_valid)), log(1 + count_valid_all))
  ft.e.S2[pos, 2:5] <- as.matrix(tmp$scores[, 2:5])
  colnames(ft.e2) <- c("1stCode", "Pk", "ChP", "1stScore", "logN")
  colnames(ft.e.S2) <- c("1stCode", "1stScore", "2ndScore", "3rdScore", "4thScore", "logN")
  rownames(ft.e.S2) <- rownames(ft.e2) <- names(PKTS2) <- valid  
  list(
    TrainFt = ft.e,
    TrainSc = ft.e.S,
    ValidFt = ft.e2,
    ValidSc = ft.e.S2,
    TrainPK = PKTS,
    ValidPK = PKTS2,
    FPCA = FPCA)
}

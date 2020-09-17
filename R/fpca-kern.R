FPC_Kern_S <- function(x, t, N, h1, h2) {
  grp <- rep(seq(N), N)
  M <- outer(t, x, "-")
  D1 <- dnorm(M, 0, h1)
  D2 <- if (h1 == h2) D1 else dnorm(M, 0, h2) 
  v <- rowsum(D2, grp)
  list(f_mu = colSums(D1), 
       Cov_G = crossprod(v) - crossprod(D2))
}

#' @title Mean Density and Covariance Functions for FPCA
#' @description Compute mean density function and covariance function for functional principal components analysis (FPCA)
#' @param x time grid between 0 and 1 
#' @param t observed event times of all the individuals, can have duplicates
#' @param N vector, contains num of observed event of each patient
#' @param h1 bandwidth for mean density function
#' @param h2 bandwidth for covariance function
#' @param bw bandwidth selection method
#' @param nsubs number of subsets for parallel computing
#' @param n_core number of cores for parallel computing
#' @export
FPC.Kern.S <- function(x, t, N, h1 = NULL, h2 = NULL, bw = "ucv", nsubs = NULL, n_core = NULL) {
  if (is.null(n_core)) n_core <- parallel::detectCores()
  if (is.null(nsubs)) nsubs <- n_core * 5
  if (is.null(h1) | is.null(h2)) {
    h <- switch(bw, 
                "nrd0" = bw.nrd0(t),
                "nrd" = bw.nrd(t),
                "ucv" = bw.ucv(t), # leave one out cross validation
                "bcv" = bw.bcv(t), # biased cross validation            
                "SJ-ste" = bw.SJ(t, method = "ste"), 
                "SJ-dpi" = bw.SJ(t, method = "dpi"))    
  }
  if (is.null(h1)) h1 <- h
  if (is.null(h2)) h2 <- h
  n <- length(N)
  a <- n %% nsubs
  b <- floor(n / nsubs)
  subsize <- rep((b + 1):b, c(a, nsubs - a))
  cumsumsub <- cumsum(c(0, subsize))
  grp <- rep(seq(subsize), subsize)
  cumsumN <- c(0, cumsum(rowsum(N, grp))) 
  registerDoParallel(cores = n_core)  
  tmplist <- foreach(i = 1:nsubs) %dopar% {
    tsub <- t[(cumsumN[i] + 1):cumsumN[i + 1]]
    Nsub <- N[(cumsumsub[i] + 1):cumsumsub[i + 1]]
    FPC_Kern_S(x, tsub, Nsub, h1, h2)
  }
  L1 <- lapply(tmplist, `[[`, 1)
  L2 <- lapply(tmplist, `[[`, 2)
  f_mu <- Reduce("+", L1) / sum(N)
  G <- Reduce("+", L2) / sum(N * (N - 1)) - tcrossprod(f_mu)
  list(f_mu = f_mu, G = G)
}

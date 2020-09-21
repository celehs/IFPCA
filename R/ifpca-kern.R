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
#' @description Compute mean density function and covariance function for 
#' functional principal components analysis (FPCA)
#' @param x time grid between 0 and 1 
#' @param t observed event times of all the individuals, can have duplicates
#' @param N vector, contains num of observed event of each patient
#' @param h1 bandwidth for mean density function
#' @param h2 bandwidth for covariance function
#' @param bw bandwidth selection method
#' @export
FPC.Kern.S <- function(x, t, N, h1 = NULL, h2 = NULL, bw = "ucv") {
  if (is.null(h1) | is.null(h2)) {
    h <- switch(
      bw, 
      "nrd0" = bw.nrd0(t),
      "nrd" = bw.nrd(t),
      "ucv" = bw.ucv(t), # leave one out cross validation
      "bcv" = bw.bcv(t), # biased cross validation            
      "SJ-ste" = bw.SJ(t, method = "ste"), 
      "SJ-dpi" = bw.SJ(t, method = "dpi"))    
  }
  if (is.null(h1)) h1 <- h
  if (is.null(h2)) h2 <- h
  grp <- rep(seq(N), N)  
  M <- outer(t, x, "-")  
  D1 <- dnorm(M, 0, h1)
  D2 <- if (h1 == h2) D1 else dnorm(M, 0, h2) 
  v <- rowsum(D2, grp)  
  f_mu <- colSums(D1) / sum(N)
  G <- (crossprod(v) - crossprod(D2)) / 
    sum(N * (N - 1)) - tcrossprod(f_mu)
  list(f_mu = f_mu, G = G)
}

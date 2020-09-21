######################################################################
### Second version of the functional principal component analysis on rare events (Wu et al.,	2013). 
######################################################################
## FPCA approach by Wu et al (2013)
## n = num of patients
## t: observed event times of all the individuals, can have duplicates
## h: bandwidth
## N: vector, contains num of observed event of each patient
## index: index of the patient; default is NULL, i.e., patient are labelled from 1 to n.
## ngrid: number of grid points in estimating covariance function g
## subdivisions: number of subdivisions used in GCV
## propvar: proportion of variation used to select number of FPCs
## PPIC: if TRUE, the propvar will be ignored and PPIC will be used to select number of FPCS
## polybinx: if use the same partition (x) for the polynomial regression
########################################################################
PP_FPCA <- function(t, h1 = NULL, h2 = NULL, N, bw = "ucv", Tend = 1, # assume it is transformed to [0,1]
                    ngrid = 101, K.select = c('PropVar','PPIC'), Kmax = 10,
                    propvar = 0.9, ## For K.select == PropVar
                    density.method = c("kernel","local linear"), ## PPIC
                    polybinx = FALSE, derivatives = FALSE) {
  ## eleminate patients with 0 observed event
  if (sum(N == 0) > 0) {
    NN <- N
    N.index <- N!=0
    N <- N[N.index]
    # cat("Note: patients with zero observed event have been eliminated from analysis!","\n")
  }
  n <- length(N) # number of patients with at least one observed event
  ## if h is null then set it to be the optimal bandwidth under GCV
  if (is.null(h1) & is.null(h2)) {
    if (bw == "nrd0") {
      h1 <- h2 <- bw.nrd0(t)
    } else if (bw == "nrd") {
      h1 <- h2 <- bw.nrd(t)
    } else if (bw == "ucv") { # leave one out cross validation
      h1 <- h2 <- bw.ucv(t)
    } else if (bw == "bcv") { # biased cross validation
      h1 <- h2 <- bw.bcv(t)
    } else if (bw == "SJ-ste") {
      h1 <- h2 <- bw.SJ(t, method = "ste")
    } else if (bw == "SJ-dpi") {
      h1 <- h2 <- bw.SJ(t, method = "dpi")
    } else {
      h1 <- h2 <- bw.ucv(t)
    }
  }
  ## get a fine grid (x,y)
  # tmp  = range(t)
  x <- y <- seq(0, Tend, length.out = ngrid)
  ## estimate the mean density f_mu and the intermediate value g that is used to 
  tmp <- FPC_Kern_S(x, t, N, h1, h2)
  f_mu <- as.numeric(tmp$f_mu) / sum(N)
  G <- tmp$Cov_G / sum(N * (N - 1)) - outer(f_mu, f_mu)
  G.eigen <- svd(G)
  delta <- sqrt(x[2] - x[1])
  baseline <- as.numeric(t(G.eigen$u[, 1:Kmax]) %*% f_mu * delta) # second term in xi
  cumsumN2 <- c(0, cumsum(N))
  # interpolate the eigenvectors to get eigenfunctions
  tmp2 <- apply(G.eigen$u[, 1:Kmax], 2, function(s) approx(x = x, y = s, xout = t)$y) / delta ### check whether we need parallel this
  indx <- rep(1:n, N)
  xi <- -baseline + t(apply(tmp2, 2, FUN = function(x) tapply(x, indx, mean)))# FPC scores, ith column for ith patient
  rm(indx)
  ## select number of FPCs
  K.select <- match.arg(K.select)
  if (K.select == "PropVar") {
    # method 1: proportion of variation >= 90%
    K <- cumsum(G.eigen$d) / sum(G.eigen$d)
    # cat("Propvar=", K, "\n")
    K <- min(which(K >= propvar))
    if (K > Kmax) K <- Kmax
  } else if (K.select == "PPIC") {
    K <- Kmax
    if (missing(density.method)) density.method <- "kernel" else density.method <- match.arg(density.method)
    if (density.method == "local linear") {
      if (polybinx) {
        f_locpoly <- den_locpoly(partition = x, t = t, N = N)
      } else {
        f_locpoly <- den_locpoly2(t = t, N = N)
      }
    } else {
      f_locpoly <- sapply(1:n, function(j) {
        tmp = density(t[{cumsumN2[j] + 1}:cumsumN2[j + 1]], bw = "nrd")
        tmp$y[sapply(t[{cumsumN2[j] + 1}:cumsumN2[j + 1]], function(s) which.min(abs(s - tmp$x)))]
      })
    }
    K <- sapply(1:K, function(k) PPIC(
      K = k, f_locpoly = f_locpoly, t=t, N = N, 
      f_mu = f_mu, G.eigen_v = G.eigen$u, xi = xi, xgrid = x, delta = delta))
    K <- which.min(K)
  }
  ## get density functions 
  if (K == 1) {
    tmp <- f_mu + outer(G.eigen$u[, 1], xi[1, ]) / delta # density functions
    scores <- data.frame(xi[1:max(K, Kmax), ])
  } else {
    tmp <- f_mu + {G.eigen$u[, 1:K] / delta} %*% xi[1:K, ] # density functions
    scores <- data.frame(t(xi[1:max(K, Kmax), ]))
  }
  ## make adjustments to get valid density functions (nonnegative+integrate to 1)
  tmp <- apply(tmp, 2, function(x) {
    x[x < 0] <- 0     # non-negative
    x <- x / sum(x) # integrate to delta^2
    return(x)
  })
  tmp <- tmp / delta^2
  names(scores) <- paste0("score_", seq(1:max(K, Kmax)))  # FPC scores
  basis <- data.frame(cbind(x, G.eigen$u[, 1:max(K, Kmax)] / delta))
  names(basis) <- c("x", paste0("basis_", seq(1:max(K, Kmax)))) # FPC eigenfunctions
  ## name the density functions 
  if (exists("N.index")) {
    colnames(tmp) <- paste0("f", which(N.index))
    scores <- data.frame(id = as.character(which(N.index)), 
                         scores, stringsAsFactors = FALSE) # add id to scores
  } else {
    colnames(tmp) <- paste0("f", seq(1, n))
    scores <- data.frame(id = as.character(seq(1, n)), 
                         scores, stringsAsFactors = FALSE)
  }
  tmp <- as.data.frame(tmp)
  tmp <- data.frame(x = x, tmp)
  # ## return K and prop of var by now
  # tmp     = data.frame(tmp,info=c(K,sum(G.eigen$d[1:K])/sum(G.eigen$d),rep(0,ngrid-2)))
  ## get derivatives if derivatives = TRUE
  if (derivatives) {
    tmp2 <- {tmp[-c(1:2), -1] - tmp[-c(1:2 + ngrid - 2), -1]} / diff(x, lag = 2) 
    tmp2 <- as.data.frame(tmp2)
    tmp2 <- data.frame(x = x[-c(1, ngrid)], tmp2)
    colnames(tmp2) <- paste0("d", colnames(tmp))
    return(list(scores = scores, 
                densities = tmp, 
                derivatives = tmp2, 
                mean = data.frame(x = x, f_mu = f_mu),
                cov = G, 
                basis = basis, 
                baseline = baseline,
                K = K))
  } else {
    return(list(scores = scores, 
                densities = tmp, 
                mean = data.frame(x = x, f_mu = f_mu), 
                cov = G, 
                basis = basis, 
                baseline = baseline,
                K = K))
  }
}

### Predict Scores for New Subjects
### make prediction for patients with or without codes
PP.FPCA.Pred <- function(t, N, mean.fun, eigen.fun, K) {
  NZ <- N == 0
  delta <- mean.fun[2, 1] - mean.fun[1, 1]
  baseline <- as.numeric(t(eigen.fun[, -1]) %*% mean.fun[, -1] * delta) # second term in xi
  tmp2 <- apply(eigen.fun[, -1], 2, function(s) approx(x = mean.fun$x, y = s, xout = t)$y) ### check whether we need parallel this
  indx <- rep(1:length(N[!NZ]), N[!NZ])
  xi <- -baseline + t(apply(tmp2, 2, FUN = function(x) tapply(x, indx, mean))) # FPC scores, ith column for ith patient
  if (K == 1) {
    tmp <- mean.fun[, -1] + outer(eigen.fun[, 2], xi[1, ]) # density functions
  } else {
    tmp <- as.numeric(mean.fun[, -1]) + as.matrix(eigen.fun[, 1:K + 1]) %*% xi[1:K, ] # density functions
  }
  tmp <- apply(tmp, 2, function(x) {
    x[x < 0] <- 0 # non-negative
    x <- x / sum(x) # integrate to delta^2
    return(x)
  })
  tmp <- tmp / delta
  tmp2 <- {tmp[-c(1:2), ] -tmp[-c(1:2 + length(mean.fun[, 1]) - 2), ]} / diff(mean.fun[, 1], lag = 2)
  derivatives <- cbind({mean.fun[-c(1:2), 1] + mean.fun[-c(1:2 + length(mean.fun[, 1]) - 2), 1]} / 2, tmp2)
  scores <- t(xi)
  densities <- cbind(mean.fun[, 1], tmp)  
  rownames(scores) <- colnames(densities)[-1] <- colnames(derivatives)[-1] <- names(N)
  list(scores = scores,
       densities = densities,
       derivatives = derivatives,
       baseline = baseline)
}

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

VTM <- function(vc, dm) matrix(vc, ncol = length(vc), nrow = dm, byrow = TRUE)

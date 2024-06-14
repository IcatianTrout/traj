#Compute Calinski-Harabasz index
CalinskiHarabasz <- function (x, clustering, cn = max(clustering)) 
{
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  cln <- rep(0, cn)
  W <- matrix(0, p, p)
  for (i in 1:cn) cln[i] <- sum(clustering == i)
  for (i in 1:cn) {
    clx <- x[clustering == i, ]
    cclx <- var(as.matrix(clx))
    if (cln[i] < 2) 
      cclx <- 0
    W <- W + ((cln[i] - 1) * cclx)
  }
  S <- (n - 1) * var(x)
  B <- S - W
  out <- (n - cn) * sum(diag(B))/((cn - 1) * sum(diag(W)))
  out
}
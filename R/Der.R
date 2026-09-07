# Function that approximates the derivative of a function as a weighted average of the left and right derivatives
Der <- function(x, y) {
  
  # Here x and y are vectors, and y is interpreted as f(x) for the function f that we're trying to approximate the derivative of
  
  if (!length(y) == length(x)) {
    stop("y and x must be vectors of the same length.")
  }
  
  m <- length(x)
  
  if (!identical(order(x), c(1:m))) {
    stop("The elements of the 'x' vector must be strictly increasing.")
  }
  DR <- rep(NA, m)
  DL <- rep(NA, m)
  D <- c()
  
  for (i in 1:(m - 1)) {
    DR[i] <- (y[i + 1] - y[i]) / (x[i + 1] - x[i]) ## Approximate the right derivative
  }
  
  for (i in 2:m) {
    DL[i] <- (y[i - 1] - y[i]) / (x[i - 1] - x[i]) ## Approximate the left derivative
  }
  
  D[1] <- DR[1]
  D[m] <- DL[m]
  
  wL <- rep(NA, m)
  wR <- rep(NA, m)
  
  for (i in 2:(m - 1)) {
    wL[i] <- (x[i + 1] - x[i]) / (x[i + 1] - x[i - 1]) ## Weight for the left derivative
    wR[i] <- 1 - wL[i] ## Weight for the right derivative
    D[i] <- wL[i] * DL[i] + wR[i] * DR[i] ## The approximation is the weighted average of the left and right derivatives
  }
  
  return(D)
}

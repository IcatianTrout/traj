# Function that approximates the mean of a function using the mean of the left and right Riemann sums

FctMean <- function (x, y) {
  
  # Here x and y are vectors, and y is interpreted as f(x) for the function f that we're trying to approximate the mean of
  
  if (!length(y) == length(x)) {
    stop("y and x must be vectors of the same length.")
  }
  
  m <- length(x)
  
  if (!identical(order(x), c(1:m))) {
    stop("The elements of the 'x' vector must be strictly increasing.")
  }
  
  RR <- 0
  RL <- 0
  for (i in 1:(m - 1)) {
    RR <- RR + y[i] * (x[i + 1] - x[i])
    RL <- RL + y[i + 1] * (x[i + 1] - x[i])
  }
  
  mu <- (1 / 2) * (RR + RL) / (x[m] - x[1])
  
  return(mu)
}

#'@title Computes the third quartile
#'
#'@param x vector of real numbers
#'
Q3 <- function(x) {
  return(quantile(x , probs = c(.75)))
}
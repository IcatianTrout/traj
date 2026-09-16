#'@title Computes the first quartile
#'
#'@param x vector of real numbers
#'
Q1 <- function(x) {
  return(quantile(x , probs = c(.25)))
}

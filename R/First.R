#'@title Function that returns the first non-NA value of a vector
#'
#'@param v vector of real numbers
#'
First <- function(v) {
  if (!(FALSE %in% is.na(v))) {
    stop("Argument must contain at least one non-NA entry.")
  }
  
  w <- v[complete.cases(v)]
  
  return(w[1])
}

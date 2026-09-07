# Function that returns the last non-NA value of a vector

Last <- function(v) {
  if (!(FALSE %in% is.na(v))) {
    stop("Argument must contain at least one non-NA entry.")
  }
  
  w <- v[complete.cases(v)]
  
  m <- length(w)
  
  return(w[m])
}

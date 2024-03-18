#returns the first non-NA coordinate of a vector
First <- function(v) {
  if (!(FALSE %in% is.na(v))) {
    stop("Argument must contain at least one non-NA entry.")
  }
  
  w <- v[complete.cases(v)]
  
  return(w[1])
}

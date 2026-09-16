#'@title Function to prevent another function from printing things as part of its execution
#'
#'@param x code 
#'
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
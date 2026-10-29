#'@title Function that approximates the derivative of a function as a weighted average of the left and right derivatives
#'
#'@param x vector of strictly increasing real numbers
#'@param y vector of real numbers corresponding to f(x)
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
  
  D[1] <- (y[2] - y[1]) / (x[2] - x[1])
  D[m] <- (y[m - 1] - y[m]) / (x[m - 1] - x[m])
  
  ## Compute derivative at the interior points from the standard three-point finite-difference formula for unequally spaced observations
  
  for (i in 2:(m - 1)) {
    hm <- x[i] - x[i - 1]
    hp <- x[i + 1] - x[i]
    
    D[i] <- (-hp * y[i - 1]) / (hm * (hm + hp)) + ((hp - hm) * y[i]) / (hm * hp) + (hm * y[i + 1]) / (hp * (hm + hp)) 
  }
  
  if(2 %in% which){
    
    ## Compute second derivative at the boundary
    A <- c()
    
    A[1] <- (D[2] - D[1]) / (x[2] - x[1])
    A[m] <- (D[m - 1] - D[m]) / (x[m - 1] - x[m])
    
    ## Compute derivative at the interior points from the standard three-point finite-difference formula for unequally spaced observations
    
    for (i in 2:(m - 1)) {
      hm <- x[i] - x[i - 1]
      hp <- x[i + 1] - x[i]
      
      A[i] <- 2 / (hm + hp) * ((y[i + 1] - y[i]) / hp - (y[i] - y[i - 1]) / hm)
    }
  }
  
  if(which == 1){
    return(D)
  }
  if(which == 2){
    return(A)
  }
  if(which == c(1,2)){
    return(rbind(D,A))
  }
}




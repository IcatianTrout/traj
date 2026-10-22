#'@title Function that approximates the first and second derivative of a function.
#'
#'@param x vector of strictly increasing real numbers
#'@param y vector of real numbers corresponding to f(x)
#'@param which vector containing 1 to signal that the first derivative is required and 2 to signal that the second derivative is required 
#'
#'@return If the \code{which} argument contains both 1 and 2, then the function returns a matrix whose first row is the first derivative and the second row is the second derivative. If the \code{which} argument contains only 1 or only 2, the function returns a vector containing the approximate values of the first and second derivative respectively.
#'
#'@details At the interior points, the standard three-point finite-difference formula for unequally spaced observations is used. At the boundary, the first derivative is approximated directly from one-sided derivatives. Likewise, the second derivative at the boundary points is approximated from the one-sided derivatives of the first derivative.
#'
Der <- function(x, y, which = c(1, 2)) {
  
  ## Compute derivative at the boundary
  m <- length(x)
  
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




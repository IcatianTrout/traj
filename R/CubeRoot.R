#'@title Compute the non negative cube root of a real number 
#'
#'@param x a non negative number
#'
CubeRoot <- function(x){
  
  sign(x)*abs(x)^(1/3)
  
}

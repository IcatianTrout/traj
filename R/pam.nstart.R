## The clusGap function takes in its 'FUNcluster' argument a function that will perform
## the clustering but it does not allow us to set the arguments of that
## function, so we have to do so beforehand by defining a new function:

pam.nstart <- function (x, k, metric = metric) {
  return(pam(
    x = x,
    k = k,
    diss = FALSE,
    metric = metric
  ))
}
## The clusGap function takes in its 'FUNcluster' argument a function that will perform
## the clustering but it does not allow us to set the arguments of that
## function, so we have to do so beforehand by defining a new function:
kmeans.nstart <- function (x, k) {
  return(kmeans(
    x = x,
    centers = k,
    nstart = nstart,
    iter.max = iter.max
  ))
}

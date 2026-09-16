#'@title Construct a similarity matrix 
#'
#'@description Define a function that construct a sparse similarity matrix (aka adjacency matrix) based on the K-nearest neighbor (KNN) principle. It takes as input a matrix X of "proximities" between points as well as the number K of neighbors to consider and outputs a sparse matrix where cell (i,j) is 0 if neither point i nor point j is among the KNN of the other, 1 if i is among the KNN of j and vice versa, and 1/2 otherwise (i.e. if i is among the KNN of j *or* vice versa, but not both)
#'
#'@param X matrix X of "proximities" between the trajectories. These are non-negative numbers with the interpretation that a small proximity between two trajectories means a large similarity.
#'@param K The number of nearest neighbors to consider for each trajectories
#'
knn_adjacency <- function(X, K){
  W <- matrix(0, nrow = nrow(X), ncol = nrow(X)) 
  for(i in seq_len(nrow(X))){
    knn <- order(X[i, ])[-which(order(X[i, ]) == i)][seq_len(K)] ## Find the K nearest neighbors to the ith data point, excluding the ith data point itself
    W[i, knn] <- 1
  }
  return((W + t(W)) / 2)
}
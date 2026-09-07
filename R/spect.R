# Function that creates a sparse similarity matrix and applies to it the spectral clustering algorithm described in Meila (2005)

spect <- function(x, k, nstart, fuzzy){ 
  
  # Define a function that construct an adjacency matrix (aka similarity matrix) based on the K-nearest neighbor (KNN) principle. It takes as input a matrix X of "distances" between points as well as the number K of neighbors to consider and outputs a sparse matrix where cell (i,j) is 0 if neither point i nor point j is among the KNN of the other, 1 if i is among the KNN of j and vice versa, and 1/2 otherwise (i.e. if i is among the KNN of j *or* vice versa, but not both)
  knn_adjacency <- function(X, K){
    W <- matrix(0, nrow = nrow(X), ncol = nrow(X)) 
    for(i in seq_len(nrow(X))){
      knn <- order(X[i, ])[-which(order(X[i, ]) == i)][seq_len(K)] ## Find the K nearest neighbors to the ith data point, excluding the ith data point itself
      W[i, knn] <- 1
    }
    return((W + t(W)) / 2)
  }
  
  n <- nrow(x)
  K <- max(4, min(8, floor(n / (2 * k)))) ## The number of nearest neighbors to consider is a number between 4 and 8 depending on the sample size n and the number of clusters k 
  S <- Matrix::Matrix(knn_adjacency(X = as.matrix(stats::dist(x)), K = K), sparse = TRUE) ## The similarity matrix
  Dsq.inv <- Matrix::Diagonal(x = 1/sqrt(Matrix::rowSums(S)))
  L <- Dsq.inv %*% S %*% Dsq.inv 
  
  eigen_result <- RSpectra::eigs_sym(L, k + 5, which = "LA") ## Requesting 5 more eigenvalues than the k that are needed improves estimation  
  eigenvalues <- eigen_result$values[seq_len(k)]
  eigenvectors <- eigen_result$vectors[, seq_len(k)]
  
  # If the first eigenvector is a multiple of (1, 1, ...,1), remove it
  if(length(unique(eigenvectors[, 1])) == 1){
    Y <- eigenvectors[, -1, drop = FALSE]
  } else{
    Y <- eigenvectors
  }
  X <- Dsq.inv %*% Y
  
  # Normalize the rows of X
  for(l in seq_len(n)){
    if(!(sum(X[l, ])^2 == 0)){
      X[l, ] <- X[l, ] / sqrt(sum(X[l, ]^2))
    }
  }
  
  # Initialize a bunch of things
  CH.v <- c() 
  partition.matrix <- matrix(NA, nrow = nstart, ncol = n)
  row.centers.matrix <- matrix(NA, nrow = nstart, ncol = k)
  fuzzy.partition.list <- list()
  sq.sum <- function(y, g, cl){return((y - cl$centers[g, ])^2)} ## Function that computes the squared difference between between y and cl$centers[g, ]
  find.row <- function(y, w, aux){return(identical(y, aux[w, ]))} ## Function that outputs TRUE if y coincides with aux[w, ]
  
  # Perform k-means (or fuzzy k-means) clustering on the rows of X but initialize the algorithm with random cluster centers that are close to being mutually orthogonal (because this is what we expect in the ideal scenario where the graph corresponding to S has k connected components). We do this by looking at the cosine of the dot product, which is close to 1 (its theoretical maximum) iff the vectors are close to orthogonal. Do this 'nstart' times, storing away the results at each iteration.
  for(s in seq_len(nstart)){
    
    # Pick the first center randomly among the rows of X
    center.matrix <- as.matrix(X[sample(seq_len(n), size = 1), , drop = F])
    
    # Pick the rest of the centers among the rows of X to be as mutually orthogonal as possible
    for(p in 2:k){
      dot.prod <- center.matrix %*% t(X)
      cos.matrix <- cos(dot.prod)
      w <- which(apply(cos.matrix, 2, min) == max(apply(cos.matrix, 2, min))) ## The rows of X with the largest minimum cos dot product with the centers 
      # In case of a tie at the cosine level, favor large angles over small ones. I.e. favor negative dot products over positive ones
      if(length(w) > 1){
        # If one of the centers is 0, remove from contention those rows of X which are 0
        if(0 %in% apply(t(apply(center.matrix, 1, abs)), 1, sum)){ 
          w0 <- which(apply(t(apply(X[w, ], 1, abs)), 1, sum) == 0) 
          w <- w[-w0] 
        }
        w <- w[which( colSums(dot.prod[, w, drop = FALSE] < 0) == max(colSums(dot.prod[, w, drop = FALSE] < 0)) )] ## Retain the rows of X that have negative dot product with the most amount of centers
        # If there are more than one rows with these properties, pick one randomly
        if(length(w) > 1){
          w <- sample(w, 1)  
        }
      }
      center.matrix <- rbind(center.matrix, X[w, ])
    }
    
    if(fuzzy == FALSE){
      cl <- stats::kmeans(X, centers = center.matrix, iter.max = 100) ## Cluster the "X representation" of the trajectories using k-means
      CH.v[s] <- ((n - k) * cl$betweenss) / ((k - 1) * cl$tot.withinss) ## The Calinski-Harabasz index of the clustering
      partition.matrix[s, ] <- cl$cluster
      fuzzy.partition.list <- NULL
      
      for(g in seq_len(k)){
        aux <- X[which(cl$cluster == g), , drop = FALSE]
        # Compute the distance (in the "X representation") between the cluster means and each element of said cluster
        if(ncol(X) > 1){
          distance <- sqrt(colSums(apply(aux, 1, FUN = sq.sum, g = g, cl = cl))) 
        } else {
          distance <- abs(aux - cl$centers[g, ])
        }
        w <- which(distance == min(distance))[1]
        row.centers.matrix[s, g] <- which(apply(X, 1, find.row, w = w, aux = aux) == TRUE)[1] ## Find the (first) row of X that's identical to the row aux which is closest to the center. This will be the centroid for cluster g
      }
    }
    if(fuzzy == TRUE){
      cl <- e1071::cmeans(X, centers = center.matrix, iter.max = 100) ## Cluster the "X representation" of the trajectories using fuzzy k-means
      CH.v[s] <- unlist(clusterCrit::intCriteria(as.matrix(X), part = cl$cluster, crit = "calinski_harabasz")$calinski_harabasz) ## The Calinski-Harabasz index of the corresponding hard clustering
      partition.matrix[s, ] <- cl$cluster
      fuzzy.partition.list[[s]] <- cl$membership
      for(g in seq_len(k)){
        aux <- X[which(cl$cluster == g), , drop = FALSE]
        # Compute the distance (in the "X representation") between the cluster means and each element of said cluster
        if(ncol(X) > 1){
          distance <- sqrt(colSums(apply(aux, 1, FUN = sq.sum, g = g, cl = cl))) 
        } else {
          distance <- abs(aux - cl$centers[g, ])
        }
        w <- which(distance == min(distance))[1]
        row.centers.matrix[s, g] <- which(apply(X, 1, find.row, w = w, aux = aux) == TRUE)[1] ## Find the (first) row of X that's identical to the row aux which is closest to the center. This will be the centroid for cluster g
      }
    }
  }
  
  w <- which(CH.v == max(CH.v))[1] ## The iteration that gave rise to the largest the Calinski-Harabasz index
  if(!is.null(fuzzy.partition.list)){fuzzy.partition <- fuzzy.partition.list[[w]]} else(fuzzy.partition <- NULL)
  output <- list(cluster = partition.matrix[w, ], 
                 row.centers = row.centers.matrix[w, ], 
                 fuzzy.partition = fuzzy.partition, 
                 X = X)
  
  return(output)
}
spect <- function(x, k, nstart, fuzzy){ 
  
  knn_adjacency <- function(X, K){
    W <- matrix(0, nrow = nrow(X), ncol = nrow(X)) 
    for(j in 1:ncol(X)){
      knn <- order(X[j, ])[-which(order(X[j, ]) == j)][1:K] # find the K nearest points to the j-th data point, excluding the j-th data point.
      W[knn, j] <- 1
    }
    return((W + t(W))/2)
  }
  
  n <- nrow(x)
  K <- max(4, min(8, floor(n/(2*k)))) #number between 4 and 8
  S <- knn_adjacency(X = as.matrix(stats::dist(x)), K = K)
  D <- diag(rowSums(S))
  L <- solve(sqrt(D)) %*% S %*% solve(sqrt(D))
  eigen_result <- eigen(L)
  eigenvalues <- eigen_result$values
  eigenvectors <- eigen_result$vectors
  
  
  if( abs(eigenvalues[2] - 1) > 1e-6 ){
    Y <- as.matrix(eigenvectors[, c(2:k)])
  } else {
    Y <- as.matrix(eigenvectors[, c(1:k)])
  }
  X <- solve(sqrt(D)) %*% Y
  for(l in 1:n){
    if(!(sum(X[l, ])^2 == 0)){
      X[l, ] <- X[l, ]/sqrt(sum(X[l, ]^2))
    }
  }
  
  CH.v <- c()
  partition.matrix <- matrix(NA, nrow = nstart, ncol = n)
  row.centers.matrix <- matrix(NA, nrow = nstart, ncol = k)
  fuzzy.partition.list <- list()
  sq.sum <- function(y, g, cl){return((y - cl$centers[g, ])^2)}
  find.row <- function(y, w, aux){return(identical(y, aux[w, ]))}
  
  for(s in 1:nstart){
    
    center.matrix <- as.matrix(X[sample(1:n, size = 1), , drop = F]) #pick the first center randomly
    
    #pick the rest of the centers to be as mutually orthogonal as possible
    for(p in 2:k){
      dot.prod <- center.matrix %*% t(X)
      cos.matrix <- cos(dot.prod)
      w <- which(apply(cos.matrix, 2, min) == max(apply(cos.matrix, 2, min)))
      if(length(w) > 1){ #in case of a tie at the cosine level, favor large angles over small ones
        if(0 %in% apply(t(apply(center.matrix, 1, abs)), 1, sum)){ #If one of the centers is 0, remove from contention those rows of X which are 0
          w0 <- which(apply(t(apply(X[w, ], 1, abs)), 1, sum) == 0) 
          w <- w[-w0] 
        }
        w <- w[which( colSums(dot.prod[, w, drop = FALSE] < 0) == max(colSums(dot.prod[, w, drop = FALSE] < 0)))]
        if(length(w) > 1){
          w <- sample(w, 1)
        }
      }
      center.matrix <- rbind(center.matrix, X[w, ])
    }
    
    if(fuzzy == FALSE){
      cl <- stats::kmeans(X, centers = center.matrix, iter.max = 100) 
      CH.v[s] <- ((n-k)*cl$betweenss)/((k-1)*cl$tot.withinss)
      partition.matrix[s, ] <- cl$cluster
      fuzzy.partition.list <- NULL
      
      for(g in 1:k){
        aux <- X[which(cl$cluster == g), , drop = FALSE]
        if(ncol(X) > 1){
          distance <- sqrt(colSums(apply(aux, 1, FUN = sq.sum, g = g, cl = cl)))
        } else {
          distance <- abs(aux - cl$centers[g, ])
        }
        w <- which(distance == min(distance))[1]
        row.centers.matrix[s, g] <- which(apply(X, 1, find.row, w = w, aux = aux) == TRUE)[1]
      }
    }
    if(fuzzy == TRUE){
      cl <- e1071::cmeans(X, centers = center.matrix, iter.max = 100)
      CH.v[s] <- unlist(clusterCrit::intCriteria(as.matrix(X), part = cl$cluster, crit = "calinski_harabasz")$calinski_harabasz)
      partition.matrix[s, ] <- cl$cluster
      fuzzy.partition.list[[s]] <- cl$membership
      for(g in 1:k){
        aux <- X[which(cl$cluster == g), , drop = FALSE]
        if(ncol(X) > 1){
          distance <- sqrt(colSums(apply(aux, 1, FUN = sq.sum, g = g, cl = cl)))
        } else {
          distance <- abs(aux - cl$centers[g, ])
        }
        w <- which(distance == min(distance))[1]
        row.centers.matrix[s, g] <- which(apply(X, 1, find.row, w = w, aux = aux) == TRUE)[1]
      }
    }
  }
  
  w <- which(CH.v == max(CH.v))[1] #the clustering that maximizes the Calinski-Harabasz index
  if(!is.null(fuzzy.partition.list)){fuzzy.partition <- fuzzy.partition.list[[w]]} else(fuzzy.partition <- NULL)
  output <- list(cluster = partition.matrix[w, ], row.centers = row.centers.matrix[w, ], fuzzy.partition = fuzzy.partition, X = X)
  
  return(output)
}
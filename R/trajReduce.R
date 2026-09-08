#'@title Select a Subset of the Measures Using a Similarity Index on the Set of Clusterings
#'
#'@description This function examines the effect of reducing the number of measures on which the trajectories are clustered. Specifically, starting from a clustering \eqn{C} in the form of an object of class \code{trajClusters} and a choice of a similarity index to compare clusterings, this function finds the subset of measures which results in the clustering most similar to \eqn{C}.
#'
#'@param Measures object of class \code{trajMeasures} as returned by \code{\link[traj]{trajMeasures}}.
#'@param Clusters object of class \code{trajClusters} as returned by \code{\link[traj]{trajClusters}}.
#'@param index The similarity index. Either "ARI" for the Adjusted Rand Index of Hubert and Arabie (1985), "nVId" for the normalized variation of information distance (eg. Meila (2007)) or "nSJd" for the normalized split/joint distance of van Dongen (2000).
#'@param keep The number of measures to keep. Defaults to 3.
#'
#'@details The Rand index ranges from 0 to 1 with 0 indicating identical clusters and 1 indicating maximally different clusters. The normalized variation of information distance (nVId) and normalized split-join distance (nSJd) and have the opposite interpretation with 0 indicating maximally different clusters and 1 indicating identical clusters. Therefor, to facilitate comparison, we plot 1 - nVId (resp. 1 - nSJd) instead of nVId (resp. nSJd).
#'
#'The function runs the \code{\link[traj]{trajClusters}} function on every combination of \code{keep} measures and identifies the reduced representation yielding the clustering most similar to the (hard) clustering in \code{Clusters}. If the (hard) clustering in \code{Clusters} derives from a soft clustering, the iterations of \code{\link[traj]{trajClusters}}  will also be ran using \code{fuzzy = TRUE} and the corresponding hard clusterings will be used for comparison.
#'
#'@importFrom igraph compare
#'
#'@references 
#'Hubert L, Arabie P. Comparing partitions. Journal of Classification 2:193-218, 1985.
#'
#'Meila M. Comparing clusterings -- an information based distance. Journal of Multivariate Analysis, 98, pp 873-895, 2007.
#'  
#'van Dongen S. Performance criteria for graph clustering and Markov cluster experiments. Technical Report INS-R0012, National Research Institute for Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000.
#'  
#'
#' @examples
#' \dontrun{
#'data("trajdata")
#'trajdata.noGrp <- trajdata[, -which(colnames(trajdata) == "Group")] #remove the Group column
#'
#'m = trajMeasures(trajdata.noGrp, ID = TRUE, measures = 1:20)
#'s2.4 <- trajClusters(m, nclusters = 4)
#'trajReduce(m, s2.4)
#'}
#'
#'@rdname trajReduce
#'
#'@export
trajReduce <-
  function (Measures, Clusters, index = "ARI", keep = 3) {
    
    if(is.null(Clusters$nclusters)){
      stop("The 'Clusters' argument must be the output of the 'trajClusters' function in which the 'nclusters' argument is not 'NULL'.")
    }
      
    # Function that inputs a list of vectors and outputs a new list containing every vector obtained by removing one element at a time from each vector in the input list
      useful.fct <- function(list_of_vectors){
        output.list <- list()
        
        for(i in seq_along(list_of_vectors)){
          ith.vector <- list_of_vectors[[i]]
          n.i <- length(ith.vector)
          
          for(j in seq_len(n.i)){
            output.list[[length(output.list) + 1]] <- ith.vector[-j]
          }
        }
        
        return(output.list)
      }
      
      traj <- as.factor(Clusters$partition[, "Cluster"]) ## The clustering in the 'Clusters' argument
      nCluster <- Clusters$nclusters ## The number of clusters
      n <- nrow(Clusters$partition) ## The sample size
      
      combin <- utils::combn(Clusters$select, keep) ## All subsets of keep elements among Clusters$select
      
      criterion.v <- c()
      
      ## Run trajClusters() on each combination of measures and compute the similarity index of the resulting clustering with the original clustering  
      for(i in seq_len(ncol(combin))){
        s3.i <- quiet(trajClusters(Measures, 
                                   select = combin[, i], 
                                   fuzzy = Clusters$fuzzy, 
                                   nclusters =  Clusters$nclusters))
        
        traj.i <- as.factor(s3.i$partition[, "Cluster"])
        
        if (index == "ARI") {
          criterion.v <- c(criterion.v, igraph::compare(traj, traj.i, method = "adjusted.rand"))
        }
        if (index == "nVId") {
          criterion.v <- c(criterion.v, 1 - igraph::compare(traj, traj.i, method = "vi") / (2 * log(nCluster)))
        }
        if (index == "nSJd") {
          criterion.v <- c(criterion.v, 1 - igraph::compare(traj, traj.i, method = "split.join") / (2 * nCluster * (n / nCluster - ceiling(n / (nCluster^2)))) )
        }
      }
      
      w <- which(criterion.v == max(criterion.v))[1] ## The simplest combination of measure for which the similarity index is maximal
      
      Clusters.red <- quiet(trajClusters(Measures, 
                                         select = combin[, w], 
                                         fuzzy = Clusters$fuzzy, 
                                         nclusters = Clusters$nclusters))
      
      traj.red <- as.factor(Clusters.red$partition[, "Cluster"])
      
      output <- list(
        reduced.list = combin[, w],
        index.value = criterion.v[w][1],
        partition.red = traj.red,
        fuzzy.partition.red = Clusters.red$fuzzy.partition
      )
      
      return(output)
  }
#'@title Classify the Longitudinal Data Based on the Measures
#'
#'@description Classifies the trajectories by applying a nonparametric clustering algorithm to the measures computed by \code{trajMeasures()}.
#'
#'@param Measures object of class \code{trajMeasures} as returned by the function
#'  \code{trajMeasures()}.
#'@param select an optional vector of positive integers corresponding to the
#'  measures to use in the clustering. Defaults to \code{NULL}, which uses all the measures contained in \code{Measures}.
#'@param fuzzy logical. If FALSE, each trajectory is assigned to a unique group. If TRUE, each trajectory is assigned a "degree of membership" to each group. Defaults to FALSE.
#'@param nclusters The desired number of clusters. If \code{NULL}, clustering is carried out for every number of clusters between 2 and (up to) 8 and the "best" number of clusters is used, as judged by the combination of three internal cluster validity indices. See section 'Value' for more details. Defaults to \code{NULL}.
#'@param subset.n A positive numerical integer smaller than the number of trajectories. If \code{nclusters} is \code{NULL}, \code{subset} is the number of trajectories, randomly sampled from the complete data set, that will be used to determine the optimal number of clusters in the interest of speeding up the process.
#'@param nstart The number of random starts. Defaults to \code{50}.
#'@param x object of class \code{trajClusters}.
#'@param object object of class \code{trajClusters}.
#'@param ... further arguments passed to or from other methods.
#'
#'@details The spectral clustering algorithm presented in Meila (2005) is implemented in which the similarity matrix \eqn{S} is built from a binary K nearest neighbors similarity function (\eqn{S=(W+W^T)/2}, where \eqn{W_{ij}=1} if data point \eqn{j} is among the nearest points to data point \eqn{i} and \eqn{W_{ij}=0} otherwise). 
#'
#'When \code{nclusters = NULL}, the function evaluates candidate clusterings with number of clusters \eqn{k} ranging from 2 to 8 using three internal validity indices: C-index, Calinski-Harabasz and Wemmert-Gançarski. These indices are normalized so that the highest value is 1 and the lowest is 0, and so that a high value is synonymous with high validity. The optimal number of clusters is determined according to a ranked voting system in which each index contributes a fractional vote according to its ranking of the candidate solutions. Specifically, each index casts a vote worth 1 in favor of \eqn{k} if it takes its greatest value when the number of groups is \eqn{k}, worth 5/6 if it takes its second greatest value when the number of groups is \eqn{k}, and so on down to a vote worth 0 if the index takes its smallest value when the number of groups is \eqn{k}. The favorability of \eqn{k} is the sum of the 3 votes. 
#'
#'@return An object of class \code{trajClusters}; a list containing the result of the clustering, as well as a curated form of the arguments.
#'
#'@import cluster fclust clusterCrit 
#'@importFrom stats quantile kmeans
#'@importFrom e1071 cmeans
#'
#'@references 
#'
#'Meila, M., Spectral Clustering. Handbook of Cluster Analysis, Chapter 7, Chapman and Hall/CRC, 2005.
#'
#' @examples
#' 
#'data(trajdata)
#'
#'dat <- trajdata[, -2] #remove the Group column
#'
#'m <- trajMeasures(dat, ID = TRUE, measures = c(1:20))
#'
#'s <- trajClusters(m, nclusters = 4)
#'plot(s, which.plots = 2, ask = FALSE)
#'
#'dat$clusters <- s$partition[match(trajdata$id, s$partition[, 1]), 2]
#'tail(dat)
#'
#'
#'@rdname trajClusters
#'
#'@export

trajClusters <-
  function (Measures,
            select = NULL,
            fuzzy = FALSE,
            nclusters = NULL,
            subset.n = NULL,
            nstart = 50
  ) {
    
    # Perform checks on the arguments
    if ((!is.null(select)) & ((!is.numeric(select)) | (!is.vector(select)))) {
      stop("Argument 'select' must be either NULL or a numerical vector.")
    } else {
      if(is.null(select)) {
        select <- Measures$measures.arg
      }
      m.select <- paste("m", select, sep = "")
      if (FALSE %in% (m.select %in% colnames(Measures$measures[, -1, drop = FALSE]))) {
        stop("The 'select' argument must only contain measures included in Measures.")
      }
    }
    
    if (!(fuzzy %in% c("TRUE", "FALSE"))) {
      stop("'fuzzy' should be either 'TRUE' or 'FALSE'.")
    }
    
    if(is.null(nclusters)){
      if ( !is.null(subset.n) && !( (length(subset.n) == 1) && (subset.n %in% seq_len(nrow(Measures$measures))) ) ){
        stop("'subset.n' should be a numerical integer smaller than the total number of admissible trajectories.")
      }
    } else if (!is.numeric(nclusters)){
      stop("The number 'nclusters' of requested clusters should be a numerical integer.")
    } else if (nclusters > nrow(Measures$measures)) {
      stop("The number 'nclusters' of requested clusters cannot exceed the number of trajectories.")
    }
    
    # Initiate a bunch of variables to be used later
    ID <- Measures$measures[, 1]
    nclusters.input <- nclusters
    partition.summary <- NULL
    clust.by.id <- Measures$data[, 1, drop = FALSE]
    rownames(clust.by.id) <- NULL
    ICV.raw <- NULL
    ICV <- NULL
    bins <- NULL
    selection <- Measures$measures[, c("ID", m.select), drop = FALSE]
    
    # Check if there are constant measures. If so, remove them from the analysis since they are non discriminating
    dat <- Measures$measures[, m.select, drop = FALSE]
    w <- which(apply(dat, 2, sd) == 0)
    if(length(w) == 1){
      meas.rmv <- colnames(dat)[w]
      warning(paste("Being constant, measure ", noquote(paste(meas.rmv, collapse = ", ")), " has been removed.", sep = ""))
      dat <- dat[, -w, drop = FALSE]
    } else if(length(w) > 1){
      meas.rmv <- colnames(dat)[w]
      warning(paste("Being constant, measures ", noquote(paste(meas.rmv, collapse = ", ")), " have been removed.", sep = ""))
      dat <- dat[, -w, drop = FALSE]
    } 
    
    # Standardize the measures to be clustered
    dat <- data.frame(apply(dat, 2, scale))
    
    # If the desired number of clusters 'nclusters' was left unspecified, define it as the winner of the ranked voting system of three internal cluster validity (ICV) criteria
    if (is.null(nclusters)) {
      
      k.max <- min(ceiling(sqrt(nrow(Measures$measures))), 8) ## The maximum number of clusters to be investigated is set to 8, or to the square root of the number n of trajectories, if n < 50
      
      if(is.null(subset.n)){
        dat0 <- dat
      } else{
        dat0 <- dat[sample(seq_len(nrow(dat)), subset.n, replace = FALSE), ]
      }
      
      crit.list <- c("C_index", "Calinski_Harabasz", "Wemmert_Gancarski")
      
      ICV <- matrix(NA, nrow = length(crit.list), ncol = (k.max - 1))
      rownames(ICV) <- crit.list
      colnames(ICV) <- paste("k=", 2:k.max, sep = "")
      
      for(p in 2:k.max){
        ICV[, p - 1] <-  unlist(clusterCrit::intCriteria(as.matrix(dat0), part = spect(x = dat0, k = p, fuzzy = FALSE, nstart = nstart)$cluster, crit = crit.list))
      }
      
      ICV.raw <- ICV
      
      # Rescaling of the ICVs
      for(i in seq_len(nrow(ICV))){
        if(row.names(ICV)[i] %in% c("Calinski_Harabasz", "Wemmert_Gancarski")){
          v <- ICV[i, ] - min(ICV[i, ], na.rm = TRUE)
          v <- v/max(v, na.rm = TRUE)
          ICV[i, ] <- v
        }
        
        if(row.names(ICV)[i] %in% c("C_index")){
          v <- -ICV[i, ]
          v <- v - min(v, na.rm = TRUE)
          v <- v/max(v, na.rm = TRUE)
          ICV[i, ] <- v
        }
      }
      
      bins <- rep(0, ncol(ICV))
      names(bins) <- paste("k=", 2:k.max, sep = "")
      
      for(i in seq_len(nrow(ICV))){
        for(j in seq_len(ncol(ICV))){
          bins[order(ICV[i, ], decreasing = T)[j]] <- bins[order(ICV[i, ], decreasing = T)[j]] + 1 - (j - 1) / (ncol(ICV) - 1)
        }
      }
      
      nclusters <- order(bins, decreasing = TRUE)[1] + 1
    }
    
    # Perform spectral clustering
    c <- spect(
      x = dat,
      k = nclusters,
      nstart = nstart,
      fuzzy = fuzzy
    )
    
    partition <- c$cluster
    row.centers <- c$row.centers ## The rows of dat corresponding to the centroids
    fuzzy.partition <- c$fuzzy.partition
    
    # Re-label the groups from largest in size to smallest
    decr.order <- rev(order(summary(factor(partition))))
    
    w <- list()
    for (g in seq_len(nclusters)) {
      w[[g]] <- which(partition == g)
    }
    
    for (g in seq_len(nclusters)) {
      partition[w[[decr.order[g]]]] <- g
    }
    
    row.centers <- row.centers[decr.order] ## Reorder the centers to match the new labeling
    if(fuzzy == TRUE){
      fuzzy.partition <- fuzzy.partition[, decr.order]
      colnames(fuzzy.partition) <- as.character(seq_len(nclusters))
    }
    
    partition.summary <- summary(factor(partition))
    
    clust.by.id <- cbind(clust.by.id, partition)
    colnames(clust.by.id)[2] <- "Cluster"
    
    
    trajClusters <-
      structure(
        list(
          data = Measures$data,
          time = Measures$time,
          select = select,
          selection = selection,
          fuzzy = fuzzy,
          standardized.data = dat,
          nclusters.input = nclusters.input,
          raw.cluster.validity.indices = ICV.raw,
          cluster.validity.indices = ICV,
          ranked.voting.results = bins,
          nclusters = nclusters,
          partition = clust.by.id,
          partition.summary = partition.summary,
          fuzzy.partition = fuzzy.partition,
          ID.centers = ID[row.centers]
        ),
        class = "trajClusters"
      )
    
    return(trajClusters)
  }
#' @rdname trajClusters
#' @method print trajClusters
#' @export
print.trajClusters <- function(x, ...) {
  
  # If the 'nclusters' argument was unspecified in trajClusters(), display the cluster validity index (CVI) values by number of clusters and print a sentence saying what the ranked voting system determined is the optimal number of clusters prior to disclosing the cluster sizes
  if(is.null(x$nclusters.input)){
    print(round(x$raw.cluster.validity.indices, 3))
    
    cat("\n")  
    
    cat(paste("Using the combined input from the C-index, Calinski-Harabasz and Wemmert-Gancarski internal cluster validity indices, it was determined that an appropriate number of clusters for this data is ", x$nclusters, ". The clusters are labeled ", paste( names(x$partition.summary), collapse = ", ", sep = ""), " and are of respective size ", paste(x$partition.summary, collapse = ", ", sep = ""), ". The exact clustering is as follows.\n\n", sep = ""))
    
    print(x$partition, row.names = FALSE)
    
    cat("\n")
    
    cat("From here, you can use \n
        - CVIplot() for a graphical representation of the internal cluster validity indices used to determine the number of groups;\n
        - trajPartition() to extract the partition;\n
        - trajFuzzyPartition() to extract the fuzzy partition (if applicable);\n
        - trajClusterSummary() for the cluster-wise summaries of the measures;\n
        - trajAnalysis() for order analysis on the cluster-wise medians;\n
        - plot() to see the centroid trajectories and a sample from each groups;\n
        - trajScatter() for scatter plots of all the pairs of measures;\n
        - trajReduce() to investigate the possibility of reducing the number of measures used in the classification (optional).","\n")
  } else{
    
    cat(paste("The clusters are labeled ", paste( names(x$partition.summary), collapse = ", ", sep = ""), " and are of respective size ", paste(x$partition.summary, collapse = ", ", sep = ""), ". The exact clustering is as follows.\n\n", sep = ""))
    
    print(x$partition, row.names = FALSE)
    
    cat("\n")
    
    cat("From here, you can use \n
        - trajPartition() to extract the partition;\n
        - trajFuzzyPartition() to extract the fuzzy partition (if applicable);\n
        - trajClusterSummary() for the cluster-wise summaries of the measures;\n
        - trajAnalysis() for order analysis on the cluster-wise medians;\n
        - plot() to see the centroid trajectories and a sample from each groups;\n
        - trajScatter() for scatter plots of all the pairs of measures;\n
        - trajReduce() to investigate the possibility of reducing the number of measures used in the classification (optional).","\n")
  }
}
#' @rdname trajClusters
#' @method summary trajClusters
#' @export
summary.trajClusters <- function(object, ...) {
    
    # Construct a table 'clust.dist' containing the cluster frequencies, both absolute and relative
    clust.dist <- data.frame(matrix(nrow = 2, ncol = (object$nclusters + 1)))
    clust.dist[1,] <- signif(c(object$partition.summary, sum(object$partition.summary)))
    clust.dist[2,] <- signif(c(object$partition.summary / sum(object$partition.summary), sum(object$partition.summary) / sum(object$partition.summary)), 2)
    rownames(clust.dist) <- c("Absolute", "Relative")
    colnames(clust.dist) <- c(seq_len(object$nclusters), "Total")

    cat("\n")
    
    cat("Cluster frequencies:\n")
    print(clust.dist)
}
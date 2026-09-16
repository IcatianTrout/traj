#' Extract Cluster Assignments
#'
#' Extracts the cluster assignment of each trajectory from a
#' \code{trajClusters} object.
#'
#' @param object An object of class \code{trajClusters}.
#'
#' @return A data frame with two columns: \code{ID}, containing the
#'   trajectory identifiers, and \code{Cluster}, containing the cluster
#'   assignment of each trajectory.
#'
#' @examples
#' data(trajdata)
#' m <- trajMeasures(trajdata[, -2], ID = TRUE, measures = 1:20)
#' c <- trajClusters(m, nclusters = 4)
#' trajPartition(c)
#'
#' @export
trajPartition <- function(object) {
  if (!inherits(object, "trajClusters")) {
    stop("'object' must be an object of class 'trajClusters'.")
  }
  
  object$partition
}
#' 
#'
#' Extract Fuzzy Cluster Memberships
#'
#' Extracts the fuzzy cluster membership matrix from a
#' \code{trajClusters} object.
#'
#' @param object An object of class \code{trajClusters}.
#'
#' @return A matrix containing the membership value of each trajectory
#'   in each cluster. Each row corresponds to a trajectory and each
#'   column to a cluster. The values in each row represent the degree
#'   of membership of the trajectory in the corresponding clusters.
#'   Returns \code{NULL} if the clustering was performed with
#'   \code{fuzzy = FALSE}.
#'
#' @examples
#' data(trajdata)
#' m <- trajMeasures(trajdata[, -2], ID = TRUE, measures = 1:20)
#' c <- trajClusters(m, nclusters = 4, fuzzy = TRUE)
#' trajFuzzyPartition(c)
#'
#' @export
trajFuzzyPartition <- function(object) {
  if (!inherits(object, "trajClusters")) {
    stop("'object' must be an object of class 'trajClusters'.")
  }
  
  object$fuzzy.partition
}
#' 
#'
#' Computes the cluster-wise summaries of the measures
#'
#' Computes the cluster-wise summaries of the measures from a
#' \code{trajClusters} object.
#'
#' @param object An object of class \code{trajClusters}.
#'
#' @return A list where each entry contains the summaries of the measures for the trajectories in the corresponding cluster
#'
#' @examples
#' data(trajdata)
#' m <- trajMeasures(trajdata[, -2], ID = TRUE, measures = 1:20)
#' c <- trajClusters(m, nclusters = 4, fuzzy = TRUE)
#' trajClusterSummary(c)
#'
#' @export
trajClusterSummary <- function(object) {
  
  if (!inherits(object, "trajClusters")) {
    stop("'object' must be an object of class 'trajClusters'.")
  }
  
  # Construct cluster-specific summary tables of measures and store them in a list called 'groupwise.summaries'. 
  Q1 <- function(x) {
    return(quantile(x, probs = .25))
  }
  
  Q3 <- function(x) {
    return(quantile(x, probs = .75))
  }
  
  groupwise.summaries <- list()
  
  for (i in seq_len(object$nclusters)) {
    measures.summary <- data.frame(matrix(nrow = 6, ncol = ncol(object$selection) - 1))
    rownames(measures.summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    colnames(measures.summary) <- colnames(object$selection)[-1]
    
    which.id.i <- object$partition[which(object$partition[, 2] == i), 1]
    
    selection.cluster.i <- object$selection[which(object$selection$ID %in% which.id.i), ]
    
    measures.summary[1, ] <- apply(selection.cluster.i[, -1], 2, min)
    measures.summary[2, ] <- apply(selection.cluster.i[, -1], 2, Q1)
    measures.summary[3, ] <- apply(selection.cluster.i[, -1], 2, median)
    measures.summary[4, ] <- apply(selection.cluster.i[, -1], 2, mean)
    measures.summary[5, ] <- apply(selection.cluster.i[, -1], 2, Q3)
    measures.summary[6, ] <- apply(selection.cluster.i[, -1], 2, max)
    
    groupwise.summaries[[i]] <- measures.summary
  }
return(groupwise.summaries)
}
#' 
#' 
#' Perform order analysis on the cluster-wise medians
#'
#' Identifies which measures are most extreme for each cluster of a \code{trajClusters} object. 
#'
#' @param object An object of class \code{trajClusters}.
#' @param top_p The \code{top_p} most discriminating measures for each cluster to be reported in the summary.
#' 
#' @details For each cluster, the median of each standardized measure is computed and the difference (\code{delta}) between this median and the median of the medians is computed. A large positive or negative value of this difference is interpreted to mean that the measure is ``extreme'' for the cluster. The \code{top_p} most extreme measures for each cluster are reported under the column \code{measure}. The \code{delta} column logs the value of the difference. The number under the column \code{rank} is 1 if the absolute value of \code{delta} is the largest among all clusters, it is 2 if the absolute value of \code{delta} is the second largest among all clusters, etc. The \code{direction} column reads \code{largest} (resp. \code{smallest}) if the cluster's median is largest (resp. smallest) among all clusters. Otherwise, it reads \code{large} (resp. \code{small}) if \code{delta} is positive (resp. negative).
#'
#' @export
trajAnalysis <- function(object, top_p = 3) {
  
  # Perform various checks on the arguments
  if (!inherits(object, "trajClusters")) {
    stop("'object' must be an object of class 'trajClusters'.")
  }
  
  if(!is.numeric(top_p)) stop(paste("top_p must be an integer greater than 1", sep = ""))
  if(!(length(top_p) == 1)) stop(paste("top_p must be an integer greater than 1", sep = ""))
  if(!((top_p > 1) & (top_p %% 1 == 0))) stop(paste("top_p must be an integer greater than 1", sep = ""))
  
  cl.medians <- data.frame(matrix(NA, nrow = object$nclusters, ncol = length(object$select)))
  colnames(cl.medians) <- colnames(object$selection)[-1]
  
  for (i in seq_len(object$nclusters)) {
    which.i <- which(object$partition[, 2] == i)
    cl.medians[i, ] <- apply(object$standardized.data[which.i, ], 2, median)
  }
  
  # Initialize tables 'ranks', 'dirs' and 'deltas' of the same dimensions as cl.medians. Here, 'dirs' stand for directions; 'deltas' stands for the difference between a group's median and median of all the group medians; 'ranks' stands for where the group's median ranks from most extreme (largest absolute value of delta) to least extreme (smallest absolute value of delta)
  ranks <- cl.medians
  ranks[seq_len(nrow(ranks)), seq_len(ncol(ranks))] <- NA
  deltas <- dirs <- ranks
  
  for(j in seq_len(ncol(ranks))){
    
    median.j <- median(cl.medians[, j]) ## The median of the group medians
    
    deltas[, j] <- round(cl.medians[, j] - median.j, 4)
    abs.deltas <- abs(deltas[, j])
    n.unique <- length(unique(abs.deltas))
    
    # There might be multiple groups whose abs.deltas are the same. In this case, we give them all the same rank
    for(k in seq_len(n.unique)){
      w <- which(abs.deltas == unique(abs.deltas)[order(unique(abs.deltas), decreasing = TRUE)[k]])
      ranks[w, j] <- k
    }
    
    for(i in seq_len(nrow(ranks))){
      if((ranks[i, j] > 1) & (deltas[i, j] > 0)){dirs[i, j] <- "large"} ## large means delta > 0 but also rank > 1 so it's not the largest
      if((ranks[i, j] > 1) & (deltas[i, j] < 0)){dirs[i, j] <- "small"} ## small means delta < 0 but also rank > 1 so it's not the smallest
      if(cl.medians[i, j] == max(cl.medians[, j])){dirs[i, j] <- "largest"} 
      if(cl.medians[i, j] == min(cl.medians[, j])){dirs[i, j] <- "smallest"} 
      if(deltas[i, j] == 0){dirs[i, j] <- " "; ranks[i, j] <- 9999}} ## If delta = 0, there's no direction to speak of so we put " ". 
  }
  
  # Construct a table 'analysis' that reports, for each group, the top_p measures assuming the highest ranks, along with those measure's directions and deltas
  analysis <- data.frame(matrix(NA, ncol = 5, nrow = top_p * object$nclusters))
  colnames(analysis) <- c("cluster", "measure", "rank", "direction", "delta")
  
  measure.names <- colnames(ranks)
  for(m in seq_len(length(measure.names))){
    if(colnames(ranks[m]) == "m1"){ measure.names[m]  <- paste(colnames(ranks[m])," (max)", sep = "")}
    if(colnames(ranks[m]) == "m2"){ measure.names[m]  <- paste(colnames(ranks[m])," (min)", sep = "")}
    if(colnames(ranks[m]) == "m3"){ measure.names[m]  <- paste(colnames(ranks[m])," (range)", sep = "")}
    if(colnames(ranks[m]) == "m4"){ measure.names[m]  <- paste(colnames(ranks[m]),": mean)", sep = "")}
    if(colnames(ranks[m]) == "m5"){ measure.names[m]  <- paste(colnames(ranks[m])," (SD)", sep = "")}
    if(colnames(ranks[m]) == "m6"){ measure.names[m]  <- paste(colnames(ranks[m])," (slope)", sep = "")}
    if(colnames(ranks[m]) == "m7"){ measure.names[m]  <- paste(colnames(ranks[m])," (intercept)", sep = "")}
    if(colnames(ranks[m]) == "m8"){ measure.names[m]  <- paste(colnames(ranks[m])," (R^2)", sep = "")}
    if(colnames(ranks[m]) == "m9"){ measure.names[m]  <- paste(colnames(ranks[m])," (int. rate)", sep = "")}
    if(colnames(ranks[m]) == "m10"){ measure.names[m]  <- paste(colnames(ranks[m])," (net vari)", sep = "")}
    if(colnames(ranks[m]) == "m11"){ measure.names[m]  <- paste(colnames(ranks[m])," (contrast)", sep = "")}
    if(colnames(ranks[m]) == "m12"){ measure.names[m]  <- paste(colnames(ranks[m])," (tot vari)", sep = "")}
    if(colnames(ranks[m]) == "m13"){ measure.names[m]  <- paste(colnames(ranks[m])," (spikiness)", sep = "")}
    if(colnames(ranks[m]) == "m14"){ measure.names[m]  <- paste(colnames(ranks[m])," (max f')", sep = "")}
    if(colnames(ranks[m]) == "m15"){ measure.names[m]  <- paste(colnames(ranks[m])," (min f')", sep = "")}
    if(colnames(ranks[m]) == "m16"){ measure.names[m]  <- paste(colnames(ranks[m])," (SD f')", sep = "")}
    if(colnames(ranks[m]) == "m17"){ measure.names[m]  <- paste(colnames(ranks[m])," (f' net vari)", sep = "")}
    if(colnames(ranks[m]) == "m18"){ measure.names[m]  <- paste(colnames(ranks[m])," (max f'')", sep = "")}
    if(colnames(ranks[m]) == "m19"){ measure.names[m]  <- paste(colnames(ranks[m])," (min f'')", sep = "")}
    if(colnames(ranks[m]) == "m20"){ measure.names[m]  <- paste(colnames(ranks[m])," (SD f'')", sep = "")}
  }
  
  for(i in seq_len(object$nclusters)){
    # For a given cluster, order the measures by increasing value of rank and, among measures of a given rank, by decreasing absolute value of delta
    w <- order(unlist(ranks[i, ]), -unlist(abs(deltas[i, ])))[1:top_p]
    
    analysis[(i-1)*top_p + c(1:top_p), 1] <- i
    analysis[(i-1)*top_p + c(1:top_p), 2] <- measure.names[w]
    analysis[(i-1)*top_p + c(1:top_p), 3] <- unlist(ranks[i, w])
    analysis[(i-1)*top_p + c(1:top_p), 4] <- unlist(dirs[i, w])
    analysis[(i-1)*top_p + c(1:top_p), 5] <- unlist(deltas[i, w])
  }
  
  return(analysis)
}
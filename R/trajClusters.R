#'@title Classify the Longitudinal Data Based on the Measures.
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
#'@param top_p The \code{top_p} most discriminating measures for each cluster to be reported in the summary.
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
#' \dontrun{
#'data("trajdata")
#'trajdata.noGrp <- trajdata[, -which(colnames(trajdata) == "Group")] # remove the Group column
#'
#'m = trajMeasures(trajdata.noGrp, ID = TRUE, measures = 1:20)
#'
#'s2.3 <- trajClusters(m, nclusters = 3)
#'plot(s2.3)
#'
#'s2.4 <- trajClusters(m, nclusters = 4)
#'plot(s2.4)
#'
#'s2.5 <- trajClusters(m, nclusters = 5)
#'plot(s2.5)
#'
#'groups <- trajClusters(m, nclusters = 4)$partition
#'}
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
      
      for(i in 1:nrow(ICV)){
        for(j in 1:ncol(ICV)){
          bins[order(ICV[i, ], decreasing = T)[j]] <- bins[order(ICV[i, ], decreasing = T)[j]] + 1 - (j - 1)/(ncol(ICV) - 1)
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
    partition <- decr.order[partition]
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
    
    cat("From here, use the plot() function to see the centroid trajectories and a sample from each groups. Use CVIplot() for a graphical representation of the internal cluster validity indices used to determine the number of groups. For a better understanding of how the measures were used to discriminate amongst the groups, use scatterplot() for scatter plots of all the pairs of measures. To investigate the possibility of reducing the number of measures used in the classification (optional), use trajReduce().")
  } else{
    
    cat(paste("The clusters are labeled ", paste( names(x$partition.summary), collapse = ", ", sep = ""), " and are of respective size ", paste(x$partition.summary, collapse = ", ", sep = ""), ". The exact clustering is as follows.\n\n", sep = ""))
    
    print(x$partition, row.names = FALSE)
    
    cat("\n")
    
    cat("From here, use the plot() function to see the centroid trajectories and a sample from each groups. For a better understanding of how the measures were used to discriminate amongst the groups, use scatterplot() for scatter plots of all the pairs of measures. To investigate the possibility of reducing the number of measures used in the classification (optional), use trajReduce().")
  }
}
#' @rdname trajClusters
#' @method summary trajClusters
#' @export
summary.trajClusters <- function(object, top_p = 3, ...) {
    
    # Perform various checks on the top_p argument
    if(!is.numeric(top_p)) stop(paste("top_p must be an integer greater than 1", sep = ""))
    if(!(length(top_p) == 1)) stop(paste("top_p must be an integer greater than 1", sep = ""))
    if(!((top_p > 1) & (top_p %% 1 == 0))) stop(paste("top_p must be an integer greater than 1", sep = ""))
    
    # Construct a table 'clust.dist' containing the cluster frequencies, both absolute and relative
    clust.dist <- data.frame(matrix(nrow = 2, ncol = (object$nclusters + 1)))
    clust.dist[1,] <- signif(c(object$partition.summary, sum(object$partition.summary)))
    clust.dist[2,] <- signif(c(object$partition.summary / sum(object$partition.summary), sum(object$partition.summary) / sum(object$partition.summary)), 2)
    rownames(clust.dist) <- c("Absolute", "Relative")
    colnames(clust.dist) <- c(1:object$nclusters, "Total")

    # Construct cluster-specific summary tables of measures and store them in a list called 'groupwise.summaries'. Separately, construct a table 'cl.medians' of cluster-wise medians for each measures.
    Q1 <- function(x) {
      return(quantile(x, probs = .25))
    }

    Q3 <- function(x) {
      return(quantile(x, probs = .75))
    }
    
    cl.medians <- data.frame(matrix(NA, nrow = object$nclusters, ncol = length(object$select)))
    colnames(cl.medians) <- colnames(object$selection)[-1]
    
    groupwise.summaries <- list()
    group.sizes <- c()
    
    for (i in seq_len(object$nclusters)) {
      measures.summary <- data.frame(matrix(nrow = 6, ncol = ncol(object$selection) - 1))
      rownames(measures.summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
      colnames(measures.summary) <- colnames(object$selection)[-1]
      
      which.id.i <- object$partition[which(object$partition[, 2] == i), 1]
      which.i <- which(object$partition[, 2] == i)
      
      selection.cluster.i <- object$selection[which(object$selection$ID %in% which.id.i), ]
      
      measures.summary[1, ] <- apply(selection.cluster.i[, -1], 2, min)
      measures.summary[2, ] <- apply(selection.cluster.i[, -1], 2, Q1)
      measures.summary[3, ] <- apply(selection.cluster.i[, -1], 2, median)
      measures.summary[4, ] <- apply(selection.cluster.i[, -1], 2, mean)
      measures.summary[5, ] <- apply(selection.cluster.i[, -1], 2, Q3)
      measures.summary[6, ] <- apply(selection.cluster.i[, -1], 2, max)
      
      cl.medians[i, ] <- apply(object$standardized.data[which.i, ], 2, median)
      
      groupwise.summaries[[i]] <- measures.summary
      group.sizes[i] <- object$partition.summary[i]
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
    
    structure(
      list(
        clust.dist = clust.dist,
        groupwise.summaries = groupwise.summaries,
        group.sizes = group.sizes,
        analysis = analysis
      ),
      class = "summary.trajClusters"
    )
}
#' @rdname trajMeasures
#' @method print summary.trajClusters
#' @export
print.summary.trajClusters <- function(x, ...) {
  
  cat("Description of the measures:\n")
  cat("m1: Maximum\n")
  cat("m2: Minimum\n")
  cat("m3: Range\n")
  cat("m4: Mean\n")
  cat("m5: Standard deviation\n")
  cat("m6: Slope of the affine approximation\n")
  cat("m7: Intercept of the affine approximation\n")
  cat("m8: Proportion of variance explained by the affine approximation\n")
  cat("m9: Rate of intersection with the best affine approximation\n")
  cat("m10: Net variation per unit of time\n")
  cat("m11: Late variation to early variation contrast\n")
  cat("m12: Total variation per unit of time\n")
  cat("m13: Spikiness\n")
  cat("m14: Maximum of the first derivative\n")
  cat("m15: Minimum of the first derivative\n")
  cat("m16: Standard deviation of the first derivative\n")
  cat("m17: First derivative's net variation per unit of time\n")
  cat("m18: Maximum  of the second derivative\n")
  cat("m19: Minimum  of the second derivative\n")
  cat("m20: Standard deviation of the second derivative\n")
  
  cat("\n")
  
  cat("Cluster frequencies:\n")
  print(x$clust.dist)
  
  cat("\n")
  cat("Summary of selected measures by cluster:\n")
  
  for(i in seq_len(length(x$groupwise.summaries))){
    cat(paste("Cluster ", i, " (size ", x$group.sizes[i], "):", sep = ""))
    cat("\n")
    print(x$groupwise.summaries[[i]])
    cat("\n")
  }
  
  cat("\n")
  print(x$analysis, rownames = FALSE)
  cat("\n")
}

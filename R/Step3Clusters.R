#'@title Classify the Longitudinal Data Based on the Selected Measures.
#'
#'@description Classifies the trajectories by applying the K-mean clustering
#'  algorithm to the measures selected by \code{Step2Selection}.
#'
#'@param trajSelection object of class \code{trajSelection} as returned by
#'  \code{Step2Selection}.
#'@param nclusters either NULL or the desired number of clusters. If NULL, the
#'  number of clustersis determined using the GAP criterion as implemented in
#'  the \code{\link[cluster]{clusGap}} function.
#'@param nstart to be passed to the \code{nstart} argument of
#'  \code{\link[stats]{kmeans}}.
#'@param iter.max to be passed to the \code{iter.max} argument of
#'  \code{\link[stats]{kmeans}}.
#'@param K.max to be passed to the \code{K.max} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param B to be passed to the \code{SE.factor} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param d.power to be passed to the \code{B} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param spaceH0 to be passed to the \code{spaceH0} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param method to be passed to the \code{method} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param SE.factor to be passed to the \code{SE.factor} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param x object of class trajClusters
#'@param object object of class trajClusters
#'@param ... further arguments passed to or from other methods.
#'
#'@return An object of class \code{trajClusters}; a list containing the result
#'  of the clustering, the output of the \code{clusGap} function, as well as a
#'  curated form of the arguments.
#'
#'@import cluster
#'@importFrom stats kmeans na.omit qt quantile
#'
#' @examples
#' \dontrun{
#'data("trajdata")
#'trajdata.noGrp <- trajdata[, which(colnames(trajdata) == "Group")] #remove the Group column
#'
#'m = Step1Measures(trajdata.noGrp, ID = TRUE, measures = 1:18)
#'s = Step2Selection(m)
#'
#'s$RC$loadings
#'
#'s2 = Step2Selection(m, select = c(10, 12, 8, 4))
#'
#'c3.part <- Step3Clusters(s2, nclusters = 3)$partition
#'c4.part <- Step3Clusters(s2, nclusters = 4)$partition
#'c5.part <- Step3Clusters(s2, nclusters = 5)$partition
#'
#'}
#'
#'@seealso \code{\link[traj]{Step2Selection}}
#'
#'@rdname Step3Clusters
#'
#'@export
Step3Clusters <-
  function (trajSelection,
            nclusters = NULL,
            nstart = 200,
            iter.max = 100,
            K.max = 8,
            B = 500,
            d.power = 2,
            spaceH0 = "scaledPCA",
            method = "Tibs2001SEmax",
            SE.factor = 1) {
    nclusters.input <- nclusters
    
    ID <- trajSelection$selection[, 1]
    
    #standardize the measures to be clustered:
    
    data <-
      data.frame(apply(data.frame(trajSelection$selection[,-1]), 2, scale))
    
    if (!is.null(nclusters) && (nclusters > nrow(data))) {
      stop(
        "The number 'nclusters' of requested clusters cannot exceed the number of trajectories."
      )
    }
    if (is.null(nclusters)) {
      ## The clusGap function takes as argument a function that will perform
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
      GAP <-
        cluster::clusGap(
          x = data,
          FUNcluster = kmeans.nstart,
          K.max = K.max,
          B = B,
          d.power = d.power,
          spaceH0 = spaceH0
        )
      nclusters <-
        maxSE(
          f = GAP$Tab[, "gap"],
          SE.f = GAP$Tab[, "SE.sim"],
          method = method,
          SE.factor = SE.factor
        )
      cluster.est2 <-
        kmeans(
          x = data,
          centers = nclusters,
          iter.max = iter.max,
          nstart = nstart
        )
      partition <- cluster.est2$cluster
    } else {
      GAP <- NULL
      cluster.est2 <-
        kmeans(
          x = data,
          centers = nclusters,
          iter.max = iter.max,
          nstart = nstart
        )
      partition <- cluster.est2$cluster
    }
    
    #relabel the groups from largest to smallest
    
    decr.order <- rev(order(summary(factor(partition))))
    
    w <- list()
    for (g in seq_len(nclusters)) {
      w[[g]] <- which(partition == g)
    }
    
    for (g in seq_len(nclusters)) {
      partition[w[[g]]] <- decr.order[g]
    }
    
    partition.summary <- summary(factor(partition))
    
    clust.by.id <-
      as.data.frame(cbind(trajSelection$data[, 1], partition))
    colnames(clust.by.id) <- c("ID", "Cluster")
    
    trajClusters <-
      structure(
        list(
          data = trajSelection$data,
          time = trajSelection$time,
          selection = trajSelection$selection,
          GAP = GAP,
          nclusters = nclusters,
          partition = clust.by.id,
          partition.summary = partition.summary
        ),
        class = "trajClusters"
      )
    
    return(trajClusters)
    
  }

#'@rdname Step3Clusters
#'
#'@export
print.trajClusters <- function(x, ...) {
  cat(
    paste(
      "The trajectories were grouped in ",
      x$nclusters,
      " clusters labeled ",
      paste(
        names(x$partition.summary),
        collapse = ", ",
        sep = ""
      ),
      " of respective size ",
      paste(
        x$partition.summary,
        collapse = ", ",
        sep = ""
      ),
      ". The exact clustering is as follows.\n\n",
      sep = ""
    )
  )
  
  print(x$partition, row.names = FALSE)
}

#'@rdname Step3Clusters
#'
#'@export
summary.trajClusters <- function(object, ...) {
  if (!is.null(object$GAP)) {
    gap.stat <-
      cbind(seq_len(nrow(object$GAP$Tab)), object$GAP$Tab[, 3:4])
    colnames(gap.stat) <- c("K", "GAP(K)", "SE")
    cat("GAP statistic as a function of the number of clusters:\n")
    print(as.data.frame(gap.stat), row.names = FALSE)
  }
  
  cat("\n")
  
  cat("Cluster frequencies:\n")
  clust.dist <-
    data.frame(matrix(nrow = 2, ncol = (object$nclusters + 1)))
  clust.dist[1,] <-
    signif(c(
      object$partition.summary,
      sum(object$partition.summary)
    ))
  clust.dist[2,] <-
    signif(c(
      object$partition.summary / sum(object$partition.summary),
      sum(object$partition.summary) / sum(object$partition.summary)
    ), 2)
  rownames(clust.dist) <- c("Absolute", "Relative")
  colnames(clust.dist) <- c(1:object$nclusters, "Total")
  print(clust.dist)
  
  cat("\n")
  cat("Summary of selected measures by cluster:\n")
  
  Q1 <- function(x) {
    return(quantile(x, probs = .25))
  }
  
  Q2 <- function(x) {
    return(quantile(x, probs = .5))
  }
  
  Q3 <- function(x) {
    return(quantile(x, probs = .75))
  }
  
  for (i in 1:object$nclusters) {
    measures.summary <-
      data.frame(matrix(
        nrow = 6,
        ncol = ncol(object$selection) - 1
      ))
    rownames(measures.summary) <-
      c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    colnames(measures.summary) <-
      colnames(object$selection)[-1]
    
    which.i <- which(object$partition[, 2] == i)
    
    selection.cluster.i <- object$selection[which.i,]
    
    measures.summary[1,] <- apply(selection.cluster.i, 2, min)[-1]
    measures.summary[2,] <- apply(selection.cluster.i, 2, Q1)[-1]
    measures.summary[3,] <- apply(selection.cluster.i, 2, Q2)[-1]
    measures.summary[4,] <- apply(selection.cluster.i, 2, mean)[-1]
    measures.summary[5,] <- apply(selection.cluster.i, 2, Q3)[-1]
    measures.summary[6,] <- apply(selection.cluster.i, 2, max)[-1]
    
    cat(paste(
      "Cluster ",
      i,
      " (size ",
      object$partition.summary[i],
      "):",
      sep = ""
    ))
    cat("\n")
    print(measures.summary)
    cat("\n")
  }
}

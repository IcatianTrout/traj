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
#'@param nstart The number of random starts. Defaults to \code{50}.
#'@param x object of class \code{trajClusters}.
#'@param object object of class \code{trajClusters}.
#'@param ... further arguments passed to or from other methods.
#'
#'@details The spectral clustering algorithm presented in Meila (2005) is implemented in which the similarity matrix \eqn{S} is built from a binary K nearest neighbors similarity function (\eqn{S=(W+W^T)/2}, where \eqn{W_{ij}=1} if data point \eqn{j} is among the nearest points to data point \eqn{i} and \eqn{W_{ij}=0} otherwise). 
#'
#'@return An object of class \code{trajClusters}; a list containing the result of the clustering, as well as a curated form of the arguments. If \code{nclusters} is set to \code{NULL}, clustering is carried out for each number \eqn{k} of clusters between 2 and (up to) 8 and a plot is produced representing the value of three internal cluster validity indices (C-index, Calinski-Harabasz, Wemmert-Gancarski) as a function of \eqn{k}. As in the 'KmL' package of Genolini et al., these validity indices are presented on a scale from 0 to 1, with 1 corresponding to the highest validity score and 0 corresponding to the lowest. From this, a "best" value of \eqn{k} is determined using a ranked voting system.
#'
#'@import cluster fclust clusterCrit 
#'@importFrom stats quantile kmeans
#'@importFrom e1071 cmeans
#'
#'@references 
#'Genolini, C. et al., kml: K-Means for Longitudinal Data, https://CRAN.R-project.org/package=kml
#'
#'Meila, M., Spectral Clustering. Handbook of Cluster Analysis, Chapter 7, Chapman and Hall/CRC, 2005.
#'
#' @examples
#' \dontrun{
#'data("trajdata")
#'trajdata.noGrp <- trajdata[, -which(colnames(trajdata) == "Group")] # remove the Group column
#'
#'m = trajMeasures(trajdata.noGrp, ID = TRUE, measures = 1:19)
#'
#'s2.3 <- trajClusters(m, nclusters = 3)
#'plot(s2.3)
#'
#'#'s2.4 <- trajClusters(m, nclusters = 4)
#'plot(s2.4)
#'
#'#'s2.5 <- trajClusters(m, nclusters = 5)
#'plot(s2.5)
#'
#'groups <- s2.4 <- trajClusters(m, nclusters = 4)$partition
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
            nstart = 50
  ) {
    
    if ((!is.null(select)) & ((!is.numeric(select)) | (!is.vector(select)))) {
      stop("Argument 'select' must be either NULL or a numerical vector.")
    }
    
    if (!(fuzzy %in% c("TRUE", "FALSE"))) {
      stop("'fuzzy' should be either 'TRUE' or 'FALSE'.")
    }
    
    k.max <- min(ceiling(sqrt(nrow(Measures$measures))), 8)
    ID <- Measures$measures[, 1]
    nclusters.input <- nclusters
    partition.summary <- NULL
    clust.by.id <- Measures$data[, 1, drop = FALSE]
    rownames(clust.by.id) <- NULL
    ICV.raw <- NULL
    ICV <- NULL
    bins <- NULL
    
    #standardize the measures to be clustered:
    dat <-
      data.frame(apply(data.frame(Measures$measures[,-1]), 2, scale))
    
    if(is.null(select)) {
      select <- Measures$measures.arg
    }
    if (!is.null(select)) {
      m.select <- paste("m", select, sep = "")
      if (FALSE %in% (m.select %in% colnames(dat))) {
        stop("The 'select' argument must only contain measures included in Measures.")
      } else {
        selection <- Measures$measures[, c("ID", m.select)]
        dat <- dat[, m.select, drop = FALSE]
      }
    }
    
    if (!is.null(nclusters) && (nclusters > nrow(dat))) {
      stop(
        "The number 'nclusters' of requested clusters cannot exceed the number of trajectories."
      )
    }
    
    w <- which(is.na(apply(dat, 2, sd)))
    if(length(w) == 1){
      meas.rmv <- colnames(dat)[w]
      warning(paste("Being constant, measure ", noquote(paste(meas.rmv, collapse = ", ")), " has been removed.", sep = ""))
      dat <- dat[, -w, drop = FALSE]
    }
    if(length(w) > 1){
      meas.rmv <- colnames(dat)[w]
      warning(paste("Being constant, measures ", noquote(paste(meas.rmv, collapse = ", ")), " have been removed.", sep = ""))
      dat <- dat[, -w, drop = FALSE]
    }
    
    if (is.null(nclusters)) {
      
      crit.list <- c("C_index", "Calinski_Harabasz", "Wemmert_Gancarski")
      
      ICV <- matrix(NA, nrow = length(crit.list), ncol = (k.max - 1))
      rownames(ICV) <- crit.list
      colnames(ICV) <- paste("k=",2:k.max,sep="")
      
      for(p in 2:k.max){
        ICV[, p-1] <-  unlist(clusterCrit::intCriteria(as.matrix(dat), part = spect(x = dat, k = p, fuzzy = FALSE, nstart = nstart)$cluster, crit = crit.list))
      }
      
      ICV.raw <- ICV
      
      # Rescaling 
      for(i in 1:nrow(ICV)){
        if(row.names(ICV)[i] %in% c("Calinski_Harabasz", "Wemmert_Gancarski")){
          v <- ICV[i, ] - min(ICV[i, ], na.rm = T)
          v <- v/max(v, na.rm = T)
          ICV[i,] <- v
        }
        
        if( row.names(ICV)[i] %in% c("C_index")){
          v <- -ICV[i, ]
          v <- v - min(v, na.rm = T)
          v <- v/max(v, na.rm = T)
          ICV[i,] <- v
        }
      }
      
      bins <- rep(0, ncol(ICV))
      names(bins) <- paste("k=", 2:k.max, sep = "")
      
      for(i in 1:nrow(ICV)){
        for(j in 1:ncol(ICV)){
          bins[order(ICV[i, ], decreasing = T)[j]] <- bins[order(ICV[i, ], decreasing = T)[j]] + 1 - (j - 1)/(ncol(ICV) - 1)
        }
      }
      
      nclusters <- order(bins, decreasing=T)[1] + 1
    }
    
    c <- spect(
      x = dat,
      k = nclusters,
      nstart = nstart,
      fuzzy = fuzzy
    )
    
    partition <- c$cluster
    row.centers <- c$row.centers
    fuzzy.partition <- c$fuzzy.partition
    
    #re-label the groups from largest to smallest
    
    decr.order <- rev(order(summary(factor(partition))))
    row.centers <- row.centers[decr.order]
    fuzzy.partition <- fuzzy.partition[, decr.order]
    
    w <- list()
    for (g in seq_len(nclusters)) {
      w[[g]] <- which(partition == g)
    }
    
    for (g in seq_len(nclusters)) {
      partition[w[[decr.order[g]]]] <- g
    }
    
    partition.summary <- summary(factor(partition))
    
    clust.by.id <-
      cbind(clust.by.id, partition)
    colnames(clust.by.id)[2] <- "Cluster"
    
    
    trajClusters <-
      structure(
        list(
          data = Measures$data,
          time = Measures$time,
          select = select,
          selection = selection,
          fuzzy = fuzzy,
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
#'@rdname trajClusters
#'
#'@export
print.trajClusters <- function(x, ...) {
  if(is.null(x$nclusters.input)){
    print(round(x$raw.cluster.validity.indices, 3))
    
    cat("\n")  
    
    cat(
      paste(
        "Using the combined input from the C-index, Calinski-Harabasz and Wemmert-Gancarski internal cluster validity indices, it was determined that an appropriate number of clusters for this data is ",
        x$nclusters,
        ". The clusters are labeled ",
        paste(
          names(x$partition.summary),
          collapse = ", ",
          sep = ""
        ),
        " and are of respective size ",
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
    
    cat("\n")
    
    cat("From here, use the plot() function to see the centroid trajectories and a sample from each groups. Use CVIplot() for a graphical representation of the internal cluster validity indices used to determine the number of groups. For a better understanding of how the measures were used to discriminate amongst the groups, use scatterplot() for scatter plots of all the pairs of measures. To investigate the possibility of reducing the number of measures used in the classification (optional), use trajReduce().")
  } else{
    cat(
      paste("The clusters are labeled ",
        paste(
          names(x$partition.summary),
          collapse = ", ",
          sep = ""
        ),
        " and are of respective size ",
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
    
    cat("\n")
    
    cat("From here, use the plot() function to see the centroid trajectories and a sample from each groups. For a better understanding of how the measures were used to discriminate amongst the groups, use scatterplot() for scatter plots of all the pairs of measures. To investigate the possibility of reducing the number of measures used in the classification (optional), use trajReduce().")
  }
}
#'@rdname trajClusters
#'
#'@export
summary.trajClusters <- function(object, ...) {
  if(!is.null(object$nclusters)){
    
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
}


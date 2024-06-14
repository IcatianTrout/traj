#'@title Classify the Longitudinal Data Based on the Selected Measures.
#'
#'@description Classifies the trajectories by applying the k-means clustering
#'  algorithm to the measures selected by \code{Step2Selection}.
#'
#'@param trajSelection object of class \code{trajSelection} as returned by
#'  \code{Step2Selection}.
#'@param algorithm either \code{"k-medoids"} or \code{"k-means"}. Determines the clustering algorithm to use. Defaults to \code{"k-medoids"}.
#'@param metric to be passed to the \code{metric} argument of
#'  \code{\link[cluster]{pam}} if \code{"k-medoids"} is the chosen \code{algorithm}. Defaults to \code{"euclidean"}.
#'@param nstart to be passed to the \code{nstart} argument of
#'  \code{\link[stats]{kmeans}} if \code{"k-means"} is the chosen \code{algorithm}. Defaults to \code{200}.
#'@param iter.max to be passed to the \code{iter.max} argument of
#'  \code{\link[stats]{kmeans}} if \code{"k-means"} is the chosen \code{algorithm}. Defaults to \code{100}.
#'@param nclusters either \code{NULL} or the desired number of clusters. If \code{NULL}, the
#'  number of clusters is determined using the criterion chosen in \code{criterion}. Defaults to \code{NULL}.
#'@param criterion criterion to determine the optimal number of clusters if \code{nclusters} is \code{NULL}. Either \code{"GAP"} or \code{"Calinski-Harabasz"}. Defaults to \code{"Calinski-Harabasz"}.
#'@param K.max maximum number of clusters to be considered if \code{nclusters} is set to \code{NULL}. Defaults to \code{15}.
#'@param boot logical. If \code{TRUE}, and if \code{"Calinski-Harabasz"} is the chosen \code{criterion}, the optimal number of clusters will be the first mode of sampling distribution of the optimal number of clusters obtained by bootstrap. Defaults to \code{TRUE}.
#'@param R the number of bootstrap replicate if \code{boot} is set to \code{TRUE}. Defaults to \code{100}.
#'@param B to be passed to the \code{B} argument of
#'  \code{\link[cluster]{clusGap}} if \code{"GAP"} is the chosen \code{criterion}.
#'@param x object of class \code{trajClusters}.
#'@param object object of class \code{trajClusters}.
#'@param ... further arguments passed to or from other methods.
#'
#'@return An object of class \code{trajClusters}; a list containing the result
#'  of the clustering, as well as a curated form of the arguments.
#'
#'@details If \code{"GAP"} is the chosen \code{criterion} for determining the optimal number of clusters, the method described by Tibshirani et al. is implemented by the \code{\link[cluster]{clusGap}} function.
#'
#'Instead, if \code{"Calinski-Harabasz"} is the chosen \code{criterion}, the Calinski-Harabasz index is computed for each possible number of clusters between 2 and \code{K.max} and the optimal number of clusters is the maximizer of the Calinski-Harabasz index. Moreover, if \code{boot} is set to \code{TRUE}, then, following the guidelines suggested by Mesidor et al., a sampling distribution of the optimal number of clusters is obtained by bootstrap and the optimal number of clusters is chosen to be the (first) mode of this sampling distribution. 
#'
#'@import cluster
#'@importFrom stats kmeans na.omit qt quantile var
#'@importFrom graphics barplot
#'
#'@references Miceline Mésidor, Caroline Sirois, Marc Simard, Denis Talbot, A Bootstrap Approach for Evaluating Uncertainty in the Number of Groups Identified by Latent Class Growth Models, American Journal of Epidemiology, Volume 192, Issue 11, November 2023, Pages 1896–1903, https://doi.org/10.1093/aje/kwad148
#'
#'Tibshirani, R., Walther, G. and Hastie, T. (2001). Estimating the number of data clusters via the Gap statistic. Journal of the Royal Statistical Society B, 63, 411–423.
#'
#'Tibshirani, R., Walther, G. and Hastie, T. (2000). Estimating the number of clusters in a dataset via the Gap statistic. Technical Report. Stanford.
#'
#' @examples
#' \dontrun{
#'data("trajdata")
#'trajdata.noGrp <- trajdata[, -which(colnames(trajdata) == "Group")] #remove the Group column
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
            algorithm = "k-medoids",
            metric = "euclidean",
            nstart = 200,
            iter.max = 100,
            nclusters = NULL,
            criterion = "Calinski-Harabasz",
            K.max = 15,
            boot = TRUE,
            R = 100,
            B = 500
  ) {
    
    if (is.null(nclusters) & !(criterion %in% c("Calinski-Harabasz", "GAP"))) { 
      stop("'criterion' should be one of 'Calinski-Harabasz' or 'GAP'.")
    }
    
    if (!(algorithm %in% c("k-medoids", "k-means"))) {
      stop("'algorithm' should be either 'k-medoids' or 'k-means'.")
    }
    
    GAP <- NULL
    CH <- NULL
    CH.boot <- NULL
    
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
      
      if (criterion == "GAP") {
        if (algorithm == "k-means") {
          
          kmeans.nstart <- function (x, k) {
            return(kmeans(
              x = x,
              centers = k,
              nstart = nstart,
              iter.max = iter.max
            ))
          }
          
          FUNcluster <- kmeans.nstart
        } else {
          
          pam.nstart <- function (x, k) {
            return(pam(
              x = x,
              k = k,
              diss = FALSE,
              metric = metric
            ))
          }
          
          FUNcluster <- pam.nstart
        }
        
        GAP <-
          cluster::clusGap(
            x = data,
            FUNcluster = FUNcluster,
            K.max = K.max,
            B = B,
            d.power = 2,
            spaceH0 = "scaledPCA"
          )
        
        nclusters <-
          cluster::maxSE(
            f = GAP$Tab[, "gap"],
            SE.f = GAP$Tab[, "SE.sim"],
            method = "Tibs2001SEmax",
            SE.factor = 1
          )
      }
      
      if (criterion == "Calinski-Harabasz") {
        if (isTRUE(boot)) {
          CH.boot <- c()
          for(i in 1:R){
            s1 = trajSelection$trajMeasures
            id = s1$measures$ID 
            indices <- sample(nrow(s1$measures), replace=TRUE)
            s1$measures <- s1$measures[indices, ]
            s1$measures$ID <- id
            s2.boot <- Step2Selection(trajMeasures = s1, num.select = s1$input$num.select, discard = s1$input$discard, select = s1$input$select)
            d = data.frame(apply(data.frame(s2.boot$selection[,-c(1), drop = FALSE]), 2, scale))
            CH.aux <- c()
            for(k in 2:K.max){
              if(algorithm == "k-medoids"){
                CH.aux[k] <- CalinskiHarabasz(x = d, clustering = cluster::pam(x = d, k = k, cluster.only = TRUE))
              }
              if(algorithm == "k-means"){
                CH.aux[k] <- CalinskiHarabasz(x = d, clustering = stats::kmeans(x = d, centers = k, iter.max = iter.max, nstart = nstart)$cluster)
              }
            }
            CH.boot[i] <- which(CH.aux == max(CH.aux, na.rm = T))
          }
          
          nclusters <- FirstMode(CH.boot)
          
        } else{
          CH <- c()
          if (algorithm == "k-medoids") {
            for(k in 2:K.max){
              CH[k] <- CalinskiHarabasz(x = data, clustering = cluster::pam(x = data, k = k, cluster.only = TRUE))
            }
            nclusters <- which(CH == max(CH, na.rm = T))
          }
          
          if (algorithm == "k-means") {
            for (k in 2:K.max) {
              CH[k] <- CalinskiHarabasz(x = data, clustering = stats::kmeans(x = data, centers = k, iter.max = iter.max, nstart = nstart)$cluster)
            }
            nclusters <- which(CH == max(CH, na.rm = T))
          }
        }
      }
    } 
    
    if (algorithm == "k-means"){
      partition <- 
        stats::kmeans(
          x = data,
          centers = nclusters,
          iter.max = iter.max,
          nstart = nstart
        )$cluster
    }
    
    if(algorithm == "k-medoids"){
      partition <- 
        cluster::pam(
          x = data,
          k = nclusters,
          diss = FALSE,
          metric = metric
        )$clustering
    }
    
    #re-label the groups from largest to smallest
    
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
    
    if(isFALSE(trajSelection$cap.outliers) & algorithm == "k-means"){
      warning("You have selected k-means as the clustering algorithm in Step3Clusters but you also have chosen to not cap the outliers in Step1Measures. As a result, it is possible that the clustering has been strongly driven by outliers. To prevent this, it is recommended to either set 'cap.outliers' to TRUE in Step1Measures or to set 'algorithm' to 'k-medoids' in Step3Clusters.")
    }
    
    trajClusters <-
      structure(
        list(
          data = trajSelection$data,
          time = trajSelection$time,
          selection = trajSelection$selection,
          GAP = GAP,
          CH = CH,
          CH.boot = CH.boot,
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
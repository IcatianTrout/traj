#' Extract Cluster Assignments
#'
#' Extracts the cluster assignment of each trajectory from a
#' \code{trajReduce} object.
#'
#' @param object An object of class \code{trajClusters}.
#'
#' @return A data frame with two columns: \code{ID}, containing the
#'   trajectory identifiers, and \code{Cluster}, containing the cluster
#'   assignment of each trajectory.
#'
#' @export
trajReducePartition <- function(object) {
  if (!inherits(object, "trajReduce")) {
    stop("'object' must be an object of class 'trajReduce'.")
  }
  
  object$partition.red
}
#' Extract Fuzzy Cluster Memberships
#'
#' Extracts the fuzzy cluster membership matrix from a
#' \code{trajReduce} object.
#'
#' @param object An object of class \code{trajReduce}.
#'
#' @export
trajFuzzyReducedPartition <- function(object) {
  if (!inherits(object, "trajReduce")) {
    stop("'object' must be an object of class 'trajReduce'.")
  }
  
  object$fuzzy.partition.red
}
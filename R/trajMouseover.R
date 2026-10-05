#' Interactive Plot of Trajectories by Cluster
#'
#' Creates an interactive \code{highcharter} plot displaying 
#' trajectories colored according to their cluster membership. 
#'
#' The plot provides interactive highlighting: when the user hovers over a
#' trajectory, only that trajectory is highlighted while the other trajectories
#' become less prominent. Clicking
#' on a cluster in the legend hides or displays all trajectories belonging to
#' that cluster. Zoom vertically by clicking and dragging.
#'
#' @param x object of class \code{trajClusters} as returned by the function
#'  \code{trajClusters()}.
#'
#' @return An interactive \code{highchart} object.
#'
#'
#' @examples
#'data(trajdata)
#'
#'m = trajMeasures(trajdata[, -2], ID = TRUE)
#'
#'c3 = trajClusters(m, nclusters = 3)
#'
#'trajMouseover(c3)
#'
#' @importFrom highcharter highchart hc_chart hc_title hc_xAxis hc_yAxis
#'   hc_tooltip hc_plotOptions hc_adataM_series JS
#' @export

trajMouseover <- function(x) {
  
  dataM <- merge(x$partition, x$data, by = "ID", all.x = TRUE)
  colnames(dataM)[-c(1:2)] <- seq_along(colnames(dataM)[-c(1:2)])
  
  # Check that highcharter is available
  if (!requireNamespace("highcharter", quietly = TRUE)) {
    stop("Package 'highcharter' is required.")
  }
  
  # Basic checks
  if (ncol(dataM) < 3) {
    stop("The data frame must contain at least 3 columns.")
  }
  
  if (!all(names(dataM)[1:2] == c("ID", "Cluster"))) {
    stop("The first two columns must be named 'ID' and 'Cluster'.")
  }
  
  # Time points
  time_columns <- names(dataM)[-(1:2)]
  time <- as.numeric(time_columns)
  
  if (anyNA(time)) {
    stop("All time-point column names must be numeric.")
  }
  
  # Cluster labels in increasing order
  cluster_labels <- sort(unique(dataM$Cluster))
  K <- length(cluster_labels)
  
  # Color palette
  color.pal <- palette.colors(
    palette = "Polychrome 36",
    alpha = 1
  )[-2]
  
  if (K > length(color.pal)) {
    stop("There are more clusters than available colors in the palette.")
  }
  
  # Sort data by cluster and ID
  dataM <- dataM[order(dataM$Cluster, dataM$ID), , drop = FALSE]
  
  # Create chart
  p <- highcharter::highchart() |>
    highcharter::hc_chart(
      type = "line",
      zoomType = "xy"
    ) |>
    highcharter::hc_title(
      text = "Trajectories by Cluster"
    ) |>
    highcharter::hc_subtitle(
      text = paste(
        "Hover over a line to highlight it,",
        "click and drag to zoom,",
        "and select a cluster in the legend to filter."
      )
    ) |>
    highcharter::hc_xAxis(
      title = list(text = "Time")
    ) |>
    highcharter::hc_yAxis(
      title = list(text = "Value")
    ) |>
    highcharter::hc_tooltip(
      pointFormat = ""
    ) |>
    highcharter::hc_plotOptions(
      series = list(
        states = list(
          hover = list(
            enabled = TRUE,
            lineWidthPlus = 2
          ),
          inactive = list(
            opacity = 0.2
          )
        )
      )
    )
  
  # AdataM trajectories
  for (i in seq_len(nrow(dataM))) {
    
    cluster <- dataM$Cluster[i]
    
    # Position of cluster in increasing cluster order
    cluster_index <- match(cluster, cluster_labels)
    
    # Trajectory data
    trajectory <- lapply(
      seq_along(time),
      function(k) {
        c(
          time[k],
          as.numeric(dataM[[time_columns[k]]][i])
        )
      }
    )
    
    # Only the first trajectory of each cluster appears in the legend
    first_in_cluster <- i == which(dataM$Cluster == cluster)[1]
    
    # JavaScript events
    series_events <- list(
      
      # Hover: highlight ONLY the trajectory being hovered over
      mouseOver = highcharter::JS(
        "function () {
           var hovered = this;
           
           this.chart.series.forEach(function (s) {
             if (s === hovered) {
               s.setState('hover');
             } else {
               s.setState('inactive');
             }
           });
         }"
      ),
      
      # Stop highlighting when the mouse leaves
      mouseOut = highcharter::JS(
        "function () {
           this.chart.series.forEach(function (s) {
             s.setState('');
           });
         }"
      )
    )
    
    # Clicking the legend entry toggles the visibility
    # of all trajectories belonging to that cluster
    if (first_in_cluster) {
      
      series_events$legendItemClick <- highcharter::JS(
        "function () {
           var cluster = this.userOptions.clusterGroup;
           var newVisibility = !this.visible;
           
           this.chart.series.forEach(function (s) {
             if (s.userOptions.clusterGroup === cluster) {
               s.setVisible(newVisibility, false);
             }
           });
           
           this.chart.redraw();
           
           return false;
         }"
      )
    }
    
    p <- p |>
      highcharter::hc_add_series(
        data = trajectory,
        type = "line",
        name = paste("Cluster", cluster),
        color = color.pal[cluster_index],
        showInLegend = first_in_cluster,
        legendIndex = cluster_index,
        clusterGroup = paste0("cluster_", cluster),
        lineWidth = 1,
        marker = list(
          enabled = FALSE
        ),
        events = series_events
      )
  }
  
  p
}
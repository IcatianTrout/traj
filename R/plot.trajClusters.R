#'@title Plot \code{trajClusters} object
#'
#'@description Plots the curves corresponding to (or closest to) the centroids of the clusters and plots a random sample from each groups.
#'
#'@param x object of class \code{trajClusters} as returned by the function
#'  \code{trajClusters()}.
#'@param sample.size the number of random trajectories to be randomly sampled
#'  from each cluster. Defaults to \code{5}.
#'@param ask logical. If \code{TRUE}, the user is asked before each plot. Defaults to
#'  \code{TRUE}.
#'@param which.plots either \code{NULL} or a vector of integers. If \code{NULL}, every
#'  available plot is displayed. If a vector is supplied, only the corresponding
#'  plots will be displayed.
#'@param ... other parameters to be passed through to plotting functions.
#'
#'@importFrom grDevices palette.colors devAskNewPage graphics.off
#'@importFrom graphics legend lines par grid
#'@importFrom stats dist
#'
#'@examples
#' \dontrun{
#'data("trajdata")
#'trajdata.noGrp <- trajdata[, -which(colnames(trajdata) == "Group")] #remove the Group column
#'
#'m = trajMeasures(trajdata.noGrp, ID = TRUE)

#'c3 = trajClusters(m, nclusters = 3)
#'
#'plot(c3)
#'
#'#The pointwise mean trajectories correspond to the third and fourth displayed plots.
#'
#'c4 = trajClusters(m, nclusters = 4)
#'
#'plot(c4, which.plots = 3:4)
#'
#'}
#'
#'
#'@rdname plot.trajClusters
#'
#'@export
plot.trajClusters <-
  function(x,
           sample.size = 5,
           ask = TRUE,
           which.plots = NULL,
           ...
           ) {
    if (!is.null(which.plots) &
        (!is.numeric(which.plots) | !is.vector(which.plots))) {
      stop(
        "The argument 'which.plots' should be a subset of the plots required, specified as a vector, so either 1, 2 or 1:2."
      )
    }
    
    if (!is.null(which.plots)) {
      which.plots <- which.plots[order(which.plots)]
    }
    
    current.ask.status <- devAskNewPage(ask = NULL)
    on.exit(devAskNewPage(ask = current.ask.status))  # Restore ask status on exit
    devAskNewPage(ask = ask)
    
    color.pal <- palette.colors(palette = "Polychrome 36", alpha = 1)[-2]
    par(mfrow = c(1, 1))
    
    if(is.null(which.plots) | 1 %in% which.plots){

      centroids.data <- x$data[x$data[, 1] %in% x$ID.centers, -1]

      plot(
        x = 0,
        y = 0,
        xlim = c(min(unlist(x$time[, -1]), na.rm = TRUE), max(unlist(x$time[, -1]), na.rm = TRUE)),
        ylim = c(min(unlist(centroids.data), na.rm = TRUE), max(unlist(centroids.data), na.rm = TRUE)),
        type = "n",
        xlab = "",
        ylab = "",
        main = "Centroids"
      )
      legend(
        "topright",
        col = color.pal[1:x$nclusters],
        legend = paste(1:x$nclusters),
        lty = rep(1, x$nclusters)
      )
      
      for (j in 1:x$nclusters) {
        lines(
          x = x$time[which(x$time[, 1] == x$ID.centers[j]), -1],
          y = x$data[which(x$data[, 1] == x$ID.centers[j]), -1],
          type = "b",
          pch = 16,
          col = color.pal[j]
        )
      }
      
    }
    
    if(is.null(which.plots) | 2 %in% which.plots){
      traj.by.clusters <- list()
      for (k in seq_len(x$nclusters)) {
        traj.by.clusters[[k]] <-
          x$data[which(x$partition[, 2] == k), -c(1), drop = FALSE]
      }
      
      time.by.clusters <- list()
      for (k in 1:x$nclusters) {
        time.by.clusters[[k]] <-
          x$time[which(x$partition[, 2] == k),-c(1), drop = FALSE]
      }
      
      # Plot (max) sample.size random trajectories from each group
      smpl.traj.by.clusters <- list()
      smpl.time.by.clusters <- list()
      
      smpl.traj <-
        matrix(nrow = 0, ncol = ncol(x$data) - 1)
      smpl.time <-
        matrix(nrow = 0, ncol = ncol(x$time) - 1)
      
      for (k in seq_len(x$nclusters)) {
        size <- min(sample.size, nrow(traj.by.clusters[[k]]))
        smpl <-
          sample(x = seq_len(nrow(traj.by.clusters[[k]])),
                 size = size,
                 replace = FALSE)
        smpl <- smpl[order(smpl)]
        
        smpl.traj.by.clusters[[k]] <- traj.by.clusters[[k]][smpl, , drop = FALSE]
        smpl.time.by.clusters[[k]] <- time.by.clusters[[k]][smpl, , drop = FALSE]
        
        smpl.traj <- rbind(smpl.traj, smpl.traj.by.clusters[[k]])
        smpl.time <- rbind(smpl.time, smpl.time.by.clusters[[k]])
      }
      
      par(mfrow = c(1, 1))
      
      plot(
        x = 0,
        y = 0,
        xlim = c(min(smpl.time, na.rm = TRUE), max(smpl.time, na.rm = TRUE)),
        ylim = c(min(smpl.traj, na.rm = TRUE), max(smpl.traj, na.rm = TRUE)),
        type = "n",
        xlab = "",
        ylab = "",
        main = "Sample trajectories"
      )
      
      for (k in seq_len(x$nclusters)) {
        for (i in seq_len(min(size, x$partition.summary[k]))) {
          lines(
            x = smpl.time.by.clusters[[k]][i,],
            y = smpl.traj.by.clusters[[k]][i,],
            type = "l",
            col = color.pal[k]
          )
        }
        legend(
          "topright",
          col = color.pal[1:k],
          legend = paste(seq_len(x$nclusters))[1:k],
          lty = rep(1, k)
        )
      }
    }
    
    print("See also 'CVIplot' for a plot of the statistic used to determined the number of clusters (if applicable) and see 'scatterplots' for scatter plots of the measures involved in the clustering.")
  }
#'@rdname plot.trajClusters
#'
#'@export
scatterplots <- function(x, ask = TRUE, ...) {
  
  current.ask.status <- devAskNewPage(ask = NULL)
  on.exit(devAskNewPage(ask = current.ask.status))  # Restore ask status on exit
  devAskNewPage(ask = ask)
  
  color.pal <- palette.colors(palette = "Polychrome 36", alpha = 1)[-2]
  
  nb.measures <- ncol(x$selection) - 1
  scatter.condition <- (nb.measures > 1)
  
  if (scatter.condition) {
    
    which.scatter <- c(1:nb.measures)
    
    selection <- x$selection[, -c(1)]
    
    # Set up the most compact grid depending on the number of selected measures
    X <- sqrt(nb.measures - 1)
    
    int.X <- floor(X)
    frac.X <- X - int.X
    
    
    if (frac.X == 0) {
      good.grid <- c(int.X, int.X)
    }
    
    if ((frac.X > 0) & (frac.X < 0.5)) {
      good.grid <- c(int.X, int.X + 1)
    }
    
    if (frac.X >= 0.5) {
      good.grid <- c(int.X + 1, int.X + 1)
    }
    
    for (m in which.scatter) {
      par(mfrow = good.grid)
      
      for (n in seq_len(nb.measures)[-m]) {
        plot(
          x = 0,
          y = 0,
          xlim = c(min(selection[, m]), max(selection[, m])),
          ylim = c(min(selection[, n]), max(selection[, n])),
          type = "n",
          xlab = paste(colnames(selection[m])),
          ylab = paste(colnames(selection[n])),
          main = paste(
            "Scatter plot of ",
            paste(colnames(selection[m])),
            " vs ",
            paste(colnames(selection)[n]),
            sep = ""
          )
        )
        
        set.seed(38550)
        S <- sample(1:nrow(selection), nrow(selection), replace = FALSE)
        
        for(s in S){
          lines(
            x = selection[s, m],
            y = selection[s, n],
            type = "p",
            pch = (x$partition[s,2] - 1),
            col = color.pal[x$partition[s,2]],
            bg = color.pal[x$partition[s,2]]
          )
        }
        
        
        legend(
          "topright",
          lty = rep(0, x$nclusters),
          pch = c(0:(x$nclusters-1)),
          col = color.pal[1:x$nclusters],
          legend = paste(seq_len(x$nclusters))[1:x$nclusters]
        )
      }
    }
  } else{
    print("Since only one measure participated in the clustering there are no scatter plots to display.")
  }
}
#'@rdname plot.trajClusters
#'@export
CVIplot <- function(x, ...) {
  
  CVI <- x$cluster.validity.indices
  
  if (!is.null(CVI)) {
    par(mfrow = c(2, 1))
    color.pal <- palette.colors(palette = "Okabe-Ito", alpha = 1)
    plot(y = 0, x = 0, xlim = c(2, 1+ncol(CVI)), ylim = c(0,1), type = "n", xlab = "k", ylab="", main = "Scaled cluster validity indices")
    grid(nx = NULL, ny = NULL,
         lty = 2,      
         col = "gray", 
         lwd = 1)      
    for(j in 1:nrow(CVI)){
      lines(y = CVI[j, ], x = 2:(1+ncol(CVI)), type = "b", pch = 16, xlab = "k", main = "Scaled cluster validity indices", col = color.pal[j])
    }
    legend(
      "topright",
      col = color.pal[1:nrow(CVI)],
      legend = c("C-index", "C-H", "W-G"),
      lty = rep(1, nrow(CVI))
    )
    
    RVR <- x$ranked.voting.results
    
    plot(y = 0, x = 0, xlim = c(2,(1+length(RVR))), ylim = c(0,3), type = "n", xlab = "k", ylab="", main = "Ranked voting results")
    grid(nx = NULL, ny = NULL,
         lty = 2,      
         col = "gray", 
         lwd = 1)      
    
    lines(y = RVR, x = 2:(1+length(RVR)), type = "b", pch = 16)
    
    
  } else {
    print("There are no cluster validity indices to plot.")
  }
}
#'@rdname plot.trajClusters
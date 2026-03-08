#'@title Plot \code{trajClusters} object
#'
#'@description Plots the curves corresponding to (or closest to) the centroids of the clusters and plots a random sample from each groups.
#'
#'@param x object of class \code{trajClusters} as returned by the function
#'  \code{trajClusters()}.
#'@param sample.size the number of trajectories to be randomly sampled
#'  from each clusters. If \code{NULL}, all the trajectories are included in a random order. Defaults to \code{5}.
#'@param ask logical. If \code{TRUE}, the user is asked before each plot. Defaults to
#'  \code{TRUE}.
#'@param which.plots either \code{NULL} or a vector of integers. If \code{NULL}, every
#'  available plot is displayed. If a vector is supplied, only the corresponding
#'  plots will be displayed.
#'@param which.scatter either \code{NULL} or a vector of integers that is a subset of the \code{measure} argument used in function \code{trajClusters} to produce object \code{x}. If \code{NULL}, every available scatter plots are displayed. If a vector is supplied, only the corresponding plots will be displayed.
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
      
      ## Plot (max) sample.size random trajectories from each group
      smpl.traj.by.clusters <- list()
      smpl.time.by.clusters <- list()
      
      smpl.traj <-
        matrix(nrow = 0, ncol = ncol(x$data))
      smpl.time <-
        matrix(nrow = 0, ncol = ncol(x$time) - 1)

      size <- c()
      
      for (k in seq_len(x$nclusters)) {
        if(!is.null(sample.size)){
        size[k] <- min(sample.size, nrow(traj.by.clusters[[k]]))
        } else{
          size[k] <- nrow(traj.by.clusters[[k]])
        }
        smpl <-
          sample(x = seq_len(nrow(traj.by.clusters[[k]])),
                 size = size[k],
                 replace = FALSE)
        smpl <- smpl[order(smpl)]
        
        smpl.traj.by.clusters[[k]] <- cbind(traj.by.clusters[[k]][smpl, , drop = FALSE],k)
        smpl.time.by.clusters[[k]] <- time.by.clusters[[k]][smpl, , drop = FALSE]
        
        smpl.traj <- rbind(smpl.traj, smpl.traj.by.clusters[[k]])
        smpl.time <- rbind(smpl.time, smpl.time.by.clusters[[k]])
      }
      
      ## Shuffle the rows
      s <- sample(seq_len(nrow(smpl.traj)), nrow(smpl.traj), replace=F)
      smpl.traj <- smpl.traj[s, ]
      smpl.time <- smpl.time[s, ]
      
      par(mfrow = c(1, 1))
      
      plot(
        x = 0,
        y = 0,
        xlim = c(min(smpl.time, na.rm = TRUE), max(smpl.time, na.rm = TRUE)),
        ylim = c(min(smpl.traj[, -ncol(smpl.traj)], na.rm = TRUE), max(smpl.traj[, -ncol(smpl.traj)], na.rm = TRUE)),
        type = "n",
        xlab = "",
        ylab = "",
        main = "Sample trajectories"
      )
      
      ## Plot the trajectories in a random order
      for(i in seq_len(nrow(smpl.traj))){
        k <- smpl.traj[i, ncol(smpl.traj)]
        
        lines(
          x = smpl.time[i,],
          y = smpl.traj[i, -ncol(smpl.traj)],
          type = "l",
          col = color.pal[k]
        )
      }
      legend(
        "topright",
        col = color.pal[seq_len(x$nclusters)],
        legend = paste(seq_len(x$nclusters)),
        lty = rep(1, x$nclusters)
      )
    }
    
    print("See also 'CVIplot' for a plot of the statistic used to determined the number of clusters (if applicable) and see 'scatterplots' for scatter plots of the measures involved in the clustering.")
  }
#'@rdname plot.trajClusters
#'
#'@export
scatterplots <- function(x, ask = TRUE, which.scatter = NULL, ...) {
  
  if( (!is.null(which.scatter)) & (sum(!(which.scatter %in% x$select)) > 0) ){stop("The argument which.scatter should be a subset of the measure argument used in function trajClusters.")}
  
  current.ask.status <- devAskNewPage(ask = NULL)
  on.exit(devAskNewPage(ask = current.ask.status))  # Restore ask status on exit
  devAskNewPage(ask = ask)
  
  color.pal <- palette.colors(palette = "Polychrome 36", alpha = 1)[-2]
  
  nb.measures <- ncol(x$selection) - 1
  scatter.condition <- (nb.measures > 1)
  
  if (scatter.condition) {
    
    if(is.null(which.scatter)){
      v <- c(1:nb.measures)
    } else{
      v <- c(1:length(which.scatter))
    }
    
    
    selection.y <- x$selection[, -c(1), drop = F]
    
    if(!is.null(which.scatter)){
      selection.x <- selection.y[, which.scatter, drop = F]
    } else{
      selection.x <- selection.y
    }
    
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
    
    for (m in v) {
      par(mfrow = good.grid)
      
      if(!is.null(which.scatter)){
        w <- which(x$select == which.scatter[m])
      } else {
        w <- m
      }
      for (n in seq_len(nb.measures)[-w]) {
        plot(
          x = 0,
          y = 0,
          xlim = c(min(selection.x[, m]), max(selection.x[, m])),
          ylim = c(min(selection.y[, n]), max(selection.y[, n])),
          type = "n",
          xlab = paste(colnames(selection.x[m])),
          ylab = paste(colnames(selection.y[n])),
          main = paste(
            "Scatter plot of ",
            paste(colnames(selection.x[m])),
            " vs ",
            paste(colnames(selection.y)[n]),
            sep = ""
          )
        )
        
        set.seed(38550)
        S <- sample(1:nrow(selection.x), nrow(selection.x), replace = FALSE)
        
        for(s in S){
          lines(
            x = selection.x[s, m],
            y = selection.y[s, n],
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
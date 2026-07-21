#'@title Plot \code{trajClusters} object
#'
#'@description \code{plot()} sequentially plots (i) the clusters centroids, (ii) a random sample of trajectories from each cluster, and (iii) the median of each standardized measure by cluster. \code{(scatterplots()} displays the two dimensional projections of the standardized data in feature space for every pair of measures. \code{CVIplot()} graphs the internal cluster validity indices as a function of the number of clusters as well as the result of the ranked voting system for determining the optimal number of clusters.
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
#'@param which.scatter either \code{NULL} or a vector of integers that is a subset of the \code{measure} argument used in the function \code{trajClusters} to produce object \code{x}. If \code{NULL}, every available scatter plots are displayed. If a vector is supplied, only the corresponding plots will be displayed.
#'@param N the maximum number of points present in each scatter plots. If a non \code{NULL} value is specified, \code{N} points are sampled randomly in a way that preserves the relative groups sizes. If \code{NULL} (the default), all the points are plotted.
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
#'
#'c3 = trajClusters(m, nclusters = 3)
#'
#'plot(c3, which.plots = 1, ask= FALSE) ## centroids
#'
#'plot(c3, which.plots = 2, ask= FALSE) ## sample trajectories
#'
#'plot(c3, which.plots = 3, ask= FALSE) ## median of standardized measures by cluster
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
        "The argument 'which.plots' should be a subset of the plots required, specified as a vector, eg. 1, 2, 2:3."
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
      s <- sample(seq_len(nrow(smpl.traj)), nrow(smpl.traj), replace = FALSE)
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
    
    if(is.null(which.plots) | 3 %in% which.plots){
      
      cl.medians <- data.frame(matrix(NA, nrow = x$nclusters, ncol = length(x$select)))
      colnames(cl.medians) <- colnames(x$selection)[-1]
      
      for (i in seq_len(x$nclusters)) {
        which.i <- which(x$partition[, 2] == i)
        cl.medians[i, ] <- apply(x$standardized.data[which.i, ], 2, median)
      }
      
      # mar=c(bottom, left, top, right)
      par(mfrow = c(1,1), mar = c(5, 8, 4, 5), xpd = TRUE)
      color.pal <- palette.colors(palette = "Polychrome 36", alpha = 1)[-2]
      
      hor.labels <- c()
      for(m in seq_len(length(x$select))){
        if(x$select[m] == 1){ hor.labels <- c(hor.labels, paste("m1 : max", sep = ""))}
        if(x$select[m] == 2){ hor.labels <- c(hor.labels, paste("m2 : min", sep = ""))}
        if(x$select[m] == 3){ hor.labels <- c(hor.labels, paste("m3 : range", sep = ""))}
        if(x$select[m] == 4){ hor.labels <- c(hor.labels, paste("m4 : mean", sep = ""))}
        if(x$select[m] == 5){ hor.labels <- c(hor.labels, paste("m5 : SD", sep = ""))}
        if(x$select[m] == 6){ hor.labels <- c(hor.labels, paste("m6 : slope", sep = ""))}
        if(x$select[m] == 7){ hor.labels <- c(hor.labels, paste("m7 : intercept", sep = ""))}
        if(x$select[m] == 8){ hor.labels <- c(hor.labels, paste("m8 : R^2", sep = ""))}
        if(x$select[m] == 9){ hor.labels <- c(hor.labels, paste("m9 : int. rate", sep = ""))}
        if(x$select[m] == 10){ hor.labels <- c(hor.labels, paste("m10 : var. rate", sep = ""))}
        if(x$select[m] == 11){ hor.labels <- c(hor.labels, paste("m11 : contrast", sep = ""))}
        if(x$select[m] == 12){ hor.labels <- c(hor.labels, paste("m12 : tot var", sep = ""))}
        if(x$select[m] == 13){ hor.labels <- c(hor.labels, paste("m13 : spikiness", sep = ""))}
        if(x$select[m] == 14){ hor.labels <- c(hor.labels, paste("m14 : max f'", sep = ""))}
        if(x$select[m] == 15){ hor.labels <- c(hor.labels, paste("m15 : min f'", sep = ""))}
        if(x$select[m] == 16){ hor.labels <- c(hor.labels, paste("m16 : SD f'", sep = ""))}
        if(x$select[m] == 17){ hor.labels <- c(hor.labels, paste("m17 : f' var. rate", sep = ""))}
        if(x$select[m] == 18){ hor.labels <- c(hor.labels, paste("m18 : max f''", sep = ""))}
        if(x$select[m] == 19){ hor.labels <- c(hor.labels, paste("m19 : min f''", sep = ""))}
        if(x$select[m] == 20){ hor.labels <- c(hor.labels, paste("m20 : SD f''", sep = ""))}
      }
      
      plot(
        x = 0,
        y = 0,
        xlim = c(min(cl.medians), max(cl.medians)),
        ylim = c(1,length(x$select)),
        type = "n",
        xlab = "",
        ylab = "",
        yaxt = "n",
        main = "Standardized feature medians by clusters"
      )
      
      graphics::axis(2,
           at = seq_len(length(x$select)),
           labels = rev(hor.labels),
           las = 1)       # horizontal labels
      
      for(m in seq_len(length(x$select))){
        lines(x = c(min(cl.medians), max(cl.medians)), y = c(m,m), col="gray")
        
        for(s in seq_len(x$nclusters)){
          lines(
            x = cl.medians[s, m],
            y = rev(1:length(x$select))[m],
            type = "p",
            pch = s - 1,
            col = color.pal[s],
            bg = color.pal[s]
          )
        }
      }
      
      
      usr <- par("usr")
      
      legend(x = usr[2],
             y = usr[4],
             legend = paste(seq_len(x$nclusters))[1:x$nclusters],
             col = color.pal[1:x$nclusters],
             lty = rep(0, x$nclusters),
             pch = c(0:(x$nclusters - 1)))
    }

print("See also 'CVIplot' for a plot of the statistic used to determined the number of clusters (if applicable) and see 'scatterplots' for scatter plots of the measures involved in the clustering.")
}
#'@rdname plot.trajClusters
#'
#'@export
scatterplots <- function(x, ask = TRUE, which.scatter = NULL, N = NULL, ...) {
  
  if( (!is.null(which.scatter)) & (sum(!(which.scatter %in% x$select)) > 0) ){stop("The argument which.scatter should be a subset of the measure argument used in function trajClusters.")}
  
  if ( !is.null(N) && !( ( is.numeric(N) && (length(N) == 1)) && (N %in% seq_len(nrow(x$selection))) ) ){
    stop("'N' should be either NULL or a numerical integer smaller than the total number of admissible trajectories.")
  }
  
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
    
    
    selection.y <- x$standardized.data
    
    if(!is.null(which.scatter)){
      selection.x <- selection.y[, which(x$select %in% which.scatter), drop = FALSE]
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
    
    
    selection.x0 <- selection.x
    selection.y0 <- selection.y
    
    grps <- x$partition[, 2]
    
    if(!is.null(N)){
      selection.x.new <- selection.x[0, , drop = FALSE]
      selection.y.new <- selection.y[0, , drop = FALSE]
      new.grp.size <- round(N*x$partition.summary/sum(x$partition.summary))
      grps <- c()
      for(k in seq_len(x$nclusters)){
        s <- sample(seq_len(x$partition.summary[k]), size = new.grp.size[k], replace = FALSE)
        aux.x <- selection.x[which(x$partition[, 2] == k), , drop = FALSE][s, , drop = FALSE]
        aux.y <- selection.y[which(x$partition[, 2] == k), , drop = FALSE][s, , drop = FALSE]
        selection.x.new <- rbind(selection.x.new, aux.x)
        selection.y.new <- rbind(selection.y.new, aux.y)
        grps <- c(grps, rep(k, new.grp.size[k]))
      }
      selection.x <- selection.x.new
      selection.y <- selection.y.new
    }
    
    for (m in v) {
      if(!is.null(which.scatter)){
        w <- which(x$select == which.scatter[m])
      } else {
        w <- m
      }
      
      # mar=c(bottom, left, top, right)
      par(mfrow = good.grid, mar = c(5, 5, 4, 5), xpd = TRUE)
      
      for (n in seq_len(nb.measures)[-w]) {
        if(colnames(selection.x0[m]) == "m1"){ main1 <- paste(colnames(selection.x0[m])," : max", sep = "")}
        if(colnames(selection.x0[m]) == "m2"){ main1 <- paste(colnames(selection.x0[m])," : min", sep = "")}
        if(colnames(selection.x0[m]) == "m3"){ main1 <- paste(colnames(selection.x0[m])," : range", sep = "")}
        if(colnames(selection.x0[m]) == "m4"){ main1 <- paste(colnames(selection.x0[m])," : mean", sep = "")}
        if(colnames(selection.x0[m]) == "m5"){ main1 <- paste(colnames(selection.x0[m])," : SD", sep = "")}
        if(colnames(selection.x0[m]) == "m6"){ main1 <- paste(colnames(selection.x0[m])," : slope", sep = "")}
        if(colnames(selection.x0[m]) == "m7"){ main1 <- paste(colnames(selection.x0[m])," : intercept", sep = "")}
        if(colnames(selection.x0[m]) == "m8"){ main1 <- paste(colnames(selection.x0[m])," : R^2", sep = "")}
        if(colnames(selection.x0[m]) == "m9"){ main1 <- paste(colnames(selection.x0[m])," : int. rate", sep = "")}
        if(colnames(selection.x0[m]) == "m10"){ main1 <- paste(colnames(selection.x0[m])," : var. rate", sep = "")}
        if(colnames(selection.x0[m]) == "m11"){ main1 <- paste(colnames(selection.x0[m])," : contrast", sep = "")}
        if(colnames(selection.x0[m]) == "m12"){ main1 <- paste(colnames(selection.x0[m])," : tot var", sep = "")}
        if(colnames(selection.x0[m]) == "m13"){ main1 <- paste(colnames(selection.x0[m])," : spikiness", sep = "")}
        if(colnames(selection.x0[m]) == "m14"){ main1 <- paste(colnames(selection.x0[m])," : max f'", sep = "")}
        if(colnames(selection.x0[m]) == "m15"){ main1 <- paste(colnames(selection.x0[m])," : min f'", sep = "")}
        if(colnames(selection.x0[m]) == "m16"){ main1 <- paste(colnames(selection.x0[m])," : SD f'", sep = "")}
        if(colnames(selection.x0[m]) == "m17"){ main1 <- paste(colnames(selection.x0[m])," : f' var. rate", sep = "")}
        if(colnames(selection.x0[m]) == "m18"){ main1 <- paste(colnames(selection.x0[m])," : max f''", sep = "")}
        if(colnames(selection.x0[m]) == "m19"){ main1 <- paste(colnames(selection.x0[m])," : min f''", sep = "")}
        if(colnames(selection.x0[m]) == "m20"){ main1 <- paste(colnames(selection.x0[m])," : SD f''", sep = "")}
        if(colnames(selection.y0[n]) == "m1"){ main2 <- paste(colnames(selection.y0[n])," : max", sep = "")}
        if(colnames(selection.y0[n]) == "m2"){ main2 <- paste(colnames(selection.y0[n])," : min", sep = "")}
        if(colnames(selection.y0[n]) == "m3"){ main2 <- paste(colnames(selection.y0[n])," : range", sep = "")}
        if(colnames(selection.y0[n]) == "m4"){ main2 <- paste(colnames(selection.y0[n])," : mean", sep = "")}
        if(colnames(selection.y0[n]) == "m5"){ main2 <- paste(colnames(selection.y0[n])," : SD", sep = "")}
        if(colnames(selection.y0[n]) == "m6"){ main2 <- paste(colnames(selection.y0[n])," : slope", sep = "")}
        if(colnames(selection.y0[n]) == "m7"){ main2 <- paste(colnames(selection.y0[n])," : intercept", sep = "")}
        if(colnames(selection.y0[n]) == "m8"){ main2 <- paste(colnames(selection.y0[n])," : R^2", sep = "")}
        if(colnames(selection.y0[n]) == "m9"){ main2 <- paste(colnames(selection.y0[n])," : int. rate", sep = "")}
        if(colnames(selection.y0[n]) == "m10"){ main2 <- paste(colnames(selection.y0[n])," : var. rate", sep = "")}
        if(colnames(selection.y0[n]) == "m11"){ main2 <- paste(colnames(selection.y0[n])," : contrast", sep = "")}
        if(colnames(selection.y0[n]) == "m12"){ main2 <- paste(colnames(selection.y0[n])," : tot var", sep = "")}
        if(colnames(selection.y0[n]) == "m13"){ main2 <- paste(colnames(selection.y0[n])," : spikiness", sep = "")}
        if(colnames(selection.y0[n]) == "m14"){ main2 <- paste(colnames(selection.y0[n])," : max f'", sep = "")}
        if(colnames(selection.y0[n]) == "m15"){ main2 <- paste(colnames(selection.y0[n])," : min f'", sep = "")}
        if(colnames(selection.y0[n]) == "m16"){ main2 <- paste(colnames(selection.y0[n])," : SD f'", sep = "")}
        if(colnames(selection.y0[n]) == "m17"){ main2 <- paste(colnames(selection.y0[n])," : f' var. rate", sep = "")}
        if(colnames(selection.y0[n]) == "m18"){ main2 <- paste(colnames(selection.y0[n])," : max f''", sep = "")}
        if(colnames(selection.y0[n]) == "m19"){ main2 <- paste(colnames(selection.y0[n])," : min f''", sep = "")}
        if(colnames(selection.y0[n]) == "m20"){ main2 <- paste(colnames(selection.y0[n])," : SD f''", sep = "")}
        
        plot(
          x = 0,
          y = 0,
          xlim = c(min(selection.x0[, m]), max(selection.x0[, m])),
          ylim = c(min(selection.y0[, n]), max(selection.y0[, n])),
          type = "n",
          xlab = paste(colnames(selection.x0[m])),
          ylab = paste(colnames(selection.y0[n])),
          main = paste(main1,main2, sep = "\n")
        )
        
        S <- sample(seq_len(nrow(selection.x)), nrow(selection.x), replace = FALSE) #randomize
        
        for(s in S){
          lines(
            x = selection.x[s, m],
            y = selection.y[s, n],
            type = "p",
            pch = (grps[s] - 1),
            col = color.pal[grps[s]],
            bg = color.pal[grps[s]]
          )
        }
        
        usr <- par("usr")
        
        legend(x = usr[2],
               y = usr[4],
               legend = paste(seq_len(x$nclusters))[1:x$nclusters],
               col = color.pal[1:x$nclusters],
               lty = rep(0, x$nclusters),
               pch = c(0:(x$nclusters-1)))
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
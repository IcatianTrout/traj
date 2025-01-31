#'@title Plots \code{trajClusters} objects
#'
#'@description Plots the
#'  cluster-specific median and mean trajectories and a random sample of
#'  trajectories from each cluster.
#'
#'@param x object of class \code{trajClusters} as returned by
#'  \code{Step3Cluster}.
#'@param sample.size the number of random trajectories to be randomly sampled
#'  from each cluster. Defaults to \code{5}.
#'@param ask logical. If \code{TRUE}, the user is asked before each plot. Defaults to
#'  \code{TRUE}.
#'@param which.plots either \code{NULL} or a vector of integers. If \code{NULL}, every
#'  available plot is displayed. If a vector is supplied, only the corresponding
#'  plots will be displayed.
#'@param spline logical. If \code{TRUE}, each trajectory will be smoothed using
#'  smoothing splines and the median and mean trajectories will be plotted from
#'  the smoothed trajectories. Defaults to \code{FALSE}
#'@param ... other parameters to be passed through to plotting functions.
#'
#'@importFrom grDevices palette.colors
#'@importFrom graphics legend lines par polygon barplot
#'@importFrom grDevices devAskNewPage graphics.off
#'@importFrom stats smooth.spline predict
#'
#'
#'@seealso \code{\link[traj]{Step3Clusters}}
#'
#'@examples
#' \dontrun{
#'data("trajdata")
#'trajdata.noGrp <- trajdata[, -which(colnames(trajdata) == "Group")] #remove the Group column
#'
#'m = Step1Measures(trajdata.noGrp, ID = TRUE)
#'s = Step2Selection(m)
#'c3 = Step3Clusters(s, nclusters = 3)
#'
#'plot(c3)
#'
#'#The pointwise mean trajectories correspond to the third and fourth displayed plots.
#'
#'c4 = Step3Clusters(s, nclusters = 4)
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
           spline = FALSE, ...) {
    if (!is.null(which.plots) &
        (!is.numeric(which.plots) | !is.vector(which.plots))) {
      stop(
        "The argument 'which.plots' should be a subset of the plots required, specified as a vector such as c(1, 3:5), or NULL if all the plots are required."
      )
    }
    
    if (!is.null(which.plots)) {
      which.plots <- which.plots[order(which.plots)]
    }
    
    current.ask.status <- devAskNewPage(ask = NULL)
    on.exit(devAskNewPage(ask = current.ask.status))  # Restore ask status on exit
    devAskNewPage(ask = ask)
    
    color.pal <- palette.colors(palette = "Polychrome 36", alpha = 1)[-2]
    
    
    ### Mean and median traj ###
    universal.time <- x$time[1,-1]
    is.different <- c()
    for (i in seq_len(nrow(x$time))) {
      is.different[i] <-
        !identical(x$time[i,-1], universal.time)
    }
    homogeneous.times <- (sum(is.different) == 0)
    
    med.mean.condition <-
      ((spline |
          homogeneous.times) &
         (is.null(which.plots) |
            sum(c(1:4) %in% which.plots) > 0))
    if (med.mean.condition) {
      if (is.null(which.plots)) {
        which.med.mean <- c(1:4)
      } else{
        which.med.mean <-
          which(c(1:4) %in% which.plots)
      }
      
      if (spline == TRUE) {
        data <- x$data[,-1]
        time <- x$time[,-1]
        
        a <- min(time, na.rm = TRUE)
        b <- max(time, na.rm = TRUE)
        
        min.delta.by.row <- c()
        for (i in seq_len(nrow(time))) {
          cc.time.i <- time[i, complete.cases(time[i,])]
          delta.i <- c()
          for (j in 1:(length(cc.time.i) - 1)) {
            delta.i[j] <- cc.time.i[j + 1] - cc.time.i[j]
          }
          min.delta.by.row[i] <- min(delta.i)
        }
        min.delta <- min(min.delta.by.row)
        N <- ceiling((b - a) / min.delta)
        
        time.new <- matrix(NA, nrow = nrow(time), ncol = N)
        data.new <- matrix(NA, nrow = nrow(time), ncol = N)
        
        for (i in seq_len(nrow(time))) {
          time.new[i,] <- seq(from = a,
                              to = b,
                              length.out = N)
          
          for (j in 1:N) {
            if (time.new[i, j] < min(time[i,], na.rm = TRUE) |
                time.new[i, j] > max(time[i,], na.rm = TRUE)) {
              time.new[i, j] <- NA
            }
          }
          spl.i <-
            smooth.spline(x = na.omit(time[i,]),
                          y = na.omit(data[i,]),
                          cv = FALSE)
          pred.i <- predict(spl.i, x = na.omit(time.new[i,]))
          data.new[i,!is.na(time.new[i,])] <- pred.i$y
          
        }
        
        time.new <- cbind(x$time[, "ID"], time.new)
        data.new <- cbind(x$data[, "ID"], data.new)
        
      }
      
      # Set up the most compact grid depending on the number of clusters
      X <- sqrt(x$nclusters)
      
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
      
      mean.na.rm <- function(x) {
        return(mean(x, na.rm = TRUE))
      }
      
      sd.na.rm <- function(x) {
        return(sd(x, na.rm = TRUE))
      }
      
      median.na.rm <- function(x) {
        return(median(x, na.rm = TRUE))
      }
      
      mean.traj <- list()
      
      median.traj <- list()
      
      time.for.plots <- list()
      
      CI.low <- list()
      
      CI.high <- list()
      
      
      for (j in seq_len(x$nclusters)) {
        which.j <- which(x$partition[, 2] == j)
        
        if (spline == TRUE) {
          data.cluster.j <- data.new[which.j, -1, drop = FALSE]
        } else{
          data.cluster.j <- x$data[which.j, -1, drop = FALSE]
        }
        
        col.keep <-
          colSums(is.na(data.cluster.j)) < nrow(data.cluster.j) #identifies the columns which are not mono NA
        
        data.cluster.j <- data.cluster.j[, col.keep, drop = FALSE]
        
        mean.traj.j <- apply(data.cluster.j, 2, mean.na.rm)
        
        mean.traj[[j]] <- mean.traj.j
        
        median.traj[[j]] <- apply(data.cluster.j, 2, median.na.rm)
        
        if (spline == TRUE) {
          time.for.plots[[j]] <-
            seq(from = a,
                to = b,
                length.out = N)[col.keep]
        } else {
          time.for.plots[[j]] <- universal.time[col.keep]
        }
        
        sd.traj.j <- apply(data.cluster.j, 2, sd.na.rm)
        
        df <- length(which.j) - 1
        
        if(df > 0){
          CI.low[[j]] <-
            mean.traj.j - qt(p = 0.025,
                             df = df,
                             lower.tail = FALSE) * sd.traj.j / sqrt(length(which.j))
          
          CI.high[[j]] <-
            mean.traj.j + qt(p = 0.025,
                             df = df,
                             lower.tail = FALSE) * sd.traj.j / sqrt(length(which.j)) 
        } else {
          CI.low[[j]] <- CI.high[[j]] <- mean.traj.j
        }
      }
      
      if (spline == TRUE) {
        universal.time <- seq(from = a,
                              to = b,
                              length.out = N)
        main.med <-
          paste("Pointwise median \n of smoothed trajectories", sep = "")
        main.mean <-
          paste("Pointwise mean \n of smoothed trajectories", sep = "")
        main.mean2 <-
          paste("Pointwise mean of smoothed \n trajectory and 95% CI" , sep =
                  "")
      } else {
        main.med <- paste("Pointwise median trajectories", sep = "")
        main.mean <- paste("Pointwise mean trajectories", sep = "")
        main.mean2 <-
          paste("Pointwise mean trajectory \n and 95% CI" , sep = "")
      }
      
      # Plot all the median trajectories together
      if (1 %in% which.med.mean) {
        par(mfrow = c(1, 1))
        
        plot(
          x = 0,
          y = 0,
          xlim = c(min(universal.time), max(universal.time)),
          ylim = c(min(unlist(median.traj)), max(unlist(median.traj))),
          type = "n",
          xlab = "",
          ylab = "",
          main = main.med
        )
        legend(
          "topright",
          col = color.pal[1:x$nclusters],
          legend = paste(1:x$nclusters),
          lty = rep(1, x$nclusters)
        )
        
        for (j in 1:x$nclusters) {
          lines(
            x = time.for.plots[[j]],
            y = median.traj[[j]],
            type = "l",
            col = color.pal[j]
          )
        }
      }
      
      # Plot median trajectories one by one
      if (2 %in% which.med.mean) {
        par(mfrow = good.grid)
        
        for (j in 1:x$nclusters) {
          plot(
            x = time.for.plots[[j]],
            y = median.traj[[j]],
            type = "l",
            xlab = paste(
              "Cluster ",
              j,
              " (n = ",
              x$partition.summary[j],
              ")",
              sep = ""
            ),
            ylab = "",
            main = main.med,
            col = color.pal[j]
          )
        }
      }
      
      # Plot all the mean trajectories together
      if (3 %in% which.med.mean) {
        par(mfrow = c(1, 1))
        
        plot(
          x = 0,
          y = 0,
          xlim = c(min(universal.time), max(universal.time)),
          ylim = c(min(unlist(mean.traj)), max(unlist(mean.traj))),
          type = "n",
          xlab = "",
          ylab = "",
          main = main.mean
        )
        legend(
          "topright",
          col = color.pal[1:x$nclusters],
          legend = paste(1:x$nclusters),
          lty = rep(1, x$nclusters)
        )
        
        
        for (j in 1:x$nclusters) {
          lines(
            x = time.for.plots[[j]],
            y = mean.traj[[j]],
            type = "l",
            col = color.pal[j]
          )
        }
      }
      
      # Plot mean curves one by one
      if (4 %in% which.med.mean) {
        par(mfrow = good.grid)
        
        for (j in 1:x$nclusters) {
          plot(
            x = 0,
            y = 0,
            xlim = c(
              min(time.for.plots[[j]], na.rm = TRUE),
              max(time.for.plots[[j]], na.rm = TRUE)
            ),
            ylim = c(min(CI.low[[j]], na.rm = TRUE), max(CI.high[[j]], na.rm = TRUE)),
            type = "n",
            xlab = paste(
              "Cluster ",
              j,
              " (n = ",
              x$partition.summary[j],
              ")",
              sep = ""
            ),
            ylab = "",
            main = main.mean2
          )
          
          w <- which(is.na(CI.high[[j]]))
          if (length(w) > 0) {
            polygon(
              c(rev(time.for.plots[[j]][-w]), time.for.plots[[j]][-w]),
              c(rev(CI.high[[j]][-w]), CI.low[[j]][-w]),
              col = 'grey80',
              border = NA
            )
          } else{
            polygon(
              c(rev(time.for.plots[[j]]), time.for.plots[[j]]),
              c(rev(CI.high[[j]]), CI.low[[j]]),
              col = 'grey80',
              border = NA
            )
          }
          
          lines(
            x = time.for.plots[[j]],
            y = mean.traj[[j]],
            col = color.pal[j],
            type = "l"
          )
        }
      }
    }
    
    if (!(spline | homogeneous.times)) {
      print(
        paste(
          "No pointwise median and pointwise mean trajectories were generated because the observation times are not the same for each trajectories. Set the 'spline' argument to TRUE to plot median and mean trajectories based on smoothing splines."
        )
      )
    }
    
    
    
    
    ### Randomly sampled trajectories from each clusters ###
    random.condition <- (is.null(which.plots) | 5 %in% which.plots)
    
    if (random.condition) {
      traj.by.clusters <- list()
      for (k in seq_len(x$nclusters)) {
        traj.by.clusters[[k]] <-
          x$data[which(x$partition[, 2] == k),-c(1), drop = FALSE]
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
            x = na.omit(smpl.time.by.clusters[[k]][i,]),
            y = na.omit(smpl.traj.by.clusters[[k]][i,]),
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
    print("See also 'plot.crit' for a plot of the statistic used to determined the number of clusters and see 'plot.scatter' for scatter plots of the measures involved in the clustering.")
  }

#'@rdname Step1Measures
#'
#'@export
plot.scatter <- function(x, ask = TRUE, ...) {
  
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
            pch = (x$partition$Cluster[s] - 1),
            col = color.pal[x$partition$Cluster[s]],
            bg = color.pal[x$partition$Cluster[s]]
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

#'@rdname Step1Measures
#'
#'@export
plot.crit <- function(x, ...) {
  
  if (!is.null(x$GAP)) {
    par(mfrow = c(1, 1))
    plot(x$GAP, main = "Gap statistic up to one SE")
  } else if (!is.null(x$CH)) {
    par(mfrow = c(1, 1))
    plot(x$CH, type = "b", xlab = "k", ylab = "Index", main = "Calinski-Harabasz index")
  } else {
    print("There are no GAP or Calinski-Harabasz statistics to plot.")
  }
}
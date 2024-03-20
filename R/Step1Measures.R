#'@title Compute Measures for Identifying Patterns of Change in Longitudinal
#'  Data
#'
#'@description \code{Step1Measures} computes up to 18 measures for each
#'  longitudinal trajectory. See Details for the list of measures.
#'@param Data a matrix or data frame in which each row contains the longitudinal
#'  data (trajectories).
#'@param Time either \code{NULL}, a vector or a matrix/data frame of the same dimension
#'  as \code{Data}. If a vector, matrix or data frame is supplied, its entries
#'  are assumed to be measured at the times of the corresponding cells in
#'  \code{Data}. When set to \code{NULL} (the default), the times are assumed
#'  equidistant.
#'@param ID logical. Set to \code{TRUE} if the first columns of \code{Data} and
#'  \code{Time} corresponds to an \code{ID} variable identifying the
#'  trajectories. Defaults to \code{FALSE}.
#'@param measures a vector containing the numerical identifiers of the measures
#'  to compute (see "Details" section below). The default, 1:17, corresponds to
#'  measures 1-17 and thus excludes the measures which require specifying a
#'  midpoint.
#'@param midpoint specifies which column of \code{Time} to use as the midpoint
#'  in measure 18. Can be \code{NULL}, an integer or a vector of integers of length
#'  the number of rows in \code{Time}. The default is \code{NULL}, in which case the
#'  midpoint is the time closest to the median of the Time vector specific to
#'  each trajectory.
#'@param cap.outliers logical. If \code{TRUE}, extreme values of the measures will be caped (see details below). Defaults to \code{TRUE}.
#'@param x object of class trajMeasures.
#'@param object object of class trajMeasures.
#'@param ... further arguments passed to or from other methods.
#'
#'@return An object of class \code{trajMeasures}; a list containing the values
#'  of the measures, a table of the outliers which have been capped, as well as
#'  a curated form of the function's arguments.
#'
#'@details Each trajectory must have a minimum of 3 observations otherwise it
#'  will be omitted from the analysis.
#'
#'  The 18 measures and their numerical identifiers are listed below. Please
#'  refer to the vignette for the specific formulas used to compute them.
#'\enumerate{
#'\item  Maximum\cr
#'\item  Range (max - min)\cr
#'\item  Mean value\cr
#'\item  Standard deviation\cr
#'\item  Slope of the linear model\cr
#'\item  \eqn{R^2}: Proportion of variance explained by the linear model\cr
#'\item  Curve length (total variation)\cr
#'\item  Rate of intersection with the mean\cr
#'\item  Proportion of time spent under the mean\cr
#'\item  Minimum  of the first derivative\cr
#'\item  Maximum of the first derivative\cr
#'\item  Mean of the first derivative\cr
#'\item  Standard deviation of the first derivative\cr
#'\item  Minimum  of the second derivative\cr
#'\item  Maximum of the second derivative\cr
#'\item  Mean of the second derivative\cr
#'\item  Standard deviation of the second derivative\cr
#'\item  Early change/Later change\cr
#'}
#'
#'  In the presence of highly correlated measures (Pearson correlation >
#'  0.98), the function selects the highest-ranking measure on the list (see
#'  \code{\link[traj]{Step1Measures}}) and discards the others. Because the
#'  K-means algorithm is sensitive to outliers, the measures are prevented from taking extreme or
#'  infinite values (caused by a possible division by 0 in m18). Nishiyama's improved Chebychev
#'  bound is used to determine extreme values for each measure, corresponding to
#'  a 0.3% probability threshold. Extreme values beyond the threshold are then capped
#'  to the 0.3% probability threshold. If applicable, the values of m18 which
#'  would be of the form 0/0 are set to 1. PCA is applied on the remaining
#'  measures using the \code{\link[psych]{principal}} function from the
#'  \code{psych} package.
#'
#'@importFrom stats complete.cases coefficients lm median sd density
#'
#'@references Leffondre K, Abrahamowicz M, Regeasse A, Hawker GA, Badley EM,
#'  McCusker J, Belzile E. Statistical measures were proposed for identifying
#'  longitudinal patterns of change in quantitative health indicators. J Clin
#'  Epidemiol. 2004 Oct;57(10):1049-62. doi: 10.1016/j.jclinepi.2004.02.012.
#'  PMID: 15528056.
#'
#'  Nishiyama T, Improved Chebyshev inequality: new probability bounds with
#'  known supremum of PDF, arXiv:1808.10770v2 stat.ME
#'  https://doi.org/10.48550/arXiv.1808.10770
#'
#'@examples
#'\dontrun{
#'data("trajdata")
#'trajdata.noGrp <- trajdata[, which(colnames(trajdata) == "Group")] #remove the Group column
#'
#'m1 = Step1Measures(trajdata.noGrp, ID = TRUE, measures = 18, midpoint = NULL)
#'m2 = Step1Measures(trajdata.noGrp, ID = TRUE, measures = 18, midpoint = 3)
#'
#'identical(m1$measures, m2$measures)
#'}
#'
#'@rdname Step1Measures
#'
#'@export
Step1Measures <-
  function (Data,
            Time = NULL,
            ID = FALSE,
            measures = c(1:17),
            midpoint = NULL,
            cap.outliers = TRUE) {
    
    ###############################################################
    ##         Perform checks and uniformization of data         ##
    ###############################################################
    
    if (is.null(Time)) {
      if (ID == TRUE) {
        Time <- c(1:(ncol(Data) - 1))
      }
      if (ID == FALSE) {
        Time <- c(1:(ncol(Data)))
      }
    }
    
    Time.is.vector <- is.vector(Time)
    
    data <- Data
    data <- data.frame(data)
    names(data) <- NULL
    
    time <- Time
    time <- data.frame(time)
    names(time) <- NULL
    
    if (ID == TRUE) {
      IDvector <- data[, 1]
      data <- data[,-1]
      if (Time.is.vector == FALSE) {
        ID.time <- time[, 1]
        time <- time[,-1]
        if (identical(IDvector, ID.time) == FALSE) {
          stop("ID vector in Data differs from ID vector in Time.")
        }
      }
    } else{
      IDvector <- seq_len(nrow(data))
    }
    
    if (identical(dim(time), dim(data))) {
      if (identical(is.na(data), is.na(time)) == FALSE) {
        stop(
          "Data has cells with no corresponding time or Time has cells with no corresponding data."
        )
      }
      
      data2 <- data
      time2 <- time
      for (i in seq_len(nrow(data))) {
        NA.str_i <- is.na(data)[i,]
        w <- unname(which(NA.str_i == FALSE))
        if (length(w) < 3) {
          data2 <- data2[-i,]
          time2 <- time2[-i,]
          warning (
            paste(
              "Row ",
              i,
              " of Data contains less than 3 observations; it has been removed.",
              sep = ""
            )
          )
          IDvector <- IDvector[-i]
        } else if (!identical(w, seq_len(length(w)))) {
          stop(
            paste(
              "Row ",
              i,
              " of Data is not formatted correctly. Rows should be of the form X Y ... Z NA ... NA.",
              sep = ""
            )
          )
        }
      }
      data <- data2
      time <- time2
      
      for (i in seq_len(nrow(data))) {
        w2 <- unname(which(!is.na(time)[i,]))
        non.NA <- unname(unlist(time[i, w2]))
        if (!length(unique(non.NA)) == length(non.NA)) {
          stop(
            paste(
              "Line ",
              i,
              " of Time contains duplicates. The rows of Time should be strictly increasing sequences.",
              sep = ""
            )
          )
        }
        if (!identical(w2, order(non.NA))) {
          stop(
            paste(
              "The elements in row ",
              i,
              " of Time are not ordered chronologically.",
              sep = ""
            )
          )
        }
      }
      
    } else if (is.vector(Time) & (length(Time) == ncol(data))) {
      time <- unname(unlist(Time))
      if (TRUE %in% is.na(time)) {
        stop("If Time is supplied as a vector, it must not contain NA.")
      }
      if (!length(unique(time)) == length(time)) {
        stop(
          "The Time vector contains duplicates. The elements of Time should form a strictly increasing sequence."
        )
      }
      if (!identical(order(time), seq_len(length(time)))) {
        stop("The elements of Time are not ordered chronologically.")
      }
      if (length(time) < 3) {
        stop("Time must contain at least 3 elements.")
      }
      
      data2 <- data
      for (i in seq_len(nrow(data))) {
        NA.str_i <- is.na(data)[i,]
        w <- unname(which(NA.str_i == FALSE))
        if (length(w) < 3) {
          data2 <- data2[-i,]
          warning(
            paste(
              "Row ",
              i,
              " of Data contains less than 3 observations; it has been removed.",
              sep = ""
            )
          )
          if (ID == TRUE) {
            IDvector <- IDvector[-i]
          }
        }
      }
      data <- data2
      
      data2 <-
        unname(data.frame(matrix(
          NA, ncol = ncol(data), nrow = nrow(data)
        )))
      time2 <-
        unname(data.frame(matrix(
          NA, ncol = ncol(data), nrow = nrow(data)
        )))
      
      for (i in seq_len(nrow(data))) {
        w <- which(!is.na(data[i,]))
        data2[i, seq_len(length(w))] <- data[i, w]
        time2[i, seq_len(length(w))] <- Time[w]
      }
      data <- data2
      time <- time2
      
    } else{
      stop("The dimension of Data and Time are incompatible.")
    }
    
    
    ######################################################
    ##         Construct vector of "mid points"         ##
    ######################################################
    
    if (TRUE %in% (18 %in% measures)) {
      mid.position <- c()
      
      if (is.null(midpoint)) {
        median.time <- (max(time, na.rm = TRUE) + min(time, na.rm = TRUE)) / 2
        flag <- c()
        for (i in seq_len(nrow(data))) {
          v <- time[i, ][!is.na(time[i, ])]
          w <- which.min(abs(v - median.time))
          mid.position[i] <- w
          if (w == length(v) | w == 1) {
            flag <-
              c(flag, i)  #  the mid point can't be the first or the last point, so if this occurs, the ith row gets "flagged"
          }
        }
        if (length(flag) > 0) {
          mid.position <- mid.position[-flag]
          data <- data[-flag, ]
          time <- time[-flag, ]
          if (ID == TRUE) {
            Lines <- which(Data[, 1] %in% IDvector[flag])
          } else{
            Lines <- IDvector[flag]
          }
          if (length(Lines) == 1) {
            warning(
              paste(
                "When left blank, the 'midpoint' argument defaults to the observation time closests to the median time 0.5*(max(Time)-min(Time)), but this can't be either the first or last observation time. As a result, row ",
                Lines,
                " has been removed. To avoid this, consider excluding measure 18 from the analysis or providing custom 'midpoint' values.",
                sep = ""
              )
            )
          }
          if (length(Lines) > 1) {
            warning(
              paste(
                "When left blank, the 'midpoint' argument defaults to the observation time closests to the median time 0.5*(max(Time)-min(Time)), but this can't be either the first or last observation time. As a result, rows ",
                noquote(paste(Lines, collapse = ", ")),
                " have been removed. To avoid this, consider excluding measure 18 from the analysis or providing custom 'midpoint' values.",
                sep = ""
              )
            )
          }
          
          IDvector <- IDvector[-flag]
          
        }
      } else if (!is.vector(midpoint)) {
        stop(
          "'midpoint' should be either NULL or an integer/vector of integers corresponding to a column of Time."
        )
      } else{
        if (length(midpoint) == nrow(data)) {
          mid.position <- midpoint
        } else if (is.vector(Time) &
                   (length(Time) == ncol(data)) & length(midpoint) == 1) {
          mid.position <- rep(midpoint, nrow(data))
        } else{
          stop("'midpoint' does not have the correct format.")
        }
      }
      # Check to see if the midpoints are all greater than 1 but less than the number of observations
      for (i in seq_len(nrow(data))) {
        v <- time[i,][!is.na(time[i, ])]
        if (mid.position[i] <= 1) {
          stop(
            paste(
              "Error in 'midpoint' for subject ",
              i,
              "; 'midpoint' must be greater than 1.",
              sep = ""
            )
          )
        } else if (mid.position[i] >= length(v)) {
          stop(
            paste(
              "Error in 'midpoint' for subject ",
              i,
              "; 'midpoint' must be less than the number of observations.",
              sep = ""
            )
          )
        }
      }
    } else{
      mid.position <- NULL
    }
    
    ##########################################################################################
    ##         Initialize main output data frame and compute the requested measures         ##
    ##########################################################################################
    
    output <-
      data.frame(matrix(ncol = 1 + length(measures), nrow = nrow(data)))
    colnames(output) <- c("ID", paste("m", measures, sep = ""))
    output$ID <- IDvector
    
    data <- as.matrix(data)
    time <- as.matrix(time)
    
    
    #max
    if (1 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m1[i] <-
          max(data[i, ], na.rm = TRUE)
      }
    }
    
    # range 
    if (2 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m2[i] <-
          max(data[i, ], na.rm = TRUE) - min(data[i, ], na.rm = TRUE)
      }
    }
    
    #mean
    if (sum(c(3, 4, 8, 9) %in% measures) > 0) {
      m3 <- c()
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i, ])]
        x <- time[i, complete.cases(time[i, ])]
        m3[i] <- FctMean(x, y)
      }
      if (3 %in% measures) {
        output$m3 <- m3
      }
    }
    
    #SD
    if (4 %in% measures) {
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i, ])]
        x <- time[i, complete.cases(time[i, ])]
        output$m4[i] <- sqrt(FctMean(x, (y - m3[i]) ^ 2))
      }
    }
    
    #slope of linear model
    if (5 %in% measures) {
      for (i in seq_len(nrow(data))) {
        b <- coefficients(lm(data[i, ] ~ time[i, ]))
        output$m5[i] <- b[2]
      }
    }
    
    # Turn off warnings coming from summary(model) when the residuals are 0, which can happen when the trajectory is a straight line.
    defaultW <- getOption("warn")
    options(warn = -1)
    
    #R^2
    if (6 %in% measures) {
      for (i in seq_len(nrow(data))) {
        model <- lm(data[i, ] ~ time[i, ])
        output$m6[i] <- summary(model)$r.squared
      }
    }
    options(warn = defaultW)
    
    #curve length:
    if (7 %in% measures) { 
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i, ])]
        x <- time[i, complete.cases(time[i, ])]
        Length <- 0
        for(k in seq_along(y)[-length(y)]){
          Length <- Length + sqrt( (x[k+1] - x[k])^2 + (y[k+1] - y[k])^2 )
        }
        output$m7[i] <- Length
      }
    }
    
    #rate of intersections with the mean
    
    if (8 %in% measures) { 
      intersection.count <- c()
      for (i in seq_len(nrow(data))) {
        norm <- data[i, ] - m3[i]
        y <- norm[complete.cases(norm)]
        x <- time[i, complete.cases(time[i, ])]
        
        intersection.count[i] <- 0
        for(k in  seq_along(y)[-length(y)]){
          w <- which(sign(y[k] * y[(k+1):length(y)]) !=0)
          if(length(w) > 0){
            if(sign(y[k] * y[(k+1):length(y)])[w[1]] == -1){
              intersection.count[i] <- intersection.count[i] + 1 
            }
          }
        }
        output$m8[i] <- intersection.count[i] /( max(x) - min(x) )
      }
    }

    #fraction of time spent under the mean (if this is large, then there are big sharp spikes)
    if (9 %in% measures) {
      for (i in seq_len(nrow(data))) {
        norm <- data[i, ] - m3[i]
        y <- norm[complete.cases(norm)]
        x <- time[i, complete.cases(time[i, ])]
        
        blue.time <- 0
        wb <- which(y > 0)
        if(length(wb > 0)){
          for (k in wb) {
            if(k == 1){
              blue.time <- blue.time + 0.5*(x[k+1] - x[k])
            }
            if (! (k %in% c(1, length(x)))) {
              blue.time <- blue.time + 0.5*(x[k+1] - x[k-1])
            }
            if(k == length(x)){
              blue.time <- blue.time + 0.5*(x[k] - x[k-1])
            }
          }
        }
        
        red.time <- 0
        wr <- which(y < 0)
        if(length(wr > 0)){
          for (k in wr) {
            if(k == 1){
              red.time <- red.time + 0.5*(x[k+1] - x[k])
            }
            if (! (k %in% c(1, length(x)))) {
              red.time <- red.time + 0.5*(x[k+1] - x[k-1])
            }
            if(k == length(x)){
              red.time <- red.time + 0.5*(x[k] - x[k-1])
            }
          }
        }
        output$m9[i] <- blue.time/(blue.time + red.time)
      }
    }
    
    
    ### Measures on y'(t) ###
    
    # If a measure involving the speed or the acceleration was requested, compute the derivative:
    if (sum(measures %in% 10:17) > 0) {
      speed.data <- data
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i, ])]
        x <- time[i, complete.cases(time[i, ])]
        speed.data[i, complete.cases(speed.data[i, ])] <- Der(x, y)
      }
    }
    
    #min
    if (10 %in% measures) {
      for (i in seq_len(nrow(speed.data))) {
        output$m10[i] <-
          min(speed.data[i, ], na.rm = TRUE)
      }
    }
    
    #max
    if (11 %in% measures) {
      for (i in seq_len(nrow(speed.data))) {
        output$m11[i] <-
          max(speed.data[i, ], na.rm = TRUE)
      }
    }
    
    #mean 
    if (sum(c(12, 13) %in% measures) > 0) {
      m12 <- c()
      for (i in seq_len(nrow(speed.data))) {
        y <- speed.data[i, complete.cases(speed.data[i, ])]
        x <- time[i, complete.cases(time[i, ])]
        m12[i] <- FctMean(x, y)
      }
      if (12 %in% measures) {
        output$m12 <- m12
      }
    }
    
    #SD
    if (13 %in% measures) {
      for (i in seq_len(nrow(speed.data))) {
        y <- speed.data[i, complete.cases(speed.data[i, ])]
        x <- time[i, complete.cases(time[i, ])]
        output$m13[i] <- sqrt(FctMean(x, (y - m12[i]) ^ 2))
      }
    }
    
    ### Measures on y''(t) ###
    
    #if a measure involving the acceleration was requested, compute the derivative of the derivative:
    
    
    
    if (sum(measures %in% 14:17) > 0) {
      accel.data <- speed.data
      for (i in seq_len(nrow(accel.data))) {
        y <- speed.data[i, complete.cases(speed.data[i, ])]
        x <- time[i, complete.cases(time[i, ])]
        accel.data[i, complete.cases(accel.data[i, ])] <- Der(x, y)
      }
    }
    
    #min
    if (14 %in% measures) {
      for (i in seq_len(nrow(accel.data))) {
        output$m14[i] <-
          min(accel.data[i, ], na.rm = TRUE)
      }
    }
    
    #max
    if (15 %in% measures) {
      for (i in seq_len(nrow(accel.data))) {
        output$m15[i] <-
          max(accel.data[i, ], na.rm = TRUE)
      }
    }
    
    #mean 
    if (sum(c(16, 17) %in% measures) > 0) {
      m16 <- c()
      for (i in seq_len(nrow(accel.data))) {
        y <- accel.data[i, complete.cases(accel.data[i, ])]
        x <- time[i, complete.cases(time[i, ])]
        m16[i] <- FctMean(x, y)
      }
      if (16 %in% measures) {
        output$m16 <- m16
      }
    }
    
    #SD
    if (17 %in% measures) {
      for (i in seq_len(nrow(accel.data))) {
        y <- accel.data[i, complete.cases(accel.data[i, ])]
        x <- time[i, complete.cases(time[i, ])]
        output$m17[i] <- sqrt(FctMean(x, (y - m16[i]) ^ 2))
      }
    }
    
    #later change/early change
    if (18 %in% measures) {
      for (i in seq_len(nrow(data))) {
        early <- Last(data[i, 1:mid.position[i]]) - First(data[i, 1:mid.position[i]])
        later <- Last(data[i, mid.position[i]:ncol(data)]) - First(data[i, mid.position[i]:ncol(data)])
        output$m18[i] <- later/early
      }
    }

    output[is.na(output)] <- 1 #  NAs correspond to quotient measures of the form 0/0
    
    
    ######################################
    ##         Cap the outliers         ##
    ######################################
    
    if (cap.outliers == FALSE) {
      outliers <- NULL
    } 
    
    if (cap.outliers == TRUE){
      outliers <-
        data.frame(matrix(nrow = nrow(output), ncol = ncol(output)))
      colnames(outliers) <- colnames(output)
      outliers[, 1] <- output$ID
      
      for (j in 2:ncol(output)) {
        y <- output[, j]
        y.TRUE <- y
        n <- length(y)
        which.inf <- which(is.infinite(y.TRUE))
        
        if (length(which.inf) > 0) {
          y.TRUE <- y.TRUE[-which.inf]
        }
        
        if (length(y.TRUE) > 2) {
          top <-
            rev(order(abs(y.TRUE - median(y.TRUE))))[1:ceiling(n * 0.01)] #  if n < 100, remove 1 element, so this is never empty
          y.TRUE <- y.TRUE[-top]
        }
        
        mu <- mean(y.TRUE)
        sigma <- sd(y.TRUE)
        
        k_Cheb <-
          sqrt(100 / 0.3) #  The classical Chebychev bound. Approximately 18.26.
        k <- seq(from = 0.1, to = 18.26, by = 0.1)
        M <- c()
        for (i in seq(length(k))) {
          max.left <-
            max(density(
              y.TRUE,
              from = (mu - 18.3 * sigma),
              to = (mu + 18.3 * sigma)
            )$y[density(y.TRUE,
                        from = (mu - 18.3 * sigma),
                        to = (mu + 18.3 * sigma))$x > mu + k[i] * sigma])
          max.right <-
            max(density(
              y.TRUE,
              from = (mu - 18.3 * sigma),
              to = (mu + 18.3 * sigma)
            )$y[density(y.TRUE,
                        from = (mu - 18.3 * sigma),
                        to = (mu + 18.3 * sigma))$x < mu - k[i] * sigma])
          M[i] <- max(max.left, max.right)
        }
        
        p <- 2 * pi * k ^ 2 * (exp(1) - 2 * pi / 3)
        q <-
          2 * (2 * pi * k) ^ 3 / 27 - (2 * pi) ^ 2 * k ^ 3 * exp(1) / 3 - 2 * pi * exp(1) /
          (sigma * M)
        
        root <-
          CubeRoot(-q / 2 + sqrt(q ^ 2 / 4 + p ^ 3 / 27)) + CubeRoot(-q / 2 - sqrt(q ^
                                                                                     2 / 4 + p ^ 3 / 27)) + 2 * pi * k / 3
        
        w <- which(root * sigma * M < 0.3 / 100)
        if (length(w) > 0) {
          k.opt <- min(k_Cheb, k[w[1]])
        } else{
          k.opt <- k_Cheb
        }
        
        cap <- which(abs(y - mu) > k.opt * sigma)
        if (length(cap) > 0) {
          outliers[cap, j] <- signif(y[cap], 3)
          
          y[cap] <- mu + sign(y[cap]) * k.opt * sigma
          output[, j] <- y
        }
      }
      
      row.rm <- which(rowSums(!is.na(outliers[, -c(1), drop = FALSE])) == 0)
      outliers <- outliers[-row.rm, , drop = FALSE]
      col.kp <- which(colSums(!is.na(outliers)) != 0)
      outliers <- outliers[, col.kp, drop = FALSE]
    }

    ID <- IDvector
    
    trajMeasures <-
      structure(
        list(
          measures = output,
          outliers = outliers,
          mid = mid.position,
          data = cbind(ID, data),
          time = cbind(ID, time)
        ),
        class = "trajMeasures"
      )
    return(trajMeasures)
    
  }

#'@rdname Step1Measures
#'
#'@export
print.trajMeasures <- function(x, ...) {
  cat("Measures and time at midpoint (if applicable):\n")
  
  if (is.null(x$mid)) {
    print(x$measures)
  } else{
    measure.plus <- cbind(x$measures$ID, x$mid)
    measure.plus <- cbind(measure.plus, x$measures[, -1, drop = FALSE])
    colnames(measure.plus)[1:2] <- c("ID", "mid time")
    print(measure.plus, row.names = FALSE)
  }
}

#'@rdname Step1Measures
#'
#'@export
summary.trajMeasures <- function(object, ...) {
  
  cat("Description of the measures:\n")
  cat("m1: Maximum\n")
  cat("m2: Range\n")
  cat("m3: Mean value\n")
  cat("m4: Standard deviation\n")
  cat("m5: Slope of the linear model\n")
  cat("m6: Proportion of variance explained by the linear model (R squared)\n")
  cat("m7: Curve length (total variation)\n")
  cat("m8: Number of times crossing the mean per unit time\n")
  cat("m9: Proportion of time spent under the mean\n")
  cat("m10: Minimum of the first derivative\n")
  cat("m11: Maximum of the first derivative\n")
  cat("m12: Mean of the first derivative\n")
  cat("m13: Standard deviation of the first derivative\n")
  cat("m14: Minimum of the second derivative\n")
  cat("m15: Maximum of the second derivative\n")
  cat("m16: Mean of the second derivative\n")
  cat("m17: Standard deviation of the second derivative\n")
  cat("m18: Early change/Later change\n")
  
  cat("\n")
  
  cat("Summary of measures:\n")
  
  Q1 <- function(x) {
    return(quantile(x , probs = c(.25)))
  }
  
  Q2 <- function(x) {
    return(quantile(x , probs = c(.5)))
  }
  
  Q3 <- function(x) {
    return(quantile(x , probs = c(.75)))
  }
  
  measures.summary <-
    data.frame(matrix(nrow = 6, ncol = ncol(object$measures) - 1))
  rownames(measures.summary) <-
    c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  colnames(measures.summary) <- colnames(object$measures)[-1]
  
  measures.summary[1,] <- apply(object$measures, 2, min)[-1]
  measures.summary[2,] <- apply(object$measures, 2, Q1)[-1]
  measures.summary[3,] <- apply(object$measures, 2, Q2)[-1]
  measures.summary[4,] <- apply(object$measures, 2, mean)[-1]
  measures.summary[5,] <- apply(object$measures, 2, Q3)[-1]
  measures.summary[6,] <- apply(object$measures, 2, max)[-1]
  
  print(measures.summary)
  
  
  if(!is.null(object$outliers)){
    cat("\n")
    
    cat("Outliers before capping:\n")
    outliers <- object$outliers
    outliers.pre <- outliers
    outliers.pre[is.na(outliers.pre)] <- ""
    print(outliers.pre, row.names = FALSE)
    
    cat("Outliers after capping:\n")
    if (nrow(outliers) == 0) {
      print(outliers, row.names = FALSE)
    } else{
      outliers.post <- outliers
      for (j in seq_len(nrow(outliers))) {
        for (k in 2:ncol(outliers)) {
          if (is.na(outliers[j, k]) == FALSE) {
            outliers.post[j, k] <-
              signif(object$measures[object$measures$ID == outliers$ID[j], colnames(outliers)[k]], 3)
          }
        }
      }
      outliers.post[is.na(outliers.post)] <- ""
      print(outliers.post, row.names = FALSE)
    }
  }
}

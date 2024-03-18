#'@title Select a Subset of the Measures Using Factor Analysis
#'
#'@description This function applies the following dimension reduction algorithm
#'  to the measures computed by \code{\link[traj]{Step1Measures}}:
#' \enumerate{
#'   \item Drop the measures whose values are constant across the trajectories;
#'   \item Use principal component analysis (PCA) on the measures to form factors summarizing the variability in the measures;
#'   \item Drop the factors whose variance is smaller than any one of the standardized measures;
#'   \item Perform a varimax rotation on the remaining factors;
#'   \item For each rotated factor, select the measure that has the highest correlation (aka factor loading) with it and that hasn't yet been selected;
#'   \item Drop the remaining measures.
#' }
#'
#'@param trajMeasures object of class \code{trajMeasures} as returned by
#'  \code{\link[traj]{Step1Measures}}.
#'@param num.select an optional positive integer indicating the number of
#'  factors to keep in the second stage of the algorithm. Defaults to \code{NULL} so
#'  that all factors with variance greater than any one of the normalized
#'  measures are selected.
#'@param discard an optional vector of positive integers corresponding to the
#'  measures to be dropped from the analysis. See
#'  \code{\link[traj]{Step1Measures}} for the list of measures. Defaults to
#'  \code{NULL}.
#'@param select an optional vector of positive integers corresponding to the
#'  measures to forcefully select. Defaults to \code{NULL}. If a vector is supplied,
#'  the five-steps selection algorithm described above is bypassed and the
#'  corresponding measures are selected instead.
#'@param x object of class \code{trajSelection}.
#'@param object object of class  \code{trajSelection}.
#'@param ... further arguments passed to or from other methods.
#'
#'@return An object of class \code{trajSelection}; a list containing the values
#'  of the selected measures, the output of the principal component analysis as
#'  well as a curated form of the arguments.
#'
#'@importFrom psych principal
#'@importFrom stats cor
#'
#'@references Leffondre K, Abrahamowicz M, Regeasse A, Hawker GA, Badley EM,
#'  McCusker J, Belzile E. Statistical measures were proposed for identifying
#'  longitudinal patterns of change in quantitative health indicators. J Clin
#'  Epidemiol. 2004 Oct;57(10):1049-62. doi: 10.1016/j.jclinepi.2004.02.012.
#'  PMID: 15528056.
#'
#' @examples
#' \dontrun{
#'data("trajdata")
#'trajdata.noGrp <- trajdata[, which(colnames(trajdata) == "Group")] #remove the Group column
#'
#'m = Step1Measures(trajdata.noGrp, measure = c(1:18), ID = TRUE)
#'s = Step2Selection(m)
#'
#'print(s)
#'
#'s2 = Step2Selection(m, select = c(13, 3, 12, 9))
#'}
#'
#'
#'@seealso \code{\link[psych]{principal}} \code{\link[traj]{Step1Measures}}
#'
#'@rdname Step2Selection
#'
#'@export
Step2Selection <-
  function (trajMeasures,
            num.select = NULL,
            discard = NULL,
            select = NULL
            ) {
    input <- list(num.select, discard, select)
    names(input) <- c("num.select", "discard", "select")
    
    if ((!is.null(select)) & (!is.numeric(select))) {
      stop("Argument 'select' must be either NULL or a numerical vector.")
    }
    
    data <- data.frame(trajMeasures$measures)
    ID <- data[, 1]
    data <- data.frame(data[, -1])
    
    if (!is.null(select)) {
      m.select <- paste("m", select, sep = "")
      if (FALSE %in% (m.select %in% colnames(data))) {
        stop("The 'select' argument must only contain measures included in Step1Measures.")
      } else {
        output <- cbind(ID, data[, m.select, drop = FALSE])
        #colnames(output) <- c("ID", paste("m", select, sep = ""))
      }
      
      trajSelection <-
        structure(
          list(
            selection = output,
            PC = NULL,
            RC = NULL,
            constant.measures = NULL,
            correlated.measures = NULL,
            measures = trajMeasures$measures,
            data = trajMeasures$data,
            time = trajMeasures$time,
            input = input
          ),
          class = "trajSelection"
        )
      
    } else {
      if ((!is.null(discard)) & (!is.numeric(discard))) {
        stop("Argument 'discard' must be either NULL or a numerical vector.")
      }
      if (!is.null(discard)) {
        mes.to.discard <- paste("m", discard, sep = "")
        if (FALSE %in% (mes.to.discard %in% colnames(data))) {
          stop("Can't discard a measure which was not included in Step1Measures.")
        }
        w <- which(colnames(data) %in% mes.to.discard)
        data <- data[, -w, drop = FALSE]
      }
      
      if (!is.null(num.select)) {
        if (!is.numeric(num.select)) {
          stop("Argument 'num.select' must be a numerical vector of length 1.")
        }
        if (!is.vector(num.select)) {
          stop("Argument 'num.select' must be a numerical vector of length 1.")
        }
        if (!(length(num.select) == 1)) {
          stop("Argument 'num.select' must be a numerical vector of length 1.")
        }
        if (!is.null(discard) & (num.select > ncol(data))) {
          stop(
            "After discarding the measures specified in 'discard', the requested number 'num.select' of measures to retain exceeds the number of measures available."
          )
        }
        if (is.null(discard) & (num.select > ncol(data))) {
          stop(
            "The requested number 'num.select' of measures to retain exceeds the number of measures included in Step1Measures."
          )
        }
      }
      
      # Remove the measures that are constant because (1) these are not useful for
      # discriminating between the trajectories and (2) they cause problem with the
      # CheckCorrelation function later because their variance is 0 so division by 0
      # occurs when computing correlation.
      flag <- c()
      cst.measures <- NULL
      for (j in seq_len(ncol(data))) {
        col <- data[, j]
        if (max(col) == min(col)) {
          flag <- c(flag, j)
        }
      }
      if (length(flag) > 0) { 
        cst.measures <- colnames(data)[flag]
        data <- data[, -flag, drop = FALSE]    
      }
      if (ncol(data) == 0) {
        stop("All the measures are constant.")
      }
      
      corr.vars <-
        CheckCorrelation(data, verbose = FALSE, is.return = TRUE)
      
      if (!is.null(corr.vars)) {
        correlated.measures <- corr.vars
        corr.vars.pos <- which(names(data) %in% corr.vars[, 1])
        data <- data[, -corr.vars.pos]
      } else{
        correlated.measures <- NULL
      }
      
      if (num.select > ncol(data) && !is.null(num.select)) {
        stop(
          "After discarding the perfectly or almost perfectly correlated measures, there are ",
          ncol(data),
          " measures left, which is less than the number 'num.select' of requested measures to retain."
        )
      }
      
      Z <- scale(x = data,
                 center = TRUE,
                 scale = TRUE)
      
      if (is.null(num.select)) {
        eigen.values <-
          psych::principal(Z, rotate = "none", nfactors = ncol(data))$values
        num.select <- max(1, length(which(eigen.values > 1)))
      }
      
        PC <- psych::principal(Z, rotate = "none", nfactors = ncol(data))
        RC <-
          psych::principal(Z, rotate = "varimax", nfactors = num.select)
        principal.factors <- RC
        
        principal.variables <- c()
        
        if (is.null(select)) {
          for (j in 1:num.select) {
            if (j == 1) {
              aux <- principal.factors$loadings
            }
            if (j > 1) {
              w <-
                which(row.names(principal.factors$loadings) %in% principal.variables)
              aux <- principal.factors$loadings[-w, ]
            }
            principal.variables[j] <- rownames(aux)[which.max(abs(aux[, j, drop = FALSE]))]
          }
          output <- cbind(ID, data[, principal.variables, drop = FALSE])
        }
      
      trajSelection <-
        structure(
          list(
            selection = output,
            PC = PC,
            RC = RC,
            constant.measures = cst.measures,
            correlated.measures = correlated.measures,
            measures = trajMeasures$measures,
            data = trajMeasures$data,
            time = trajMeasures$time,
            input = input
          ),
          class = "trajSelection"
        )
    }
    
    return(trajSelection)
    
  }


#'@rdname Step2Selection
#'
#'@export
print.trajSelection <- function(x, ...) {
  
  if(length(x$constant.measures) == 1){
    cat(
      paste(
        paste(x$constant.measures, collapse = ", "),
        " has been discarded due to being constant.\n",
        sep = ""
      )
    )
    cat("\n")
    }
  
  if(length(x$constant.measures) > 1){
    cat(
      paste(
        paste(x$constant.measures, collapse = ", "),
        " have been discarded due to being constant.\n",
        sep = ""
      )
    )
    cat("\n")
  }
  
  if (length(x$correlated.measures) > 0) {
    for (i in 1:nrow(x$correlated.measures)) {
      cat(
        paste(
          x$correlated.measures[i, 1],
          " has been discarded due to being perfectly or almost perfectly correlated with ",
          x$correlated.measures[i, 2],
          ".\n",
          sep = ""
        )
      ) 
    }
    cat("\n")
  }
  
  if (!is.null(x$input$select)) {
    if(length(colnames(x$selection)[-1]) > 1){
      cat(
        paste(
          "The selected measures are ",
          paste(
            colnames(x$selection)[-1],
            collapse = ', ',
            sep = ""
          ),
          ".",
          sep = ""
        )
      )
    }
    if(length(colnames(x$selection)[-1]) == 1){
      cat(
        paste(
          "The selected measure is ",
          paste(
            colnames(x$selection)[-1],
            collapse = ', ',
            sep = ""
          ),
          ".",
          sep = ""
        )
      )
    }

    
  } else {
    
    if (length(colnames(x$selection)[-1]) > 1) {
      cat(
        paste(
          "In decreasing order of variance explained, the selected measures are ",
          paste(
            colnames(x$selection)[-1],
            collapse = ', ',
            sep = ""
          ),
          ".",
          sep = ""
        )
      )
      cat("\n")
    } else {
      cat(
        paste(
          "The selected measure is ",
          paste(
            colnames(x$selection)[-1],
            collapse = ', ',
            sep = ""
          ),
          ".",
          sep = ""
        )
      )
      cat("\n")
    }
  }

  
  if (!is.null(x$RC)) {
    print(x$RC$loadings)
  }
}

#'@rdname Step2Selection
#'
#'@export
summary.trajSelection <- function(object, ...) {
  if (!is.null(object$input$select)) {
    if(length(object$input$select) == 1){
      cat(paste(
        "The measure ",
        paste(
          colnames(object$selection)[-1],
          collapse = ', ',
          sep = ""
        ),
        " was selected.",
        sep = ""
      ))
    }
    if(length(object$input$select) > 1){
      cat(paste(
        "The measures ",
        paste(
          colnames(object$selection)[-1],
          collapse = ', ',
          sep = ""
        ),
        " were selected.",
        sep = ""
      ))
    }
} else{
    if (!is.null(object$input$num.select)) {
      if (!is.null(object$correlated.measures)) {
        dropped <- unique(object$correlated.measures[, 1])
        cat(
          paste(
            "The measures ",
            paste(dropped, collapse = ", "),
            " were discarded because they were perfectly or almost perfectly correlated with another measure. Upon forming the principal components from the remaining measures, the ",
            object$input$num.select,
            " factors that held the most variance were retained. Together, they explained ",
            round(
              100 * sum(object$PC$values[1:(length(colnames(object$selection)) - 1)]) /
                length(object$PC$values),
              1
            ) ,
            "% of the total variance. A varimax rotation was performed to maximize the correlation with the original measures without affecting the proportion of explained variance. For each rotated factor, the measure that had the highest correlation (loading) with it was selected. As a result of this procedure, the selected measures are, in decreasing order of variance explained, ",
            paste(
              colnames(object$selection)[-1],
              collapse = ', ',
              sep = ""
            ),
            ". Use print() to see more detailed informations.",
            sep = ""
          )
        )
      } else{
        cat(
          paste(
            "Upon forming the principal components from the measures, the ",
            object$input$num.select,
            " factors that held the most variance were retained. Together, they explained ",
            round(
              100 * sum(object$PC$values[1:(length(colnames(object$selection)) - 1)]) /
                length(object$PC$values),
              1
            ) ,
            "% of the total variance. A varimax rotation was performed to maximize the correlation with the original measures without affecting the proportion of explained variance. For each rotated factor, the measure that had the highest correlation (loading) with it was selected. As a result of this procedure, the selected measures are, in decreasing order of variance explained, ",
            paste(
              colnames(object$selection)[-1],
              collapse = ', ',
              sep = ""
            ),
            ". Use print() to see more detailed informations.",
            sep = ""
          )
        )
      }
    } else{
      if (!is.null(object$correlated.measures)) {
        dropped <- unique(object$correlated.measures[, 1])
        cat(
          paste(
            "The measures ",
            paste(dropped, collapse = ", "),
            " were discarded because they were perfectly or almost perfectly correlated with another measure. Upon forming the principal components from the remaining measures, ",
            length(colnames(object$selection)) - 1,
            " of them had a variance greater than any one of the normalized measures. Together, they explained ",
            round(
              100 * sum(object$PC$values[1:(length(colnames(object$selection)) - 1)]) /
                length(object$PC$values),
              1
            ) ,
            "% of the total variance. A varimax rotation was performed to maximize the correlation with the original measures without affecting the proportion of explained variance. For each rotated factor, the measure that had the highest correlation (loading) with it was selected. As a result of this procedure, the selected measures are, in decreasing order of variance explained, ",
            paste(
              colnames(object$selection)[-1],
              collapse = ', ',
              sep = ""
            ),
            ". Use print() to see more detailed informations.",
            sep = ""
          )
        )
      } else{
        cat(
          paste(
            "Upon forming the principal components from the measures, ",
            length(colnames(object$selection)) - 1,
            " of them had a variance greater than any one of the normalized measures. Together, they explained ",
            round(
              100 * sum(object$PC$values[1:(length(colnames(object$selection)) - 1)]) /
                length(object$PC$values),
              1
            ) ,
            "% of the total variance. A varimax rotation was performed to maximize the correlation with the original measures without affecting the proportion of explained variance. For each rotated factor, the measure that had the highest correlation (loading) with it was selected. As a result of this procedure, the selected measures are, in decreasing order of variance explained, ",
            paste(
              colnames(object$selection)[-1],
              collapse = ', ',
              sep = ""
            ),
            ". Use print() to see more detailed informations.",
            sep = ""
          )
        )
      }
    }
  }
}

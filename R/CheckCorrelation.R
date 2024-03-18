CheckCorrelation <-
  function(output,
           verbose = TRUE,
           is.return = FALSE,
           tresh = 0.98) {
    cor.mat <- cor(output)
    
    mes.names <- names(output)
    
    is.corr <- FALSE
    
    corr.var <- NULL
    
    for (i_row in mes.names[-length(mes.names)]) {
      i_pos <- which(mes.names == i_row)
      res.names <- mes.names[(i_pos + 1):length(mes.names)]
      
      for (i_col in res.names) {
        if (abs(cor.mat[i_row, i_col]) > tresh) {
          corr.var <- rbind(corr.var, c(i_col, i_row))
          
          if (verbose) {
            print(
              paste(
                "Variables ",
                i_row,
                " and ",
                i_col,
                " are perfectly or almost perfectly correlated.",
                sep = ""
              )
            )
            is.corr <- TRUE
          }
        }
      }
    }
    if (!is.corr && verbose)
      print("None of the measures are perfectly or almost perfectly correlated.")
    if (is.return)
      return(corr.var)
  }

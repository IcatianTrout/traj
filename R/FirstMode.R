#'@title Function that returns the mode of a vector x. In case of a tie, it returns the first value of x that is a mode
#'
#'@param x vector of real numbers
#'
FirstMode <- function(x) {
y <- unique(x)
y[which.max(tabulate(match(x, y)))]
}
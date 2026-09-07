# Function that returns the mode of a vector x. In case of a tie, it returns the first value of x that is a mode

FirstMode <- function(x) {
y <- unique(x)
y[which.max(tabulate(match(x, y)))]
}
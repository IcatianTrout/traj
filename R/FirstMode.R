#returns the smallest mode
FirstMode <- function(x) {
y <- unique(x)
y[which.max(tabulate(match(x, y)))]
}
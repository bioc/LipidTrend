# Generalized logarithmic transformation with base 10
# Implementation based on MKmisc::glog10
.glog10 <- function(x) {
    log(x + sqrt(x^2 + 1), base=10) - log(2, base=10)
}

.countDistance <- function(dist_input) {
    dist_input <- as.matrix(dist_input)
    dim <- ncol(dist_input)
    x_distance <- min(dist(unique(dist_input[,1])))
    dist_input[,1] <- dist_input[,1]/x_distance
    if (dim == 2) {
        y_distance <- min(dist(unique(dist_input[,2])))
        dist_input[,2] <- dist_input[,2]/y_distance
    } else {
        y_distance <- NULL
    }
    return(list(
        normalized_matrix=dist_input, x_distance=x_distance,
        y_distance=y_distance, dimensions=dim))
}

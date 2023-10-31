#' Modified Z-score Function
#'
#' This function calculates and returns standardized (z-score like) values based on the
#' Median Absolute Deviation (MAD) method. It is a robust method to detect outliers in
#' data that might not be normally distributed. The function is a replication of the outliers()
#' that was once available in the spatialEco package.
#'
#' @param x A numeric vector for which outliers need to be detected.
#' @param s A constant used in the MAD calculation, default is 1.4826.
#' @return A numeric vector of standardized values. Outliers typically have values far from zero.
#' @examples
#' data <- c(1, 2, 3, 4, 5, 100)
#' outliers(data)
#' @export
modified_z <- function(x, s = 1.4826) {
  e <- (length(x) - 1) / sqrt(length(x))

  mad <- function (x, center=stats::median(x), constant=s,
                   low=FALSE, high=FALSE) {
    n <- length(x)
    constant * if ((low || high) && n%%2 == 0) {
      if (low && high)
        stop("'low' and 'high' cannot be both TRUE")
      n2 <- n%/%2 + as.integer(high)
      sort(abs(x - center), partial = n2)[n2]
    }
    else stats::median(abs(x - center))
  }

  return((0.6745 * (x - stats::median(x))) / mad(x))
}

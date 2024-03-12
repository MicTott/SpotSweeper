#' Plot Outliers for a Single Sample in a SingleCellExperiment
#'
#' This function generates a plot for a specified sample within a SingleCellExperiment
#' object, highlighting outliers based on a specified metric. The plot visualizes
#' the metric of interest and indicates outliers with a distinct color.
#'
#' @param spe A SpatialExperiment object containing the data to be plotted.
#' @param sample_id A character string specifying the column name in `colData(spe)`
#'   that contains unique sample identifiers. Default is "sample_id".
#' @param sample A character string or numeric value specifying the sample to be
#'   plotted. By default, it plots the first unique sample found in `spe$sample_id`.
#' @param metric A character string specifying the metric to be visualized
#'   in the plot. This metric should be a column name in `colData(spe)`.
#' @param outliers A character string specifying the column name in `colData(spe)`
#'   that indicates whether a data point is considered an outlier. Default is
#'   "local_outliers".
#' @param low_color A character string indicating the color to be used for the low end
#'   of the gradient scale. Default is "white".
#' @param high_color A character string indicating the color to be used for the high end
#'   of the gradient scale. Default is "black".
#' @param stroke A numeric value specifying the border thickness for outlier
#'   points. Default is 1.
#'
#' @return The function returns a plot object created by `make_escheR` and modified
#'   with additional layers for visualizing the specified metric and outliers. The
#'   plot is not explicitly printed by the function and should be printed by the caller.
#'
#' @examples
#' plot <- plotOutliers(spe)
#' print(plot)
#'
#' @export

plotOutliers <- function(spe, sample_id = "sample_id",
                        sample=unique(spe$sample_id)[1], metric="detected",
                        outliers="local_outliers",
                        low_color="white", high_color="black", stroke=1) {


  spe.subset <- spe[ ,colData(spe)[[sample_id]] == sample]

  make_escheR(spe.subset) |>
    add_fill(var = metric) |>
    add_ground(var = outliers, stroke = stroke) +
    scale_color_manual(
      name = "", # turn off legend name for ground_truth
      values = c(
        "TRUE" = "red2",
        "FALSE" = "transparent")
    ) +
    scale_fill_gradient(low =low_color,high =  high_color)

}

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
#'   NULL.
#' @param colors A character vector specifying the colors to be used for the
#'   gradient scale. If length is 2, the gradient will be a single color gradient.
#' @param stroke A numeric value specifying the border thickness for outlier
#'   points. Default is 1.
#'
#' @return The function returns a plot object created by `make_escheR` and modified
#'   with additional layers for visualizing the specified metric and outliers. The
#'   plot is not explicitly printed by the function and should be printed by the caller.
#'
#' @importFrom escheR make_escheR add_fill add_ground
#'
#' @examples
#' plot <- plotOutliers(spe)
#' plot
#'
#' @export

plotOutliers <- function(spe, sample_id = "sample_id",
                        sample=unique(spe$sample_id)[1], metric="detected",
                        outliers=NULL, point_size=2,
                        colors=c("white","black"), stroke=1) {


  # Subset the data to the specified sample
  spe.subset <- spe[ ,colData(spe)[[sample_id]] == sample]

  # Start building the plot
  p <- make_escheR(spe.subset) |>
    add_fill(var = metric, point_size = point_size)

  # Conditionally add outliers if they are not NULL
  if (!is.null(outliers)) {
    p <- p |> add_ground(var = outliers, stroke = stroke, point_size = point_size)
  }

  # Add title to the plot
  p <- p + ggtitle(paste0("Sample: ", sample))

  # Remove and replace scales to avoid warnings for re-coloring
  p$scales$scales <- list()

  # Conditionally add color scales based on the length of 'colors'
  if (length(colors) == 2) {
    p <- p + scale_fill_gradient(low = colors[1], high = colors[2]) +
      scale_y_reverse()

    # If outliers are not NULL, add a manual color scale for them
    if (!is.null(outliers)) {
      p <- p + scale_color_manual(
        name = "", # turn off legend name for ground
        values = c("TRUE" = "red", "FALSE" = "transparent")
      )
    }
  } else if (length(colors) > 2) {
    p <- p + scale_fill_gradient2(low = colors[1], mid = colors[2], high = colors[3]) +
      scale_y_reverse()

    # If outliers are not NULL, add a manual color scale for them
    if (!is.null(outliers)) {
      p <- p + scale_color_manual(
        name = "", # turn off legend name for ground
        values = c("TRUE" = "red", "FALSE" = "transparent")
      )
    }
  }

  return(p)
}

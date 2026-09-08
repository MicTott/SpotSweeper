#' Plot Outlier Metrics to PDF
#'
#' This function generates a PDF file containing plots for each sample in the
#' SpatialExperiment object, highlighting outliers based on specified metrics.
#' Each plot visualizes outlier metrics for a single sample, allowing for
#' easy comparison and analysis across samples.
#'
#' @param spe A SpatialExperiment object containing the data to be plotted.
#' @param sample_id A character string specifying the column name in
#'   `colData(spe)` that contains unique sample identifiers. Default is
#'   'sample_id'.
#' @param metric A character string specifying the metric to be visualized
#'   in the plot. This metric should be a column name in `colData(spe)`.
#' @param outliers A character string specifying the column name in
#'   `colData(spe)` that indicates whether a data point is considered an
#'   outlier. Default is local_outliers'.
#' @param colors A character vector specifying the colors to be used for the
#'  gradient scale. If length is 2, the gradient will be a single color gradient
#' @param stroke A numeric value specifying the border thickness for outlier
#'   points. Default is 1.
#' @param width A numeric value indicating the width of the plot. Default
#'   is 5.
#' @param height A numeric value indicating the height of the plot. Default
#'   is 5.
#' @param fname A character string specifying the path and name of the output
#'  PDF file.
#' @param point_size A numeric value specifying the size of the points in the
#'  plot. Default is 2.
#'
#' @return ggplot object if specified. Generates a plot otherwise.
#'
#' @importFrom escheR make_escheR add_fill add_ground
#' @importFrom grDevices dev.off pdf
#'
#' @examples
#' library(SpotSweeper)
#' library(SpatialExperiment)
#' library(escheR)
#'
#' # load example data
#' spe <- STexampleData::Visium_humanDLPFC()
#'
#' tempFilePath <- file.path(tempdir(), "examplePlot.pdf")
#'
#' # change from gene id to gene names
#' rownames(spe) <- rowData(spe)$gene_name
#'
#' # drop out-of-tissue spots
#' spe <- spe[, spe$in_tissue == 1]
#' spe <- spe[, !is.na(spe$ground_truth)]
#'
#' # Identifying the mitochondrial transcripts in our SpatialExperiment.
#' is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]
#'
#' # Calculating QC metrics for each spot using scuttle
#' spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
#' colnames(colData(spe))
#'
#' # Identifying local outliers using SpotSweeper
#' spe <- localOutliers(spe,
#'                      metric = "sum",
#'                      direction = "lower",
#'                      log = TRUE
#' )
#'
#' plotQCpdf(spe,
#'           metric="sum",
#'           outliers="sum_outliers",
#'           fname=tempFilePath)
#'
#' @export
plotQCpdf <- function(
        spe, sample_id = "sample_id", metric = "detected",
        outliers = "local_outliers", colors = c("white", "black"), stroke = 1,
        point_size = 2, width = 5, height = 5, fname) {
    # Get a list of unique sample IDs
    unique_sample_ids <- unique(colData(spe)[[sample_id]])

    # initialize PDF and loop through plots
    pdf(width = width, height = height, fname)
    for (sample in unique_sample_ids) {
        p <- plotQCmetrics(spe,
            sample_id = sample_id, sample = sample,
            metric = metric, outliers = outliers, colors = colors,
            stroke = stroke, point_size = point_size
        )

        # print
        print(p)
    }
    dev.off()
}

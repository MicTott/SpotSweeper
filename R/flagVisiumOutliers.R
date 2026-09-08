#' Flag systematic Visium outlier positions
#'
#' `flagVisiumOutliers()` marks six Visium array positions that have shown
#' systematically low library sizes across multiple datasets. It accepts either
#' a `SpatialExperiment` or `Seurat` object and adds a logical
#' `systematic_outliers` column to its metadata.
#'
#' @param spe A `SpatialExperiment` or `Seurat` object containing Visium array
#'   coordinates. `SpatialExperiment` inputs must contain `array_row` and
#'   `array_col` metadata columns. Seurat inputs may contain the same columns or
#'   expose them through a Visium image.
#' @param image_id For Seurat inputs, the image containing Visium array
#'   coordinates. The first image is used when `NULL`.
#'
#' @return `spe` with a logical `systematic_outliers` metadata column.
#' @export
#'
#' @examples
#' spe <- STexampleData::Visium_humanDLPFC()
#' spe <- flagVisiumOutliers(spe)
#' spe <- spe[, !spe$systematic_outliers]
flagVisiumOutliers <- function(spe, image_id = NULL) {
  .check_supported_object(spe)
  metadata <- getMetadata(spe)

  if (all(c("array_row", "array_col") %in% colnames(metadata))) {
    array_coordinates <- metadata[, c("array_row", "array_col"), drop = FALSE]
  } else if (is_seurat(spe)) {
    image_coordinates <- .get_seurat_array_coordinates(spe, image_id)
    array_coordinates <- data.frame(
      array_row = rep(NA_real_, nrow(metadata)),
      array_col = rep(NA_real_, nrow(metadata)),
      row.names = rownames(metadata)
    )
    shared_spots <- intersect(
      rownames(array_coordinates),
      rownames(image_coordinates)
    )
    if (length(shared_spots) == 0L) {
      stop(
        "The image coordinates do not match any spots in the Seurat object.",
        call. = FALSE
      )
    }
    array_coordinates[shared_spots, ] <- image_coordinates[shared_spots, ]
  } else {
    stop(
      "SpatialExperiment inputs must contain 'array_row' and 'array_col' ",
      "metadata columns."
    )
  }

  data_environment <- new.env(parent = emptyenv())
  utils::data(
    "biased_spots",
    package = "SpotSweeper",
    envir = data_environment
  )
  biased_spots <- get("biased_spots", envir = data_environment)

  spot_keys <- paste(
    array_coordinates$array_row,
    array_coordinates$array_col,
    sep = ":"
  )
  biased_keys <- paste(biased_spots$row, biased_spots$col, sep = ":")
  metadata$systematic_outliers <- spot_keys %in% biased_keys

  setMetadata(spe, metadata)
}

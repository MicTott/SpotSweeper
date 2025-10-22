#' Flag Visium Outliers in Spatial Objects
#'
#' The `flagVisiumOutliers` function identifies and flags Visium systematic outlier spots in a
#' `SpatialExperiment` or `Seurat` object based on array coordinates. These outliers are marked
#' in the metadata, allowing users to exclude them from downstream analyses to enhance data
#' quality and reliability.
#'
#' @param spe A `SpatialExperiment` or `Seurat` object containing Visium spatial transcriptomics data.
#'   For SpatialExperiment: must include `array_row` and `array_col` columns in `colData`.
#'   For Seurat: must include spatial data with an image containing array coordinates.
#' @param image_id For Seurat objects, which image to use (default: first image).
#'
#' @return The input object with an additional logical column `systematic_outliers` in its metadata.
#'   This column indicates whether each spot is flagged as a technical outlier (`TRUE`) or not (`FALSE`).
#'
#' @importFrom utils data
#' @import SpatialExperiment
#'
#' @export
#'
#' @examples
#' library(SpotSweeper)
#' library(SpatialExperiment)
#'
#' # Example with SpatialExperiment
#' spe <- STexampleData::Visium_humanDLPFC()
#' spe <- flagVisiumOutliers(spe)
#' # drop outlier spots
#' spe <- spe[, !colData(spe)$systematic_outliers]
#'
#' # Example with Seurat object (if Seurat is available):
#' # library(Seurat)
#' # seurat_obj <- flagVisiumOutliers(seurat_obj)
#' # # drop outlier spots
#' # seurat_obj <- subset(seurat_obj, subset = systematic_outliers == FALSE)
#'
flagVisiumOutliers <- function(spe, image_id = NULL) {

  # Check if 'spe' is a supported object type
  if (!is_seurat(spe) && !is_spatial_experiment(spe)) {
    stop("Input must be a SpatialExperiment or Seurat object.")
  }

  # Load the biased_spots dataset
  data("biased_spots", package = "SpotSweeper", envir = environment())

  # Get metadata using compatibility layer
  metadata <- getMetadata(spe)

  # Create a logical mask for cells to drop
  drop_mask <- rep(FALSE, nrow(metadata))

  # Get array coordinates based on object type
  if (is_seurat(spe)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required but not available.\n",
           "Install with: install.packages('Seurat')")
    }

    # Get available images
    available_images <- names(spe@images)
    if (length(available_images) == 0) {
      stop("No spatial images found in Seurat object")
    }

    # Use specified image or first available
    if (is.null(image_id)) {
      image_id <- available_images[1]
      message("Using image: ", image_id)
    }

    if (!image_id %in% available_images) {
      stop("Image '", image_id, "' not found. Available images: ",
           paste(available_images, collapse = ", "))
    }

    # Get tissue coordinates which include array positions
    coords <- Seurat::GetTissueCoordinates(spe, image = image_id)

    # For Visium data, array coordinates should be in the data
    # Check for common column names
    if (all(c("row", "col") %in% colnames(coords))) {
      array_row <- coords$row
      array_col <- coords$col
    } else if (all(c("array_row", "array_col") %in% colnames(coords))) {
      array_row <- coords$array_row
      array_col <- coords$array_col
    } else {
      stop("Could not find array row/col coordinates in Seurat object. ",
           "Available columns: ", paste(colnames(coords), collapse = ", "))
    }

  } else if (is_spatial_experiment(spe)) {
    # SpatialExperiment case
    if (!all(c("array_row", "array_col") %in% colnames(metadata))) {
      stop("SpatialExperiment object must have 'array_row' and 'array_col' in colData")
    }
    array_row <- metadata$array_row
    array_col <- metadata$array_col
  }

  # For each pair of (row, col) in biased_spots, mark the matching cells for removal
  for (i in 1:nrow(biased_spots)) {
    row_match <- array_row == biased_spots$row[i]
    col_match <- array_col == biased_spots$col[i]

    # Mark cells for removal where both conditions are true
    drop_mask <- drop_mask | (row_match & col_match)
  }

  # Add the drop_mask to metadata using compatibility layer
  metadata$systematic_outliers <- drop_mask
  spe <- setMetadata(spe, metadata)

  return(spe)
}

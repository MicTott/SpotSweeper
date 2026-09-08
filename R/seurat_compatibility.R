# Object adapters ---------------------------------------------------------

is_seurat <- function(x) {
  inherits(x, "Seurat")
}

is_spatial_experiment <- function(x) {
  inherits(x, "SpatialExperiment")
}

.check_supported_object <- function(x) {
  if (!is_spatial_experiment(x) && !is_seurat(x)) {
    stop("Object must be a SpatialExperiment or Seurat object.")
  }
  invisible(TRUE)
}

.require_seurat_object <- function() {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "Package 'SeuratObject' is required for Seurat inputs. ",
      "Install it with install.packages('SeuratObject')."
    )
  }
}

#' Get column metadata from a spatial object
#'
#' @param x A `SpatialExperiment` or `Seurat` object.
#'
#' @return A `data.frame` with one row per spot, in object order.
#' @export
#'
#' @examples
#' spe <- STexampleData::Visium_humanDLPFC()
#' metadata <- getMetadata(spe)
#' head(metadata)
getMetadata <- function(x) {
  .check_supported_object(x)

  if (is_spatial_experiment(x)) {
    return(as.data.frame(SummarizedExperiment::colData(x)))
  }

  .require_seurat_object()
  x[[]]
}

#' Replace column metadata in a spatial object
#'
#' @param x A `SpatialExperiment` or `Seurat` object.
#' @param metadata A data-frame-like object with one row per spot. Named rows
#'   are reordered to match the object before assignment.
#'
#' @return `x` with replaced column metadata.
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#' spe <- STexampleData::Visium_humanDLPFC()
#' metadata <- getMetadata(spe)
#' metadata$example_column <- seq_len(nrow(metadata))
#' spe <- setMetadata(spe, metadata)
setMetadata <- function(x, metadata) {
  .check_supported_object(x)
  metadata <- as.data.frame(metadata)
  spot_names <- colnames(x)

  if (nrow(metadata) != length(spot_names)) {
    stop("'metadata' must have one row for every spot in 'x'.")
  }

  metadata_names <- rownames(metadata)
  default_names <- as.character(seq_len(nrow(metadata)))
  if (is.null(metadata_names) || identical(metadata_names, default_names)) {
    rownames(metadata) <- spot_names
  } else if (!identical(metadata_names, spot_names)) {
    if (anyDuplicated(metadata_names) ||
        !setequal(metadata_names, spot_names)) {
      stop("The row names of 'metadata' must match the spot names in 'x'.")
    }
    metadata <- metadata[spot_names, , drop = FALSE]
  }

  if (is_spatial_experiment(x)) {
    SummarizedExperiment::colData(x) <- S4Vectors::DataFrame(metadata)
    return(x)
  }

  .require_seurat_object()
  x[[]] <- metadata
  x
}

#' Subset a spatial object by spots
#'
#' @param x A `SpatialExperiment` or `Seurat` object.
#' @param indices Logical, integer, or character spot indices. If
#'   `column_name` is supplied, these are values to retain from that column.
#' @param column_name Optional metadata column used for value-based subsetting.
#'
#' @return A subset of `x` containing the selected spots.
#' @export
#'
#' @examples
#' spe <- STexampleData::Visium_humanDLPFC()
#' spe_subset <- subsetSpatialObject(spe, seq_len(100))
subsetSpatialObject <- function(x, indices, column_name = NULL) {
  .check_supported_object(x)

  if (!is.null(column_name)) {
    validateMetadataColumns(x, column_name)
    indices <- getMetadata(x)[[column_name]] %in% indices
  }

  x[, indices]
}

#' Validate columns in spatial-object metadata
#'
#' @param x A `SpatialExperiment` or `Seurat` object.
#' @param required_columns Character vector of required metadata column names.
#'
#' @return `TRUE`, invisibly. An error is raised if columns are missing.
#' @export
#'
#' @examples
#' spe <- STexampleData::Visium_humanDLPFC()
#' validateMetadataColumns(spe, c("sample_id", "in_tissue"))
validateMetadataColumns <- function(x, required_columns) {
  metadata <- getMetadata(x)
  missing_columns <- setdiff(required_columns, colnames(metadata))

  if (length(missing_columns) > 0L) {
    stop(
      "Required columns missing from metadata: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  invisible(TRUE)
}

# Coordinate adapters -----------------------------------------------------

#' Get spatial coordinates from a spatial object
#'
#' @param x A `SpatialExperiment` or `Seurat` object.
#' @param image_id For a Seurat object, an optional image name. If omitted,
#'   coordinates from all images are combined and returned in object order.
#'
#' @return A numeric matrix with spots as rows and spatial dimensions as
#'   columns.
#' @export
#'
#' @examples
#' spe <- STexampleData::Visium_humanDLPFC()
#' coordinates <- getSpatialCoords(spe)
#' head(coordinates)
getSpatialCoords <- function(x, image_id = NULL) {
  .check_supported_object(x)

  if (is_spatial_experiment(x)) {
    return(SpatialExperiment::spatialCoords(x))
  }

  .require_seurat_object()
  image_names <- SeuratObject::Images(x)
  if (length(image_names) == 0L) {
    stop(
      "The Seurat object has no spatial images. Supply coordinates with ",
      "the 'coords' argument."
    )
  }

  if (!is.null(image_id)) {
    if (length(image_id) != 1L || !image_id %in% image_names) {
      stop(
        "Unknown Seurat image. Available images: ",
        paste(image_names, collapse = ", ")
      )
    }
    image_names <- image_id
  }

  coordinate_list <- lapply(image_names, function(image_name) {
    coordinates <- SeuratObject::GetTissueCoordinates(
      x,
      image = image_name
    )
    .coordinate_matrix(coordinates, image_name)
  })
  coordinates <- do.call(rbind, coordinate_list)

  if (anyDuplicated(rownames(coordinates))) {
    stop("A spot occurs in more than one Seurat spatial image.")
  }

  if (is.null(image_id)) {
    spot_names <- colnames(x)
    missing_spots <- setdiff(spot_names, rownames(coordinates))
    if (length(missing_spots) > 0L) {
      stop(
        "Spatial images do not contain every spot in the Seurat object. ",
        "Supply a complete matrix with the 'coords' argument."
      )
    }
    coordinates <- coordinates[spot_names, , drop = FALSE]
  }

  coordinates
}

.coordinate_matrix <- function(coordinates, image_name) {
  coordinates <- as.data.frame(coordinates)

  if (is.null(rownames(coordinates)) && "cell" %in% colnames(coordinates)) {
    rownames(coordinates) <- coordinates$cell
  }
  if (is.null(rownames(coordinates))) {
    stop("Coordinates for Seurat image '", image_name, "' are not named.")
  }

  coordinate_pairs <- list(
    c("x", "y"),
    c("imagecol", "imagerow"),
    c("col", "row")
  )
  selected <- NULL
  for (pair in coordinate_pairs) {
    if (all(pair %in% colnames(coordinates))) {
      selected <- pair
      break
    }
  }
  if (is.null(selected)) {
    stop(
      "Could not identify spatial coordinate columns for Seurat image '",
      image_name, "'."
    )
  }

  result <- as.matrix(coordinates[, selected, drop = FALSE])
  storage.mode(result) <- "double"
  colnames(result) <- c("x", "y")
  result
}

.get_seurat_array_coordinates <- function(x, image_id = NULL) {
  .require_seurat_object()
  image_names <- SeuratObject::Images(x)
  if (length(image_names) == 0L) {
    stop("The Seurat object has no spatial images with array coordinates.")
  }

  if (is.null(image_id)) {
    image_id <- image_names[[1L]]
  } else if (length(image_id) != 1L || !image_id %in% image_names) {
    stop(
      "Unknown Seurat image. Available images: ",
      paste(image_names, collapse = ", ")
    )
  }

  coordinates <- tryCatch(
    SeuratObject::GetTissueCoordinates(
      x,
      image = image_id,
      scale = NULL,
      cols = c("row", "col")
    ),
    error = function(error) NULL
  )
  coordinates <- as.data.frame(coordinates)

  coordinate_names <- if (all(c("row", "col") %in% colnames(coordinates))) {
    c("row", "col")
  } else if (all(c("array_row", "array_col") %in%
                 colnames(coordinates))) {
    c("array_row", "array_col")
  } else {
    stop(
      "Seurat image '", image_id, "' does not expose Visium array row and ",
      "column coordinates. Add 'array_row' and 'array_col' to the object's ",
      "metadata or select a Visium image."
    )
  }

  result <- coordinates[, coordinate_names, drop = FALSE]
  colnames(result) <- c("array_row", "array_col")
  result
}

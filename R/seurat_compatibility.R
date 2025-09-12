#' SpotSweeper Seurat Compatibility Layer
#' 
#' Functions to enable SpotSweeper compatibility with Seurat spatial objects
#' alongside SpatialExperiment objects.

#' Check if object is a Seurat object
#' @param x Object to check
#' @return Logical indicating if object is Seurat
is_seurat <- function(x) {
  inherits(x, "Seurat")
}

#' Check if object is a SpatialExperiment object
#' @param x Object to check  
#' @return Logical indicating if object is SpatialExperiment
is_spatial_experiment <- function(x) {
  inherits(x, "SpatialExperiment")
}

#' Get spatial coordinates from either SpatialExperiment or Seurat object
#' 
#' @param x SpatialExperiment or Seurat object
#' @param image_id For Seurat objects, which image to use (default: first image)
#' @return Matrix of spatial coordinates with rows as spots and columns as x,y coordinates
#' @export
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' coords <- getSpatialCoords(spe)
#' head(coords)
getSpatialCoords <- function(x, image_id = NULL) {
  if (is_seurat(x)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required but not available.\n",
           "Install with: install.packages('Seurat')")
    }
    
    # Get available images
    available_images <- names(x@images)
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
    
    # Extract coordinates using Seurat function
    coords <- Seurat::GetTissueCoordinates(x, image = image_id)
    
    # Convert to matrix format matching spatialCoords output
    coords_matrix <- as.matrix(coords[, c("imagerow", "imagecol")])
    colnames(coords_matrix) <- c("x", "y")
    
    return(coords_matrix)
    
  } else if (is_spatial_experiment(x)) {
    if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
      stop("Package 'SpatialExperiment' is required but not available.\n",
           "Install with: BiocManager::install('SpatialExperiment')")
    }
    
    return(SpatialExperiment::spatialCoords(x))
    
  } else {
    stop("Object must be either a Seurat or SpatialExperiment object")
  }
}

#' Get metadata from either SpatialExperiment or Seurat object
#' 
#' @param x SpatialExperiment or Seurat object
#' @return Data.frame of metadata
#' @export
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' metadata <- getMetadata(spe)
#' colnames(metadata)
getMetadata <- function(x) {
  if (is_seurat(x)) {
    return(x@meta.data)
    
  } else if (is_spatial_experiment(x)) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Package 'SummarizedExperiment' is required but not available.\n",
           "Install with: BiocManager::install('SummarizedExperiment')")
    }
    
    # Convert S4 DataFrame to regular data.frame for consistency
    return(as.data.frame(SummarizedExperiment::colData(x)))
    
  } else {
    stop("Object must be either a Seurat or SpatialExperiment object")
  }
}

#' Set metadata for either SpatialExperiment or Seurat object
#' 
#' @param x SpatialExperiment or Seurat object
#' @param metadata Data.frame of metadata to set
#' @return Modified object with updated metadata
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' metadata <- getMetadata(spe)
#' metadata$new_column <- 1
#' spe_updated <- setMetadata(spe, metadata)
#' "new_column" %in% colnames(getMetadata(spe_updated))
setMetadata <- function(x, metadata) {
  if (is_seurat(x)) {
    x@meta.data <- metadata
    return(x)
    
  } else if (is_spatial_experiment(x)) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Package 'SummarizedExperiment' is required but not available.\n",
           "Install with: BiocManager::install('SummarizedExperiment')")
    }
    
    # Convert data.frame to DataFrame for SpatialExperiment compatibility
    if (is.data.frame(metadata)) {
      metadata <- DataFrame(metadata)
    }
    SummarizedExperiment::colData(x) <- metadata
    return(x)
    
  } else {
    stop("Object must be either a Seurat or SpatialExperiment object")
  }
}

#' Subset spatial object by cell/spot indices
#' 
#' @param x SpatialExperiment or Seurat object
#' @param indices Logical or integer vector for subsetting
#' @param column_name Column name for logical subsetting (for Seurat objects)
#' @return Subsetted object
#' @export
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' # Subset first 100 spots
#' spe_subset <- subsetSpatialObject(spe, 1:100)
#' ncol(spe_subset)
subsetSpatialObject <- function(x, indices, column_name = NULL) {
  if (is_seurat(x)) {
    if (is.logical(indices) && !is.null(column_name)) {
      # Subset by logical column
      cells_to_keep <- rownames(x@meta.data)[x@meta.data[[column_name]] %in% indices]
      return(subset(x, cells = cells_to_keep))
    } else {
      # Direct subsetting by indices
      if (is.logical(indices)) {
        cells_to_keep <- rownames(x@meta.data)[indices]
      } else {
        cells_to_keep <- rownames(x@meta.data)[indices]
      }
      return(subset(x, cells = cells_to_keep))
    }
    
  } else if (is_spatial_experiment(x)) {
    return(x[, indices])
    
  } else {
    stop("Object must be either a Seurat or SpatialExperiment object")
  }
}

#' Validate that required columns exist in metadata
#' 
#' @param x SpatialExperiment or Seurat object
#' @param required_columns Character vector of required column names
#' @return Logical indicating if all columns exist
#' @export
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' # Check if required columns exist
#' validateMetadataColumns(spe, c("sample_id", "in_tissue"))
validateMetadataColumns <- function(x, required_columns) {
  metadata <- getMetadata(x)
  missing_cols <- setdiff(required_columns, colnames(metadata))
  
  if (length(missing_cols) > 0) {
    stop("Required columns missing from metadata: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  return(TRUE)
}
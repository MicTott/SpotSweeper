#' Detect local outliers in spatial transcriptomics data
#'
#' `localOutliers()` compares each spot's quality-control metric with the same
#' metric in its nearest neighbors. It supports `SpatialExperiment` and
#' `Seurat` objects and stores the resulting robust z-scores and outlier calls
#' in the object's column metadata.
#'
#' The robust z-score for spot i is
#' `qnorm(0.75) * (x_i - median(x_neighbors)) / MAD_unscaled(x_neighbors)`.
#' The focal spot is not included when estimating the neighborhood median or
#' median absolute deviation. A score of zero is returned when the neighborhood
#' MAD is zero or non-finite.
#'
#' @param spe A `SpatialExperiment` or `Seurat` object.
#' @param metric A single column name in the object's column metadata containing
#'   a numeric QC metric.
#' @param direction Direction of outlier detection: `"higher"`, `"lower"`, or
#'   `"both"`.
#' @param n_neighbors Number of nearest neighbors used as the local reference.
#' @param samples A single metadata column name containing sample identifiers.
#' @param log Whether to apply `log1p()` to `metric` before scoring.
#' @param cutoff Non-negative absolute robust z-score cutoff.
#' @param workers Number of workers passed to `BiocNeighbors::findKNN()`.
#' @param coords Optional numeric coordinate matrix with spots as rows. When
#'   `NULL`, spatial coordinates are obtained from `spe`.
#'
#' @return The input object with `<metric>_z` and `<metric>_outliers` metadata
#'   columns. When `log = TRUE`, `<metric>_log` is also added.
#'
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocParallel MulticoreParam
#' @export
#'
#' @examples
#' library(SpotSweeper)
#' library(SpatialExperiment)
#'
#' spe <- STexampleData::Visium_humanDLPFC()
#' rownames(spe) <- rowData(spe)$gene_name
#' spe <- spe[, spe$in_tissue == 1]
#' spe <- spe[, !is.na(spe$ground_truth)]
#'
#' is.mito <- grepl("^MT-", rownames(spe))
#' spe <- scuttle::addPerCellQCMetrics(
#'     spe,
#'     subsets = list(Mito = is.mito)
#' )
#'
#' spe <- localOutliers(
#'     spe,
#'     metric = "sum",
#'     direction = "lower",
#'     log = TRUE
#' )
#'
#' # Seurat objects use the same interface. Custom coordinates can also be
#' # supplied for either object type.
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'     counts <- matrix(
#'         seq_len(24),
#'         nrow = 4,
#'         dimnames = list(paste0("gene", 1:4), paste0("spot", 1:6))
#'     )
#'     seurat <- SeuratObject::CreateSeuratObject(counts)
#'     seurat$sample_id <- "sample"
#'     seurat$detected <- c(100, 10, 11, 9, 12, 10)
#'     example_coords <- cbind(x = seq_len(6), y = 0)
#'     seurat <- localOutliers(
#'         seurat,
#'         n_neighbors = 4,
#'         coords = example_coords
#'     )
#' }
localOutliers <- function(
    spe, metric = "detected", direction = "lower", n_neighbors = 36,
    samples = "sample_id", log = TRUE, cutoff = 3, workers = 1,
    coords = NULL) {
  if (!is_spatial_experiment(spe) && !is_seurat(spe)) {
    stop("'spe' must be a SpatialExperiment or Seurat object.")
  }

  if (!is.character(metric) || length(metric) != 1L || is.na(metric)) {
    stop("'metric' must be a single column name.")
  }
  if (!is.character(samples) || length(samples) != 1L || is.na(samples)) {
    stop("'samples' must be a single column name.")
  }
  if (!is.character(direction) || length(direction) != 1L ||
      !direction %in% c("lower", "higher", "both")) {
    stop("'direction' must be one of 'lower', 'higher', or 'both'.")
  }
  if (!.is_positive_integer(n_neighbors)) {
    stop("'n_neighbors' must be a positive integer.")
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L ||
      !is.finite(cutoff) || cutoff < 0) {
    stop("'cutoff' must be a single non-negative numeric value.")
  }
  if (!is.logical(log) || length(log) != 1L || is.na(log)) {
    stop("'log' must be TRUE or FALSE.")
  }
  if (!.is_positive_integer(workers)) {
    stop("'workers' must be a positive integer.")
  }

  metadata <- getMetadata(spe)
  validateMetadataColumns(spe, c(metric, samples))

  if (!is.numeric(metadata[[metric]])) {
    stop("'metric' must identify a numeric metadata column.")
  }
  if (anyNA(metadata[[samples]])) {
    stop("The sample identifier column must not contain missing values.")
  }

  spot_names <- colnames(spe)
  coords <- .prepare_coordinates(spe, coords, spot_names)

  metric_to_use <- metric
  if (log) {
    metric_to_use <- paste0(metric, "_log")
    metadata[[metric_to_use]] <- log1p(metadata[[metric]])
  }
  if (any(!is.finite(metadata[[metric_to_use]]))) {
    stop("The metric contains non-finite values after transformation.")
  }

  z_scores <- numeric(nrow(metadata))
  sample_ids <- unique(metadata[[samples]])

  for (sample_id in sample_ids) {
    sample_indices <- which(metadata[[samples]] == sample_id)
    if (length(sample_indices) <= n_neighbors) {
      stop(
        "Each sample must contain more spots than 'n_neighbors'. ",
        "Sample '", sample_id, "' contains ", length(sample_indices),
        " spots."
      )
    }

    sample_coords <- coords[sample_indices, , drop = FALSE]
    sample_values <- metadata[[metric_to_use]][sample_indices]
    parallel_param <- if (workers == 1L) {
      BiocParallel::SerialParam()
    } else {
      BiocParallel::MulticoreParam(workers = workers)
    }
    neighbors <- BiocNeighbors::findKNN(
      sample_coords,
      k = n_neighbors,
      warn.ties = FALSE,
      BPPARAM = parallel_param
    )$index

    sample_z <- vapply(seq_along(sample_values), function(i) {
      .local_modified_z(sample_values[[i]], sample_values[neighbors[i, ]])
    }, numeric(1))
    z_scores[sample_indices] <- sample_z
  }

  metadata[[paste0(metric, "_z")]] <- z_scores
  metadata[[paste0(metric, "_outliers")]] <- switch(
    direction,
    higher = z_scores > cutoff,
    lower = z_scores < -cutoff,
    both = abs(z_scores) > cutoff
  )

  setMetadata(spe, metadata)
}

.local_modified_z <- function(focal, neighbors) {
  center <- stats::median(neighbors)
  mad_unscaled <- stats::median(abs(neighbors - center))

  if (!is.finite(mad_unscaled) || mad_unscaled == 0) {
    return(0)
  }

  stats::qnorm(0.75) * (focal - center) / mad_unscaled
}

.is_positive_integer <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 &&
    x == as.integer(x)
}

.prepare_coordinates <- function(spe, coords, spot_names) {
  if (is.null(coords)) {
    coords <- getSpatialCoords(spe)
  }
  if (!is.matrix(coords) || !is.numeric(coords)) {
    stop("'coords' must be a numeric matrix.")
  }
  if (nrow(coords) != length(spot_names)) {
    stop("'coords' must have one row for every spot in 'spe'.")
  }
  if (ncol(coords) < 1L || any(!is.finite(coords))) {
    stop("'coords' must contain at least one finite numeric column.")
  }

  coord_names <- rownames(coords)
  if (!is.null(coord_names)) {
    if (anyDuplicated(coord_names) || !setequal(coord_names, spot_names)) {
      stop("The row names of 'coords' must match the spot names in 'spe'.")
    }
    coords <- coords[spot_names, , drop = FALSE]
  } else {
    rownames(coords) <- spot_names
  }

  coords
}

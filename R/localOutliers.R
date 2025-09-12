#' localOutliers Function
#'
#' This function detects local outliers in spatial transcriptomics data based on
#' standard quality control metrics, such as library size, unique genes, and
#' mitochondrial ratio. Local outliers are defined as spots with low/high
#' quality metrics compared to their surrounding neighbors, based on a modified
#' z-score statistic.
#'
#' @param spe SpatialExperiment, SingleCellExperiment, or Seurat object
#' @param metric Metadata QC metric to use for outlier detection (colData for SpatialExperiment, meta.data for Seurat)
#' @param direction Direction of outlier detection (higher, lower, or both)
#' @param n_neighbors Number of nearest neighbors to use for outlier detection
#' @param samples Column name in metadata to use for sample IDs
#' @param log Logical indicating whether to log1p transform the features
#' (default is TRUE)
#' @param cutoff Cutoff for outlier detection (default is 3)
#' @param workers Number of workers for parallel processing (default is 1)
#'
#' @return SpatialExperiment, SingleCellExperiment, or Seurat object with updated metadata containing outputs
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom spatialEco outliers
#' @importFrom BiocParallel MulticoreParam
#'
#' @export localOutliers
#'
#' @examples
#' library(SpotSweeper)
#' library(SpatialExperiment)
#'
#' # load example data
#' spe <- STexampleData::Visium_humanDLPFC()
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
#' # Example with Seurat object (if Seurat is available):
#' # library(Seurat) 
#' # # Assuming 'seurat_obj' is a Seurat object with spatial data
#' # seurat_obj <- localOutliers(seurat_obj,
#' #                             metric = "nCount_RNA",
#' #                             direction = "lower", 
#' #                             samples = "orig.ident")
#'
localOutliers <- function(
    spe, metric = "detected",
    direction = "lower", n_neighbors = 36, samples = "sample_id",
    log = TRUE, cutoff = 3, workers=1) {

  # ===== Validity checks =====
  # Check if 'spe' is a supported object type
  if (!is_seurat(spe) && !is_spatial_experiment(spe)) {
    stop("Input must be a SpatialExperiment, SingleCellExperiment, or Seurat object.")
  }

  # Validate 'direction'
  if (!direction %in% c("lower", "higher", "both")) {
    stop("'direction' must be one of 'lower', 'higher', or 'both'.")
  }

  # Check 'n_neighbors' is a positive integer
  if (!is.numeric(n_neighbors) ||
      n_neighbors <= 0 ||
      n_neighbors != round(n_neighbors)) {
    stop("'n_neighbors' must be a positive integer.")
  }

  # Check 'cutoff' is a numeric value
  if (!is.numeric(cutoff)) {
    stop("'cutoff' must be a numeric value.")
  }

  # ===== Start function =====
  # Get metadata using compatibility layer
  metadata <- getMetadata(spe)
  
  # Validate required columns exist
  validateMetadataColumns(spe, c(metric, samples))
  
  # log transform specified metric
  if (log) {
    metric_log <- paste0(metric, "_log")
    metadata[[metric_log]] <- log1p(metadata[[metric]])
    metric_to_use <- metric_log
  } else {
    metric_to_use <- metric
  }

  # Get a list of unique sample IDs
  unique_sample_ids <- unique(metadata[[samples]])

  # Initialize list to store each metadata
  metadata_list <- sapply(unique_sample_ids, FUN = function(x) NULL)

  # Loop through each unique sample ID
  for (sample in unique_sample_ids) {
    # Subset the data for the current sample
    sample_indices <- metadata[[samples]] == sample
    spe_subset <- subsetSpatialObject(spe, sample_indices)

    # Get metadata and spatial coordinates for subset
    subset_metadata <- getMetadata(spe_subset)
    spatial_coords <- getSpatialCoords(spe_subset)

    # Find nearest neighbors
    dnn <- BiocNeighbors::findKNN(spatial_coords,
                                  k = n_neighbors, warn.ties = FALSE,
                                  BPPARAM = BiocParallel::MulticoreParam(workers=workers))$index

    # get neighborhood metrics
    neighborhoods <- lapply(seq_len(nrow(dnn)), function(i) {
      indices <- dnn[i, ]
      indices <- indices[indices != 0]
      indices <- c(i, indices)

      subset_metadata[indices, metric_to_use]
    })

    # Compute modified-z and return the middle spot
    mod_z_matrix <- vapply(neighborhoods, function(x) {
      spatialEco::outliers(x)[1]
    }, numeric(1))

    # Handle non-finite values
    mod_z_matrix[!is.finite(mod_z_matrix)] <- 0

    # find outliers based on cutoff, store in metadata
    metric_outliers <- paste0(metric, "_outliers")
    subset_metadata[[metric_outliers]] <- switch(direction,
                                          higher = sapply(
                                            mod_z_matrix,
                                            function(x) x > cutoff
                                          ),
                                          lower = sapply(
                                            mod_z_matrix,
                                            function(x) x < -cutoff
                                          ),
                                          both = sapply(
                                            mod_z_matrix,
                                            function(x) {
                                              x > cutoff | x <
                                                -cutoff
                                            }
                                          )
    )

    # add z-scores to metadata
    metric_z <- paste0(metric, "_z")
    subset_metadata[[metric_z]] <- mod_z_matrix

    # Store the modified metadata dataframe in the list
    metadata_list[[sample]] <- subset_metadata
  }

  # rbind the list of dataframes
  metadata_aggregated <- do.call(rbind, metadata_list)

  # replace metadata using compatibility layer
  spe <- setMetadata(spe, metadata_aggregated)

  return(spe)
}

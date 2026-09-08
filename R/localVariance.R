#' localVariance Function
#'
#' This function calculates the local variance based on kNN.
#'
#' @param spe SpatialExperiment object with the following columns in colData:
#'        sample_id, sum_umi, sum_gene
#' @param n_neighbors Number of nearest neighbors to use for variance
#'        calculation
#' @param metric Metric to use for variance calculation
#' @param samples Column in colData to use for sample ID
#' @param log Whether to log1p transform the metric
#' @param name Name of the new column to add to colData
#' @param workers Number of workers to use for parallel computation
#'
#' @return SpatialExperiment object with metric variance added to colData
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocParallel MulticoreParam
#' @importFrom MASS rlm
#'
#' @export localVariance
#'
#' @examples
#'
#' # for more details see extended example in vignettes
#' library(SpotSweeper)
#' library(SpatialExperiment)
#' library(escheR)
#'
#' # load example data
#' spe <- STexampleData::Visium_humanDLPFC()
#'
#' # change from gene id to gene names
#' rownames(spe) <- rowData(spe)$gene_name
#'
#' # show column data before SpotSweepR
#' colnames(colData(spe))
#'
#' # drop out-of-tissue spots
#' spe <- spe[, spe$in_tissue == 1]
#' spe <- spe[, !is.na(spe$ground_truth)]
#'
#' # Identifying the mitochondrial transcripts in our SpatialExperiment.
#' is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]
#'
#' # Calculating QC metric for each spot using scuttle
#' spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
#' colnames(colData(spe))
#'
#' spe <- localVariance(spe,
#'     metric = "subsets_Mito_percent",
#'     n_neighbors = 36,
#'     name = "local_mito_variance_k36",
#'     workers = 1
#'     )
#'
#'
localVariance <- function(spe, n_neighbors = 36,
                          metric = c("expr_chrM_ratio"),
                          samples = "sample_id", log = FALSE, name = NULL,
                          workers = 1) {

  # ===== Validity checks =====
  # Check if 'spe' is (or inherits from) SpatialExperiment
  if (!inherits(spe, "SpatialExperiment")) {
    stop("Input data must be a SpatialExperiment or inherit from SpatialExperiment.")
  }

  # Validate 'metric' is a character vector
  if (!is.character(metric)) {
    stop("'metric' must be a character vector.")
  }

  # validate 'metric' are valid colData
  if (!all(metric %in% colnames(colData(spe)))) {
    stop("Metric must be present in colData.")
  }

  # validate samples is valid colData
  if (!samples %in% colnames(colData(spe))) {
    stop("Samples column must be present in colData.")
  }

  # Check 'n_neighbors' is a positive integer
  if (!is.numeric(n_neighbors) ||
      n_neighbors <= 0 ||
      n_neighbors != round(n_neighbors)) {
    stop("'n_neighbors' must be a positive integer.")
  }

  # Check 'workers' is a positive integer
  if (!is.numeric(workers) ||
      workers <= 0 ||
      workers != round(workers)) {
    stop("'workers' must be a positive integer.")
  }

  # ===== start function =====
  # log1p transform specified metric
  if (log) {
    metric_log <- paste0(metric, "_log")
    colData(spe)[metric_log] <- log1p(colData(spe)[[metric]])
    metric_to_use <- metric_log
  } else {
    metric_to_use <- metric
  }

  # Get a list of unique sample IDs
  unique_sample_ids <- unique(colData(spe)[[samples]])

  # Initialize list to store each columnData dataframe
  columnData_list <- vector("list", length(unique_sample_ids))

  # Loop through each unique sample ID
  for (sample_id in seq_along(unique_sample_ids)) {
    # Subset the data for the current sample
    sample <- unique_sample_ids[sample_id]
    spe_subset <- subset(spe, , sample_id == sample)

    # Create a list of spatial coordinates and qc metric
    columnData <- colData(spe_subset)
    columnData$coords <- spatialCoords(spe_subset)

    # Find nearest neighbors
    parallel_param <- if (workers == 1L) {
      BiocParallel::SerialParam()
    } else {
      BiocParallel::MulticoreParam(workers = workers)
    }
    dnn <- BiocNeighbors::findKNN(spatialCoords(spe_subset),
                                  k = n_neighbors,
                                  BPPARAM = parallel_param,
                                  warn.ties = FALSE
    )$index

    #  === Compute local variance ===
    # get neighborhood metric
    neighborhoods <- lapply(seq_len(nrow(dnn)), function(i) {
      indices <- dnn[i, ]
      indices <- indices[indices != 0]
      indices <- c(i, indices)

      columnData[indices, metric_to_use]
    })

    # Compute variance and mean
    stats_matrix <- t(vapply(neighborhoods, function(x) {
      c(var = var(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE))
    }, numeric(2)))

    # Handle non-finite values
    stats_matrix[!is.finite(stats_matrix)] <- 0

    # Perform robust linear regression to regress out mean-var bias
    fit.irls <- MASS::rlm(log2(var) ~ mean,
                          data = as.data.frame(stats_matrix))

    var_resid <- resid(fit.irls) # Get residuals

    # add local variance to columnData dataframe
    if (!is.null(name)) {
      columnData[name] <- var_resid
    } else {
      metric_var <- paste0(metric, "_var")
      columnData[metric_var] <- var_resid
    }

    # Store the modified columnData dataframe in the list
    columnData_list[[sample_id]] <- columnData
  }

  # rbind the list of dataframes
  columnData_aggregated <- do.call(rbind, columnData_list)

  # replace SPE column data with aggregated data
  colData(spe) <- columnData_aggregated

  return(spe)
}

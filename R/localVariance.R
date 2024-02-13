#' localVariance Function
#'
#' This function does calculates the local variance based on kNN.
#'
#' @param spe SpatialExperiment object with the following columns in colData: sample_id, sum_umi, sum_gene
#' @param n_neighbors Number of nearest neighbors to use for variance calculation
#' @param features Features to use for variance calculation
#' @param samples Column in colData to use for sample ID
#' @param log2 Whether to log2 transform the features
#' @param n_cores Number of cores to use for parallelization
#' @param name Name of the new column to add to colData
#'
#' @return SpatialExperiment object with feature variance added to colData
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocParallel MulticoreParam
#'
#' @export localVariance
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
#' # Calculating QC metrics for each spot using scuttle
#' spe<- scuttle::addPerCellQCMetrics(spe, subsets=list(Mito=is.mito))
#' colnames(colData(spe))
#'
#' spe <- localVariance(spe,
#'                      features = "subsets_Mito_percent",
#'                      n_neighbors=36,
#'                      name="local_mito_variance_k36"
#'                      )
localVariance <- function(spe, n_neighbors = 36, features = c("expr_chrM_ratio"),
                          samples = "sample_id", log2 = FALSE, n_cores = 1, name=NULL) {

  # log2 transform specified features
  features_to_use <- character()
  if (log2) {
    for (feature in features) {
      feature_log2 <- paste0(feature, "_log2")
      colData(spe)[feature_log2] <- log2(colData(spe)[[feature]])
      features_to_use <- c(features_to_use, feature_log2)
    }
  } else {
    features_to_use <- features
  }

  # Get a list of unique sample IDs
  unique_sample_ids <- unique(colData(spe)[[samples]])

  # Initialize list to store each spaQC dataframe
  spaQC_list <- vector("list", length(unique_sample_ids))

  # Loop through each unique sample ID
  for (sample_id in seq_along(unique_sample_ids)) {
    # Subset the data for the current sample
    sample <- unique_sample_ids[sample_id]
    spe_subset <- subset(spe, , sample_id == sample)

    # Create a list of spatial coordinates and qc features
    spaQC <- colData(spe_subset)
    spaQC$coords <- spatialCoords(spe_subset)

    # Find nearest neighbors
    suppressWarnings(
    dnn <- BiocNeighbors::findKNN(spatialCoords(spe_subset),
                                  k = n_neighbors,
                                  BPPARAM=MulticoreParam(n_cores))$index
    )

    #  === Get local variance ===
    # Initialize a matrix to store variance for each feature
    var_matrix <- matrix(NA, nrow(spaQC), length(features_to_use))
    colnames(var_matrix) <- features_to_use

    mean_matrix <- matrix(NA, nrow(spaQC), length(features_to_use))
    colnames(mean_matrix) <- features_to_use

    # Loop through each row in the nearest neighbor index matrix
    for (i in 1:nrow(dnn)) {
      dnn.idx <- dnn[i, ]
      for (j in seq_along(features_to_use)) {
        var_matrix[i, j] <- var(spaQC[c(i, dnn.idx[dnn.idx != 0]), ][[features_to_use[j]]], na.rm=TRUE)[1]
        mean_matrix[i, j] <- mean(spaQC[c(i, dnn.idx[dnn.idx != 0]), ][[features_to_use[j]]], na.rm=TRUE)[1]
      }
    }

    # Handle non-finite values
    var_matrix[!is.finite(var_matrix)] <- 0
    mean_matrix[!is.finite(mean_matrix)] <- 0

    # ==== If remove.bias == TRUE, regress out mean-variance bias for each feature ====
    for (feature_idx in seq_along(features_to_use)) {

      # Prepare data.frame for current feature
      mito_var_df <- data.frame(
        mito_var = log2(var_matrix[, feature_idx]),  # log2 variance of the current feature
        mito_mean = log2(mean_matrix[, feature_idx])  # log2 mean of the current feature
      )

      # Perform robust linear regression (IRLS) of variance vs mean for the current feature
      fit.irls <- MASS::rlm(mito_var ~ mito_mean, data = mito_var_df)

      # Get residuals and update the variance matrix for the current feature
      resid.irls <- resid(fit.irls)

      # Replace original variance values with residuals for the current feature
      # Note: You might need to back-transform the residuals if you want to maintain the original scale
      var_matrix[, feature_idx] <- resid.irls  # Back-transform if necessary

    }

    # add local variance to spaQC dataframe
    if (!is.null(name)) {
      spaQC[name] <- var_matrix[, j]

    } else {
      feature_var <- paste0(features[j], "_var")
      spaQC[feature_var] <- var_matrix[, j]
    }

    # Store the modified spaQC dataframe in the list
    spaQC_list[[sample_id]] <- spaQC
  }

  # rbind the list of dataframes
  spaQC_aggregated <- do.call(rbind, spaQC_list)

  # replace SPE column data with aggregated data
  colData(spe) <- spaQC_aggregated

  return(spe)
}

#' localOutliers Function
#'
#' This function detects local outliers based on k-nearest neighbors based on either a univariate
#' (z-score thresholds per QC metrics) or multivariate approach (Local Outlier Factor).
#'
#' @param spe SpatialExperiment object with the following columns in colData: sample_id, sum_umi, sum_gene
#' @param n_neighbors Number of nearest neighbors to use for outlier detection
#' @param threshold Threshold for outlier detection
#' @param features Vector of features to use for outlier detection
#' @param method Method to use for outlier detection (univariate or multivariate)
#' @param samples Column name in colData to use for sample IDs
#' @param log2 Logical indicating whether to log2 transform the features
#' @param cutoff Cutoff for outlier detection
#' @param scale Logical indicating whether to scale the features for LOF calculation (recommended)
#' @param minPts Minimum number of points (nearest neighbors) to use for LOF calculation
#' @param data_output Logical indicating whether to output the z-scores for each feature
#' @param n_cores Number of cores to use for parallelization in the findKNN function
#'
#' @return SpatialExperiment object with updated colData
#'
#' @importFrom dbscan lof
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocParallel MulticoreParam
#'
#' @export localOutliers
#'
#' @examples
#' library(SpotSweeper)
#' library(spatialLIBD)
#'
#' spe <- fetch_data(type = "spe")
#'
#' spe <- localOutliers(spe, n_neighbors=36)
localOutliers <- function(spe, n_neighbors = 36, features = c("sum_umi","sum_gene", "expr_chrM_ratio"), method = "multivariate", samples = "sample_id", log2 = TRUE, cutoff = 2.58,
    scale = TRUE, minPts = 20, data_output = FALSE, n_cores = 1) {
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

        # Initialize a matrix to store z-scores for each feature
        mod_z_matrix <- matrix(NA, nrow(spaQC), length(features_to_use))
        colnames(mod_z_matrix) <- features_to_use

        # Initialize a matrix to store variance for each feature
        var_matrix <- matrix(NA, nrow(spaQC), length(features_to_use))
        colnames(var_matrix) <- features_to_use

        # Loop through each row in the nearest neighbor index matrix
        for (i in 1:nrow(dnn)) {
            dnn.idx <- dnn[i, ]
            for (j in seq_along(features_to_use)) {
                mod_z_matrix[i, j] <- modified_z(spaQC[c(i, dnn.idx[dnn.idx != 0]), ][[features_to_use[j]]])[1]
            }
        }

        # Handle non-finite values
        mod_z_matrix[!is.finite(mod_z_matrix)] <- 0

        # Determine local outliers
        if (method == "univariate") {
            # Determine local outliers across all features using z threshold
            spaQC$local_outliers <- as.factor(apply(mod_z_matrix, 1, function(x) any(x > cutoff | x < -cutoff)))
        } else if (method == "multivariate") {
            # convert matrix to df and set column names
            df <- as.data.frame(mod_z_matrix)
            colnames(df) <- features_to_use

            # calculate local outlier factor
            if (scale) {
                LOF_outs <- lof(scale(df), minPts = minPts)
            } else {
                LOF_outs <- lof(df, minPts = minPts)
            }

            spaQC$local_outliers <- as.factor(ifelse(LOF_outs > cutoff, TRUE, FALSE))
        }

        # output z features if desired
        if (data_output) {
            for (j in seq_along(features_to_use)) {
                feature_z <- paste0(features[j], "_z")
                spaQC[feature_z] <- mod_z_matrix[, j]
            }
            if (method == "multivariate") {
                spaQC$LOF <- LOF_outs
            }
        }

        # Store the modified spaQC dataframe in the list
        spaQC_list[[sample_id]] <- spaQC
    }

    # rbind the list of dataframes
    spaQC_aggregated <- do.call(rbind, spaQC_list)

    # replace column data
    colData(spe) <- spaQC_aggregated

    return(spe)
}

#' localOutliers Function
#'
#' This function does detects local outliers based on k-nearest neighbors.
#'
#' @param spe SpatialExperiment object with the following columns in colData: sample_id, sum_umi, sum_gene
#' @param k Number of nearest neighbors to use for outlier detection
#' @param threshold Threshold for outlier detection
#' @return SpatialExperiment object with update z_umi_outlier columns
#' @export localOutliers
#' @examples
#' localOutliers(spe_example, k = 15, threshold = 2)
localOutliers <- function(spe, k = 36, features = c("sum_umi"), samples = "sample_id", log2 = TRUE, z_threshold = 3, output_z = FALSE) {
    # log2 transform specified features
    features_log2 <- character()
    if (log2) {
        for (feature in features) {
            feature_log2 <- paste0(feature, "_log2")
            colData(spe)[feature_log2] <- log2(colData(spe)[[feature]])
            features_log2 <- c(features_log2, feature_log2)
        }
    } else {
        features_log2 <- features
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
        dnn <- findKNN(spatialCoords(spe_subset), k = k)$index

        # Initialize a matrix to store z-scores for each feature
        mod_z_matrix <- matrix(NA, nrow(spaQC), length(features_log2))
        colnames(mod_z_matrix) <- features_log2

        # Loop through each row in the nearest neighbor index matrix
        for (i in 1:nrow(dnn)) {
            dnn.idx <- dnn[i, ]
            for (j in seq_along(features_log2)) {
                mod_z_matrix[i, j] <- modified_z(spaQC[c(i, dnn.idx[dnn.idx != 0]), ][[features_log2[j]]])[1]
            }
        }

        # Handle non-finite values
        mod_z_matrix[!is.finite(mod_z_matrix)] <- 0

        # Determine local outliers across all features
        spaQC$local_outliers <- as.factor(apply(mod_z_matrix, 1, function(x) any(x > z_threshold | x < -z_threshold)))

        # output z features if desired
        if (output_z) {
            for (j in seq_along(features_log2)) {
                feature_z <- paste0(features[j], "_z")
                spaQC[feature_z] <- mod_z_matrix[, j]
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

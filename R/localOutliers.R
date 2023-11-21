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
localOutliers <- function(spe, k = 36, feature='sum_umi', samples='sample_id', log2=TRUE, z_threshold = 3, output_z=FALSE) {

    # log2 transform the sum_umi and sum_gene features
    if (log2) {
      feature_log2 <- paste0(feature, '_log2')
      colData(spe)[feature_log2] <- log2(colData(spe)[[feature]])
      feature2use <- feature_log2
    } else {
      feature2use <- feature
    }

    # Get a list of unique sample IDs
    unique_sample_ids <- unique(colData(spe)[[samples]])

    # Initialize list variables to store the results
    mod_z <- vector("list", length(unique_sample_ids))

    # Initialize a list to store each spaQC dataframe
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

        # Initialize variables for the current sample
        mod_z[[sample_id]] <- rep(NA, nrow(spaQC))

        # Loop through each row in the nearest neighbor index matrix
        for(i in 1:nrow(dnn)) {
          dnn.idx <- dnn[i,]
          mod_z[[sample_id]][i] <- modified_z(spaQC[c(i, dnn.idx[dnn.idx != 0]),][[feature_log2]])[1]
        }

        # Handle non-finite values
        mod_z[[sample_id]][!is.finite(mod_z[[sample_id]])] <- 0

        # Save stats to the spaQC dataframe
        spaQC$local_outliers <- as.factor(ifelse(mod_z[[sample_id]] > z_threshold | mod_z[[sample_id]] < -z_threshold, TRUE, FALSE))

        # output z features if desired
        if (output_z) {
          feature_z <- paste0(feature, '_z')
          spaQC[feature_z] <- mod_z[[sample_id]]
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

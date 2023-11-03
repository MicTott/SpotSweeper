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
localOutliers <- function(spe, k = 15, threshold = 3) {

    # log2 transform the sum_umi and sum_gene features
    colData(spe)$sum_umi_log2 <- log2(spe$sum_umi)
    colData(spe)$sum_gene_log2 <- log2(spe$sum_gene)

    # Get a list of unique sample IDs
    unique_sample_ids <- unique(spe$sample_id)

    # Initialize list variables to store the results
    var.umi <- vector("list", length(unique_sample_ids))
    z.umi <- vector("list", length(unique_sample_ids))

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
        var.umi[[sample_id]] <- rep(NA, nrow(spaQC))
        z.umi[[sample_id]] <- rep(NA, nrow(spaQC))

        # Loop through each row in the nearest neighbor index matrix
        for(i in 1:nrow(dnn)) {
          dnn.idx <- dnn[i,]
          var.umi[[sample_id]][i] <- var(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_umi_log2, na.rm=TRUE)
          z.umi[[sample_id]][i] <- modified_z(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_umi_log2)[1]
        }

        # Handle non-finite values
        z.umi[[sample_id]][!is.finite(z.umi[[sample_id]])] <- 0

        # Save stats to the spaQC dataframe
        spaQC$var.umi <- var.umi[[sample_id]]
        spaQC$z.umi <- z.umi[[sample_id]]
        spaQC$z_umi_outlier <- ifelse(spaQC$z.umi > threshold | spaQC$z.umi < -threshold, TRUE, FALSE)

        # Store the modified spaQC dataframe in the list
        spaQC_list[[sample_id]] <- spaQC
    }

    # rbind the list of dataframes
    spaQC_aggregated <- do.call(rbind, spaQC_list)

    # replace column data
    colData(spe) <- spaQC_aggregated

    return(spe)
}

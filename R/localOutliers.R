#' localOutliers Function
#'
#' This function detects local outliers in spatial transcriptomics data based on standard
#' quality control metrics, such as library size, unique genes, and mitochondrial ratio.
#' Local outliers are defined as spots with low/high quality metrics compared to their
#' surrounding neighbors, based on a modified z-score statistic.
#'
#' @param spe SpatialExperiment object
#' @param metric colData QC metric to use for outlier detection
#' @param direction Direction of outlier detection (higher, lower, or both)
#' @param n_neighbors Number of nearest neighbors to use for outlier detection
#' @param samples Column name in colData to use for sample IDs
#' @param log Logical indicating whether to log2 transform the features (default is TRUE)
#' @param cutoff Cutoff for outlier detection (default is 3)
#'
#' @return SpatialExperiment object with updated colData containing outputs
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocNeighbors findKNN
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
#' # Calculating QC features for each spot using scuttle
#' spe<- scuttle::addPerCellQC(spe, subsets=list(Mito=is.mito))
#'
# Identifying local outliers suing SpotSweepR
#' spe<- localOutliers(spe,
#'                     metric="detected",
#'                     direction="lower"
#'                     )
localOutliers <- function(spe,
                          metric="detected",
                          direction="lower",
                          n_neighbors = 36,
                          samples = "sample_id",
                          log = TRUE,
                          cutoff = 3) {

    # log transform specified metric
    if (log) {
        metric_log <- paste0(metric, "_log2")
        colData(spe)[metric_log] <- log2(colData(spe)[[metric]])
        metric_to_use <- metric_log

    } else {
       metric_to_use <- metric
    }

    # Get a list of unique sample IDs
    unique_sample_ids <- unique(colData(spe)[[samples]])

    # Initialize list to store each spaQC dataframe
    spaQC_list <- vector("list", length(unique_sample_ids))

    # Loop through each unique sample ID
    for (sample in unique_sample_ids) {
        # Subset the data for the current sample
        spe_subset <- spe[ ,colData(spe)[[samples]] == sample]

        # Create a list of spatial coordinates and qc features
        spaQC <- colData(spe_subset)

        # Find nearest neighbors
        dnn <- BiocNeighbors::findKNN(spatialCoords(spe_subset),
                                      k = n_neighbors,
                                      warn.ties = FALSE)$index

        # Initialize a matrix to store z-scores
        mod_z_matrix <- array(NA, nrow(dnn))

        # Loop through each row in the nearest neighbor index matrix
        for (i in seq_len(nrow(dnn))) {
            dnn.idx <- dnn[i, ]
            mod_z_matrix[i] <- modified_z(spaQC[c(i, dnn.idx[dnn.idx != 0]), ][[metric_to_use]])[1]
        }

        # Handle non-finite values
        mod_z_matrix[!is.finite(mod_z_matrix)] <- 0

        # find outliers based on cutoff, store in colData
        metric_outliers <- paste0(metric, "_outliers")
        spaQC[metric_outliers]  <- switch(direction,
                              "higher" = as.factor(apply(mod_z_matrix, 1, function(x) x > cutoff)),
                              "lower" = as.factor(apply(mod_z_matrix, 1, function(x) x < -cutoff)),
                              "both" = as.factor(apply(mod_z_matrix, 1, function(x) x > cutoff | x < -cutoff))
                              )

        # add z-scores to colData
        metric_z <- paste0(metric, "_z")
        spaQC[metric_z] <- mod_z_matrix[]

        # Store the modified spaQC dataframe in the list
        spaQC_list[[sample]] <- spaQC
    }

    # rbind the list of dataframes
    spaQC_aggregated <- do.call(rbind, spaQC_list)

    # replace column data
    colData(spe) <- spaQC_aggregated

    return(spe)
}

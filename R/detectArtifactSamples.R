#' Detect Artifact Samples in Spatial Transcriptomics Data
#'
#' This function detects artifact samples in spatial transcriptomics data based on the bimodality of local variance
#' distributions. It leverages local variance calculations and the `is.bimodal` function from the LaplacesDemon package
#' to assess bimodality in the local variance of specified features across samples.
#'
#' @param spe A `SpatialExperiment` object containing the spatial transcriptomics data.
#' @param n_neighbors The number of neighbors to consider for local variance calculation. Default is 18.
#' @param features A character vector of column names in `colData` representing the features for which to calculate local variance.
#'        Default is `c("expr_chrM_ratio")`.
#' @param samples The name of the column in `colData` that contains sample identifiers. Default is "sample_id".
#' @param log2 A logical value indicating whether to apply a log2 transformation to the features before calculating local variance.
#'        Default is FALSE.
#' @param n_cores The number of cores to use for parallel processing. Default is 1.
#' @param name An optional name for the column in `colData` where the local variance results will be stored. If NULL,
#'        a default name will be generated based on the feature name and number of neighbors.
#'
#' @return A named list where each element corresponds to a sample and indicates whether the distribution
#'         of local variance for the specified features is bimodal (TRUE) or not (FALSE). The names of the list
#'         elements are the sample identifiers.
#'
#' @examples
#' bimodal_results <- detectArtifactSamples(spe,
#'                                          n_neighbors = 18,
#'                                          features = c("expr_chrM_ratio"),
#'                                          samples = "sample_id",
#'                                          log2 = FALSE,
#'                                          n_cores = 1,
#'                                          name = "local_variance")
#'
#' @importFrom LaplacesDemon is.bimodal
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocParallel MulticoreParam
#'
#'
#' @export
detectArtifactSamples <- function(spe, n_neighbors = 18, features = c("expr_chrM_ratio"),
                                  samples = "sample_id", log2 = FALSE, n_cores = 1, name=NULL) {

  # get sample names
  samples <- unique(colData(spe)[[samples]])

  # initialize list
  bimodal_logical <- list()

  # loop through samples
  for (sample in samples) {

    # make tmp.name for reference
    if (is.null(name)) {
      tmp.name <- paste0(features, "_k", n_neighbors)
    } else {
      tmp.name <- name
    }

    # subset and get local variance
    spe.subset <- subset(spe, ,sample_id == sample)
    spe.subset <- localVariance(spe.subset,
                                n_neighbors=n_neighbors,
                                n_cores=n_cores,
                                name=name,
                                log2=log2)

    # plot histogram for each sample
    hist(colData(spe.subset)[[name]], col =alpha('black',0.4), breaks=100)

    # is data bimodal?
    bimodal_logical[[sample]] <- is.bimodal(colData(spe.subset)[[name]])
  }


  # name bimodal list with sample IDs
  #names(bimodal_logical) <- samples

  return(bimodal_logical)
}

#' Identify and annotate artifacts in spatial transcriptomics data
#'
#' This function identifies and annotates potential artifacts in spatial transcriptomics data.
#' Artifacts are detected based on local mito variance, and the results are added to the original
#' SpatialExperiment (sce) object.
#'
#' @param spe A SingleCellExperiment object.
#' @param mito_percent The column name representing the mitochondrial percent. Default is "expr_chrM_ratio".
#' @param mito_sum The column name representing sum mitochondrial expression. Default is "expr_chrM".
#' @param samples The column name representing sample IDs. Default is "sample_id".
#' @param n_rings The number of rings for local mito variance calculation. Default is 5.
#' @param n_cores Number of cores to use for parallel processing. Default is 1.
#' @param log2 Logical, indicating whether to log2 transform specified features. Default is TRUE.
#' @param name Prefix for the local variance column names. Default is "artifact".
#' @param var_output Logical, indicating whether to include local variances in the output. Default is TRUE.
#'
#' @return Returns the modified SingleCellExperiment object with artifact annotations.
#'
#' @seealso
#' \code{\link{localVariance}}
#'
#' @import MASS
#' @import SingleCellExperiment
#' @import SpatialExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom stats prcomp kmeans resid var
#'
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
#' # find rtifacts
#' spe <- findArtifacts(spe,
#'                      mito_percent="subsets_Mito_percent",
#'                      mito_sum="subsets_Mito_sum",
#'                      n_rings=5,
#'                      name="artifact"
#'                    )
#'
#' @export
findArtifacts <- function(spe, mito_percent="expr_chrM_ratio",
                          mito_sum="expr_chrM", samples="sample_id",
                          n_rings=5, n_cores=1,
                          log2=TRUE, name="artifact", var_output=TRUE) {

  # log2 transform specified features
  features <- c(mito_percent, mito_sum)
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

  # get unique sample IDs
  sample_ids <- unique(colData(spe)[[samples]])

  # Initialize a list to store spe for each sample
  colData_list <- list()

  for (sample in sample_ids) {

    # subset by sample
    spe.temp <- spe[, colData(spe)[[samples]] == sample]

    # ======= Calculate local mito variance ========
    var_matrix <- c()
    for (i in 1:n_rings) {

      # ==== local mito variance ====
      # get n neighbors for i rings
      n_neighbors <- 3*i*(i+1)
      tmp.name <- paste0("k", n_neighbors)

      # get local variance of mito ratio
      spe.temp <- localVariance(spe.temp,
                                features=mito_percent,
                                n_neighbors=n_neighbors,
                                n_cores=n_cores,
                                name=tmp.name)

      # add to var_matrix
      var_matrix<- cbind(var_matrix, colData(spe.temp)[[tmp.name]])
    }


    # ========== PCA and clustering ==========
    # add mito and mito_percent to var dataframe
    var_matrix <- cbind(var_matrix,
                          colData(spe.temp)[[mito_percent]],
                          colData(spe.temp)[[mito_sum]]
    )

    var_df <- data.frame(var_matrix)

    # Run PCA and add to reduced dims
    pc <- prcomp(var_df, center = TRUE, scale. = TRUE)
    reducedDim(spe.temp, "PCA_artifacts") <- pc$x

    # Cluster using kmeans and add to temp sce
    clus <- kmeans(pc$x, centers = 2, nstart = 25)
    spe.temp$Kmeans <- clus$cluster

    # =========== Artifact annotation ===========

    # calculate average local variance of the two clusters
    clus1_mean <- mean(colData(spe.temp)[[paste0("k", 18)]][spe.temp$Kmeans==1])
    clus2_mean <- mean(colData(spe.temp)[[paste0("k", 18)]][spe.temp$Kmeans==2])

    artifact_clus <- which.min(c(clus1_mean, clus2_mean))

    # create a new $artifact column; if clus1 mean < clus2 - rename Kmeans 1 to "TRUE"
    spe.temp$artifact <- FALSE
    spe.temp$artifact[spe.temp$Kmeans==artifact_clus] <- TRUE

    # add to sample list
    colData_list[[sample]] <- colData(spe.temp)
  }

  # rbind the list of dataframes
  colData_aggregated <- do.call(rbind, colData_list)

  # replace SPE column data with aggregated data
  colData(spe) <- colData_aggregated

  # return sce
  return(spe)
}


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

  # make spe.temp
  spe.temp <- spe

  # ======= Calculate local mito variance ========
  resid_matrix <- c()
  for (i in 1:n_rings) {

    # ==== local mito variance ====
    # get n neighbors for i rings
    n_neighbors <- 3*i*(i+1)
    name <- paste0("k", n_neighbors)

    # get local variance of mito ratio
    spe.temp <- localVariance(spe.temp,
                         features=mito_percent,
                         n_neighbors=n_neighbors,
                         n_cores=n_cores,
                         name=name)

    # ==== Regress out mean-variance bias ====
    # data.frame containing all local var k
    mito_var_df <- data.frame(mito_var = log2(colData(spe.temp)[[name]]),
                              mito_mean = log2(colData(spe.temp)[[paste0(name,"_mean")]])
    )

    #  (Robust) Linear regression of mito varaince vs mean using IRLS
    fit.irls <- MASS::rlm(mito_var~mito_mean, mito_var_df)

    # get residuals
    resid.irls <- resid(fit.irls)
    colData(spe.temp)[[name]] <- resid.irls

    # add to resid_matrix
    resid_matrix<- cbind(resid_matrix, resid.irls)
  }


  # ========== PCA and clustering ==========
  # add mito and mito_percent to resid dataframe
  resid_matrix <- cbind(resid_matrix,
                     colData(spe)[[mito_percent]],
                     colData(spe)[[mito_sum]]
                     )

  resid_df <- data.frame(resid_matrix)

  # Run PCA and add to reduced dims
  pc <- prcomp(resid_df, center = TRUE, scale. = TRUE)
  reducedDim(spe.temp, "PCA_artifacts") <- pc$x

  # Cluster using kmeans and add to temp sce
  clus <- kmeans(pc$x, centers = 2, nstart = 25)
  spe.temp$Kmeans <- clus$cluster

  # =========== Artifact annotation ===========

  # calculate average local variance of the two clusters
  clus1_mean <- mean(colData(spe.temp)[[paste0("k", 18)]][spe.temp$Kmeans==1])
  print(clus1_mean)

  clus2_mean <- mean(colData(spe.temp)[[paste0("k", 18)]][spe.temp$Kmeans==2])
  clus2_mean

  artifact_clus <- which.min(c(clus1_mean, clus2_mean))

  # create a new $artifact column; if clus1 mean < clus2 - rename Kmeans 1 to "TRUE"
  spe.temp$artifact <- FALSE
  spe.temp$artifact[spe.temp$Kmeans==artifact_clus] <- TRUE

  # =========== Add to sce ===========
  # add artifact to sce
  colData(spe)$artifact <- spe.temp$artifact

  # add residuals to sce if desired
  if (var_output) {
    for (i in 1:n_rings) {
      name <- paste0("k", 3*i*(i+1))
      colData(spe)[[name]] <- colData(spe.temp)[[name]]
      reducedDim(spe, "PCA_artifact") <- pc$x
    }
  }

  # return sce
  return(spe)
}





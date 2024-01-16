#' findArtifacts Function
#'
#' This function uses the local variance of mitochondrial ratio to detect large artifacts
#' such as tissue hangnails.
#'
#' @return SpatialExperiment object with feature variance added to colData
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom MASS rlm
#'
#'
#' @export findArtifacts
#'
#' @examples
#' spe <- findArtifacts(spe,
#'                      mito_ratio="expr_chrM_ratio",
#'                      mito_expr="expr_chrM",
#'                      name="artifact"
#'                    )
#'


findArtifacts <- function(spe, mito_ratio="expr_chrM_ratio",
                          mito_expr="expr_chrM", samples="sample_id",
                          n_rings=5, n_cores=1,
                          log2=TRUE, name="artifact", var_output=TRUE) {

  # log2 transform specified features
  features <- c(mito_ratio, mito_expr)
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
  # add mito and mito_ratio to resid dataframe
  resid_matrix <- cbind(resid_matrix,
                     colData(spe)[[mito_ratio]],
                     colData(spe)[[mito_expr]]
                     )

  resid_df <- data.frame(resid_matrix)

  # Run PCA and add to reduced dims
  pc <- prcomp(resid_df, center = TRUE, scale. = TRUE)
  reducedDim(spe.temp, "PCA_artifacts") <- pc$x

  # Cluster using kmeans and add to temp sce
  clus <- kmeans(pc$x, centers = 2, nstart = 25)
  spe.temp$Kmeans <- clus$cluster

  # =========== Artifact annotation ===========
  print(unique(colnames(colData(spe.temp))))
  print(unique(spe.temp$Kmeans))

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
    }
  }

  # return sce
  return(spe)
}





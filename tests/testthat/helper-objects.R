make_test_spe <- function(
    values = c(100, 1, 2, 3, 4, 5),
    samples = rep("sample", length(values)),
    coordinates = cbind(x = seq_along(values), y = 0)) {
  spot_names <- paste0("spot", seq_along(values))
  counts <- matrix(
    0,
    nrow = 1,
    ncol = length(values),
    dimnames = list("gene", spot_names)
  )
  rownames(coordinates) <- spot_names
  metadata <- S4Vectors::DataFrame(
    sample_id = samples,
    metric = values,
    row.names = spot_names
  )

  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = metadata,
    spatialCoords = coordinates
  )
}

make_test_seurat <- function(with_image = FALSE) {
  spot_names <- paste0("spot", seq_len(6))
  counts <- matrix(
    seq_len(24),
    nrow = 4,
    dimnames = list(paste0("gene", seq_len(4)), spot_names)
  )
  metadata <- data.frame(
    sample_id = "sample",
    metric = c(100, 1, 2, 3, 4, 5),
    row.names = spot_names
  )
  object <- suppressWarnings(
    SeuratObject::CreateSeuratObject(
      counts,
      assay = if (with_image) "Spatial" else "RNA",
      meta.data = metadata
    )
  )

  if (with_image) {
    coordinates <- data.frame(
      tissue = 1,
      row = c(38, 1, 2, 3, 4, 5),
      col = c(88, 2, 4, 6, 8, 10),
      imagerow = seq(0, 50, by = 10),
      imagecol = seq(0, 100, by = 20),
      row.names = spot_names
    )
    image <- methods::new(
      "VisiumV1",
      image = array(0, dim = c(2, 2, 3)),
      scale.factors = Seurat::scalefactors(
        spot = 1,
        fiducial = 1,
        hires = 1,
        lowres = 1
      ),
      coordinates = coordinates,
      spot.radius = 1,
      assay = "Spatial",
      key = "slice1_"
    )
    object[["slice1"]] <- image
  }

  object
}

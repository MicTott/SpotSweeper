test_that("metadata adapters preserve SpatialExperiment ordering", {
  spe <- make_test_spe()
  metadata <- getMetadata(spe)
  metadata$new_column <- seq_len(nrow(metadata))
  metadata <- metadata[rev(rownames(metadata)), , drop = FALSE]

  result <- setMetadata(spe, metadata)

  expect_s4_class(result, "SpatialExperiment")
  expect_identical(rownames(getMetadata(result)), colnames(spe))
  expect_identical(result$new_column, seq_len(ncol(spe)))
})

test_that("Seurat metadata adapters use public object accessors", {
  skip_if_not_installed("SeuratObject")
  object <- make_test_seurat()
  metadata <- getMetadata(object)
  metadata$new_column <- seq_len(nrow(metadata))
  metadata <- metadata[rev(rownames(metadata)), , drop = FALSE]

  result <- setMetadata(object, metadata)
  subset <- subsetSpatialObject(result, c("spot2", "spot4"))

  expect_s4_class(result, "Seurat")
  expect_identical(rownames(getMetadata(result)), colnames(object))
  expect_identical(unname(result$new_column), seq_len(ncol(object)))
  expect_identical(colnames(subset), c("spot2", "spot4"))
})

test_that("localOutliers works on a real Seurat object", {
  skip_if_not_installed("SeuratObject")
  object <- make_test_seurat()
  coordinates <- cbind(x = seq_len(ncol(object)), y = 0)
  rownames(coordinates) <- colnames(object)

  result <- localOutliers(
    object,
    metric = "metric",
    direction = "higher",
    n_neighbors = 4,
    log = FALSE,
    coords = coordinates
  )

  expect_s4_class(result, "Seurat")
  expect_true(result$metric_outliers[[1]])
  expect_true(all(c("metric_z", "metric_outliers") %in%
                  colnames(getMetadata(result))))
})

test_that("Seurat spatial-image coordinates are extracted and aligned", {
  skip_if_not_installed("Seurat")
  object <- make_test_seurat(with_image = TRUE)

  coordinates <- getSpatialCoords(object)
  result <- localOutliers(
    object,
    metric = "metric",
    direction = "higher",
    n_neighbors = 4,
    log = FALSE
  )

  expect_identical(rownames(coordinates), colnames(object))
  expect_identical(colnames(coordinates), c("x", "y"))
  expect_true(result$metric_outliers[[1]])
})

test_that("flagVisiumOutliers supports Seurat Visium images", {
  skip_if_not_installed("Seurat")
  object <- make_test_seurat(with_image = TRUE)

  result <- flagVisiumOutliers(object, image_id = "slice1")

  expect_s4_class(result, "Seurat")
  expect_true(result$systematic_outliers[[1]])
  expect_false(any(result$systematic_outliers[-1]))
})

test_that("flagVisiumOutliers preserves SpatialExperiment support", {
  spe <- make_test_spe()
  spe$array_row <- c(38, 1, 2, 3, 4, 5)
  spe$array_col <- c(88, 2, 4, 6, 8, 10)

  result <- flagVisiumOutliers(spe)

  expect_true(result$systematic_outliers[[1]])
  expect_false(any(result$systematic_outliers[-1]))
})

test_that("compatibility adapters reject unsupported objects", {
  invalid <- list()

  expect_false(SpotSweeper:::is_seurat(invalid))
  expect_false(SpotSweeper:::is_spatial_experiment(invalid))
  expect_error(getSpatialCoords(invalid), "SpatialExperiment or Seurat")
  expect_error(getMetadata(invalid), "SpatialExperiment or Seurat")
  expect_error(setMetadata(invalid, data.frame()), "SpatialExperiment or Seurat")
})

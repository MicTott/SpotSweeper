test_that("custom coordinate rows are aligned by spot name", {
  spe <- make_test_spe()
  custom_coordinates <- SpatialExperiment::spatialCoords(spe)
  custom_coordinates <- custom_coordinates[rev(rownames(custom_coordinates)), ]

  expected <- localOutliers(
    spe,
    metric = "metric",
    n_neighbors = 4,
    log = FALSE
  )
  observed <- localOutliers(
    spe,
    metric = "metric",
    n_neighbors = 4,
    log = FALSE,
    coords = custom_coordinates
  )

  expect_equal(observed$metric_z, expected$metric_z)
})

test_that("custom coordinates are validated", {
  spe <- make_test_spe()

  expect_error(
    localOutliers(spe, metric = "metric", coords = seq_len(6)),
    "numeric matrix"
  )
  expect_error(
    localOutliers(
      spe,
      metric = "metric",
      coords = matrix(seq_len(4), ncol = 2)
    ),
    "one row for every spot"
  )

  character_coordinates <- matrix(letters[seq_len(12)], ncol = 2)
  expect_error(
    localOutliers(spe, metric = "metric", coords = character_coordinates),
    "numeric matrix"
  )

  mismatched_coordinates <- matrix(seq_len(12), ncol = 2)
  rownames(mismatched_coordinates) <- paste0("other", seq_len(6))
  expect_error(
    localOutliers(spe, metric = "metric", coords = mismatched_coordinates),
    "row names"
  )
})

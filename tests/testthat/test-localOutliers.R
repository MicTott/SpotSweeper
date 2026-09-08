test_that("modified z-score uses one unscaled MAD", {
  neighbors <- c(1, 2, 3, 4, 5)
  expected <- stats::qnorm(0.75) * (100 - 3) / 1

  expect_equal(
    SpotSweeper:::.local_modified_z(100, neighbors),
    expected
  )
})

test_that("modified z-score handles zero neighborhood MAD", {
  expect_equal(SpotSweeper:::.local_modified_z(10, rep(1, 4)), 0)
})

test_that("localOutliers scores the focal spot against its neighbors", {
  spe <- make_test_spe()
  result <- localOutliers(
    spe,
    metric = "metric",
    direction = "higher",
    n_neighbors = 4,
    log = FALSE
  )

  expected <- stats::qnorm(0.75) * (100 - 2.5)
  expect_s4_class(result, "SpatialExperiment")
  expect_equal(result$metric_z[[1]], expected)
  expect_true(result$metric_outliers[[1]])
  expect_identical(colnames(result), colnames(spe))
})

test_that("outlier directions use the requested cutoff", {
  spe <- make_test_spe(values = c(100, 1, 2, 3, 4, -100))

  higher <- localOutliers(
    spe,
    metric = "metric",
    direction = "higher",
    n_neighbors = 4,
    log = FALSE
  )
  lower <- localOutliers(
    spe,
    metric = "metric",
    direction = "lower",
    n_neighbors = 4,
    log = FALSE
  )
  both <- localOutliers(
    spe,
    metric = "metric",
    direction = "both",
    n_neighbors = 4,
    log = FALSE
  )

  expect_true(higher$metric_outliers[[1]])
  expect_false(higher$metric_outliers[[6]])
  expect_false(lower$metric_outliers[[1]])
  expect_true(lower$metric_outliers[[6]])
  expect_true(all(both$metric_outliers[c(1, 6)]))
})

test_that("log transformation is retained in metadata", {
  spe <- make_test_spe()
  result <- localOutliers(
    spe,
    metric = "metric",
    n_neighbors = 4,
    log = TRUE
  )

  expect_equal(result$metric_log, log1p(spe$metric))
  expect_true(all(is.finite(result$metric_z)))
})

test_that("interleaved sample order is preserved", {
  values <- c(100, 100, 1, 1, 2, 2, 3, 3)
  samples <- rep(c("a", "b"), 4)
  spe <- make_test_spe(values = values, samples = samples)
  result <- localOutliers(
    spe,
    metric = "metric",
    n_neighbors = 3,
    log = FALSE
  )

  expect_identical(colnames(result), colnames(spe))
  expect_identical(as.character(result$sample_id), samples)
})

test_that("localOutliers validates inputs", {
  spe <- make_test_spe()

  expect_error(localOutliers(list()), "SpatialExperiment or Seurat")
  expect_error(localOutliers(spe, metric = "missing"), "Required columns")
  expect_error(localOutliers(spe, n_neighbors = 0), "positive integer")
  expect_error(localOutliers(spe, cutoff = -1), "non-negative")
  expect_error(
    localOutliers(spe, metric = "metric", n_neighbors = ncol(spe)),
    "more spots than 'n_neighbors'"
  )

  spe$metric <- as.character(spe$metric)
  expect_error(
    localOutliers(spe, metric = "metric"),
    "numeric metadata column"
  )
})

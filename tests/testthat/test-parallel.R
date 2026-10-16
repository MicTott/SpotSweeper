test_that("worker backends are selected automatically", {
  for (os_type in c("windows", "unix")) {
    expect_s4_class(.parallel_param(1, os_type), "SerialParam")
  }
  sockets <- .parallel_param(2, "windows")
  expect_s4_class(sockets, "SnowParam")
  expect_equal(BiocParallel::bpworkers(sockets), 2L)
  if (.Platform$OS.type != "windows") {
    expect_s4_class(.parallel_param(2, "unix"), "MulticoreParam")
  }
})

test_that("both functions agree between serial and parallel search", {
  set.seed(27)
  spe <- make_test_spe(
    values = stats::runif(40, 1, 100),
    coordinates = matrix(stats::runif(80), ncol = 2)
  )
  serial_outliers <- localOutliers(
    spe, metric = "metric", n_neighbors = 6, workers = 1
  )
  serial_variance <- localVariance(
    spe, metric = "metric", n_neighbors = 6, workers = 1
  )

  check_parallel <- function() {
    parallel_outliers <- localOutliers(
      spe, metric = "metric", n_neighbors = 6, workers = 2
    )
    parallel_variance <- localVariance(
      spe, metric = "metric", n_neighbors = 6, workers = 2
    )
    expect_equal(parallel_outliers$metric_z, serial_outliers$metric_z)
    expect_identical(
      parallel_outliers$metric_outliers, serial_outliers$metric_outliers
    )
    expect_equal(parallel_variance$metric_var, serial_variance$metric_var)
    expect_identical(colnames(parallel_outliers), colnames(spe))
    expect_identical(colnames(parallel_variance), colnames(spe))
  }

  check_parallel()
  # Exercise actual socket processes on Unix as well as in Windows CI.
  original_param <- .parallel_param
  local_mocked_bindings(
    .parallel_param = function(workers) original_param(workers, "windows")
  )
  check_parallel()
})

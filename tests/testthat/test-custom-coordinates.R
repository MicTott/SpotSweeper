test_that("localOutliers works with custom coordinates", {
  
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("STexampleData")
  
  # Load test data
  spe <- STexampleData::Visium_humanDLPFC()
  spe <- spe[1:100, 1:50]  # Small subset for testing
  
  # Test default behavior (spatial coordinates)
  expect_no_error(
    spe_default <- localOutliers(spe, 
                                metric = "in_tissue",
                                direction = "lower",
                                n_neighbors = 10)
  )
  
  # Test with custom coordinates (random for testing)
  set.seed(42)
  custom_coords <- matrix(rnorm(ncol(spe) * 2), ncol = 2)
  rownames(custom_coords) <- colnames(spe)
  
  expect_no_error(
    spe_custom <- localOutliers(spe,
                               metric = "in_tissue", 
                               direction = "lower",
                               n_neighbors = 10,
                               coords = custom_coords)
  )
  
  # Check that results are different (since neighborhoods are different)
  default_z <- getMetadata(spe_default)$in_tissue_z
  custom_z <- getMetadata(spe_custom)$in_tissue_z
  
  # Results should be different due to different neighborhood definitions
  expect_false(identical(default_z, custom_z))
  
  # Both should have the same structure
  expect_equal(length(default_z), length(custom_z))
  expect_equal(length(default_z), ncol(spe))
})

test_that("localOutliers coords validation works", {
  
  skip_if_not_installed("SpatialExperiment") 
  skip_if_not_installed("STexampleData")
  
  spe <- STexampleData::Visium_humanDLPFC()
  spe <- spe[1:50, 1:20]
  
  # Test invalid coords - not a matrix
  expect_error(
    localOutliers(spe, coords = c(1, 2, 3)),
    "'coords' must be a numeric matrix"
  )
  
  # Test wrong dimensions
  wrong_coords <- matrix(rnorm(10), ncol = 2)
  expect_error(
    localOutliers(spe, coords = wrong_coords),
    "same number of rows as spots"
  )
  
  # Test non-numeric coords
  char_coords <- matrix(letters[1:40], ncol = 2)
  expect_error(
    localOutliers(spe, coords = char_coords),
    "'coords' must be a numeric matrix"
  )
})
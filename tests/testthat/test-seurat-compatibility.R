test_that("Seurat compatibility layer works", {
  
  # Skip if Seurat not available
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("STexampleData")
  
  # Test with SpatialExperiment first
  spe <- STexampleData::Visium_humanDLPFC()
  spe <- spe[1:100, 1:50]  # Small subset for testing
  
  # Test compatibility layer functions with SpatialExperiment
  expect_true(is_spatial_experiment(spe))
  expect_false(is_seurat(spe))
  
  coords_spe <- getSpatialCoords(spe)
  expect_true(is.matrix(coords_spe))
  expect_equal(ncol(coords_spe), 2)
  
  metadata_spe <- getMetadata(spe)
  expect_true(is.data.frame(metadata_spe))
  expect_equal(nrow(metadata_spe), ncol(spe))
  
  # Test validateMetadataColumns
  expect_true(validateMetadataColumns(spe, "in_tissue"))
  expect_error(validateMetadataColumns(spe, "nonexistent_column"))
})

test_that("Compatibility layer error handling works", {
  
  # Test with invalid object
  invalid_obj <- list(test = "data")
  
  expect_false(is_seurat(invalid_obj))
  expect_false(is_spatial_experiment(invalid_obj))
  
  expect_error(getSpatialCoords(invalid_obj))
  expect_error(getMetadata(invalid_obj))
  expect_error(setMetadata(invalid_obj, data.frame()))
})
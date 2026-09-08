# Changelog

## Version 1.9.1

### Bug fixes

- Corrected the robust local z-score calculation in
  [`localOutliers()`](https://mictott.github.io/SpotSweeper/reference/localOutliers.md).
  Scores now use the focal spot value and the unscaled median absolute
  deviation of neighboring spots, matching the published method.
  Previous versions applied the normal-consistency scaling twice and
  could associate a neighbor’s score with the focal spot.
- Removed the `spatialEco` dependency, which is no longer needed for
  local z-score calculation.

### Seurat compatibility

- Reworked Seurat support around small, tested object adapters that
  preserve spot names and ordering and use public `SeuratObject`
  accessors.
- Added integration tests for
  [`localOutliers()`](https://mictott.github.io/SpotSweeper/reference/localOutliers.md),
  spatial-coordinate extraction, metadata replacement, subsetting, and
  [`flagVisiumOutliers()`](https://mictott.github.io/SpotSweeper/reference/flagVisiumOutliers.md)
  with real Seurat objects.
- Added explicit guidance when
  [`findArtifacts()`](https://mictott.github.io/SpotSweeper/reference/findArtifacts.md)
  is used with data lacking mitochondrial signal, including Visium FFPE
  probe-panel data.

## Version 1.5.0

### Major Features

- **Seurat Compatibility**: Added comprehensive compatibility layer
  enabling SpotSweeper functions to work seamlessly with Seurat spatial
  objects alongside existing SpatialExperiment support. This major
  enhancement expands SpotSweeper’s usability across the spatial
  transcriptomics ecosystem.
- **Flexible Coordinate Systems**: Enhanced
  [`localOutliers()`](https://mictott.github.io/SpotSweeper/reference/localOutliers.md)
  to support custom coordinate systems for neighborhood detection,
  enabling outlier detection in PCA, UMAP, or any reduced-dimension
  space.

### New Functions

- **[`getSpatialCoords()`](https://mictott.github.io/SpotSweeper/reference/getSpatialCoords.md)**:
  Universal function to extract spatial coordinates from both
  SpatialExperiment and Seurat objects
- **[`getMetadata()`](https://mictott.github.io/SpotSweeper/reference/getMetadata.md)**:
  Universal function to access metadata (colData/meta.data) from both
  object types
- **[`setMetadata()`](https://mictott.github.io/SpotSweeper/reference/setMetadata.md)**:
  Universal function to update metadata in both object types
- **[`subsetSpatialObject()`](https://mictott.github.io/SpotSweeper/reference/subsetSpatialObject.md)**:
  Universal subsetting function for both object types
- **[`validateMetadataColumns()`](https://mictott.github.io/SpotSweeper/reference/validateMetadataColumns.md)**:
  Validation function to ensure required columns exist
- **`is_seurat()`** and **`is_spatial_experiment()`**: Object type
  detection functions

### Enhanced Functions

- **[`localOutliers()`](https://mictott.github.io/SpotSweeper/reference/localOutliers.md)**:
  - Now accepts both SpatialExperiment and Seurat objects with automatic
    detection
  - New `coords` parameter allows custom coordinate matrices for
    neighborhood detection
  - Default behavior unchanged (uses spatial coordinates)
  - Supports PCA, UMAP, t-SNE, or any coordinate system for neighbor
    finding

### Bug Fixes

- **Fixed focal spot inclusion**: Removed focal spot from neighborhood
  when calculating modified z-scores in
  [`localOutliers()`](https://mictott.github.io/SpotSweeper/reference/localOutliers.md)
  to prevent circular logic and bias, resulting in more accurate outlier
  detection.

### Documentation

- Updated function documentation to reflect dual compatibility and new
  coordinate options
- Added Seurat usage examples alongside existing SpatialExperiment
  examples
- Added examples for custom coordinate systems (PCA, UMAP)
- Comprehensive test suite for compatibility layer and coordinate
  functionality

## Version 1.3.2

### New Features

- **Added** the ‘flagVisiumOutliers()’ function to identify and flag
  systematic outlier spots in Visium datasets. This feature enhances
  data quality by allowing users to efficiently detect and exclude
  problematic spots from downstream analyses.

### Minor Changes

- **Broadened Compatibility**: Updated all functions to use
  `inherits(spe, "SpatialExperiment")` instead of checking `class(spe)`
  directly. This change ensures that derived classes (e.g.,
  `SpatialFeatureExperiment`) are also supported, improving flexibility
  and ease of use.

## Version 1.3.1

### Major Changes

- **Function Renaming**: The function `plotQC` has been renamed to
  `plotQCmetrics` to better reflect its purpose. The new function
  `plotQCmetrics` should be used moving forward. This change improves
  clarity in the package’s API by specifying that this function is
  designed for plotting QC metrics.

### New Features and Enhancements

- **`shape` argument**: Added a `shape` argument to `findArtifacts`,
  allowing users to specify the neighborhood shape as either
  `"hexagonal"` or `"square"` for local variance calculations. This
  enhancement provides flexibility for different spatial arrangements in
  spatial transcriptomics data.

- **Updated `n_order` parameter**: Renamed the `n_rings` parameter to
  `n_order` in the `findArtifacts` function to better describe its
  purpose of specifying the N-order neighbors for local mitochondrial
  variance calculations.

- **Parallelization**: Added a `workers` argument for parallel
  processing using `BiocParallel` in both `localOutlier` and
  `localVariance` functions. This allows for faster computation,
  particularly on larger datasets.

### Deprecations

- **`plotQC` Function Deprecated**: The `plotQC` function is now
  deprecated. While it remains available for backward compatibility,
  users are encouraged to transition to `plotQCmetrics`. Calling
  `plotQC` will display a warning, reminding users of the deprecation.

This change is backward compatible; existing code using `plotQC` will
still work but will show a warning. We recommend updating your code to
use `plotQCmetrics` to avoid any issues in future versions where
`plotQC` may be removed.

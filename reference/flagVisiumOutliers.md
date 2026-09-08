# Flag systematic Visium outlier positions

`flagVisiumOutliers()` marks six Visium array positions that have shown
systematically low library sizes across multiple datasets. It accepts
either a `SpatialExperiment` or `Seurat` object and adds a logical
`systematic_outliers` column to its metadata.

## Usage

``` r
flagVisiumOutliers(spe, image_id = NULL)
```

## Arguments

- spe:

  A `SpatialExperiment` or `Seurat` object containing Visium array
  coordinates. `SpatialExperiment` inputs must contain `array_row` and
  `array_col` metadata columns. Seurat inputs may contain the same
  columns or expose them through a Visium image.

- image_id:

  For Seurat inputs, the image containing Visium array coordinates. The
  first image is used when `NULL`.

## Value

`spe` with a logical `systematic_outliers` metadata column.

## Examples

``` r
spe <- STexampleData::Visium_humanDLPFC()
#> 
#> see ?STexampleData and browseVignettes('STexampleData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
spe <- flagVisiumOutliers(spe)
spe <- spe[, !spe$systematic_outliers]
```

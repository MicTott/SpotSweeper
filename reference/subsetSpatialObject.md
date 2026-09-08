# Subset a spatial object by spots

Subset a spatial object by spots

## Usage

``` r
subsetSpatialObject(x, indices, column_name = NULL)
```

## Arguments

- x:

  A `SpatialExperiment` or `Seurat` object.

- indices:

  Logical, integer, or character spot indices. If `column_name` is
  supplied, these are values to retain from that column.

- column_name:

  Optional metadata column used for value-based subsetting.

## Value

A subset of `x` containing the selected spots.

## Examples

``` r
spe <- STexampleData::Visium_humanDLPFC()
#> see ?STexampleData and browseVignettes('STexampleData') for documentation
#> loading from cache
spe_subset <- subsetSpatialObject(spe, seq_len(100))
```

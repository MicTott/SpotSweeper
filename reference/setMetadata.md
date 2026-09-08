# Replace column metadata in a spatial object

Replace column metadata in a spatial object

## Usage

``` r
setMetadata(x, metadata)
```

## Arguments

- x:

  A `SpatialExperiment` or `Seurat` object.

- metadata:

  A data-frame-like object with one row per spot. Named rows are
  reordered to match the object before assignment.

## Value

`x` with replaced column metadata.

## Examples

``` r
spe <- STexampleData::Visium_humanDLPFC()
#> see ?STexampleData and browseVignettes('STexampleData') for documentation
#> loading from cache
metadata <- getMetadata(spe)
metadata$example_column <- seq_len(nrow(metadata))
spe <- setMetadata(spe, metadata)
```

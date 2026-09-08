# Validate columns in spatial-object metadata

Validate columns in spatial-object metadata

## Usage

``` r
validateMetadataColumns(x, required_columns)
```

## Arguments

- x:

  A `SpatialExperiment` or `Seurat` object.

- required_columns:

  Character vector of required metadata column names.

## Value

`TRUE`, invisibly. An error is raised if columns are missing.

## Examples

``` r
spe <- STexampleData::Visium_humanDLPFC()
#> see ?STexampleData and browseVignettes('STexampleData') for documentation
#> loading from cache
validateMetadataColumns(spe, c("sample_id", "in_tissue"))
```

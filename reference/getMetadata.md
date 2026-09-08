# Get column metadata from a spatial object

Get column metadata from a spatial object

## Usage

``` r
getMetadata(x)
```

## Arguments

- x:

  A `SpatialExperiment` or `Seurat` object.

## Value

A `data.frame` with one row per spot, in object order.

## Examples

``` r
spe <- STexampleData::Visium_humanDLPFC()
#> see ?STexampleData and browseVignettes('STexampleData') for documentation
#> loading from cache
metadata <- getMetadata(spe)
head(metadata)
#>                            barcode_id     sample_id in_tissue array_row
#> AAACAACGAATAGTTC-1 AAACAACGAATAGTTC-1 sample_151673         0         0
#> AAACAAGTATCTCCCA-1 AAACAAGTATCTCCCA-1 sample_151673         1        50
#> AAACAATCTACTAGCA-1 AAACAATCTACTAGCA-1 sample_151673         1         3
#> AAACACCAATAACTGC-1 AAACACCAATAACTGC-1 sample_151673         1        59
#> AAACAGAGCGACTCCT-1 AAACAGAGCGACTCCT-1 sample_151673         1        14
#> AAACAGCTTTCAGAAG-1 AAACAGCTTTCAGAAG-1 sample_151673         1        43
#>                    array_col ground_truth reference cell_count
#> AAACAACGAATAGTTC-1        16         <NA>      <NA>         NA
#> AAACAAGTATCTCCCA-1       102       Layer3    Layer3          6
#> AAACAATCTACTAGCA-1        43       Layer1    Layer1         16
#> AAACACCAATAACTGC-1        19           WM        WM          5
#> AAACAGAGCGACTCCT-1        94       Layer3    Layer3          2
#> AAACAGCTTTCAGAAG-1         9       Layer5    Layer5          4
```

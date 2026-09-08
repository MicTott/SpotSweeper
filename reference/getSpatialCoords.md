# Get spatial coordinates from a spatial object

Get spatial coordinates from a spatial object

## Usage

``` r
getSpatialCoords(x, image_id = NULL)
```

## Arguments

- x:

  A `SpatialExperiment` or `Seurat` object.

- image_id:

  For a Seurat object, an optional image name. If omitted, coordinates
  from all images are combined and returned in object order.

## Value

A numeric matrix with spots as rows and spatial dimensions as columns.

## Examples

``` r
spe <- STexampleData::Visium_humanDLPFC()
#> see ?STexampleData and browseVignettes('STexampleData') for documentation
#> loading from cache
coordinates <- getSpatialCoords(spe)
head(coordinates)
#>                    pxl_col_in_fullres pxl_row_in_fullres
#> AAACAACGAATAGTTC-1               3913               2435
#> AAACAAGTATCTCCCA-1               9791               8468
#> AAACAATCTACTAGCA-1               5769               2807
#> AAACACCAATAACTGC-1               4068               9505
#> AAACAGAGCGACTCCT-1               9271               4151
#> AAACAGCTTTCAGAAG-1               3393               7583
```

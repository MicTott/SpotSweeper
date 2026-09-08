# Biased Spots Data

The `biased_spots` dataset is a `data.frame` containing information
about specific spatial spots identified as technical outliers in spatial
transcriptomics experiments. Each entry represents a biased spot
characterized by its spatial coordinates (row and column) and a unique
barcode. This dataset is utilized by the `flagVisiumOutliers` function
to flag and exclude these outlier spots from downstream analyses,
thereby enhancing data quality and reliability.

## Usage

``` r
data(biased_spots)
```

## Format

A `data.frame` with the following columns:

- row:

  Integer. The row position of a biased spot within the spatial grid.

- col:

  Integer. The column position of a biased spot within the spatial grid.

- barcode:

  Character. A unique identifier corresponding to the spatial
  transcriptomics barcode of the biased spot.

## Source

The `biased_spots.rds` file was generated in the analysis of local
outliers. See
https://github.com/boyiguo1/Manuscript-SpotSweeper/blob/main/code/03_Visium/figure_3.R
for more details.

## Examples

``` r
data(biased_spots)
```

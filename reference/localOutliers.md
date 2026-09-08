# Detect local outliers in spatial transcriptomics data

`localOutliers()` compares each spot's quality-control metric with the
same metric in its nearest neighbors. It supports `SpatialExperiment`
and `Seurat` objects and stores the resulting robust z-scores and
outlier calls in the object's column metadata.

## Usage

``` r
localOutliers(
  spe,
  metric = "detected",
  direction = "lower",
  n_neighbors = 36,
  samples = "sample_id",
  log = TRUE,
  cutoff = 3,
  workers = 1,
  coords = NULL
)
```

## Arguments

- spe:

  A `SpatialExperiment` or `Seurat` object.

- metric:

  A single column name in the object's column metadata containing a
  numeric QC metric.

- direction:

  Direction of outlier detection: `"higher"`, `"lower"`, or `"both"`.

- n_neighbors:

  Number of nearest neighbors used as the local reference.

- samples:

  A single metadata column name containing sample identifiers.

- log:

  Whether to apply [`log1p()`](https://rdrr.io/r/base/Log.html) to
  `metric` before scoring.

- cutoff:

  Non-negative absolute robust z-score cutoff.

- workers:

  Number of workers passed to
  [`BiocNeighbors::findKNN()`](https://rdrr.io/pkg/BiocNeighbors/man/findKNN.html).

- coords:

  Optional numeric coordinate matrix with spots as rows. When `NULL`,
  spatial coordinates are obtained from `spe`.

## Value

The input object with `<metric>_z` and `<metric>_outliers` metadata
columns. When `log = TRUE`, `<metric>_log` is also added.

## Details

The robust z-score for spot i is
`qnorm(0.75) * (x_i - median(x_neighbors)) / MAD_unscaled(x_neighbors)`.
The focal spot is not included when estimating the neighborhood median
or median absolute deviation. A score of zero is returned when the
neighborhood MAD is zero or non-finite.

## Examples

``` r
library(SpotSweeper)
library(SpatialExperiment)

spe <- STexampleData::Visium_humanDLPFC()
#> see ?STexampleData and browseVignettes('STexampleData') for documentation
#> loading from cache
rownames(spe) <- rowData(spe)$gene_name
spe <- spe[, spe$in_tissue == 1]
spe <- spe[, !is.na(spe$ground_truth)]

is.mito <- grepl("^MT-", rownames(spe))
spe <- scuttle::addPerCellQCMetrics(
    spe,
    subsets = list(Mito = is.mito)
)
#> Warning: 'scuttle::addPerCellQCMetrics' is deprecated.
#> Use 'scrapper::quickRnaQc.se' instead.
#> See help("Deprecated")
#> Warning: 'perCellQCMetrics' is deprecated.
#> Use 'scrapper::computeRnaQcMetrics' instead.
#> See help("Deprecated")
#> using unknown matrix fallback for 'dgTMatrix'

spe <- localOutliers(
    spe,
    metric = "sum",
    direction = "lower",
    log = TRUE
)

# Seurat objects use the same interface. Custom coordinates can also be
# supplied for either object type.
if (requireNamespace("SeuratObject", quietly = TRUE)) {
    counts <- matrix(
        seq_len(24),
        nrow = 4,
        dimnames = list(paste0("gene", 1:4), paste0("spot", 1:6))
    )
    seurat <- SeuratObject::CreateSeuratObject(counts)
    seurat$sample_id <- "sample"
    seurat$detected <- c(100, 10, 11, 9, 12, 10)
    example_coords <- cbind(x = seq_len(6), y = 0)
    seurat <- localOutliers(
        seurat,
        n_neighbors = 4,
        coords = example_coords
    )
}
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
```

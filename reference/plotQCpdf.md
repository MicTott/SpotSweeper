# Plot Outlier Metrics to PDF

This function generates a PDF file containing plots for each sample in
the SpatialExperiment object, highlighting outliers based on specified
metrics. Each plot visualizes outlier metrics for a single sample,
allowing for easy comparison and analysis across samples.

## Usage

``` r
plotQCpdf(
  spe,
  sample_id = "sample_id",
  metric = "detected",
  outliers = "local_outliers",
  colors = c("white", "black"),
  stroke = 1,
  point_size = 2,
  width = 5,
  height = 5,
  fname
)
```

## Arguments

- spe:

  A SpatialExperiment object containing the data to be plotted.

- sample_id:

  A character string specifying the column name in `colData(spe)` that
  contains unique sample identifiers. Default is 'sample_id'.

- metric:

  A character string specifying the metric to be visualized in the plot.
  This metric should be a column name in `colData(spe)`.

- outliers:

  A character string specifying the column name in `colData(spe)` that
  indicates whether a data point is considered an outlier. Default is
  local_outliers'.

- colors:

  A character vector specifying the colors to be used for the gradient
  scale. If length is 2, the gradient will be a single color gradient

- stroke:

  A numeric value specifying the border thickness for outlier points.
  Default is 1.

- point_size:

  A numeric value specifying the size of the points in the plot. Default
  is 2.

- width:

  A numeric value indicating the width of the plot. Default is 5.

- height:

  A numeric value indicating the height of the plot. Default is 5.

- fname:

  A character string specifying the path and name of the output PDF
  file.

## Value

ggplot object if specified. Generates a plot otherwise.

## Examples

``` r
library(SpotSweeper)
library(SpatialExperiment)
library(escheR)

# load example data
spe <- STexampleData::Visium_humanDLPFC()
#> see ?STexampleData and browseVignettes('STexampleData') for documentation
#> loading from cache

tempFilePath <- file.path(tempdir(), "examplePlot.pdf")

# change from gene id to gene names
rownames(spe) <- rowData(spe)$gene_name

# drop out-of-tissue spots
spe <- spe[, spe$in_tissue == 1]
spe <- spe[, !is.na(spe$ground_truth)]

# Identifying the mitochondrial transcripts in our SpatialExperiment.
is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]

# Calculating QC metrics for each spot using scuttle
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
#> Warning: 'scuttle::addPerCellQCMetrics' is deprecated.
#> Use 'scrapper::quickRnaQc.se' instead.
#> See help("Deprecated")
#> Warning: 'perCellQCMetrics' is deprecated.
#> Use 'scrapper::computeRnaQcMetrics' instead.
#> See help("Deprecated")
#> using unknown matrix fallback for 'dgTMatrix'
colnames(colData(spe))
#>  [1] "barcode_id"            "sample_id"             "in_tissue"            
#>  [4] "array_row"             "array_col"             "ground_truth"         
#>  [7] "reference"             "cell_count"            "sum"                  
#> [10] "detected"              "subsets_Mito_sum"      "subsets_Mito_detected"
#> [13] "subsets_Mito_percent"  "total"                

# Identifying local outliers using SpotSweeper
spe <- localOutliers(spe,
                     metric = "sum",
                     direction = "lower",
                     log = TRUE
)

plotQCpdf(spe,
          metric="sum",
          outliers="sum_outliers",
          fname=tempFilePath)
#> agg_record_2904a7bb26d 
#>                      2 
```

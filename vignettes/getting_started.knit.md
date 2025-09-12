---
title: "Getting started with 'SpotSweeper'"
author: 
  - name: Michael Totty
    affiliation: &id1 "Johns Hopkins Bloomberg School of Public Health, 
    Baltimore, MD, USA"
  - name: Boyi Guo
    affiliation: *id1
  - name: Stephanie Hicks
    affiliation: *id1
date: "2025-09-12"
output: BiocStyle::html_document
vignette: >
  %\VignetteIndexEntry{Getting Started with `SpotSweeper`}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
link-citations: true
---


<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>




## Introduction

`SpotSweeper` is an R package for spatial transcriptomics data quality control 
(QC). It provides functions for detecting and visualizing spot-level local 
outliers and artifacts using spatially-aware methods. The package is designed 
to work with [SpatialExperiment](https://github.com/drighelli/SpatialExperiment)
objects, and is compatible with data from 10X Genomics Visium and other spatial 
transcriptomics platforms.

## Installation

Currently, the only way to install `SpotSweeper` is by downloading the 
development version which can be installed from 
[GitHub](https://github.com/MicTott/SpotSweeper) using the following: 


``` r
if (!require("devtools")) install.packages("devtools")
remotes::install_github("MicTott/SpotSweeper")
```

Once accepted in [Bioconductor](http://bioconductor.org/), `SpotSweeper` will be
installable using:


``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("SpotSweeper")
```

## Spot-level local outlier detection

### Loading example data

Here we'll walk you through the standard workflow for using 'SpotSweeper' to 
detect and visualize local outliers in spatial transcriptomics data. We'll use 
the `Visium_humanDLPFC` dataset from the `STexampleData` package, which is a 
`SpatialExperiment` object.

Because local outliers will be saved in the `colData` of the `SpatialExperiment`
object, we'll first view the `colData` and drop out-of-tissue spots before 
calculating quality control (QC) metrics and running `SpotSweeper`. 



``` r
library(SpotSweeper)

# load  Maynard et al DLPFC daatset
spe <- STexampleData::Visium_humanDLPFC()
```

```
## see ?STexampleData and browseVignettes('STexampleData') for documentation
```

```
## loading from cache
```

``` r
# show column data before SpotSweeper
colnames(colData(spe))
```

```
## [1] "barcode_id"   "sample_id"    "in_tissue"    "array_row"    "array_col"   
## [6] "ground_truth" "reference"    "cell_count"
```

``` r
# drop out-of-tissue spots
spe <- spe[, spe$in_tissue == 1]
```

### Calculating QC metrics using `scuttle`

We'll use the `scuttle` package to calculate QC metrics. To do this, we'll need 
to first change the `rownames` from gene id to gene names. We'll then get the 
mitochondrial transcripts and calculate QC metrics for each spot using 
`scuttle::addPerCellQCMetrics`.


``` r
# change from gene id to gene names
rownames(spe) <- rowData(spe)$gene_name

# identifying the mitochondrial transcripts
is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]

# calculating QC metrics for each spot using scuttle
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
colnames(colData(spe))
```

```
##  [1] "barcode_id"            "sample_id"             "in_tissue"            
##  [4] "array_row"             "array_col"             "ground_truth"         
##  [7] "reference"             "cell_count"            "sum"                  
## [10] "detected"              "subsets_Mito_sum"      "subsets_Mito_detected"
## [13] "subsets_Mito_percent"  "total"
```


### Identifying local outliers using `SpotSweeper`


We can now use `SpotSweeper` to identify local outliers in the spatial 
transcriptomics data. We'll use the `localOutliers` function to detect local 
outliers based on the unique detected genes, total library size, and percent of 
the total reads that are mitochondrial. These methods assume a normal 
distribution, so we'll use the log-transformed sum of the counts and the 
log-transformed number of detected genes. For mitochondrial percent, we'll use 
the raw mitochondrial percentage. 


``` r
# library size
spe <- localOutliers(spe,
    metric = "sum",
    direction = "lower",
    log = TRUE
)

# unique genes
spe <- localOutliers(spe,
    metric = "detected",
    direction = "lower",
    log = TRUE
)

# mitochondrial percent
spe <- localOutliers(spe,
    metric = "subsets_Mito_percent",
    direction = "higher",
    log = FALSE
)
```

The `localOutlier` function automatically outputs the results to the `colData` 
with the naming convention `X_outliers`, where `X` is the name of the input 
`colData`. We can then combine all outliers into a single column called 
`local_outliers` in the `colData` of the `SpatialExperiment` object.


``` r
# combine all outliers into "local_outliers" column
spe$local_outliers <- as.logical(spe$sum_outliers) |
    as.logical(spe$detected_outliers) |
    as.logical(spe$subsets_Mito_percent_outliers)
```

### Visualizing local outliers

We can visualize the local outliers using the `plotQCmetrics` function. This 
function creates a scatter plot of the specified metric and highlights the 
local outliers in red using the `escheR` package. Here, we'll visualize local 
outliers of library size, unique genes, mitochondrial percent, and finally, all
local outliers. We'll then arrange these plots in a grid using 
`ggpubr::arrange`.


``` r
library(escheR)
```

```
## Loading required package: ggplot2
```

``` r
# all local outliers
plotQCmetrics(spe, metric = "sum_log", outliers = "local_outliers", point_size = 1.1, 
       stroke = 0.75) +
      ggtitle("All Local Outliers")
```

<img src="/Users/michael.totty/Documents/R/SpotSweeper/test_output3/getting_started_files/figure-html/local_outlier_plot-1.png" width="100%" />


## Removing technical artifacts using `SpotSweeper`

### Loading example data


``` r
# load in DLPFC sample with hangnail artifact
data(DLPFC_artifact)
spe <- DLPFC_artifact

# inspect colData before artifact detection
colnames(colData(spe))
```

```
##  [1] "sample_id"          "in_tissue"          "array_row"         
##  [4] "array_col"          "key"                "sum_umi"           
##  [7] "sum_gene"           "expr_chrM"          "expr_chrM_ratio"   
## [10] "ManualAnnotation"   "subject"            "region"            
## [13] "sex"                "age"                "diagnosis"         
## [16] "sample_id_complete" "count"              "sizeFactor"
```


### Visualizing technical artifacts

Technical artifacts can commonly be visualized by standard QC metrics, including
library size, unique genes, or mitochondrial percentage. We can first visualize
the technical artifacts using the `plotQCmetrics` function. In this sample, we can 
clearly see a hangnail artifact on the right side of the tissue section
in the mitochondrial ratio plot.


``` r
plotQCmetrics(spe,
    metric = "expr_chrM_ratio",
    outliers = NULL, point_size = 1.1
) +
    ggtitle("Mitochondrial Percent")
```

<img src="/Users/michael.totty/Documents/R/SpotSweeper/test_output3/getting_started_files/figure-html/artifact_QC_plots-1.png" width="100%" />

### Identifying artifacts using `SpotSweeper`

We can then use the `findArtifacts` function to identify artifacts in the 
spatial transcriptomics (data. This function identifies technical artifacts 
based on the first principle component of the local variance of the specified QC
metric (`mito_percent`) at numerous neighorhood sizes (`n_order=5`). Currently,
`kmeans` clustering is used to cluster the technical artifact vs high-quality 
Visium spots. Similar to `localOutliers`, the `findArtifacts` function then
outputs the results to the `colData`.


``` r
# find artifacts using SpotSweeper
spe <- findArtifacts(spe,
    mito_percent = "expr_chrM_ratio",
    mito_sum = "expr_chrM",
    n_order = 5,
    name = "artifact"
)

# check that "artifact" is now in colData
colnames(colData(spe))
```

```
##  [1] "sample_id"           "in_tissue"           "array_row"          
##  [4] "array_col"           "key"                 "sum_umi"            
##  [7] "sum_gene"            "expr_chrM"           "expr_chrM_ratio"    
## [10] "ManualAnnotation"    "subject"             "region"             
## [13] "sex"                 "age"                 "diagnosis"          
## [16] "sample_id_complete"  "count"               "sizeFactor"         
## [19] "expr_chrM_ratio_log" "coords"              "k6"                 
## [22] "k18"                 "k36"                 "k60"                
## [25] "k90"                 "artifact"
```


### Visualizing artifacts

We can visualize the artifacts using the `escheR` package. Here, we'll visualize
the artifacts using the `plotQCmetrics` function and arrange these plots using 
`ggpubr::arrange`.


``` r
plotQCmetrics(spe,
    metric = "expr_chrM_ratio",
    outliers = "artifact", point_size = 1.1
) +
    ggtitle("Hangnail artifact")
```

<img src="/Users/michael.totty/Documents/R/SpotSweeper/test_output3/getting_started_files/figure-html/artifact_visualization-1.png" width="100%" />
# Session information


``` r
utils::sessionInfo()
```

```
## R version 4.4.3 Patched (2025-02-28 r87922)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sequoia 15.6.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/New_York
## tzcode source: internal
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] escheR_1.6.0                ggplot2_3.5.2              
##  [3] STexampleData_1.14.1        SpatialExperiment_1.16.0   
##  [5] SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
##  [7] Biobase_2.66.0              GenomicRanges_1.58.0       
##  [9] GenomeInfoDb_1.42.3         IRanges_2.40.1             
## [11] S4Vectors_0.44.0            MatrixGenerics_1.18.1      
## [13] matrixStats_1.5.0           ExperimentHub_2.14.0       
## [15] AnnotationHub_3.14.0        BiocFileCache_2.14.0       
## [17] dbplyr_2.5.0                BiocGenerics_0.52.0        
## [19] SpotSweeper_1.3.3           BiocStyle_2.34.0           
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.2.1        viridisLite_0.4.2       dplyr_1.1.4            
##  [4] farver_2.1.2            blob_1.2.4              Biostrings_2.74.1      
##  [7] filelock_1.0.3          fastmap_1.2.0           digest_0.6.37          
## [10] mime_0.12               lifecycle_1.0.4         KEGGREST_1.46.0        
## [13] terra_1.8-29            RSQLite_2.3.9           magrittr_2.0.3         
## [16] compiler_4.4.3          rlang_1.1.6             sass_0.4.9             
## [19] tools_4.4.3             yaml_2.3.10             knitr_1.49             
## [22] labeling_0.4.3          S4Arrays_1.6.0          bit_4.6.0              
## [25] curl_6.2.1              DelayedArray_0.32.0     RColorBrewer_1.1-3     
## [28] abind_1.4-8             BiocParallel_1.40.0     withr_3.0.2            
## [31] purrr_1.0.4             grid_4.4.3              beachmat_2.22.0        
## [34] scales_1.4.0            MASS_7.3-65             tinytex_0.56           
## [37] cli_3.6.5               rmarkdown_2.29          crayon_1.5.3           
## [40] generics_0.1.3          httr_1.4.7              rjson_0.2.23           
## [43] scuttle_1.16.0          DBI_1.2.3               cachem_1.1.0           
## [46] zlibbioc_1.52.0         parallel_4.4.3          AnnotationDbi_1.68.0   
## [49] BiocManager_1.30.25     XVector_0.46.0          vctrs_0.6.5            
## [52] Matrix_1.7-3            jsonlite_1.9.1          bookdown_0.42          
## [55] BiocNeighbors_2.0.1     bit64_4.6.0-1           magick_2.8.5           
## [58] jquerylib_0.1.4         glue_1.8.0              codetools_0.2-20       
## [61] gtable_0.3.6            BiocVersion_3.20.0      UCSC.utils_1.2.0       
## [64] tibble_3.2.1            pillar_1.10.1           rappdirs_0.3.3         
## [67] htmltools_0.5.8.1       GenomeInfoDbData_1.2.13 R6_2.6.1               
## [70] evaluate_1.0.3          lattice_0.22-6          png_0.1-8              
## [73] memoise_2.0.1           bslib_0.9.0             Rcpp_1.1.0             
## [76] SparseArray_1.6.2       xfun_0.51               spatialEco_2.0-2       
## [79] pkgconfig_2.0.3
```

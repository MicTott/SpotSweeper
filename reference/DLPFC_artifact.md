# human DLPFC dataset with a technical artifact (hangnail).

The `DLPFC_artifact` dataset is a `SpatialExperiment` object containing
a single-sample subset of the human dorsolateral prefrontal cortex
(DLPFC) dataset from Hukki-Myers et al. 2023. This particular sample
('Br2743_ant') is included to demonstrate the identification and removal
of technical artifacts within spatial transcriptomics data. The dataset
serves as an example for artifact detection using the 'SpotSweeper'
workflow.

## Usage

``` r
data(DLPFC_artifact)
```

## Format

An SpatialExperiment object.

## Source

[spatialLIBD](https://github.com/LieberInstitute/spatialDLPFC)

## References

Hukki-Myers et al. (2023) bioRxiv
([bioRxiv](https://www.biorxiv.org/content/10.1101/2023.02.15.528722v1))

## Examples

``` r
data(DLPFC_artifact)
```

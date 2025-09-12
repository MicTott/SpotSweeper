# run examples from localOutliers() function documentation
library(SpotSweeper)
library(SpatialExperiment)

# load example data
spe <- STexampleData::Visium_humanDLPFC()

# change from gene id to gene names
rownames(spe) <- rowData(spe)$gene_name

# drop out-of-tissue spots
spe <- spe[, spe$in_tissue == 1]
spe <- spe[, !is.na(spe$ground_truth)]

# Identifying the mitochondrial transcripts in our SpatialExperiment.
is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]

# Calculating QC metrics for each spot using scuttle
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
colnames(colData(spe))

# Identifying local outliers using SpotSweeper
spe <- localOutliers(spe,
                     metric = "sum",
                     direction = "lower",
                     log = TRUE
)

# === Tests ===
test_that("example objects have correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("outlier detection functionality works", {
  # Check that outlier columns were created
  expect_true("sum_outliers" %in% colnames(colData(spe)))
  expect_true("sum_z" %in% colnames(colData(spe)))
  expect_true("sum_log" %in% colnames(colData(spe)))
  
  # Check that some outliers were found (but don't specify exact number)
  outlier_count <- sum(as.logical(spe$sum_outliers))
  expect_gt(outlier_count, 0)
  expect_lt(outlier_count, ncol(spe)) # Less than total spots
  
  # Check that z-scores are numeric and finite
  expect_true(is.numeric(spe$sum_z))
  expect_true(all(is.finite(spe$sum_z)))
  
  # Check that log column was created and is numeric
  expect_true(is.numeric(spe$sum_log))
  expect_true(all(is.finite(spe$sum_log)))
})

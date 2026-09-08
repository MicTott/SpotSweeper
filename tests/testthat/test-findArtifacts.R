# run examples from findArtifacts() function documentation
data(DLPFC_artifact)
spe <- DLPFC_artifact

# find artifacts using
set.seed(123)
spe <- findArtifacts(spe,
                     mito_percent = "expr_chrM_ratio",
                     mito_sum = "expr_chrM",
                     n_order = 2,
                     name = "artifact"
)

# === Tests ===
test_that("example objects have correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("example object contains artifacts colDta", {
  expect_equal(dim(spe), c(5000, 3529))
})

test_that("findArtifacts explains missing mitochondrial signal", {
  ffpe <- DLPFC_artifact[, seq_len(100)]
  ffpe$expr_chrM_ratio <- 0
  ffpe$expr_chrM <- 0

  expect_error(
    findArtifacts(
      ffpe,
      mito_percent = "expr_chrM_ratio",
      mito_sum = "expr_chrM",
      n_order = 1
    ),
    "No mitochondrial signal"
  )
})

# test_that("examples gives correct number of artifact spots", {
#   expect_equal(sum(spe$artifact), 517)
# })

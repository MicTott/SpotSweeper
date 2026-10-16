.parallel_param <- function(workers, os_type = .Platform$OS.type) {
  if (workers == 1L) {
    BiocParallel::SerialParam()
  } else if (os_type == "windows") {
    BiocParallel::SnowParam(workers = workers, type = "SOCK")
  } else {
    BiocParallel::MulticoreParam(workers = workers)
  }
}

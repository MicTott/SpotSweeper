#' Spatial Experiment Data Fetcher
#'
#' This function fetches spatial experiment data. It currently supports fetching data for
#' the DLPFC (Dorsolateral Prefrontal Cortex) brain region. The data is retrieved using
#' the `spatialLIBD` package.
#'
#' @param data A character string specifying the type of data to fetch.
#'             Default is "DLPFC", which fetches the spatialExperiment object
#'             for Dorsolateral Prefrontal Cortex.
#'
#' @return A `SpatialExperiment` object containing the fetched data.
#'
#' @examples
#' spe <- example_SPE()
#' spe <- example_SPE(data = "DLPFC")
#'
#' @export
#' @importFrom spatialLIBD fetch_data
example_SPE <- function(data = "DLPFC") {
  if (data == "DLPFC") {
    # Fetch the data
    spe <- spatialLIBD::fetch_data(type = "spe")

    # Return the spatialExperiment object
    return(spe)
  }
}

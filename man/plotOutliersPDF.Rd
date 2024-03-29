% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotOutliersPDF.R
\name{plotOutliersPDF}
\alias{plotOutliersPDF}
\title{Plot Outlier Metrics to PDF}
\usage{
plotOutliersPDF(
  spe,
  sample_id = "sample_id",
  metric = "detected",
  outliers = "local_outliers",
  colors = c("white", "black"),
  stroke = 1,
  width = 5,
  height = 5,
  fname
)
}
\arguments{
\item{spe}{A SpatialExperiment object containing the data to be plotted.}

\item{sample_id}{A character string specifying the column name in \code{colData(spe)}
that contains unique sample identifiers. Default is "sample_id".}

\item{metric}{A character string specifying the metric to be visualized
in the plot. This metric should be a column name in \code{colData(spe)}.}

\item{outliers}{A character string specifying the column name in \code{colData(spe)}
that indicates whether a data point is considered an outlier. Default is
"local_outliers".}

\item{colors}{A character vector specifying the colors to be used for the
gradient scale. If length is 2, the gradient will be a single color gradient.}

\item{stroke}{A numeric value specifying the border thickness for outlier
points. Default is 1.}

\item{width}{A numeric value indicating the width of the plot. Default
is 5.}

\item{height}{A numeric value indicating the height of the plot. Default
is 5.}

\item{fname}{A character string specifying the path and name of the output PDF file.}
}
\value{
Invisible NULL. The function's primary effect is the creation of a PDF
file at the specified location.
}
\description{
This function generates a PDF file containing plots for each sample in the
SpatialExperiment object, highlighting outliers based on specified metrics.
Each plot visualizes outlier metrics for a single sample, allowing for
easy comparison and analysis across samples.
}
\examples{
plotOutliersPDF(spe, fname="outliers_visualization.pdf")

}

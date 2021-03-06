% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metric.R
\name{dists_to_mds}
\alias{dists_to_mds}
\title{Multidimensional scaling for pairwise sample distances}
\usage{
dists_to_mds(dists)
}
\arguments{
\item{data}{distance table produced by '*_dist' functions}
}
\value{
a data table holding sample ids and their 2D coordinates
}
\description{
Performs a multidimensional scaling for pairwise sample distance
tables produced by '*_dist' functions
}

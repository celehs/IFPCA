% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ifpca-kern.R
\name{FPC.Kern.S}
\alias{FPC.Kern.S}
\title{Mean Density and Covariance Functions for FPCA}
\usage{
FPC.Kern.S(x, t, N, h1 = NULL, h2 = NULL, bw = "ucv")
}
\arguments{
\item{x}{time grid between 0 and 1}

\item{t}{observed event times of all the individuals, can have duplicates}

\item{N}{vector, contains num of observed event of each patient}

\item{h1}{bandwidth for mean density function}

\item{h2}{bandwidth for covariance function}

\item{bw}{bandwidth selection method}
}
\description{
Compute mean density function and covariance function for 
functional principal components analysis (FPCA)
}

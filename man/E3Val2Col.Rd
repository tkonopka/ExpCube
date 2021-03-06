% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ExpCube-plothelp.R
\name{E3Val2Col}
\alias{E3Val2Col}
\title{Convert a number into a color using  transparency}
\usage{
E3Val2Col(x, col = c("#0000ff", "#ff0000"), maxval = 1)
}
\arguments{
\item{x}{- numeric matrix.}

\item{col}{- vector of two colors in #XXXXXX format. First element determines
color associated with negative values. Second element determines color associated
with positive values}

\item{maxval}{- numeric. Value for which color reaches saturation}
}
\description{
See also E3ValMat2ColMat. Uses val2hex from package Rpipelines
}


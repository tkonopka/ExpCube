% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ExpCube-plot.R
\name{E3PlotBoxes}
\alias{E3PlotBoxes}
\title{Draw a boxplot (using Rcssplot)}
\usage{
E3PlotBoxes(vlist, quantiles = c(0, 0.25, 0.5, 0.75, 1), xlab = "",
  ylab = "", main = "", Rcss = "default", Rcssclass = "plotbox")
}
\arguments{
\item{vlist}{- names list of numeric vectors. Analogous to input to boxplot()}

\item{quantiles}{- numeric vector of five elements. The five elements will
be treated as quantiles for lower-whisker, lowe-box, middle, upper-box, upper-whisker}

\item{xlab}{- character string. Text to display below x axis}

\item{ylab}{- character string. Text to display below y axis.}

\item{main}{- character string. Text to display as title, above heatmap.}

\item{Rcss}{- Rcss object. Used to style the heatmap with Rcssplot.}

\item{Rcssclass}{- character vector. Classes to tune Rcssplot formatting.}
}
\description{
Draw a boxplot (using Rcssplot)
}


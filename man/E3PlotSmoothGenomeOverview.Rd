% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ExpCube-plot.R
\name{E3PlotSmoothGenomeOverview}
\alias{E3PlotSmoothGenomeOverview}
\title{Draws a heatmap representation of regional genomic patterns (using Rcssplot)}
\usage{
E3PlotSmoothGenomeOverview(fcdata, gene.positions = NULL,
  categories = list(All = colnames(fcdata)), label.chr = c(),
  chrlens = NULL, heatcols = c("#00ffff", "#ff00ff"), seg.rescale = 2,
  seg.resolution = 5e+06, seg.hlab = 0.3, seg.vlab = 1e+06,
  legend = NULL, axes.labels = NULL, Rcss = "default",
  Rcssclass = "diamond")
}
\arguments{
\item{fcdata}{- numeric matrix with gene/feature names in rownames and samples in colnames.
This data will be log2 transformed before the plotting takes place}

\item{gene.positions}{- data frame linking gene names to gene positions}

\item{categories}{- named list of vectors, splitting samples into groups
(e.g. Plates or batches). This object also determines the order of the samples on
the y axis}

\item{label.chr}{- names of chromosomes to label on horizontal axis (}

\item{chrlens}{- named vector with lengths of chromosomes (chromosomes will appear
in this order)}

\item{heatcols}{- numeric vector of lengths two. Determines the color scale of the plot.}

\item{seg.rescale}{- numeric. Used to modulate color intensity}

\item{seg.resolution}{- numeric. Width of bins to use for smoothing}

\item{seg.hlab}{- numeric. Determines spacing between labels near x axis.}

\item{seg.vlab}{- numeric. Determines spacing between labels near y axis.}

\item{legend}{- numeric. Values to use on the legend box.}

\item{axes.labels}{- character vector of length two. Text to write near
y axis to the left and to the right of the heatmap.}

\item{Rcss}{- Rcss object. Used to style the heatmap with Rcssplot.}

\item{Rcssclass}{- character vector. Classes to tune Rcssplot formatting.}
}
\description{
The function takes numeric data associated with genes, then performs a smoothing using
bins and plots the resulting smoothed signal.
}


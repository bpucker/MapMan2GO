% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funks.R
\name{plotDistAsHistAndBox}
\alias{plotDistAsHistAndBox}
\title{Generates a two row plot with the first one being a Histogram and the second
row a horizontal Boxplot.}
\usage{
plotDistAsHistAndBox(x, main,
  summary.as.title = getOption("MapMan2GO.plot.dist.summary.as.title", FALSE))
}
\arguments{
\item{x}{The values passed into \code{hist} and \code{boxplot}}

\item{summary.as.title}{boolean indicating whether to add the output of
\code{base::summary(x)} to the title of the plot. Default is
\code{getOption('MapMan2GO.plot.dist.summary.as.title', FALSE)}.}

\item{The}{main title of the resulting plot}
}
\value{
TRUE if and only if no error has occurred
}
\description{
Generates a two row plot with the first one being a Histogram and the second
row a horizontal Boxplot.
}

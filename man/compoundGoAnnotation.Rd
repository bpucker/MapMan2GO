% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funks.R
\name{compoundGoAnnotation}
\alias{compoundGoAnnotation}
\title{Extracts the unique and sorted Gene Ontology Term Annotations for a given
gene.}
\usage{
compoundGoAnnotation(gene.id, goa.tbl = getOption("MapMan2GO.goa.tbl",
  ukb.goa.hits), gene.col = getOption("MapMan2GO.goa.tbl.gene.col", 3),
  go.col = getOption("MapMan2GO.goa.tbl.go.col", 2),
  extend.goas.with.ancestors = getOption("MapMan2GO.extend.goa.with.ancestors",
  TRUE))
}
\arguments{
\item{gene.id}{The unique identifier of the gene to lookup GOA for}

\item{goa.tbl}{An instance of \code{data.frame} holding GOAs. Default is
\code{getOption('MapMan2GO.goa.tbl', ukb.goa.hits)}}

\item{gene.col}{The column of \code{goa.tbl} number or name in which to
lookup the genes. Default is \code{getOption('MapMan2GO.goa.tbl.gene.col',
3)}.}

\item{go.col}{The column of \code{goa.tbl} number or name in which to lookup
the GO terms. Default is \code{getOption('MapMan2GO.goa.tbl.go.col',2)}}

\item{extend.goas.with.ancestors}{boolean indicating whether to extend each
proteins' GOA with the ancestors of the respective GO Terms. Default is
\code{getOption('MapMan2GO.extend.goa.with.ancestors', TRUE)}.}
}
\value{
A character holding the GO terms for \code{gene.id}
}
\description{
Extracts the unique and sorted Gene Ontology Term Annotations for a given
gene.
}

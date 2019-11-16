\name{tna.graph}
\alias{tna.graph}

\title{
Compute a graph from TNA objects (deprecated function).
}

\description{
This function is provided for compatibility with older versions of \sQuote{RTN} only, and will be defunct at the next release (see details).
}

\usage{
tna.graph(object, tnet = "dpi", gtype="rmap", minRegulonSize=15, tfs=NULL, 
amapFilter="quantile", amapCutoff=NULL, ...)
}

\arguments{
\item{object}{
an object of class 'TNA' \code{\link[RTN:TNA-class]{TNA-class}}.
}
\item{tnet}{
a single character value specifying which network information should be used to compute the graph. 
Options: "ref" and "dpi".
}
\item{gtype}{
a single character value specifying the graph type.
Options: "rmap" and "amap". The "rmap" option returns regulatory maps represented by TFs and targets 
(regulons) and "amap" computes association maps among regulons (estimates the overlap using the Jaccard Coefficient).
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of 
elements in a regulon. Regulons with fewer than this number are 
removed from the graph.
}
\item{tfs}{
a vector with transcription factor identifiers.
}
\item{amapFilter}{
a single character value specifying which method should be used to filter association maps 
(only when gtype="amap"). Options: "phyper","quantile" and "custom".
}
\item{amapCutoff}{
a single numeric value (>=0 and <=1) specifying the cutoff for the association map filter. When 
amapFilter="phyper", amapCutoff corresponds to a pvalue cutoff; when amapFilter="quantile", 
amapCutoff corresponds to a quantile threshold; and when amapFilter="custom", amapCutoff is 
a JC threshold.
}
\item{...}{ additional arguments passed to tna.graph function. } 
}

\details{

  The following function is deprecated and will be made defunct; use
  the replacement indicated below:
  \itemize{

    \item{tna.graph: \code{\link{tni.graph}}}

  }
}

\value{
a graph object.
}

\author{
Mauro Castro
}

\keyword{methods}

\name{tni.graph}
\alias{tni.graph}

\title{
Compute a graph from TNI objects.
}

\description{
Extract results from a TNI object and compute a graph.
}

\usage{
tni.graph(object, tnet = "dpi", gtype="rmap", minRegulonSize=15, 
regulatoryElements=NULL, amapFilter="quantile", 
amapCutoff=NULL, ntop=NULL, ...)
}

\arguments{
\item{object}{
an object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{tnet}{
a single character value specifying which network information should be used to compute the graph. 
Options: "ref" and "dpi".
}
\item{gtype}{
a single character value specifying the graph type.
Options: "rmap", "amap", "mmap" and "mmapDetailed". The "rmap" option returns regulatory maps represented by 
TFs and targets (regulons); "amap" computes association maps among regulons (estimates the overlap using the Jaccard Coefficient); "mmap" and "mmapDetailed" returns modulated maps derived from the \code{\link[RTN:tni.conditional]{tni.conditional}} function.
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of 
elements in a regulon. Regulons with fewer than this number are 
removed from the graph.
}
\item{regulatoryElements}{
an optional vector with transcription factor identifiers.
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
\item{ntop}{
when gtype="mmapDetailed", ntop is an optional single integer value (>=1) specifying the number of TF's targets that should be used to compute the modulated map. The n targets is derived from the top ranked TF-target interations, as defined in the mutual information analysis used to construct the regulon set.
}
\item{...}{ additional arguments passed to tni.graph function. } 
}

\value{
an igraph object.
}

\author{
Mauro Castro
}

\examples{

data(tniData)

\dontrun{

rtni <- tni.constructor(expData=tniData$expData, 
        regulatoryElements=c("PTTG1","E2F2","FOXM1","E2F3","RUNX2"), 
        rowAnnotation=tniData$rowAnnotation)
rtni <- tni.permutation(rtni)
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni, eps=0)

# get a regulatory map
g <- tni.graph(rtni, tnet="dpi", gtype="rmap", 
               regulatoryElements=c("PTTG1","E2F2","FOXM1"))

# option: plot the igraph object using RedeR
library(RedeR)
rdp <- RedPort()
calld(rdp)
addGraph(rdp,g)
addLegend.shape(rdp, g)
addLegend.color(rdp, g, type="edge")
relax(rdp, p1=20, p3=50, p4=50, p5=10, p8=50, ps=TRUE)
#...it should take some time to relax!

# get an association map
resetd(rdp)
g <- tni.graph(rtni, tnet="ref", gtype="amap")
addGraph(rdp,g)
addLegend.size(rdp, g)
addLegend.size(rdp, g, type="edge")

}
}

\keyword{methods}

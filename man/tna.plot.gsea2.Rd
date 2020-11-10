\name{tna.plot.gsea2}
\alias{tna.plot.gsea2}
\title{
Plot enrichment analyses from TNA objects.
}

\description{
This function takes a TNA object and plots the two-tailed GSEA results for 
individual regulons.
}

\usage{
tna.plot.gsea2(object, labPheno="", file="tna_gsea2", filepath=".", 
regulon.order="size", ntop=NULL, tfs=NULL, ylimPanels=c(-3.0,3.0,-0.5,0.5), 
heightPanels=c(2.0,0.8,5.0), width=2.8, height=3.0, 
ylabPanels=c("Phenotype","Regulon","Enrichment score"), 
xlab="Position in the ranked list of genes", alpha=1.0, 
sparsity=10, autoformat=TRUE, plotpdf=TRUE, ...)
}

\arguments{
\item{object}{
an object of class 'TNA' \code{\link[RTN:TNA-class]{TNA-class}}.
}
\item{file}{
a character string naming a file.
}
\item{filepath}{
a single character value specifying where to store GSEA2 figures.
}
\item{regulon.order}{
a single character value specifying whether regulons should be ordered by 
'size', 'score', 'pvalue', 'adj.pvalue' and 'name' (or 'none' to keep the 
input ordering).
}
\item{ntop}{
a single integer value specifying how many regulons of top significance will be 
plotted.
}
\item{tfs}{
an optional vector with transcription factor identifiers (this option overrides 
the 'ntop' argument).
}
\item{ylimPanels}{
a numeric vector of length=4 specifying y coordinates ranges of the 1st and 3th 
plots (i.e. ylim for 'Phenotypes' and 'Running enrichment score').
}
\item{heightPanels}{
a numeric vector of length=3 specifying the relative height of each panel in 
the plot.
}
\item{width}{
a single numeric value specifying the width of the graphics region in inches.
}
\item{height}{
a single numeric value specifying the height of the graphics region in inches.
}
\item{ylabPanels}{
a character vector of length=3 specifying the the title for the y axes.
}
\item{xlab}{
a single character specifying the the title for the x axis.
}
\item{labPheno}{
a single character specifying a label for the phenotype 
(will also be used as the name of output file).
}
\item{alpha}{
a single numeric value in [0,1] specifying the transparency of the hits in 
the ranked list.
}
\item{sparsity}{
a single integer value (>1) specifying the density of the dots representing 
the running score.
}
\item{autoformat}{
a single logical value specifying to set the graph format using predefined 
themes. This option overrides the "ylimPanels" argument.
}
\item{plotpdf}{
a single logical value specifying to whether to plot a PDF file or 
directly to Viewer.
}
\item{\dots}{
other arguments used by the function pdf.
}
}

\value{
A plot showing results from the 'tna.gsea2' method.
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:tna.gsea2]{tna.gsea2}}
}

\examples{

data(tniData)
data(tnaData)

\dontrun{

rtni <- tni.constructor(expData=tniData$expData, 
        regulatoryElements=c("PTTG1","E2F2","FOXM1","E2F3","RUNX2"), 
        rowAnnotation=tniData$rowAnnotation)
rtni <- tni.permutation(rtni)
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni)
rtna <- tni2tna.preprocess(rtni, phenotype=tnaData$phenotype, 
        hits=tnaData$hits, phenoIDs=tnaData$phenoIDs)

# run GSEA2 analysis pipeline
rtna <- tna.gsea2(rtna)

# plot available GSEA2 results
tna.plot.gsea2(rtna, labPheno="test")
}
}

\keyword{GSEA2}

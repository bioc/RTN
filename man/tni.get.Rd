\name{tni.get}
\alias{tni.get}

\title{
Get information from individual slots in a TNI object.
}

\description{
Get available results from individual slots in a TNI object.
}

\usage{
tni.get(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE, 
idkey=NULL)
}

\arguments{
\item{object}{
an object of class \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{what}{
a single character value specifying which information should be retrieved from 
the slots. Options: 'summary', 'status', 'para', 'gexp','regulatoryElements', 
'targetElements', 'modulators', 'tnet', 'refnet', 'regulons', 'refregulons', 
'regulons.and.mode', 'refregulons.and.mode', 'rowAnnotation', 'colAnnotation', 
'cdt.list', 'cdt.table', 'regulonSize','regulonActivity'.
}
\item{order}{
a single logical value specifying whether or not the output data should be 
ordered by significance. Valid only for 'cdt' option.
}
\item{ntop}{
a single integer value specifying to select how many results of top 
significance from 'cdt' option.
}
\item{reportNames}{
a single logical value specifying to report regulators with 'names' 
(when reportNames=TRUE) or not (when reportNames=FALSE). This option takes 
effect on 'cdt' option if regulators are named with alternative identifiers.
}
\item{idkey}{
an optional single character value specifying an ID name from the available 
'TNI' annotation to be used as alias for data query outputs 
(obs. it has no effect on consolidated tables).
}
}

\details{

Options for the 'what' argument:

\describe{

\item{summary}{
A list summarizing parameters and results available in the TNI object 
(see \code{\link[RTN:tni.regulon.summary]{tni.regulon.summary}} for a 
summary of the network and regulons).
}
\item{status}{
A vector indicating the status of each available method in the pipeline.
}
\item{para}{
A list with the parameters used by each available method in the pipeline.
}
\item{gexp}{
A gene expression matrix.
}
\item{regulatoryElements}{
A vector of regulatory elements (e.g. transcription factors).
}
\item{targetElements}{
A vector of target elements (e.g. TF targets).
}
\item{modulators}{
A vector of modulators (e.g. TF modulators).
}
\item{tnet}{
A data matrix with MI values, evaluated by the DPI filter. MI values are 
computed between regulators and targets, with regulators on cols and targets 
on rows. Note that signals (+/-) are assigned to the inferred associations in 
order to represent the 'mode of action', which is derived from Pearson's 
correlation between regulators and targets.
}
\item{refnet}{
A data matrix with MI values (not evaluated by the DPI filter). MI values are 
computed between regulators and targets, with regulators on cols and targets 
on rows. Note that signals (+/-) are assigned to the inferred associations in 
order to represent the 'mode of action', which is derived from Pearson's 
correlation between regulators and targets.
}
\item{regulons}{
A list with regulons extracted from the 'tnet' data matrix.
}
\item{refregulons}{
A list with regulons extracted from the 'refnet' data matrix.
}
\item{regulons.and.mode}{
A list with regulons extracted from the 'tnet' data matrix, including the 
assiged 'mode of action'.
}
\item{refregulons.and.mode}{
A list with regulons extracted from the 'refnet' data matrix, including the 
assiged 'mode of action'.
}
\item{rowAnnotation}{
A data frame with probe-to-gene annotation.
}
\item{colAnnotation}{
A data frame with sample annotation.
}
\item{cdt.table}{
A data frame with results from the 
\code{\link[RTN:tni.conditional]{tni.conditional}} analysis pipeline.
}
\item{cdt.list}{
A list with results from the 
\code{\link[RTN:tni.conditional]{tni.conditional}} analysis pipeline.
}
\item{regulonSize}{
A data frame with the number of targets annotated in each regulon.
}
\item{regulonActivity}{
A list with results from the 
\code{\link[RTN:tni.gsea2]{tni.gsea2}} analysis pipeline.
}

}

}

\value{
Get the slot content from a \code{\link[RTN:TNI-class]{TNI-class}} object.
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
rtni <- tni.dpi.filter(rtni)

# check summary
tni.get(rtni, what="summary")

# get regulons
regulons <- tni.get(rtni, what = "regulons")

# get status of the pipeline
tni.get(rtni, what="status")
}
}

\keyword{methods}

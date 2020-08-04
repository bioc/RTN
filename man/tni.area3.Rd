\name{tni.area3}
\alias{tni.area3}

\title{
Compute regulon activity by calling aREA (analytic Rank-based Enrichment Analysis) algorithm
}

\description{
Uses \code{\link{aREA}} 3-tail algorithm to compute regulon activity for 
\code{\link[RTN:TNI-class]{TNI-class}} objects.
}

\usage{
tni.area3(object, minRegulonSize=15, sizeFilterMethod="posORneg", scale=FALSE, tnet="dpi", 
regulatoryElements=NULL, samples=NULL, features=NULL, refsamp=NULL, log=FALSE, verbose=TRUE, 
doSizeFilter=NULL)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of elements in a regulon. 
Regulons smaller than this number are removed from the analysis.
}
\item{sizeFilterMethod}{
a single character value specifying the use of the 'minRegulonSize' argument, which is applyed to regulon's positive and negative targets. Options: "posANDneg", "posORneg", "posPLUSneg". For "posANDneg", the number of both positive and negative targets should be > 'minRegulonSize'; for "posORneg", the number of either positive or negative targets should be > 'minRegulonSize'; and for "posPLUSneg", the number of all targets should be > 'minRegulonSize'.
}
\item{scale}{A logical value specifying if expression values should be 
centered and scaled across samples (when verbose=TRUE) or 
not (when verbose=FALSE).
}
\item{tnet}{
can take values of 'refnet', 'dpi' or 'cdt'. It refers to the version of the regulatory
network that will be used for GSEA analysis.
}
\item{regulatoryElements}{
an optional vector with transcription factor identifiers.
}
\item{samples}{
an optional string vector containing the sample names for which will be computed 
regulon activity.
}
\item{features}{
a string vector containing features for feature selection.
}
\item{refsamp}{
an optional string vector containing the names of the reference samples for differential 
expression calculations. If not provided, then the average of all samples will be 
used as reference.
}
\item{log}{
a logical value. If TRUE, differential expression calculations will be computed in log
space.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
\item{doSizeFilter}{
a logical value. If TRUE, negative and positive targets are independently verified by the 
'minRegulonSize' argument.
}
}

\value{
a list with enrichment scores for all samples in the TNI.
}

\references{
Alvarez et al. Functional characterization of somatic mutations in cancer using 
network-based inference of protein activity. Nature Genetics, 48(8):838-847, 2016.
}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
\code{\link{aREA}}
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

#run aREA algorithm
rtni <- tni.area3(rtni)

#get results
regulonActivity <- tni.get(rtni, what = "regulonActivity")


}
}

\keyword{aREA}
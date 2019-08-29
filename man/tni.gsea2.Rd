\name{tni.gsea2}
\alias{tni.gsea2}

\title{
Compute regulon activity by calling GSEA2 (two-tailed Gene Set Enrichment Analysis) 
algorithm
}

\description{
Uses GSEA2 algorithm to compute regulon activity for 
\code{\link[RTN:TNI-class]{TNI-class}} objects.
}

\usage{
tni.gsea2(object, minRegulonSize=15, doSizeFilter=FALSE, scale=FALSE, exponent=1, 
tnet="dpi", tfs=NULL, samples=NULL, features=NULL, refsamp=NULL, log=FALSE, 
alternative=c("two.sided", "less", "greater"), targetContribution=FALSE, 
additionalData=FALSE, verbose=TRUE)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of elements in 
a regulon. Regulons smaller than this number are removed from the analysis.
}
\item{doSizeFilter}{
a logical value. If TRUE, negative and positive targets are independently 
verified by the 'minRegulonSize' argument.
}
\item{scale}{A logical value specifying if expression values should be 
centered and scaled across samples (when verbose=TRUE) or 
not (when verbose=FALSE).
}
\item{exponent}{
a single integer or numeric value used in weighting phenotypes in GSEA 
(see 'gseaScores' function at HTSanalyzeR).
}
\item{tnet}{
can take values of 'refnet', 'dpi' or 'cdt'. It refers to the version of the 
regulatory network that will be used for GSEA analysis.
}
\item{tfs}{
an optional vector with transcription factor identifiers.
}
\item{samples}{
an optional string vector containing the sample names for which 
will be computed the GSEA2.
}
\item{features}{
a string vector containing features for feature selection.
}
\item{refsamp}{
an optional string vector containing the names of the reference samples for 
differential expression calculations. If not provided, then the average of 
all samples will be used as reference.
}
\item{log}{
a logical value. If TRUE, differential expression calculations will 
be computed in log space.
}
\item{alternative}{
a character string specifying the alternative hypothesis, must be one 
of "two.sided" (default), "greater" or "less".
}
\item{targetContribution}{
This argument is used for internal calls. A single logical value specifying 
to return the contribution of each target in enrichment scores 
(when verbose=TRUE) or not (when verbose=FALSE).
}
\item{additionalData}{
This argument is used for internal calls. A single logical value specifying 
to return the additional data objects (when verbose=TRUE) or not 
(when verbose=FALSE).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
}

\value{
a list with enrichment scores for all samples in the TNI.
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
\code{\link[RTN:tna.gsea2]{tna.gsea2}}
\code{\link[RTN:tna.plot.gsea2]{tna.plot.gsea2}}
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

#run GSEA2 analysis pipeline
EScores <- tni.gsea2(rtni)

#parallel version with SNOW package!
library(snow)
options(cluster=snow::makeCluster(3, "SOCK"))
EScores <- tni.gsea2(rtni)
stopCluster(getOption("cluster"))

}
}

\keyword{GSEA}
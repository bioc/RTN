\name{tni.constructor}
\alias{tni.constructor}

\title{A constructor for objects of class TNI.}

\description{
This function is the main entry point of the TNI pipeline.
}

\usage{
tni.constructor(expData, regulatoryElements, rowAnnotation=NULL, 
colAnnotation=NULL, cvfilter=FALSE, verbose=TRUE)
}

\arguments{
\item{expData}{
a gene expression matrix or 'SummarizedExperiment' object.
}
\item{regulatoryElements}{
a vector of regulatory elements (e.g. transcription factors).
}
\item{rowAnnotation}{
an optional data frame with gene annotation. Column 1 must provide all ids 
listed in the gene expression matrix. Ideally, col1 = <ID>, col2 = <GENEID>, 
and col3 = <SYMBOL>. Additional annotation can be included in the data frame 
and will be passed to the resulting TNI object. Furthermore, in order to 
eventually use the TNI object in \code{\link[RTN:AVS-class]{AVS-class}} 
methods, it should also include chromosome coordinates: columns 
<CHROM>, <START> and <END>. Values in <CHROM> should be listed in 
[chr1, chr2, chr3, ..., chrX], while <START> and <END> correspond to 
chromosome positions (see \code{\link[RTN:avs.evse]{avs.evse}}).
}
\item{colAnnotation}{
an optional data frame with sample annotation.
}
\item{cvfilter}{
a single logical value specifying to remove duplicated genes in the gene 
expression matrix using the probe-to-gene annotation. In this case, 
'rowAnnotation' must be provided, with col1 = <ID> and col2 = <GENEID>. 
Genes duplicated in col2 will be collapsed; the decision is made based on the 
maximum dinamic range (i.e. keeping the gene with max coefficient of variation 
across all samples).
}
\item{verbose}{
a single logical value specifying to display detailed messages 
(when verbose=TRUE) or not (when verbose=FALSE).
}
}

\value{
A pre-processed TNI-class object.
}

\author{Mauro Castro}

\seealso{
  \code{\link[RTN:TNI-class]{TNI-class}}
}

\examples{

data(tniData)

#--- run constructor
rtni <- tni.constructor(expData=tniData$expData, 
        regulatoryElements=c("PTTG1","E2F2","FOXM1","E2F3","RUNX2"), 
        rowAnnotation=tniData$rowAnnotation)
}

\keyword{methods}


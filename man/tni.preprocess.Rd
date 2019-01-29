\name{tni.preprocess}
\alias{tni.preprocess}

\title{A preprocessing function for objects of class TNI.}

\description{
This is a generic function, provides all preprocessing methods for the 'tni.constructor' function.
}

\usage{
tni.preprocess(object, rowAnnotation=NULL, colAnnotation=NULL, cvfilter=TRUE, 
verbose=TRUE)
}

\arguments{
\item{object}{
this argument is an object of class \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{rowAnnotation}{
an optional data frame or matrix with probe-to-gene annotation (e.g. probe-to-gene information). Column 1 must provide all probe ids listed in the 'gexp' matrix. Ideally, col1 = <PROBEID>, col2 = <GENEID>, and col3 = <SYMBOL>. Additional annotation can be included in the data frame and will be passed to the resulting TNI object. Furthermore, in order to eventually use the TNI object in \code{\link[RTN:AVS-class]{AVS-class}} methods, it should also include chromosome coordinates: columns <CHROM>, <START> and <END>. Values in <CHROM> should be listed in [chr1, chr2, chr3, ..., chrX], while <START> and <END> correspond to chromosome positions (see \code{\link[RTN:avs.evse]{avs.evse}}).
}
\item{colAnnotation}{
This option has not yet been implemented!
}
\item{cvfilter}{
a single logical value specifying to remove duplicated genes in the gene expression matrix using the probe-to-gene annotation. In this case, 'gexpIDs' must be provided, with col1 = <PROBEID> and col2 = <GENEID>. The decision is made based on the maximum dinamic range (i.e. keeping the probes with max coefficient of variation across all samples).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
}

\author{Mauro Castro}

\seealso{
  \code{\link[RTN:TNI-class]{TNI-class}}
}

\examples{

## see 'tni.constructor'!
}

\keyword{methods}


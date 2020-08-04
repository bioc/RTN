\name{tni.replace.samples}
\alias{tni.replace.samples}

\title{Entry point to assess new samples with previously calculated regulons.}

\description{
This function replaces samples of an existing TNI-class objects.
}

\usage{
tni.replace.samples(object, expData, rowAnnotation=NULL, colAnnotation=NULL, 
removeRegNotAnnotated=TRUE, verbose=TRUE)
}

\arguments{
\item{object}{
an object of class \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{expData}{
a gene expression matrix or 'SummarizedExperiment' object.
}
\item{rowAnnotation}{
an optional data frame with gene annotation. Column 1 must provide all ids 
listed in the gene expression matrix.
}
\item{colAnnotation}{
an optional data frame with sample annotation.
}
\item{removeRegNotAnnotated}{
a single logical value specifying to remove 'regulatoryElements' not annotated
in 'expData' (when removeRegNotAnnotated=TRUE) or not 
(when removeRegNotAnnotated=FALSE).
}
\item{verbose}{
a single logical value specifying to display detailed messages 
(when verbose=TRUE) or not (when verbose=FALSE).
}
}

\value{
A TNI-class object.
}

\author{Mauro Castro}

\examples{

## please the package's vignette
}

\keyword{methods}


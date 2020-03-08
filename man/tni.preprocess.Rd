\name{tni.preprocess}
\alias{tni.preprocess}

\title{A preprocessing function for objects of class TNI.}

\description{
This is a generic function, provides all preprocessing methods for 
the 'tni.constructor' function.
}

\usage{
tni.preprocess(object, rowAnnotation=NULL, colAnnotation=NULL, cvfilter=FALSE, 
    verbose=TRUE)
}

\arguments{
\item{object}{
this argument is an object of class \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{rowAnnotation}{
an optional data frame with gene annotation.
}
\item{colAnnotation}{
an optional data frame with sample annotation.
}
\item{cvfilter}{
a single logical value specifying to remove duplicated genes in the 
gene expression matrix using the gene annotation.
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
  \code{\link[RTN:tni.constructor]{tni.constructor}}
}

\examples{

## see 'tni.constructor'!
}

\keyword{methods}


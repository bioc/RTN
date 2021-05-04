\name{avs.get}
\alias{avs.get}

\title{
Get information from individual slots in an AVS object.
}

\description{
Get information from individual slots in an AVS object.
}

\usage{
avs.get(object, what="summary", report=FALSE, pValueCutoff=NULL)
}

\arguments{
\item{object}{
an object of class 'AVS' \code{\link[RTN:AVS-class]{AVS-class}}.
}
\item{what}{
a single character value specifying which information should be retrieved from the slots. 
Options: 'markers', 'validatedMarkers', 'variantSet', 'randomSet', 'linkedMarkers', 
'randomMarkers', 'vse', 'evse', 'rvse', 'pevse', 'annotation.vse', 'annotation.evse', 
'annotation.rvse', 'annotation.pevse', 'summary' and 'status'.
}
\item{report}{
a single logical value indicating whether to return results from 'vse', 'evse' as 
a consolidated table (if TRUE), or as they are (if FALSE).
}
\item{pValueCutoff}{
an optional single numeric value specifying the cutoff to retrive results for
p-values considered significant.
}

}

\value{
get the slot content from an object of class 'AVS' \code{\link[RTN:AVS-class]{AVS-class}}.
}

\author{
Mauro Castro
}

\examples{

\dontrun{
#This example requires the RTNdata package! (currently available under request)
library(RTNdata.LDHapMapRel27)
data(bcarisk, package="RTNdata.LDHapMapRel27")
bcavs <- avs.preprocess.LDHapMapRel27(bcarisk, nrand=1000)
avs.get(avs)
}
}

\keyword{methods}


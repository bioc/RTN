\name{tni.regulon.summary}
\alias{tni.regulon.summary}

\title{
Summary of regulon characteristics.
}

\description{
This function takes a TNI object and optionally a list of regulatory elements and returns a summary of the network (if no regulatory elements are given) or of the chosen regulon or regulons.
}

\usage{
tni.regulon.summary(object, regulatoryElements = NULL, verbose = TRUE)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} already evaluated by the 
\code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}} method.
}
\item{regulatoryElements}{
a vector of valid regulatory elements (e.g. transcription factors).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE).
}
}

\value{
It returns a print-out of the network summary (if verbose is TRUE) and invisibly returns a data.frame of network characteristics such as regulon size and regulon balance.
}

\author{Clarice Groeneveld}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
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

# Summary of the network
tni.regulon.summary(rtni)

# Summary of a regulon
tni.regulon.summary(rtni, regulatoryElements = "PTTG1")

}
}

\keyword{summary}

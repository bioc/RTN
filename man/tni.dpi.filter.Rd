\name{tni.dpi.filter}
\alias{tni.dpi.filter}

\title{
Data Processing Inequality (DPI) filter.
}

\description{
This function takes a TNI object and returns the transcriptional network filtered by the data processing inequality algorithm.
}

\usage{
tni.dpi.filter(object, eps=0, verbose=TRUE)
}

\arguments{

\item{object}{
a processed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} evaluated by the methods 
\code{\link[RTN:tni.permutation]{tni.permutation}} and \code{\link[RTN:tni.bootstrap]{tni.bootstrap}}.
}
\item{eps}{
a single numeric value specifying the threshold under which ARACNe algorithm should apply the dpi filter (eps>=0). If not available (i.e. eps=NA), then the threshold is estimated from the empirical null distribution computed in the permutation and bootstrap steps. For additional details see \code{\link[minet:aracne]{aracne}}.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE)
}

}

\value{
a mutual information matrix in the slot "results" containing a dpi-filtered transcriptional network, 
see 'tn.dpi' option in \code{\link[RTN:tni.get]{tni.get}}.
}

\author{Mauro Castro}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
}

\examples{

data(dt4rtn)

# select 5 regulatoryElements for a quick demonstration!
tfs4test <- dt4rtn$tfs[c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")]

\dontrun{

# preprocessing
rtni <- tni.constructor(expData=dt4rtn$gexp, regulatoryElements=tfs4test, 
        rowAnnotation=dt4rtn$gexpIDs)

# permutation analysis (infers the reference/relevance network)
rtni<-tni.permutation(rtni, nPermutations=100)

# dpi filter (infers the transcriptional network)
rtni<-tni.dpi.filter(rtni)

}
}

\keyword{methods}


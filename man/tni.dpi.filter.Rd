\name{tni.dpi.filter}
\alias{tni.dpi.filter}

\title{
Data Processing Inequality (DPI) filter.
}

\description{
This function takes a TNI object and returns the transcriptional network filtered by the data processing inequality algorithm.
}

\usage{
tni.dpi.filter(object, eps = 0, sizeThreshold = TRUE, minRegulonSize = 15, verbose = TRUE)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} evaluated by the methods 
\code{\link[RTN:tni.permutation]{tni.permutation}} and \code{\link[RTN:tni.bootstrap]{tni.bootstrap}}.
}
\item{eps}{
a single numeric value (>= 0) specifying the threshold under which ARACNe algorithm should apply the dpi filter. If not available (i.e. 'eps = NA'), then the threshold is estimated from the empirical null distribution computed in the permutation and bootstrap steps. For additional details see \code{\link[minet:aracne]{aracne}}.
}
\item{sizeThreshold}{
a logical value specifying if the 'minRegulonSize' argument should be used
(when 'sizeThreshold = TRUE') or not (when 'sizeThreshold = FALSE'). It will 
have no effect when 'eps = NA'.
}
\item{minRegulonSize}{
a single integer or numeric value. This argument prevents the DPI algorithm from removing additional targets from large unbalanced regulons, when the subset of either positive or negative targets is below the 'minRegulonSize'.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
'verbose = TRUE') or not (when 'verbose = FALSE').
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

data(tniData)

\dontrun{

# preprocessing
rtni <- tni.constructor(expData=tniData$expData, 
        regulatoryElements=c("PTTG1","E2F2","FOXM1","E2F3","RUNX2"), 
        rowAnnotation=tniData$rowAnnotation)
        
# permutation analysis (infers the reference/relevance network)
rtni <- tni.permutation(rtni)

# dpi filter (infers the transcriptional network)
rtni <- tni.dpi.filter(rtni)

}
}

\keyword{methods}


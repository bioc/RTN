\name{tni.prune}
\alias{tni.prune}

\title{
Prune regulons to remove redundant targets for regulon activity analysis
}

\description{
Uses network pruning methods to compute a 'core' regulon that retains good correlation with original regulon activity.
}

\usage{
tni.prune(object, regulatoryElements = NULL, minRegCor = 0.95, 
tarPriorityMethod = "EC", minPrunedSize = 30, verbose = TRUE, ...)
}

\arguments{

\item{object}{
a preprocessed object of \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{regulatoryElements}{
an optional vector with regulatoryElements identifiers. If NULL, all regulatoryElements
are used.
}
\item{minRegCor}{
an numeric value between 0 and 1. The minimum correlation between the original activity
values for a regulon and the activity after pruning.
}
\item{tarPriorityMethod}{
method for prioritizing targets for the target backwards elimination. One of "EC" (expression correlation), "MI" (mutual information) or "TC" (target contribution).
}
\item{minPrunedSize}{
a single integer or numeric value specifying the minimum number of elements in a regulon 
after pruning. 
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
\item{...}{
arguments passed to \code{\link[RTN:tni.gsea2]{tni.gsea2}}
}
}

\value{
a TNI-class object, with the pruned regulons.
}

\author{
Clarice Groeneveld
}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
\code{\link[RTN:tni.gsea2]{tni.gsea2}}
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

# prune the PTTG1 regulon
rtni_pruned <- tni.prune(rtni, "PTTG1", tarPriorityMethod = "TC")

#parallel version with SNOW package!
#library(snow)
#options(cluster=makeCluster(3, "SOCK"))
#rtni_pruned <- tni.prune(rtni, c("PTTG1", "E2F2"))
#stopCluster(getOption("cluster"))

}
}

\keyword{Prune}

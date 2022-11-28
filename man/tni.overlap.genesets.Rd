\name{tni.overlap.genesets}
\alias{tni.overlap.genesets}

\title{
Associate regulons with external gene set collections.
}

\description{
This function tests the overlap between gene sets and regulons.
}

\usage{
tni.overlap.genesets(object, geneSetList, regulatoryElements = NULL, 
    minSetSize = 15, sizeFilterMethod="posORneg",
    method = c("HT","JC"), pValueCutoff = 0.05, 
    pAdjustMethod = "BH", verbose = TRUE)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} already evaluated by the 
\code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}} method.
}
\item{geneSetList}{
a list with gene sets.
}
\item{regulatoryElements}{
a vector of valid regulatory elements (e.g. transcription factors).
}
\item{minSetSize}{
a single integer or numeric value specifying the minimum number of 
elements in a gene set that must map to elements of the gene universe. 
Gene sets with fewer than this number are removed from the analysis.
}
\item{sizeFilterMethod}{
a single character value specifying the use of the 'minSetSize' argument, which is applyed to regulon's positive and negative targets. Options: "posANDneg", "posORneg", "posPLUSneg". For "posANDneg", the number of both positive and negative targets should be > 'minSetSize'; for "posORneg", the number of either positive or negative targets should be > 'minRegulonSize'; and for "posPLUSneg", the number of all targets should be > 'minSetSize'.
}
\item{method}{
a string specifying the method used to assess the association between gene sets and regulons (see 'Details').
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered significant
(this parameter only affects the 'HT' option).
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be 
used (see 'p.adjust' for details) (this parameter only affects the 'HT' option).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE).
}
}

\details{
The 'HT' option assesses the overlap between gene sets and regulons using a hypergeometric test, and returns a data frame with the overlap statistics. The 'JC' option assesses the overlap between gene sets and regulons using the Jaccard Coefficient (JC), and retuns a matrix with JC values.
}

\value{
Either a data frame or a numeric matrix with association statistics between gene sets vs. regulons.
}

\author{Mauro Castro}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
}

\examples{

data(tniData)

\dontrun{

#compute regulons
rtni <- tni.constructor(expData=tniData$expData, 
        regulatoryElements=c("PTTG1","E2F2","FOXM1","E2F3","RUNX2"), 
        rowAnnotation=tniData$rowAnnotation)
rtni <- tni.permutation(rtni)
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni)

#load a gene set collection 
#here, we build three random gene sets for demonstration
geneset1 <- sample(tniData$rowAnnotation$SYMBOL,50)
geneset2 <- sample(tniData$rowAnnotation$SYMBOL,50)
geneset3 <- sample(tniData$rowAnnotation$SYMBOL,50)
geneSetList <- list(geneset1=geneset1,
                    geneset2=geneset2,
                    geneset3=geneset3)

#run the overlap analysis
ovstats <- tni.overlap.genesets(rtni, geneSetList, pValueCutoff = 1)

}
}

\keyword{annotate}

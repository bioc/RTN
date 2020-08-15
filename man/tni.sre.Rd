\name{tni.sre}
\alias{tni.sre}

\title{Subgroup Regulon Enrichment for TNI-class objects.}

\usage{

tni.sre(object, sampleGroups, regulatoryElements = NULL, 
        pValueCutoff = 0.05, pAdjustMethod = "BH")
  
}

\arguments{

\item{object}{A \linkS4class{TNI} object.}

\item{sampleGroups}{either a list featuring sample groups or a string 
indicating a group varible available in the TNI object.}

\item{regulatoryElements}{an optional string vector specifying regulons 
to use for the analysis.}

\item{pValueCutoff}{a single numeric value specifying the cutoff for 
p-values considered significant.}

\item{pAdjustMethod}{a single character value specifying the p-value 
adjustment method to be used (see 'p.adjust' for details).}

}


\value{
A TNI-class object with the results of the subgroup regulon enrichment
added to the results slot. To recover the results, use tni.get(object, "regulonEnrichment")
}

\description{
This method evaluates which regulons are enriched in sample groups, given a
grouping variable. It performs Fisher's Exact Test whether a regulon is
positively or negatively enriched in a subgroup using regulon activity.
}

\seealso{
\code{\link[RTN:tni.plot.sre]{tni.plot.sre}}
}

\examples{

# load tniData
data(tniData)

\dontrun{

# compute regulons
rtni <- tni.constructor(expData=tniData$expData, 
        regulatoryElements=c("PTTG1","E2F2","FOXM1","E2F3","RUNX2"), 
        rowAnnotation=tniData$rowAnnotation)
rtni <- tni.permutation(rtni)
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni)

#run GSEA2 analysis pipeline
rtni <- tni.gsea2(rtni)

# set sample groups
colAnnotation <- tni.get(rtni, "colAnnotation")
sampleGroups <- list(G1=colAnnotation$ID[1:60], 
                G2=colAnnotation$ID[61:90], 
                G3=colAnnotation$ID[91:120])

# run subgroup regulon enrichment analysis
rtni <- tni.sre(rtni, sampleGroups)

# get results
tni.get(rtni, "subgroupEnrichment")

# for a heatmap representation, see the tni.plot.sre() function.
}

}

\keyword{subgroups}

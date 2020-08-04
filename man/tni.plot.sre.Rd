\name{tni.plot.sre}
\alias{tni.plot.sre}

\title{Plot Subgroup Regulon Enrichment for TNI-class objects.}

\usage{
tni.plot.sre(object, nGroupsEnriched = NULL, nTopEnriched = NULL, 
             colors = c("blue","white","red"), breaks = seq(-1.5, 1.5, 0.1), 
             markEnriched = TRUE, ...)
}

\arguments{
\item{object}{A \linkS4class{TNI}-class object.}

\item{nGroupsEnriched}{a filter to keep 'nGroupsEnriched' regulons; 
a single integer specifying how many subgroups a regulon has to be 
enriched for it to appear in the rows of the heatmap 
(it must be use either 'nGroupsEnriched' or 'nTopEnriched').}

\item{nTopEnriched}{a filter to keep 'nTopEnriched' regulons; 
a single integer specifying how many regulons will be shown for 
each group. The top regulons are chosen by significance
(it must be use either 'nTopEnriched' or 'nGroupsEnriched').}

\item{colors}{a vector of color for the 'pheatmap'.}

\item{breaks}{a numerical vector of breaks for the 'pheatmap'.}

\item{markEnriched}{a single logical value. If TRUE, asterisks are added to
cells of heatmap that were found to be significant.}

\item{...}{parameters passed to 'pheatmap' for customization.}
}

\value{
A heatmap of the subgroup regulon enrichment results.
}

\description{
This method plots the results of the subgroup regulon enrichment analysis in
a heatmap. The rows of the heatmap represent enriched regulons, while the
columns show the subgroups. The plotted values correspond to average regulon
activity for a regulon in a subgroup. Enriched values can be marked.
}

\seealso{
\code{\link[RTN:tni.sre]{tni.sre}}
}

\examples{
# load tniData
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

#run GSEA2 analysis pipeline
rtni <- tni.gsea2(rtni)

# set sample groups
colAnnotation <- tni.get(rtni, "colAnnotation")
sampleGroups <- list(G1=colAnnotation$ID[1:60], 
                G2=colAnnotation$ID[61:90], 
                G3=colAnnotation$ID[91:120])

# run subgroup regulon enrichment analysis
rtni <- tni.sre(rtni, sampleGroups)

# plot results
tni.plot.sre(rtni)
}
}

\name{tna.get}
\alias{tna.get}

\title{
Get information from individual slots in a TNA object.
}

\description{
Get information from individual slots in a TNA object. Available results from a previous 
analysis can be selected either by pvalue cutoff (default) or top significance.
}

\usage{
tna.get(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE, 
idkey=NULL)
}

\arguments{
\item{object}{
an object of class \code{\link[RTN:TNA-class]{TNA-class}}.
}
\item{what}{
a single character value specifying which information should be retrieved from the slots. 
Options: 'summary', 'status', 'para', 'pheno', 'hits', 'regulatoryElements', 'tnet', 
'refnet', 'regulons', 'refregulons', 'regulons.and.mode', 'refregulons.and.mode', 
'rowAnnotation', 'colAnnotation', 'mra', 'gsea1', 'gsea2', 'gsea2summary'.
}
\item{order}{
a single logical value specifying whether or not the output data should be ordered by 
significance. Valid only for 'mra', 'gsea1' and 'gsea2' options.
}
\item{ntop}{
a single integer value specifying to select how many results of top significance from 
'mra', 'gsea1' and 'gsea2' options.
}
\item{reportNames}{
a single logical value specifying to report regulons with 'names' (when reportNames=TRUE) 
or not (when reportNames=FALSE). This option is effective only if transcription factors 
were named with alternative identifiers in the pre-processing analysis. It takes effect 
on 'mra', 'gsea1' and 'gsea2' options.
}
\item{idkey}{
an optional single character value specifying an ID name from the available 'TNA' 
annotation to be used as alias for data query outputs (obs. it has no effect on 
consolidated tables).
}
}

\details{

Options for the 'what' argument retrieve the following types of information:

\describe{

\item{summary}{
A list summarizing parameters and results available in the TNA object.
}
\item{status}{
A vector indicating the status of each available method in the pipeline.
}
\item{para}{
A list with the parameters used by each available method in the pipeline.
}
\item{pheno}{
A numeric vector of phenotypes named by gene identifiers 
(see \code{\link[RTN:tni2tna.preprocess]{tni2tna.preprocess}}).
}
\item{hits}{
A character vector of gene identifiers for those considered as hits 
(see \code{\link[RTN:tni2tna.preprocess]{tni2tna.preprocess}}).
}
\item{regulatoryElements}{
A vector of regulatory elements (e.g. transcription factors).
}
\item{tnet}{
A data matrix with MI values, evaluated by the DPI filter. MI values are 
computed between regulators and targets, with regulators on cols and targets 
on rows. Note that signals (+/-) are assigned to the inferred associations in 
order to represent the 'mode of action', which is derived from Pearson's 
correlation between regulators and targets.
}
\item{refnet}{
A data matrix with MI values (not evaluated by the DPI filter). MI values are 
computed between regulators and targets, with regulators on cols and targets 
on rows. Note that signals (+/-) are assigned to the inferred associations in 
order to represent the 'mode of action', which is derived from Pearson's 
correlation between regulators and targets.
}
\item{regulons}{
A list with regulons extracted from the 'tnet' data matrix.
}
\item{refregulons}{
A list with regulons extracted from the 'refnet' data matrix.
}
\item{regulons.and.mode}{
A list with regulons extracted from the 'tnet' data matrix, including the 
assiged 'mode of action'.
}
\item{refregulons.and.mode}{
A list with regulons extracted from the 'refnet' data matrix, including the 
assiged 'mode of action'.
}
\item{rowAnnotation}{
A data frame with probe-to-gene annotation.
}
\item{colAnnotation}{
A data frame with sample annotation.
}
\item{mra}{
A data frame with results from the \code{\link[RTN:tna.mra]{tna.mra}} 
analysis pipeline.
}
\item{gsea1}{
A data frame with results from the \code{\link[RTN:tna.gsea1]{tna.gsea1}} 
analysis pipeline.
}
\item{gsea2}{
A data frame with results from the \code{\link[RTN:tna.gsea2]{tna.gsea2}} 
analysis pipeline.
}
}

}

\value{
Get the slot content from a \code{\link[RTN:TNA-class]{TNA-class}} object.
}

\author{
Mauro Castro
}

\examples{

data(tniData)
data(tnaData)

\dontrun{

rtni <- tni.constructor(expData=tniData$expData, 
        regulatoryElements=c("PTTG1","E2F2","FOXM1","E2F3","RUNX2"), 
        rowAnnotation=tniData$rowAnnotation)
rtni <- tni.permutation(rtni)
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni)
rtna <- tni2tna.preprocess(rtni, phenotype=tnaData$phenotype, 
        hits=tnaData$hits, phenoIDs=tnaData$phenoIDs)

# run MRA analysis pipeline
rtna <- tna.mra(rtna)

# check summary
tna.get(rtna,what="summary")

# get results, e.g., from the MRA analysis
tna.get(rtna,what="mra")

}
}

\keyword{methods}

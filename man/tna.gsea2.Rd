\name{tna.gsea2}
\alias{tna.gsea2}

\title{
Two-tailed Gene Set Enrichment Analysis (GSEA) over a list of regulons.
}

\description{
This function takes a TNA object and returns a CMAP-like analysis obtained by two-tailed GSEA over a list of regulons in a transcriptional network (with multiple hypothesis testing corrections).
}

\usage{
tna.gsea2(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, 
sizeFilterMethod="posORneg", nPermutations=1000, exponent=1, tnet="dpi", 
signature=c("phenotype","hits"), tfs=NULL, verbose=TRUE, doSizeFilter=NULL)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNA' \code{\link[RTN:TNA-class]{TNA-class}}.
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered 
significant.
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be 
used (see 'p.adjust' for details).
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of 
elements in a regulon that must map to elements of the gene universe. 
Gene sets with fewer than this number are removed from the analysis.
}
\item{sizeFilterMethod}{
a single character value specifying the use of the 'minRegulonSize' argument, which is applyed to regulon's positive and negative targets. Options: "posANDneg", "posORneg", "posPLUSneg". For "posANDneg", the number of both positive and negative targets should be > 'minRegulonSize'; for "posORneg", the number of either positive or negative targets should be > 'minRegulonSize'; and for "posPLUSneg", the number of all targets should be > 'minRegulonSize'.
}
\item{nPermutations}{
a single integer or numeric value specifying the number of permutations 
for deriving p-values in GSEA.
}
\item{exponent}{
a single integer or numeric value used in weighting phenotypes in GSEA 
(see 'gseaScores' function at HTSanalyzeR).
}
\item{tnet}{
a single character value specifying which transcriptional network should to used to 
compute the GSEA analysis. Options: "dpi" and "ref".
}
\item{signature}{
a single character value specifying which signature to use in the GSEA method. This could be the 'phenotype' already provided by the user (usually log2 differential expression values) or a 'phenotype' derived from the 'hits'; see \code{\link[RTN:tni2tna.preprocess]{tni2tna.preprocess}} for details on 'phenotype' and 'hits' parameters.
}
\item{tfs}{
an optional vector with transcription factor identifiers.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE).
}
\item{doSizeFilter}{
'doSizeFilter' is deprecated, please use the 'filterSize' parameter.
}
}

\value{
a data frame in the slot "results", see 'gsea2' option in \code{\link[RTN:tna.get]{tna.get}}.
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:TNA-class]{TNA-class}}
\code{\link[RTN:tna.plot.gsea2]{tna.plot.gsea2}}
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

#run GSEA2 analysis pipeline
rtna <- tna.gsea2(rtna)

#get results
tna.get(rtna, what="gsea2")

# run parallel version with SNOW package!
library(snow)
options(cluster=snow::makeCluster(3, "SOCK"))
rtna <- tna.gsea2(rtna)
stopCluster(getOption("cluster"))
}
}

\keyword{GSEA}

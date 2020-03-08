\name{tni.plot.checks}
\alias{tni.plot.checks}

\title{Plot regulon target counts.}

\description{
This funtion can help to check whether the numbers of positive and negative targets are reasonably well balanced in the regulons.
}

\usage{
tni.plot.checks(object, minRegulonSize = 15, option = c("barplot","edf","points"))
}

\arguments{
\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}} already evaluated by the 
\code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}} method.
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of elements in a regulon (only affects the 'barplot' option).
}
\item{option}{
plot option.
}
}

\note{
We have observed that transcription factor (TF) regulons reconstructed from RTN exhibit different proportions of positive and negative targets. While the proportion can vary between different regulons, we have observed a consistent higher proportion of positive targets, especially when using RNA-seq data. RTN uses mutual information (MI) to assess TF-target associations, assigning the direction of the inferred associations by Spearman's correlations. Dam et al. (2018) have acknowledged that different RNA-seq normalization methods introduce different biases in co-expression analysis, usually towards positive correlation, possibly affected by read-depth differences between samples and the large abundance of 0 values present in RNA-seq-derived expression matrices. This funtion can help to check whether the numbers of positive and negative target genes are reasonably well balanced in the regulons.
}

\value{
A plot showing the distribution of regulons' positive and negative targets.
}

\author{Mauro Castro, Gordon Robertson}

\references{
Dam et al. Gene co-expression analysis for functional classification and gene-disease predictions. Brief Bioinform. 2018 Jul 20;19(4):575-592. doi: 10.1093/bib/bbw139.
}

\examples{

data(tniData)

\dontrun{

# preprocessing
rtni <- tni.constructor(expData=tniData$expData, 
        regulatoryElements=c("PTTG1","E2F2","FOXM1","E2F3","RUNX2"), 
        rowAnnotation=tniData$rowAnnotation)

# compute regulons
rtni <- tni.permutation(rtni, nPermutations = 1000)
rtni <- tni.permutation(rtni)
rtni <- tni.dpi.filter(rtni)

# check target distribution
tni.plot.checks(rtni)

}
}

\keyword{methods}


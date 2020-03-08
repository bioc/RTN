\name{tni.permutation}
\alias{tni.permutation}

\title{
Inference of transcriptional networks.
}

\description{
This function takes a TNI object and returns a transcriptional network inferred by mutual information
(with multiple hypothesis testing corrections).
}

\usage{
tni.permutation(object, pValueCutoff=0.01, pAdjustMethod="BH", globalAdjustment=TRUE,
        estimator="spearman", nPermutations=1000, pooledNullDistribution=TRUE, 
        boxcox=TRUE, parChunks=NULL, verbose=TRUE)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered 
significant.
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be 
used (see 'p.adjust' for details).
}
\item{globalAdjustment}{
a single logical value specifying to run global p.value adjustments (when 
globalAdjustment=TRUE) or not (when globalAdjustment=FALSE).
}
\item{estimator}{
a character string specifying the mutual information estimator. 
One of "pearson", "kendall", or "spearman" (default).
}
\item{nPermutations}{
a single integer value specifying the number of permutations for deriving TF-target 
p-values in the mutual information analysis. If running in parallel, nPermutations 
should be greater and multiple of parChunks.
}
\item{pooledNullDistribution}{
a single logical value specifying to run the permutation analysis with pooled regulons 
(when pooledNullDistribution=TRUE) or not (when pooledNullDistribution=FALSE).
}
\item{boxcox}{
a single logical value specifying to use Box-Cox procedure to find a transformation of inferred associations that approaches normality (when boxcox=TRUE) or not (when boxcox=FALSE). Dam et al. (2018) have acknowledged that different RNA-seq normalization methods introduce different biases in co-expression analysis, usually towards positive correlation, possibly affected by read-depth differences between samples and the large abundance of 0 values present in RNA-seq-derived expression matrices. In order to correct this positive correlation bias we suggest using this box-cox correction strategy. See \code{\link[car:powerTransform]{powerTransform}} and \code{\link[car:bcPower]{bcPower}}.
}
\item{parChunks}{
an optional single integer value specifying the number of permutation chunks to 
be used in the parallel analysis (effective only for "pooledNullDistribution = TRUE").
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE)
}

}

\value{
a mutual information matrix in the slot "results" containing a reference transcriptional network, 
see 'tn.ref' option in \code{\link[RTN:tni.get]{tni.get}}.
}

\author{Mauro Castro}

\references{
Dam et al. Gene co-expression analysis for functional classification and gene-disease predictions. Brief Bioinform. 2018 Jul 20;19(4):575-592. doi: 10.1093/bib/bbw139.
}

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

# linear version (set nPermutations >= 1000)
rtni <- tni.permutation(rtni, nPermutations = 100)

## parallel version with SNOW package!
#library(snow)
#options(cluster=snow::makeCluster(3, "SOCK"))
#rtni<-tni.permutation(rtni)
#stopCluster(getOption("cluster"))
}
}


\keyword{methods}



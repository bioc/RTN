\name{tni.gsea2}
\alias{tni.gsea2}

\title{
Compute regulon activity by calling GSEA2 (two-tailed Gene Set Enrichment Analysis) 
algorithm
}

\description{
Uses GSEA2 algorithm to compute regulon activity for 
\code{\link[RTN:TNI-class]{TNI-class}} objects.
}

\usage{
tni.gsea2(object, minRegulonSize=15, sizeFilterMethod="posORneg", scale=FALSE, 
exponent=1, tnet="dpi", regulatoryElements=NULL, features=NULL, samples=NULL, 
refsamp=samples, log=TRUE, alternative=c("two.sided", "less", "greater"), 
targetContribution=FALSE, additionalData=FALSE, verbose=TRUE, doSizeFilter=NULL)
}

\arguments{

\item{object}{
a preprocessed object of class 'TNI' \code{\link[RTN:TNI-class]{TNI-class}}.
}
\item{minRegulonSize}{
a single integer or numeric value specifying the minimum number of elements in 
a regulon. Regulons smaller than this number are removed from the analysis.
}
\item{sizeFilterMethod}{
a single character value specifying the use of the 'minRegulonSize' argument, which is applyed to regulon's positive and negative targets. Options: "posANDneg", "posORneg", "posPLUSneg". For "posANDneg", the number of both positive and negative targets should be > 'minRegulonSize'; for "posORneg", the number of either positive or negative targets should be > 'minRegulonSize'; and for "posPLUSneg", the number of all targets should be > 'minRegulonSize'.
}
\item{scale}{A logical value specifying if expression values should be 
centered and scaled across samples (when verbose=TRUE) or 
not (when verbose=FALSE).
}
\item{exponent}{
a single integer or numeric value used in weighting phenotypes in GSEA.
}
\item{tnet}{
can take values of 'ref', 'dpi' or 'cdt'. It refers to the version of the 
regulatory network that will be used for GSEA analysis.
}
\item{regulatoryElements}{
an optional vector with transcription factor identifiers.
}
\item{features}{
a string vector listing features for feature selection.
}
\item{samples}{
an optional string vector listing the sample names for which 
will be computed the GSEA2.
}
\item{refsamp}{
an optional string vector listing the names of the reference samples for 
differential expression calculations. If not provided, then the average of 
all samples will be used as reference.
}
\item{log}{
a logical value. If TRUE, it will check whether the expression values are provided as logged data; if not, it will performe a log2 transformation on expression values before the differential expression calculations.
}
\item{alternative}{
a character string specifying the alternative hypothesis, must be one 
of "two.sided" (default), "greater" or "less".
}
\item{targetContribution}{
This argument is used for internal calls. A single logical value specifying 
to return the contribution of each target in enrichment scores 
(when verbose=TRUE) or not (when verbose=FALSE).
}
\item{additionalData}{
This argument is used for internal calls. A single logical value specifying 
to return the additional data objects (when verbose=TRUE) or not 
(when verbose=FALSE).
}
\item{verbose}{
a single logical value specifying to display detailed messages (when verbose=TRUE) 
or not (when verbose=FALSE).
}
\item{doSizeFilter}{
'doSizeFilter' is deprecated, please use the 'sizeFilterMethod' parameter.
}
}

\value{
a list with enrichment scores for all samples in the TNI. The list contains the following elements:
  \describe{
    \item{differential:}{
      A numeric \code{"matrix"} with differential enrichment scores (dES).
    }
    \item{positive:}{
      A numeric \code{"matrix"} with enrichment scores from positive targets.
    }
    \item{negative:}{
      A numeric \code{"matrix"} with enrichment scores from negative targets.
    }
    \item{status:}{
      A numeric \code{"matrix"} with discretized scores derived from the dES values.
    } 
    \item{regulatoryElements:}{
      A character vector listing the regulatory elements assessed by the GSEA-2T algorithm.
    }     
    \item{sections:}{
      A single numeric value used in internal plots.
    }    
  }
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}
\code{\link[RTN:tna.gsea2]{tna.gsea2}}
\code{\link[RTN:tna.plot.gsea2]{tna.plot.gsea2}}
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

#run GSEA2 analysis pipeline
rtni <- tni.gsea2(rtni)

#get results
regulonActivity <- tni.get(rtni, what = "regulonActivity")

#parallel version with SNOW package!
library(snow)
options(cluster=snow::makeCluster(3, "SOCK"))
rtni <- tni.gsea2(rtni)
stopCluster(getOption("cluster"))

}
}

\keyword{GSEA}
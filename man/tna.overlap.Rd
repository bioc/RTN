\name{tna.overlap}
\alias{tna.overlap}

\title{
Overlap analysis over a list of regulons (deprecated function).
}

\description{
This function is provided for compatibility with older versions of \sQuote{RTN} only, and will be defunct at the next release (see details).
}

\usage{
tna.overlap(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, tnet="ref", 
            tfs=NULL, verbose=TRUE)
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
\item{tnet}{
a single character value specifying which transcriptional network should to used to 
compute the overlap analysis. Options: "dpi" and "ref".
}
\item{tfs}{
an optional vector with transcription factor identifiers.
}
\item{verbose}{
a single logical value specifying to display detailed messages (when 
verbose=TRUE) or not (when verbose=FALSE).
}
}

\details{
The 'tna.overlap' function is deprecated and will be made defunct. For pairwise regulon analysis, please consider using the 'RTNduals' package.
}

\value{
a data frame in the slot "results".
}

\keyword{overlap}

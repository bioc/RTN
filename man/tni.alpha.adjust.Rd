\name{tni.alpha.adjust}
\alias{tni.alpha.adjust}

\title{Adjust the significance level for two datasets.}

\description{
When analyzing two datasets that have different numbers of samples, this function can be used to assist the choice of a p value threshold for the tni.permutation() function, in order that RTN will return regulon results that have been generated with similar tradeoffs between Type I and Type II errors for both datasets. Doing this should help ensure that it is reasonable to compare the regulons in the two datasets. 
}

\usage{
tni.alpha.adjust(nB, nA, alphaA, betaA = 0.2)
}

\arguments{
\item{nB}{
a single integer specifying the number of samples in dataset 'B'.
}
\item{nA}{
a single integer specifying the number of samples in dataset 'A'.
}
\item{alphaA}{
a single numeric value specifying alpha for dataset 'A' (Type I error probability).
}
\item{betaA}{
a single numeric value specifying beta for dataset 'A' (Type II error probability).
}
}

\note{
The 'tni.alpha.adjust' function calls the \code{\link[pwr:pwr.r.test]{pwr.r.test}} function, which uses the 'uniroot' function to solve a power equation. The 'uniroot' function aims to find a root in a given interval, searching from lower to upper end-points. As the upper end-point must be strictly larger than the lower end-point, in order to avoid an error when searching the root, 'nA' must be greater than or equal to 'nB' (i.e. 'nB' is expected to be the smallest data set). Also, please note that 'uniroot' eventually will not find the root for the input arguments, especially when searching thresholds of low stringency (we suggest to avoid setting 'alphaA' > 0.01).
}

\value{
Significance level for 'nB', given 'nA', 'alphaA' and 'betaA'.
}

\author{Mauro Castro, Gordon Robertson}

\seealso{
\code{\link[pwr:pwr.r.test]{pwr.r.test}}
}

\examples{
# estimate 'alphaB' for 'nB', given 'nA', 'alphaA' and 'betaA'
alphaB <- tni.alpha.adjust(nB = 100, nA = 300, alphaA = 1e-5, betaA = 0.2)

}

\keyword{methods}


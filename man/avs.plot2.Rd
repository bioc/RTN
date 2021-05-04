\name{avs.plot2}
\alias{avs.plot2}

\title{
Plot results from AVS methods, multiple plots.
}

\description{
This function takes an AVS object and plots results from the VSE and EVSE methods.
}

\usage{
avs.plot2(object, what="evse", fname=what, width=14, height=2.5, 
width.panels=c(1,3), rmargin=1, at.x=seq(-4,8,2), decreasing=TRUE, 
ylab="Annotation", xlab="Clusters of risk-associated and linked SNPs", 
tfs=NULL)
}

\arguments{
\item{object}{
an object of class 'AVS' \code{\link[RTN:AVS-class]{AVS-class}}.
}
\item{what}{
a single character value specifying which analysis should be used. Options: "vse" and "evse".
}
\item{fname}{
a character value specifying the name of output file.
}
\item{height}{
a numeric value specifying the height of the graphics region in inches.
}
\item{width}{
a numeric value specifying the width of the graphics region in inches.
}
\item{width.panels}{
a vector of the form c(width1, width2) specifying the proportional width 
of the 1st and 2nd panels of the plot, respectively.
}
\item{rmargin}{
a numeric value specifying the right margin in inches.
}
\item{at.x}{
a numeric vector specifying which x-axis tickpoints are to be drawn.
}
\item{decreasing}{
a logical value, used to sort by EVSE scores.
}
\item{ylab}{
a character value specifying the y-axis label.
}
\item{xlab}{
a character value specifying the x-axis label (on the top of the grid image).
}
\item{tfs}{
an optional vector with annotation identifiers (e.g. transcription factor).
}
}

\value{
A plot showing results from the VSE and EVSE methods.
}

\author{
Mauro Castro
}

\examples{

# see 'avs.vse' and 'avs.evse' methods.

}

\keyword{VSE}

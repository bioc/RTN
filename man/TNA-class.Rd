\name{TNA-class}
\Rdversion{1.1}
\docType{class}
\alias{TNA-class}
\alias{tna.mra,TNA-method}
\alias{tna.gsea1,TNA-method}
\alias{tna.gsea2,TNA-method}
\alias{tna.get,TNA-method}

\title{Class \code{"TNA"}: an S4 class for Transcriptional Network Analysis.}

\description{
This S4 class includes a series of methods to do enrichment analyses in transcriptional networks.
}

\section{Objects from the Class}{
Objects can be created by calls of the form 
\code{new("TNA", referenceNetwork, transcriptionalNetwork, regulatoryElements, phenotype, hits)}.
}

\section{Slots}{
  \describe{
    \item{\code{referenceNetwork}:}{Object of class \code{"matrix"}, 
    an optional partial co-expression matrix.
    }
    \item{\code{transcriptionalNetwork}:}{Object of class \code{"matrix"}, 
    a partial co-expression matrix.
    }
    \item{\code{regulatoryElements}:}{Object of class \code{"char_Or_NULL"}, 
    a vector of regulatory elements (e.g. transcription factors).
    }
    \item{\code{phenotype}:}{Object of class \code{"num_Or_int"}, a numeric or integer 
    vector of phenotypes named by gene identifiers. 
    }
    \item{\code{hits}:}{Object of class \code{"character"}, a character vector of 
    gene identifiers for those considered as hits. 
    }
    \item{\code{gexp}:}{Object of class \code{"matrix"},
    a gene expression matrix.
    }    
    \item{\code{rowAnnotation}:}{Object of class \code{"data.frame"},
    a data frame with row annotation (e.g. probe-to-gene information).
    }
    \item{\code{colAnnotation}:}{Object of class \code{"data.frame"},
    a data frame with column annotation (e.g. sample information).
    }
    \item{\code{listOfReferenceRegulons}:}{Object of class \code{"list"}, a list of regulons 
    derived from the referenceNetwork.
    }
    \item{\code{listOfRegulons}:}{Object of class \code{"list"}, a list of regulons 
    derived from the transcriptionalNetwork (a 'regulon' is a vector of genes or potential 
    transcription factor targets). 
    }    
    \item{\code{listOfModulators}:}{Object of class \code{"list"}, a list of modulators 
    derived from the \code{\link[RTN:tni.conditional]{tni.conditional}} analysis.
    }
    \item{\code{para}:}{Object of class \code{"list"}, a list of parameters for 
    transcriptional network analysis. These parameters are those listed in the functions 
    \code{\link[RTN:tna.mra]{tna.mra}}, \code{\link[RTN:tna.gsea1]{tna.gsea1}}, 
    and \code{\link[RTN:tna.gsea2]{tna.gsea2}}.
    }
    \item{\code{results}:}{Object of class \code{"list"},
    a list of results (see return values in the functions 
		\code{\link[RTN:tna.mra]{tna.mra}},\code{\link[RTN:tna.gsea1]{tna.gsea1}},
    and \code{\link[RTN:tna.gsea2]{tna.gsea2}})
    }
    \item{\code{summary}:}{Object of class \code{"list"},
  	a list of summary information for \code{transcriptionalNetwork}, 
    \code{regulatoryElements}, \code{phenotype},\code{listOfRegulons},
    \code{para}, and \code{results}.    
    }
    \item{\code{status}:}{Object of class \code{"character"},
    a character value specifying the status of the TNI object 
    based on the available methods.
    }    
  }
}

\section{Methods}{
  \describe{
    \item{tna.mra}{\code{signature(object = "TNA")}: 
    see \code{\link[RTN:tna.mra]{tna.mra}} }
    \item{tna.gsea1}{\code{signature(object = "TNA")}: 
    see \code{\link[RTN:tna.gsea1]{tna.gsea1}} }
    \item{tna.gsea2}{\code{signature(object = "TNA")}: 
    see \code{\link[RTN:tna.gsea2]{tna.gsea2}} }
    \item{tna.get}{\code{signature(object = "TNA")}: 
    see \code{\link[RTN:tna.get]{tna.get}} }
	 }
}

\author{
Mauro Castro
}

\seealso{
\code{\link[RTN:TNI-class]{TNI-class}}.
\code{\link[RTN:tni2tna.preprocess]{tni2tna.preprocess}}.
}

\examples{

## see 'tni2tna.preprocess' method!

}

\keyword{classes}

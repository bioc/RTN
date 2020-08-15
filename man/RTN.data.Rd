\name{RTN.data}
\alias{tniData}
\alias{tnaData}
\alias{tfsData}
\alias{pksData}
\docType{data}

\title{Pre-processed datasets for the RTN package.}

\description{
Datasets used to demonstrate RTN main functions.
}

\format{
\strong{tniData}, \strong{tnaData}, and \strong{tfsData}
}

\details{

The \strong{tniData} and \strong{tnaData} datasets were extracted, pre-processed and size-reduced from Fletcher et al. (2013) and Curtis et al. (2012). They consist of two lists used in the RTN vignettes for demonstration purposes only. The \strong{tniData} list contains the 'expData', 'rowAnnotation' and 'colAnnotation' R objects, while the \strong{tnaData} list contains the 'phenotype', 'phenoIDs' and 'hits' R objects.

The \strong{tfsData} consists of a list with gene annotation for human transcription factors (TFs),  compiled from 6 resources (Lambert et. al 2018; Carro et al. 2010; Vaquerizas et al. 2009; D. L. Fulton et al. 2009; Yusuf et al. 2012; and Ravasi et al. 2010), and provided here in the form of 'data frame' objects that include ENTREZ and HGNC gene symbol annotation.

The \strong{pksData} consists of a list with gene annotation for human protein kinases, retrived from the Swiss-Prot Protein Knowledgebase (https://www.uniprot.org/docs/pkinfam), release 2020_02 of 22-Apr-2020. 

\describe{
  \item{tniData$expData}{a named gene expression matrix with 120 samples (a subset from the Fletcher2013b).}
  \item{tniData$rowAnnotation}{a data.frame of characters with gene annotation (a subset from the Fletcher2013b).}
  \item{tniData$colAnnotation}{a data.frame of characters with sample annotation (a subset from the Fletcher2013b).}  
  \item{tnaData$phenotype}{a named numeric vector with differential gene expression data.}
  \item{tnaData$phenoIDs}{a data.frame of characters with probe ids matching a secundary annotation source (e.g. Probe-to-ENTREZ).}
  \item{tnaData$hits}{a character vector with genes differentially expressed.}
  \item{tfsData$Lambert2018}{a data.frame listing TFs from Lambert et. al (2018).}
  \item{tfsData$Yusuf2012}{a data.frame listing TFs from Yusuf et al. (2012).}  
  \item{tfsData$Carro2010}{a data.frame listing TFs from Carro et al. (2010).} 
  \item{tfsData$Ravasi2010}{a data.frame listing TFs from Ravasi et al. (2010).} 
  \item{tfsData$Fulton2009}{a data.frame listing TFs from Fulton et al. (2009).} 
  \item{tfsData$Vaquerizas2009}{a data.frame listing TFs from Vaquerizas et al. (2009).} 
  \item{tfsData$all}{a data.frame listing all TFs.}
  \item{pksData$pkinfam2020}{a data.frame listing human protein kinases.}
}
}

\references{

Carro, M.S. et al., \emph{The transcriptional network for mesenchymal transformation of brain tumors.} Nature, 463(7279):318-325, 2010.

Curtis C. et al., \emph{The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups.} Nature, 486(7403):346-352, 2012.

Fletcher, M.N.C. et al., \emph{Master regulators of FGFR2 signalling and breast cancer risk.} Nature Communications, 4:2464, 2013.

Fulton, D.L. et al., \emph{TFCat: the curated catalog of mouse and human transcription factors.} Genome Biology, 10(3):R29, 2009.

Lambert, S.A. et al., \emph{The Human Transcription Factors.} Cell, 172(4):650-665, 2018.

Ravasi, T. et al., \emph{An atlas of combinatorial transcriptional regulation in mouse and man.} Cell, 140(5):744-752, 2010.

Vaquerizas, J.M. et al., \emph{A census of human transcription factors: function, expression and evolution.} Nature Reviews Genetics, 10(4):252-263, 2009.

Yusuf, D. et al., \emph{The Transcription Factor Encyclopedia.} Genome Biology, 13(3):R24, 2012.

}

\examples{
    data(tniData)
    data(tnaData)
    data(tfsData)
    data(pksData)
}

\keyword{dataset}


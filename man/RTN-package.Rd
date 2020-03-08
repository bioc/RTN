\name{RTN-package}
\alias{RTN-package}
\alias{RTN}
\docType{package}

\title{
Reconstruction and Analysis of Transcriptional Networks.
}

\description{
A transcriptional regulatory network (TRN) consists of a collection of transcription factors (TFs) and the regulated target genes. TFs are regulators that recognize specific DNA sequences and guide the expression of the genome, either activating or repressing the expression the target genes. The set of genes controlled by the same TF forms a regulon. This package provides classes and methods for the reconstruction of TRNs and analysis of regulons.
}

\section{Index}{

\tabular{ll}{
\link[RTN:TNI-class]{TNI-class}: \tab an S4 class for Transcriptional Network Inference.\cr
\link[RTN:tni.preprocess]{tni.preprocess}: \tab a preprocessing method for objects of class TNI.\cr
\link[RTN:tni.permutation]{tni.permutation}: \tab inference of transcriptional networks.\cr
\link[RTN:tni.bootstrap]{tni.bootstrap}: \tab inference of transcriptional networks.\cr
\link[RTN:tni.dpi.filter]{tni.dpi.filter}: \tab data processing inequality (DPI) filter.\cr
\link[RTN:tni.conditional]{tni.conditional}: \tab conditional mutual information analysis.\cr
\link[RTN:tni.get]{tni.get}: \tab get information from individual slots in a TNI object.\cr
\link[RTN:tni.graph]{tni.graph}: \tab compute a graph from TNI objects.\cr
\link[RTN:tni.gsea2]{tni.gsea2}: \tab compute regulon activity.\cr
\link[RTN:tni.prune]{tni.prune}: \tab prune regulons to remove redundant targets for regulon activity analysis.\cr
\link[RTN:tni.sre]{tni.sre}: \tab subgroup regulon difference analysis.\cr
\link[RTN:tni.plot.sre]{tni.plot.sre}: \tab plot subgroup regulon rnrichment .\cr
\link[RTN:tni.regulon.summary]{tni.regulon.summary}: \tab return a summary of network and regulons. \cr
\link[RTN:tni.plot.checks]{tni.plot.checks}: \tab plot regulon target counts.\cr
\link[RTN:tni.alpha.adjust]{tni.alpha.adjust}: \tab adjust the significance level for two datasets.\cr
\link[RTN:tni.replace.samples]{tni.replace.samples}: \tab replace samples of an existing TNI-class objects. \cr
\link[RTN:tni2tna.preprocess]{tni2tna.preprocess}: \tab a preprocessing method for objects of class TNI.\cr
\link[RTN:TNA-class]{TNA-class}: \tab an S4 class for Transcriptional Network Analysis.\cr
\link[RTN:tna.mra]{tna.mra}: \tab master regulator analysis (MRA) over a list of regulons.\cr
\link[RTN:tna.gsea1]{tna.gsea1}: \tab one-tailed gene set enrichment analysis (GSEA) over a list of regulons.\cr
\link[RTN:tna.gsea2]{tna.gsea2}: \tab two-tailed gene set enrichment analysis (GSEA) over a list of regulons.\cr
\link[RTN:tna.get]{tna.get}: \tab get information from individual slots in a TNA object.\cr
\link[RTN:tna.plot.gsea1]{tna.plot.gsea1}: \tab plot results from the one-tailed GSEA.\cr
\link[RTN:tna.plot.gsea2]{tna.plot.gsea2}: \tab plot results from the two-tailed GSEA.\cr
\link[RTN:AVS-class]{AVS-class}: \tab an S4 class to do enrichment analyses in associated variant sets (AVSs).\cr
\link[RTN:avs.vse]{avs.vse}: \tab variant set enrichment analysis.\cr
\link[RTN:avs.evse]{avs.evse}: \tab an eQTL/VSE pipeline for variant set enrichment analysis.\cr
\link[RTN:avs.pevse]{avs.pevse}: \tab an EVSE pipeline using precomputed eQTLs.\cr
\link[RTN:avs.get]{avs.get}: \tab get information from individual slots in an AVS object.\cr
\link[RTN:avs.plot1]{avs.plot1}: \tab plot results from AVS methods, single plots.\cr
\link[RTN:avs.plot2]{avs.plot2}: \tab plot results from AVS methods, multiple plots.\cr
}
                   
Further information is available in the vignettes by typing \code{vignette("RTN")}. Documented 
topics are also available by typing \code{help.start()} and selecting the RTN package 
from the menu.

}

\author{
Maintainer: Mauro Castro <mauro.a.castro@gmail.com>
}

\references{
Fletcher M.N.C. et al., \emph{Master regulators of FGFR2 signalling and breast cancer risk.} Nature Communications, 4:2464, 2013.

Castro M.A.A. et al., \emph{Regulators of genetic risk of breast cancer identified by integrative network analysis.} Nature Genetics, 48:12-21, 2016.
}

\keyword{ package }

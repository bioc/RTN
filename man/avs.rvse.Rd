\name{avs.rvse}
\alias{avs.rvse}

\title{An rQTL/VSE pipeline for variant set enrichment analysis.}

\description{
The VSE method (\code{\link[RTN:avs.vse]{avs.vse}}) provides a robust framework 
to cope with the heterogeneous structure of haplotype blocks, and has been 
designed to test enrichment in cistromes and epigenomes. In order to extend 
the variant set enrichment to genes this pipeline implements an additional 
step using regulon quantitative trait loci (rQTLs).
}

\usage{
avs.rvse(object, annotation, regdata, snpdata, glist, 
maxgap=250, minSize=100, pValueCutoff=0.05, 
pAdjustMethod="bonferroni", boxcox=TRUE, verbose=TRUE)
}

\arguments{
\item{object}{
an object of class \code{\link[RTN:AVS-class]{AVS-class}}.
}
\item{annotation}{
a data frame with genomic annotations listing chromosome coordinates to which 
a particular property or function has been attributed. It should include the 
following columns: <CHROM>, <START>, <END> and <ID>. The <ID> column can be 
any genomic identifier, while values in <CHROM> should be listed in 
['chr1', 'chr2', 'chr3' ..., 'chrX']. Both <START> and <END> columns 
correspond to chromosome positions mapped to the human genome assembly 
used to build the AVS object.
}
\item{regdata}{
object of class "matrix", a regulon activity matrix.
}
\item{snpdata}{
either an object of class "matrix" or "ff", a single nucleotide polymorphism 
(SNP) matrix.
}
\item{glist}{
a list with character vectors mapped to the 'annotation' data via 
<ID> column. This list is used to run a batch mode for gene sets 
and regulons.
}
\item{maxgap}{
a single integer value specifying the max distant (kb) between the AVS and the 
annotation used to compute the eQTL analysis.
}
\item{minSize}{
if 'glist' is provided, this argument is a single integer or numeric value 
specifying the minimum number of elements for each gene set in the 'glist'. 
Gene sets with fewer than this number are removed from the analysis.
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values 
considered significant.
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be 
used (see 'p.adjust' for details).
}
\item{boxcox}{
a single logical value specifying to use Box-Cox procedure to find a 
transformation of the null that approaches normality (when boxcox=TRUE) or 
not (when boxcox=FALSE). See \code{\link[car:powerTransform]{powerTransform}} 
and \code{\link[car:bcPower]{bcPower}}.
}
\item{verbose}{
a single logical value specifying to display detailed messages 
(when verbose=TRUE) or not (when verbose=FALSE).
}
}

\value{
a data frame in the slot "results", see 'what' options in 
\code{\link[RTN:avs.get]{avs.get}}.
}

\author{Mauro Castro}

\seealso{
  \code{\link[RTN:AVS-class]{AVS-class}}
}

\examples{

\dontrun{
# This example requires the RTNdata package! (currently available under request)
library(RTNdata.LDHapMapRel27.hg18)
library(Fletcher2013b)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)

##################################################
### Build AVS and random AVSs (mapped to hg18)
##################################################

#--- step 1: load 'risk SNPs' data (e.g. BCa risk SNPs from the GWAS catalog)
data(bcarisk, package = "RTNdata.LDHapMapRel27.hg18")

#--- step 2: build an AVS and 1000 matched random AVSs for the input 'risk SNPs'
bcavs <- avs.preprocess.LDHapMapRel27.hg18(bcarisk, nrand=1000)

##################################################
### Example of RVSE analysis
##################################################

#--- step 1: load a precomputed AVS (same 'bcavs' object as above!)
data(bcavs, package="RTNdata.LDHapMapRel27.hg18")

#--- step 2: load genomic annotation for all genes
genemap <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg18.knownGene))
genemap <- genemap[,c("seqnames","start","end","gene_id")]
colnames(genemap) <- c("CHROM","START","END","ID")

#--- step 3: load a TNI object and get a list with regulons
#--- (ids should be the same as in 'genemap')
data("rtni1st")
glist <- tni.get(rtni1st,what="refregulons",idkey="ENTREZ")
glist <- glist[ c("FOXA1","GATA3","ESR1") ] #reduce the list for demonstration!

#--- step 4: compute single-sample regulon activity
rtni1st <- tni.gsea2(rtni1st, regulatoryElements = c("FOXA1","GATA3","ESR1"))
regdata <- tni.get(rtni1st, what = "regulonActivity")$diff

#--- step 5: load a variation dataset matched with the regulon activity dataset
#--- here we use a "toy" dataset for demonstration purposes only.
data(toy_snpdata, package="RTNdata.LDHapMapRel27")

#--- step 6: check "snpdata"" and "regdata"" alignment (samples on cols)
regdata <- t(regdata)
all(colnames(toy_snpdata)==colnames(regdata))
#TRUE

#--- step 7: run the avs.rvse pipeline
bcavs <- avs.rvse(bcavs, 
                  annotation=genemap,
                  regdata=regdata, 
                  snpdata=toy_snpdata, 
                  glist=glist,
                  pValueCutoff=0.01)

#--- step 8: generate the RVSE plot
avs.plot2(bcavs, "rvse", height=2.5)

}
}

\keyword{methods}


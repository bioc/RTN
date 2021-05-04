\name{avs.pevse}
\alias{avs.pevse}

\title{An EVSE pipeline using precomputed eQTLs.}

\description{
The VSE method (\code{\link[RTN:avs.vse]{avs.vse}}) provides a robust framework 
to cope with the heterogeneous structure of haplotype blocks, and has been 
designed to test enrichment in cistromes and epigenomes. In order to extend the 
variant set enrichment to genes this pipeline implements an additional step 
using precomputed expression quantitative trait loci (eQTLs).
}

\usage{
avs.pevse(object, annotation, eqtls, glist, maxgap=250, minSize=100, 
pValueCutoff=0.05, pAdjustMethod="bonferroni", boxcox=TRUE,
verbose=TRUE)
}

\arguments{
\item{object}{
an object of class 'AVS' (see \code{\link[RTN:AVS-class]{AVS-class}}).
}
\item{annotation}{
a data frame with genomic annotations listing chromosome coordinates to which a 
particular property or function has been attributed. It should include the 
following columns: <CHROM>, <START>, <END> and <ID>. The <ID> column can be any 
genomic identifier, while values in <CHROM> should be listed in 
['chr1', 'chr2', 'chr3' ..., 'chrX']. Both <START> and <END> columns correspond 
to chromosome positions mapped to the human genome assembly used to build the 
AVS object.
}
\item{eqtls}{
object of class 'data.frame' with at least two columns, including the following 
column names: <RSID> and <GENEID>. 
}
\item{glist}{
a list with character vectors mapped to the 'annotation' data via 
<ID> column. This list is used to run a batch mode for gene sets and 
regulons.
}
\item{maxgap}{
a single integer value specifying the max distant (kb) used to assign the eQTLs 
in the 'eqtls' object.
}
\item{minSize}{
if 'glist' is provided, this argument is a single integer or numeric value 
specifying the minimum number of elements for each gene set in the 'glist'. 
Gene sets with fewer than this number are removed from the analysis. 
if 'fineMapping=FALSE', an alternative min size value can be provided as a 
vector of the form c(minSize1, minSize2) used to space the null distributions 
(see 'fineMapping').
}
\item{pValueCutoff}{
a single numeric value specifying the cutoff for p-values considered 
significant.
}
\item{pAdjustMethod}{
a single character value specifying the p-value adjustment method to be used 
(see 'p.adjust' for details).
}
\item{boxcox}{
a single logical value specifying to use Box-Cox procedure to find a 
transformation of the null that approaches normality (when boxcox=TRUE) or not 
(when boxcox=FALSE). See \code{\link[car:powerTransform]{powerTransform}} 
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

\author{Mauro Castro, Steve Booth}

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
### Example of EVSE analysis for sets of genomic 
### annotations (e.g. regulons, gene sets, etc.)
##################################################

#--- step 1: load a precomputed AVS
data(bcavs, package="RTNdata.LDHapMapRel27.hg18")

#--- step 2: load genomic annotation for all genes
genemap <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg18.knownGene))
genemap <- genemap[,c("seqnames","start","end","gene_id")]
colnames(genemap) <- c("CHROM","START","END","ID")

#--- step 3: load a TNI object (or any other source of regulons)
#--- and prepare a gene set list.
#--- Note: gene ids should be the same as in the 'genemap' object.
data("rtni1st")
glist <- tni.get(rtni1st,what="refregulons",idkey="ENTREZ")
glist <- glist[ c("FOXA1","GATA3","ESR1") ] #reduce the list for demonstration!

#--- step 4: load precomputed eQTLs
#--- Please note that the input data should represent eQTLs from genome-wide
#--- calls, that is, the universe size should cover 'all SNPs' vs. 'all genes'.
#--- The correct represetation of universe size is essential to build the
#--- null distributions. Or, to put it another way, the eQTL analysis 
#--- should unbiasedly test linked and random markers from the AVS.
#--- In this example it represents 2029 + 966240 SNPs:

lkMarkers <- avs.get(bcavs,what="linkedMarkers")
length(lkMarkers) # i.e. 2029 risk associated and linked SNPs 

rdMarkers <- avs.get(bcavs,what="randomMarkers")
length(rdMarkers) # i.e. 966240 random SNPs

#--- Now we prepare a 'toy' dataset for demonstration purposes only
#--- by picking (naively) SNPs within 250 kb window around the 
#--- genomic annotation.

## load HapMap SNPs (release 27) mapped to hg18
data("popsnp2")

## map SNPs to the genomic annotation
query <- with(popsnp2, 
              GRanges(chrom, IRanges(position, position), 
              id=rsid)
              )
subject <- with(genemap, 
                GRanges(CHROM, IRanges(START, END), 
                id=ID)
                )
hits <- findOverlaps(query,subject, maxgap = 250000)

## reduce 'hits' just for demonstration
hits <- hits[sort(sample(length(hits), 50000))]

## build a 'toy_eqtls' data frame
toy_eqtls <- data.frame(rsid = popsnp2$rsid[from(hits)], 
                        geneid = genemap$ID[to(hits)])

#--- step 5: run the 'avs.pevse' pipeline
#--- important: set 'maxgap' to the same searching window used 
#--- in the dQTL analysis (e.g. 250kb)
bcavs <- avs.pevse(bcavs, annotation=genemap, glist=glist, 
                   eqtls=toy_eqtls, maxgap = 250)

#--- step 6: generate the pEVSE plots
avs.plot2(bcavs,"pevse",height=2.5)


#####################################################
#--- parallel version for 'step 5' with SNOW package!
# library(snow)
# options(cluster=snow::makeCluster(3, "SOCK"))
# bcavs <- avs.pevse(bcavs, annotation=genemap, glist=glist, 
#          eqtls=toy_eqtls, maxgap = 250)
# stopCluster(getOption("cluster"))

## ps. as a technical note, the parallel version uses a
## slightly different overlap-based operation, which might
## bring slightly different counts depending on the data 
## input organization

}
}

\keyword{methods}


##-------------------------------------------------------------------------
setGeneric("tni.preprocess",
           function(object, rowAnnotation=NULL, colAnnotation=NULL, 
                    cvfilter=FALSE, verbose=TRUE)
             standardGeneric("tni.preprocess"), package="RTN")
setGeneric("tni.permutation",
           function(object, pValueCutoff=0.01, pAdjustMethod="BH", 
                    globalAdjustment=TRUE, 
                    estimator="spearman", nPermutations=1000, 
                    pooledNullDistribution=TRUE, 
                    boxcox=TRUE, parChunks=NULL, verbose=TRUE) 
             standardGeneric("tni.permutation"), package="RTN")
setGeneric("tni.bootstrap",
           function(object, nBootstraps=100, consensus=95, 
                    parChunks=NULL, verbose=TRUE)
             standardGeneric("tni.bootstrap"), package="RTN")
setGeneric("tni.dpi.filter",
           function(object, eps=0, sizeThreshold=TRUE, minRegulonSize=15, 
                    verbose=TRUE)
             standardGeneric("tni.dpi.filter"), package="RTN")
setGeneric("tni.get",
           function(object, what="summary", order=TRUE, ntop=NULL, 
                    reportNames=TRUE, 
                    idkey=NULL) 
             standardGeneric("tni.get"), package="RTN")
##-------------------------------------------------------------------------
setGeneric("tni.conditional",
           function(object, modulators, tfs=NULL, sampling=35, 
                    pValueCutoff=0.01, 
                    pAdjustMethod="bonferroni", minRegulonSize=15, 
                    minIntersectSize=5, 
                    miThreshold="md", prob=0.99, medianEffect=FALSE, 
                    iConstraint=TRUE, verbose=TRUE, ...)
             standardGeneric("tni.conditional"), package="RTN")
setGeneric("tni.gsea2",
           function(object, minRegulonSize=15, sizeFilterMethod="posORneg", 
                    scale=FALSE, exponent=1, tnet="dpi",
                    regulatoryElements=NULL, features=NULL, 
                    samples=NULL, refsamp=samples, log=TRUE,
                    alternative=c("two.sided", "less", "greater"),
                    targetContribution=FALSE, additionalData=FALSE, 
                    verbose=TRUE, doSizeFilter=NULL)
             standardGeneric("tni.gsea2"), package="RTN")
setGeneric("tni.area3",
           function(object, minRegulonSize=15, sizeFilterMethod="posORneg", 
                    scale=FALSE, tnet="dpi", regulatoryElements=NULL, 
                    samples=NULL, features=NULL, 
                    refsamp=NULL, log=FALSE, verbose=TRUE, doSizeFilter=NULL)
             standardGeneric("tni.area3"), package="RTN")
setGeneric("tni.graph",
           function(object, tnet="dpi", gtype="rmap", 
                    minRegulonSize=15, regulatoryElements=NULL, 
                    amapFilter="quantile", amapCutoff=NULL, ntop=NULL, ...)
             standardGeneric("tni.graph"), package="RTN")
setGeneric("tni.regulon.summary",
           function(object, regulatoryElements = NULL, verbose = TRUE)
               standardGeneric("tni.regulon.summary"), package="RTN")
##-------------------------------------------------------------------------
setGeneric("tni.overlap.genesets",
           function(object, geneSetList, regulatoryElements = NULL, 
                    minSetSize = 15, sizeFilterMethod="posORneg",
                    method = c("HT","JC"), pValueCutoff = 0.05, 
                    pAdjustMethod = "BH", verbose = TRUE)
             standardGeneric("tni.overlap.genesets"), package="RTN")
setGeneric("tni.annotate.regulons",
           function(object, geneSetList, sampleSetList = NULL, 
                    regulatoryElements = NULL, minSetSize = 15, 
                    sizeFilterMethod="posORneg", exponent = 1, verbose = TRUE)
             standardGeneric("tni.annotate.regulons"), package="RTN")
setGeneric("tni.annotate.samples",
           function(object, geneSetList, minSetSize = 15, exponent = 1, 
                    samples=NULL, verbose = TRUE)
             standardGeneric("tni.annotate.samples"), package="RTN")
##-------------------------------------------------------------------------
setGeneric("tni.sre",
           function(object, sampleGroups, regulatoryElements = NULL, 
                    pValueCutoff = 0.05, pAdjustMethod = "BH")
             standardGeneric("tni.sre"), package = "RTN")
setGeneric("tni.prune",
           function(object, regulatoryElements = NULL, minRegCor = 0.95, 
                    tarPriorityMethod = "EC", minPrunedSize = 30, 
                    verbose = TRUE, ...) 
               standardGeneric("tni.prune"), package="RTN")
setGeneric("tni.replace.samples",
           function(object, expData, rowAnnotation=NULL, colAnnotation=NULL,
                    removeRegNotAnnotated=TRUE, verbose=TRUE)
             standardGeneric("tni.replace.samples"), package="RTN")
setGeneric("tni2tna.preprocess",
           function(object, phenotype=NULL, hits=NULL, phenoIDs=NULL, 
                    duplicateRemoverMethod="max", verbose=TRUE)
             standardGeneric("tni2tna.preprocess"), package="RTN")
##-------------------------------------------------------------------------
setGeneric("tna.mra",
           function(object, pValueCutoff=0.05, pAdjustMethod="BH", 
                    minRegulonSize=15, tnet="dpi", tfs=NULL, verbose=TRUE) 
             standardGeneric("tna.mra"), package="RTN")
setGeneric("tna.gsea1",
           function(object, pValueCutoff=0.05, pAdjustMethod="BH",  
                    minRegulonSize=15, sizeFilterMethod="posORneg",
                    nPermutations=1000, exponent=1, 
                    tnet="dpi", signature=c("phenotype","hits"),
                    orderAbsValue=TRUE, tfs=NULL, 
                    verbose=TRUE) 
             standardGeneric("tna.gsea1"), package="RTN")
setGeneric("tna.gsea2",
           function(object, pValueCutoff=0.05, pAdjustMethod="BH",  
                    minRegulonSize=15, sizeFilterMethod="posORneg", 
                    nPermutations=1000, exponent=1, tnet="dpi", 
                    signature=c("phenotype","hits"), tfs=NULL,  
                    verbose=TRUE, doSizeFilter=NULL) 
             standardGeneric("tna.gsea2"), package="RTN")
setGeneric("tna.get",
           function(object, what="summary", order=TRUE, ntop=NULL, 
                    reportNames=TRUE, idkey=NULL) 
             standardGeneric("tna.get"), package="RTN")
##-------------------------------------------------------------------------
setGeneric("avs.vse",
           function(object, annotation, glist=NULL, maxgap=0, minSize=100, 
                    pValueCutoff=0.05, pAdjustMethod="bonferroni", 
                    boxcox=TRUE, verbose=TRUE)
             standardGeneric("avs.vse"), package="RTN")
setGeneric("avs.evse",
           function(object, annotation, gxdata, snpdata, glist=NULL, 
                    maxgap=250, minSize=100, pValueCutoff=0.05, 
                    pAdjustMethod="bonferroni", 
                    boxcox=TRUE, fineMapping=TRUE, verbose=TRUE)
             standardGeneric("avs.evse"), package="RTN")
setGeneric("avs.pevse",
           function(object, annotation, eqtls, glist, maxgap=250, minSize=100, 
                    pValueCutoff=0.05, pAdjustMethod="bonferroni", boxcox=TRUE,
                    verbose=TRUE)
             standardGeneric("avs.pevse"), package="RTN")
setGeneric("avs.rvse",
           function(object, annotation, regdata, snpdata, glist, 
                    maxgap=250, minSize=100, pValueCutoff=0.05, 
                    pAdjustMethod="bonferroni", 
                    boxcox=TRUE, verbose=TRUE)
             standardGeneric("avs.rvse"), package="RTN")
setGeneric("avs.get",
           function(object, what="summary", report=FALSE, pValueCutoff=NULL) 
             standardGeneric("avs.get"), package="RTN")


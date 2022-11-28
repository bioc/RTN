################################################################################
##########################    TNI Class Methods     ############################
################################################################################

##------------------------------------------------------------------------------
## Constructor of TNI Class objects
## Entry point for all TNI/TNA pipelines, including pre-processing
tni.constructor <- function(expData, regulatoryElements, rowAnnotation=NULL, 
                            colAnnotation=NULL, cvfilter=FALSE, verbose=TRUE){
  if(missing(expData))
    stop("'expData' is missing!",call.=FALSE)    
  if(missing(regulatoryElements))
    stop("'regulatoryElements' is missing!",call.=FALSE) 
  if(is(expData,"SummarizedExperiment") || 
     is(expData,"RangedSummarizedExperiment")){
    if( !is.null(rowAnnotation) || !is.null(colAnnotation) ){
      tp1 <- "NOTE: when using a 'SummarizedExperiment' container, "
      tp2 <- "row and col annotations are used from that object!"
      warning(tp1,tp2,call.=FALSE)
    }
    if (length(assays(expData)) > 1) {
      stop("please input a SummarizedExperiment with only one assay")
    }
    rowAnnotation <- as.data.frame(rowData(expData))
    colAnnotation <- as.data.frame(colData(expData))
    expData <- assays(expData)[[1]]
  }
  object <- new("TNI", expData=expData, regulatoryElements=regulatoryElements)
  object <- tni.preprocess(object, rowAnnotation, colAnnotation, 
                           cvfilter, verbose)
  return(object)
}

##------------------------------------------------------------------------------
## Estimate 'alphaB' for 'nB', given 'nA', 'alphaA' and 'betaA'
tni.alpha.adjust <- function(nB, nA, alphaA, betaA = 0.2){
  if(missing(nB))
    stop("'nB' is missing!",call.=FALSE) 
  if(missing(nA))
    stop("'nA' is missing!",call.=FALSE) 
  if(nA<nB)
    stop("'nA' must be greater than or equal to 'nB'!",call.=FALSE)
  if(missing(alphaA))
    stop("'alphaA' is missing!",call.=FALSE)
  tnai.checks(name="nB", para=nB)
  tnai.checks(name="nA", para=nA)
  tnai.checks(name="alphaA", para=alphaA)
  tnai.checks(name="betaA", para=betaA)
  if(nA==nB){
    alphaB <- alphaA
  } else {
    alphaB <- .alpha.adjust(nB, nA, alphaA, betaA)
  }
  return(alphaB)
}

##------------------------------------------------------------------------------
##Main initialization method
setMethod("initialize",
          "TNI",
          function(.Object, expData, regulatoryElements) {
            ##-----initialization
            .Object@gexp <- expData
            .Object@regulatoryElements <- regulatoryElements
            .Object@modulators <- character()
            ##-----result slot
            .Object@results<-list()
            ##-----status matrix
            .Object@status <- rep("[ ]", 1, 6)
            names(.Object@status) <- c("Preprocess", "Permutation", 
                                       "Bootstrap", "DPI.filter", 
                                       "Conditional","Activity")
            ##-----summary info
            ##-----regulatoryElements
            sum.info.regElements<-matrix(,1,2)
            rownames(sum.info.regElements)<-"regulatoryElements"
            colnames(sum.info.regElements)<-c("input","valid")   
            ##-----targetElements
            sum.info.targetElements<-matrix(,1,2)
            rownames(sum.info.targetElements)<-"targetElements"
            colnames(sum.info.targetElements)<-c("input","valid")  
            ##-----parameters
            sum.info.para <- list()
            sum.info.para$perm<-matrix(,1,7)
            colnames(sum.info.para$perm)<-c("pValueCutoff","pAdjustMethod", 
                                            "globalAdjustment", "estimator", 
                                            "nPermutations",
                                            "pooledNullDistribution", "boxcox")
            rownames(sum.info.para$perm)<-"Parameter"
            sum.info.para$boot<-matrix(,1,4)
            colnames(sum.info.para$boot)<-c("estimator", "nBootstraps", 
                                            "consensus", "boxcox")        
            rownames(sum.info.para$boot)<-"Parameter"
            sum.info.para$dpi<-matrix(,1,1)
            colnames(sum.info.para$dpi)<-c("eps")       
            rownames(sum.info.para$dpi)<-"Parameter"
            sum.info.para$cdt<-matrix(,1,7)
            colnames(sum.info.para$cdt)<-c("sampling","pValueCutoff",
                                           "pAdjustMethod","minRegulonSize",
                                           "minIntersectSize","miThreshold",
                                           "prob")
            rownames(sum.info.para$cdt)<-"Parameter"
            ##-----results
            sum.info.results<-list()
            sum.info.results$tnet<-matrix(,2,3)
            colnames(sum.info.results$tnet)<-c("regulatoryElements","Targets",
                                               "Edges")
            rownames(sum.info.results$tnet)<-c("tnet.ref","tnet.dpi")
            .Object@summary<-list(regulatoryElements=sum.info.regElements,
                                  targetElements=sum.info.targetElements,
                                  para=sum.info.para,results=sum.info.results)			
            .Object
          }
)

##------------------------------------------------------------------------------
setMethod(
  "tni.preprocess",
  "TNI",
  function(object, rowAnnotation=NULL, colAnnotation=NULL, 
           cvfilter=FALSE, verbose=TRUE){
    
    ##----check compatibility
    object <- upgradeTNI(object)
    
    ##----start preprocessing
    if(verbose)cat("-Preprocessing for input data...\n")
    
    ##----check arguments
    tnai.checks(name="cvfilter",para=cvfilter)
    tnai.checks(name="verbose",para=verbose)
    
    ##---- Count initial input
    object@summary$regulatoryElements[,"input"] <- length(object@regulatoryElements)
    object@summary$targetElements[,"input"] <- nrow(object@gexp)
    
    ##----initial checks for the main arguments
    tnai.checks(name="regulatoryElements",para=object@regulatoryElements)
    tmp <- .expDataChecks(object@gexp, rowAnnotation, colAnnotation)
    object@gexp <- tmp$expData
    rowAnnotation <- tmp$rowAnnotation
    colAnnotation <- tmp$colAnnotation
    
    ##----update rowAnnotation
    if(cvfilter && ncol(rowAnnotation)>1){
      #apply cvfilter if rowAnnotation is available
      if(verbose)
        cat("--Removing duplicated genes (keep max coefficient of variation!)...\n")
      #e.g. col1=id, col2=symbol (collapse cv by col2)
      cvres <- cv.filter(object@gexp, rowAnnotation)
      object@gexp <- cvres$gexp
      object@rowAnnotation <- cvres$ids
    } else {
      object@rowAnnotation <- rowAnnotation[rownames(object@gexp),,drop=FALSE]
    }
    
    ##----add colAnnotation
    object@colAnnotation <- colAnnotation[colnames(object@gexp),,drop=FALSE]
    
    ##----update regulatoryElements
    if(verbose) cat("--Checking 'regulatoryElements' in 'rowAnnotation'...\n")
    col1 <- sapply(1:ncol(object@rowAnnotation),function(i){
      sum(object@regulatoryElements%in%object@rowAnnotation[,i],na.rm=TRUE)
    })
    col1 <- which(col1==max(col1))[1]
    idx <- object@rowAnnotation[[col1]] %in% object@regulatoryElements
    object@regulatoryElements <- rownames(object@rowAnnotation)[idx]
    
    ##----check sd in expData
    if(verbose) cat("--Checking 'expData'...\n")
    sd.check <- apply(object@gexp,1,sd)
    if(any(is.na(sd.check))){
      stop("unpredicted exception found in the input data matrix! 
           ...a possible cause is the presence of 'Inf' values. ")
    }
    sd.check <- sd.check==0
    if(any(sd.check)){
      if(verbose)cat("--Removing inconsistent data: standard deviation is zero for", 
                     sum(sd.check),"gene(s)! \n")
      object@gexp <- object@gexp[!sd.check,,drop=FALSE]
      object@rowAnnotation <- object@rowAnnotation[!sd.check,,drop=FALSE]
      object@regulatoryElements <- object@regulatoryElements[
        object@regulatoryElements %in% rownames(object@rowAnnotation)]
    }
    
    ##----add targetElements
    object@targetElements <- rownames(object@rowAnnotation)
    object@summary$targetElements[,"valid"] <- nrow(object@gexp)
    
    ##----final checks for regulatoryElements
    if(length(object@regulatoryElements)==0)
      stop("input 'regulatoryElements' contains no useful data!\n")
    object@summary$regulatoryElements[,"valid"] <- length(object@regulatoryElements)
    #..make sure regulatoryElements is named
    if(!is.null(object@rowAnnotation$SYMBOL)){
      names(object@regulatoryElements)<-object@rowAnnotation[
        object@regulatoryElements,"SYMBOL"]
    } else {
      names(object@regulatoryElements)<-object@regulatoryElements
    }
    #..remove any empty space or NA from 'regulatoryElements' names
    rnames<-names(object@regulatoryElements)
    idx<-rnames==""|rnames=="NA"
    names(object@regulatoryElements)[idx]<-object@regulatoryElements[idx]
    #..sort regulatoryElements by names
    idx <- sort.list(names(object@regulatoryElements))
    object@regulatoryElements <- object@regulatoryElements[idx]
    
    ##-----updade status and return
    object@status["Preprocess"] <- "[x]"
    object@status["Permutation"] <- "[ ]"
    object@status["Bootstrap"] <- "[ ]"
    object@status["DPI.filter"] <- "[ ]"
    if(verbose)cat("-Preprocessing complete!\n\n")
    return(object)
  }
)
##------------------------------------------------------------------------------
##infer MI network
setMethod(
  "tni.permutation",
  "TNI",
  function(object, pValueCutoff=0.01, pAdjustMethod="BH", 
           globalAdjustment=TRUE, estimator="spearman", nPermutations=1000, 
           pooledNullDistribution=TRUE, boxcox=TRUE, parChunks=NULL, 
           verbose=TRUE){
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    if(object@status["Preprocess"]!="[x]")
      stop("input 'object' needs preprocessing!")
    
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="globalAdjustment",para=globalAdjustment)
    tnai.checks(name="estimator",para=estimator)  
    tnai.checks(name="nPermutations",para=nPermutations)
    tnai.checks(name="pooledNullDistribution",para=pooledNullDistribution)
    tnai.checks(name="boxcox",para=boxcox)
    tnai.checks(name="parChunks",para=parChunks)
    tnai.checks(name="verbose",para=verbose)
    object@para$perm<-list(pValueCutoff=pValueCutoff,
                           pAdjustMethod=pAdjustMethod,
                           globalAdjustment=globalAdjustment, 
                           estimator=estimator,nPermutations=nPermutations,
                           pooledNullDistribution=pooledNullDistribution,
                           boxcox=boxcox)
    object@summary$para$perm[1,]<-unlist(object@para$perm)
    ###compute reference network###
    ##---permutation analysis
    if(object@para$perm$pooledNullDistribution){
      res<-tni.perm.pooled(object, parChunks, verbose)
    } else {
      res<-tni.perm.separate(object,verbose)
    }
    # object@results$mipval <- res$mipval
    # object@results$miadjpv <- res$miadjpv
    object@results$tn.ref <- res$tn.ref * 
      tni.cor(object@gexp, res$tn.ref, estimator=object@para$perm$estimator)
    object@status["Permutation"] <- "[x]"
    if(verbose)cat("-Permutation analysis complete! \n\n")
    ##update summary and return results
    bin<-object@results$tn.ref
    bin[bin!=0]<-1
    object@summary$results$tnet[1,]<-c(ncol(bin),sum(rowSums(bin)>0),sum(bin))
    return(object)
  }
)

##------------------------------------------------------------------------------
##infer MI network
setMethod(
  "tni.bootstrap",
  "TNI",
  function(object, nBootstraps=100, consensus=95, 
           parChunks=NULL, verbose=TRUE){
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    if(object@status["Preprocess"]!="[x]")
      stop("input 'object' needs preprocessing and permutation analysis!")
    if(object@status["Permutation"]!="[x]")
      stop("input 'object' needs permutation analysis!")
    
    #----check parameters
    tnai.checks(name="nBootstraps",para=nBootstraps)    
    tnai.checks(name="consensus",para=consensus)
    tnai.checks(name="parChunks",para=parChunks)
    tnai.checks(name="verbose",para=verbose)
    
    #--- assign same estimator used in the permutation step
    estimator <- tni.get(object, "para")$perm$estimator
    boxcox <- tni.get(object, "para")$perm$boxcox
    object@para$boot<-list(estimator=estimator, nBootstraps=nBootstraps, 
                           consensus=consensus, boxcox=boxcox)
    object@summary$para$boot[1,]<-unlist(object@para$boot)
    
    #--- run bootstrap analysis
    object@results$tn.ref<-tni.boot(object,parChunks,verbose)
    object@status["Bootstrap"] <- "[x]"
    if(verbose)cat("-Bootstrap analysis complete! \n\n")
    
    #--- update summary and return results
    bin<-object@results$tn.ref
    bin[bin!=0]<-1
    object@summary$results$tnet[1,]<-c(ncol(bin),sum(rowSums(bin)>0),sum(bin))
    return(object)
  }
)

##------------------------------------------------------------------------------
##infer MI network
setMethod(
  "tni.dpi.filter",
  "TNI",
  function(object, eps=0, sizeThreshold=TRUE, minRegulonSize=15, verbose=TRUE){
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    if(object@status["Permutation"]!="[x]")
      stop("input 'object' needs permutation/bootstrep analysis!")
    
    ##---check and assign parameters
    tnai.checks(name="eps",para=eps)
    tnai.checks(name="sizeThreshold",para=sizeThreshold)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="verbose",para=verbose)
    
    ##---if not provided, estimate eps from tn.ref
    if(is.na(eps)){
      eps <- abs(object@results$tn.ref)
      eps <- min(eps[eps!=0])/2
      sizeThreshold=FALSE
    }
    object@para$dpi <- list(eps=eps)
    object@summary$para$dpi[1,]<-unlist(object@para$dpi)
    
    ##---apply dpi filter
    if(verbose)cat("-Applying dpi filter...\n")
    object@results$tn.dpi <- tni.dpi(abs(object@results$tn.ref), 
                                     eps=object@para$dpi$eps)
    object@results$tn.dpi <- object@results$tn.dpi * 
      tni.cor(object@gexp,object@results$tn.dpi, 
              estimator=object@para$perm$estimator)
    
    ##---apply sizeThreshold on small/unbalanced regulons
    if(sizeThreshold){
      rgcounts <- tni.get(object, what="regulonSize")
      rgcounts <- rgcounts[,c("Positive","Negative")]>minRegulonSize
      regs <- which(rowSums(rgcounts)==1)
      regs <- rownames(rgcounts)[regs]
      epsz <- abs(object@results$tn.ref)
      epsz <- min(epsz[epsz!=0])/2
      if(length(regs)>0 & epsz>eps){
        tnet1 <- object@results$tn.dpi
        tnet2 <- tni.dpi(abs(object@results$tn.ref), eps=epsz)
        tnet2 <- tnet2 * tni.cor(object@gexp[object@targetElements,], 
                                 tnet2, estimator=object@para$perm$estimator)
        for(reg in regs){
          idx1 <- which(!rgcounts[reg,])
          if(idx1==1){
            tp1 <- tnet1[,reg]
            tp2 <- tnet2[,reg];tp2[tp2<=0] <- NA
            tp2 <- sort(tp2, decreasing = TRUE, na.last=NA)
            if(length(tp2)>0){
              tp2 <- tp2[1:min(minRegulonSize,length(tp2))]
              tp1[names(tp2)] <- tp2
              tnet1[,reg] <- tp1
            }
          } else {
            tp1 <- tnet1[,reg]
            tp2 <- tnet2[,reg];tp2[tp2>=0] <- NA
            tp2 <- sort(tp2, decreasing = FALSE, na.last=NA)
            if(length(tp2)>0){
              tp2 <- tp2[1:min(minRegulonSize,length(tp2))]
              tp1[names(tp2)] <- tp2
              tnet1[,reg] <- tp1
            }
          }
        }
        object@results$tn.dpi <- tnet1
      }
    }
    
    ##update and return results
    bin<-object@results$tn.dpi
    bin[bin!=0]<-1
    object@summary$results$tnet[2,]<-c(ncol(bin),sum(rowSums(bin)>0),sum(bin))
    if(verbose)cat("-DPI filter complete! \n\n")
    object@status["DPI.filter"] <- "[x]"
    
    return(object)
  }
)

##------------------------------------------------------------------------------
##GSEA2 for TNI
setMethod(
  "tni.gsea2",
  "TNI",function(object, minRegulonSize=15, sizeFilterMethod="posORneg", 
                 scale=FALSE, exponent=1, tnet="dpi", regulatoryElements=NULL, 
                 features=NULL, samples=NULL, refsamp=samples, log=TRUE, 
                 alternative=c("two.sided", "less", "greater"), 
                 targetContribution=FALSE, additionalData=FALSE, verbose=TRUE, 
                 doSizeFilter=NULL){
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    if(object@status["Preprocess"]!="[x]")
      stop("TNI object requires preprocessing!")
    if(object@status["Permutation"]!="[x]")
      stop("TNI object requires permutation/bootstrap and DPI filter!")  
    if(object@status["DPI.filter"]!="[x]")
      stop("TNI object requires DPI filter!")
    
    ##-----check and assign parameters
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="scale",para=scale)
    tnai.checks(name="sizeFilterMethod",para=sizeFilterMethod)
    tnai.checks(name="exponent",para=exponent)
    tnai.checks(name="gsea.tnet",para=tnet)
    tnai.checks(name="regulatoryElements",para=regulatoryElements)
    tnai.checks(name="samples",para=samples)
    tnai.checks(name="features",para=features)
    tnai.checks(name="refsamp",para=refsamp)
    tnai.checks(name="log",para=log) 
    alternative <- match.arg(alternative)
    tnai.checks(name="targetContribution",para=targetContribution)
    tnai.checks(name="additionalData",para=additionalData)
    tnai.checks(name="verbose",para=verbose) 
    object@para$gsea2<-list(minRegulonSize=minRegulonSize, exponent=exponent,
                            tnet=tnet, sizeFilterMethod=sizeFilterMethod, 
                            alternative=alternative, scale=scale, log=log)
    
    if(!is.null(doSizeFilter)){
      warning("'doSizeFilter' is deprecated, please use the 'sizeFilterMethod' parameter.")
      tnai.checks(name="doSizeFilter",para=doSizeFilter)
      if(doSizeFilter){
        sizeFilterMethod="posANDneg"
      } else {
        sizeFilterMethod="posORneg"
      }
    }
    
    ##------ get gexp
    gexp <- object@gexp[object@targetElements,,drop=FALSE]
    if(log){
      if(verbose) cat("-Checking log space... ")
      if(.isUnloggedData(gexp)){
        if(verbose) cat("applying log2 transformation!\n")
        gexp <- .log2transform(gexp)
      } else {
        if(verbose)cat("OK!\n")
      }
    } else {
      if(.isUnloggedData(gexp)){
        tp1 <- "The 'tni.gsea2' expects expression values in log space.\n"
        tp2 <- "Please, either set 'log = TRUE' or adjust the expression data\n"
        tp3 <- "available in the TNI object."
        warning(tp1, tp2, tp3)
      }
    }
    if(scale){
      if(verbose) cat("-Applying 'scale' option...\n")
      gexp <- t(scale(t(gexp)))
    }
    
    ##------ compute reference gx vec
    if(is.null(refsamp)){
      gxref <- apply(gexp,1,mean)
    } else {
      idx <- refsamp %in% colnames(gexp)
      if(!all(idx)){
        stop("'refsamp' should list only valid names!")
      }
      gxref <- apply(gexp[,refsamp],1,mean)
    }
    
    ##----- set samples
    if(!is.null(samples)){
      idx <- samples %in% colnames(gexp)
      if(!all(idx)){
        stop("'samples' should list only valid names!")
      }
      samples<-colnames(gexp)[colnames(gexp) %in% samples]
    } else {
      samples<-colnames(gexp)
    }
    
    ##----- set features
    if(!is.null(features)){
      rowAnnotation <- object@rowAnnotation[object@targetElements,,drop=FALSE]
      col1<-sapply(1:ncol(rowAnnotation),function(i){
        sum(features%in%rowAnnotation[,i],na.rm=TRUE)
      })
      col1<-which(col1==max(col1))[1]
      idx<-rowAnnotation[[col1]]%in%features
      object@results$tn.ref[!idx,]<-0
      object@results$tn.dpi[!idx,]<-0
    }
    
    ##-----get regulons
    if(tnet=="ref"){
      listOfRegulonsAndMode<-tni.get(object,what="refregulons.and.mode")
    } else {
      listOfRegulonsAndMode<-tni.get(object,what="regulons.and.mode")
    }
    
    ##-----set regs
    if(!is.null(regulatoryElements)){
      if(sum(regulatoryElements%in%object@regulatoryElements) > 
         sum(regulatoryElements%in%names(object@regulatoryElements) ) ){
        regulatoryElements <- object@regulatoryElements[
          object@regulatoryElements%in%regulatoryElements]
      } else {
        regulatoryElements<-object@regulatoryElements[
          names(object@regulatoryElements)%in%regulatoryElements]
      }
      if(length(regulatoryElements)==0)
        stop("'regulatoryElements' argument has no valid names!")
    } else {
      regulatoryElements<-object@regulatoryElements
    }
    listOfRegulonsAndMode<-listOfRegulonsAndMode[regulatoryElements]
    
    ##-----check regulon size
    regcounts <- .regulonCounts(listOfRegulonsAndMode)
    if(sizeFilterMethod=="posANDneg"){
      idx <- regcounts$Positive >= minRegulonSize & 
        regcounts$Negative >= minRegulonSize
    } else if(sizeFilterMethod=="posORneg"){
      idx <- regcounts$Positive >= minRegulonSize | 
        regcounts$Negative >= minRegulonSize
    } else {
      idx <- regcounts$Size >= minRegulonSize
    }
    regulatoryElements <- regulatoryElements[
      regulatoryElements%in%rownames(regcounts)[idx]]
    listOfRegulonsAndMode <- listOfRegulonsAndMode[regulatoryElements]
    
    ##-----stop when no regulon passes the size requirement
    if(length(listOfRegulonsAndMode)==0){
      stop("no regulon passed the 'minRegulonSize' requirement!")
    }
    
    #-----get phenotypes
    phenotypes <- gexp-gxref
    
    #-----reset names to integer values
    listOfRegulons <- lapply(listOfRegulonsAndMode, names)
    for(i in names(listOfRegulonsAndMode)){
      reg <- listOfRegulonsAndMode[[i]]
      names(listOfRegulonsAndMode[[i]]) <- match(names(reg), rownames(phenotypes))
    }
    rnames_phenotypes <- rownames(phenotypes)
    rownames(phenotypes)<-1:nrow(phenotypes)
    
    ##-----get ranked phenotypes
    phenoranks <- apply(-phenotypes, 2, rank)
    colnames(phenoranks) <- colnames(phenotypes)
    rownames(phenoranks) <- rownames(phenotypes)

    #-----run 2t-gsea
    if(isParallel() && length(samples)>1){
      if(verbose)
        cat("-Performing two-tailed GSEA (parallel version - ProgressBar disabled)...\n")
      if(verbose)cat("--For", length(listOfRegulonsAndMode), "regulon(s) and",
                     length(samples),'sample(s)...\n')
      cl<-getOption("cluster")
      snow::clusterExport(cl, list(".run.tni.gsea2.alternative",
                                   ".fgseaScores4TNI"), 
                          envir=environment())
      regulonActivity <- list()
      res <- snow::parLapply(cl, samples, function(samp){
        .run.tni.gsea2.alternative(
          listOfRegulonsAndMode=listOfRegulonsAndMode,
          phenotype=phenotypes[, samp],
          phenorank=phenoranks[, samp],
          exponent=exponent,
          alternative=alternative
        )
      })
      regulonActivity$differential <- t(sapply(res, function(r) r$differential))
      regulonActivity$positive <- t(sapply(res, function(r) r$positive))
      regulonActivity$negative <- t(sapply(res, function(r) r$negative))
    } else {
      if(verbose)cat("-Performing two-tailed GSEA...\n")
      if(verbose)cat("--For", length(listOfRegulonsAndMode), "regulon(s) and",
                     length(samples),'sample(s)...\n')
      if(verbose)pb <- txtProgressBar(style=3)
      regulonActivity<-list()
      for(i in 1:length(samples)){
        res <- .run.tni.gsea2.alternative(
          listOfRegulonsAndMode=listOfRegulonsAndMode,
          phenotype=phenotypes[, samples[i]],
          phenorank=phenoranks[, samples[i]],
          exponent=exponent,
          alternative=alternative
        )
        regulonActivity$differential<-rbind(regulonActivity$differential,
                                            res$differential[regulatoryElements])
        regulonActivity$positive<-rbind(regulonActivity$positive,
                                        res$positive[regulatoryElements])
        regulonActivity$negative<-rbind(regulonActivity$negative,
                                        res$negative[regulatoryElements])
        if(verbose) setTxtProgressBar(pb, i/length(samples))
      }
      if(verbose) close(pb)
    }
    rownames(regulonActivity$differential)<-samples
    rownames(regulonActivity$positive)<-samples
    rownames(regulonActivity$negative)<-samples
    regulonActivity <- .tni.stratification.gsea2(regulonActivity)
    if(targetContribution){
      tc <- .target.contribution(listOfRegulonsAndMode, regulonActivity, 
                                 phenoranks, phenotypes, exponent, 
                                 alternative, verbose)
      regulonActivity$data$listOfTargetContribution <- tc
    }
    if(additionalData || targetContribution){
      regulonActivity$data$listOfRegulons <-listOfRegulons
      regulonActivity$data$listOfRegulonsAndMode <- listOfRegulonsAndMode
      regulonActivity$data$phenoranks <- phenoranks
      regulonActivity$data$phenotypes <- phenotypes
      regulonActivity$data$rnames_phenotypes <- rnames_phenotypes
      regulonActivity$data$exponent <- exponent
      regulonActivity$data$alternative <- alternative
    } else {
      colnames(regulonActivity$differential)<-names(regulatoryElements)
      colnames(regulonActivity$positive)<-names(regulatoryElements)
      colnames(regulonActivity$negative)<-names(regulatoryElements)
      colnames(regulonActivity$status)<-names(regulatoryElements)
      regulonActivity$regulatoryElements <- regulatoryElements
    }
    object@results$regulonActivity <- regulonActivity
    
    #---
    if(verbose)cat("-GSEA2 complete! \n\n")
    object@status["Activity"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##aREA-3T for TNI
setMethod(
  "tni.area3",
  "TNI",function(object, minRegulonSize=15, sizeFilterMethod="posORneg",
                 scale=FALSE, tnet="dpi", regulatoryElements=NULL, 
                 samples=NULL, features=NULL, refsamp=NULL, log=FALSE, 
                 verbose=TRUE, doSizeFilter=NULL){
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    if(object@status["Preprocess"]!="[x]")
      stop("TNI object requires preprocessing!")
    if(object@status["Permutation"]!="[x]")
      stop("TNI object requires permutation/bootstrap and DPI filter!")  
    if(object@status["DPI.filter"]!="[x]")
      stop("TNI object requires DPI filter!")
    
    ##-----check and assign parameters
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="sizeFilterMethod",para=sizeFilterMethod)
    tnai.checks(name="scale",para=scale)
    tnai.checks(name="area.tnet",para=tnet)
    tnai.checks(name="regulatoryElements",para=regulatoryElements)
    tnai.checks(name="samples",para=samples)
    tnai.checks(name="features",para=features)
    tnai.checks(name="refsamp",para=refsamp)
    tnai.checks(name="log",para=log) 
    tnai.checks(name="verbose",para=verbose) 
    object@para$area3 <- list(minRegulonSize=minRegulonSize, 
                              sizeFilterMethod=sizeFilterMethod,
                              scale=scale, tnet=tnet, log=log)
    
    if(!is.null(doSizeFilter)){
      warning("'doSizeFilter' is deprecated, please use the 'sizeFilterMethod' parameter.")
      tnai.checks(name="doSizeFilter",para=doSizeFilter)
      if(doSizeFilter){
        sizeFilterMethod="posANDneg"
      } else {
        sizeFilterMethod="posORneg"
      }
    }
    
    ##------ compute reference gx vec
    gexp <- object@gexp[object@targetElements,,drop=FALSE]
    if(scale) gexp <- t(scale(t(gexp)))
    if(is.null(refsamp)){
      gxref <- apply(gexp,1,mean)
    } else {
      idx <- refsamp %in% colnames(gexp)
      if(!all(idx)){
        stop("'refsamp' should list only valid names!")
      }
      gxref <- apply(gexp[,refsamp],1,mean)
    }
    ##----- set samples
    if(!is.null(samples)){
      idx <- samples %in% colnames(gexp)
      if(!all(idx)){
        stop("'samples' should list only valid names!")
      }
      samples<-colnames(gexp)[colnames(gexp) %in% samples]
    } else {
      samples<-colnames(gexp)
    }
    ##----- set features
    if(!is.null(features)){
      rowAnnotation <- object@rowAnnotation[object@targetElements,,drop=FALSE]
      col1<-sapply(1:ncol(rowAnnotation),function(i){
        sum(features%in%rowAnnotation[,i],na.rm=TRUE)
      })
      col1<-which(col1==max(col1))[1]
      idx<-rowAnnotation[[col1]]%in%features
      object@results$tn.ref[!idx,] <- 0
      object@results$tn.dpi[!idx,] <- 0
    }
    
    ##-----get regulons
    if(tnet=="ref"){
      listOfRegulonsAndMode <- tni.get(object,what="refregulons.and.mode")
    } else {
      listOfRegulonsAndMode <- tni.get(object,what="regulons.and.mode")
    }
    
    ##-----set regs
    if(!is.null(regulatoryElements)){
      if(sum(regulatoryElements%in%object@regulatoryElements) > 
         sum(regulatoryElements%in%names(object@regulatoryElements) ) ){
        regulatoryElements <- object@regulatoryElements[
          object@regulatoryElements%in%regulatoryElements]
      } else {
        regulatoryElements <- object@regulatoryElements[
          names(object@regulatoryElements)%in%regulatoryElements]
      }
      if(length(regulatoryElements)==0)
        stop("'regulatoryElements' argument has no valid names!")
    } else {
      regulatoryElements<-object@regulatoryElements
    }
    listOfRegulonsAndMode <- listOfRegulonsAndMode[regulatoryElements]
    
    ##-----check regulon size
    regcounts <- .regulonCounts(listOfRegulonsAndMode)
    if(sizeFilterMethod=="posANDneg"){
      idx <- regcounts$Positive >= minRegulonSize & 
        regcounts$Negative >= minRegulonSize
    } else if(sizeFilterMethod=="posORneg"){
      idx <- regcounts$Positive >= minRegulonSize | 
        regcounts$Negative >= minRegulonSize
    } else {
      idx <- regcounts$Size >= minRegulonSize
    }
    regulatoryElements <- regulatoryElements[
      regulatoryElements%in%rownames(regcounts)[idx]]
    listOfRegulonsAndMode <- listOfRegulonsAndMode[regulatoryElements]
    
    ##-----stop when no regulon passes the size requirement
    if(length(listOfRegulonsAndMode)==0){
      stop("no regulon passed the 'minRegulonSize' requirement!")
    }
    
    #--- get phenotypes
    if(log){
      phenotypes <- log2(1+gexp)-log2(1+gxref)
    } else {
      phenotypes <- gexp-gxref
    }
    
    #--- get regulons evaluated by EM algorithm
    if (verbose) {
      cat("Running EM algorithm... ")
    }
    if(tnet=="ref"){
      listOfRegulonsAndModeGmm <- tni.get(object,
                                          what="refregulons.and.mode.gmm")
    } else {
      listOfRegulonsAndModeGmm <- tni.get(object,what="regulons.and.mode.gmm")
    }
    listOfRegulonsAndModeGmm <- listOfRegulonsAndModeGmm[regulatoryElements]
    
    #--- set regulons for aREA
    arearegs <- list()
    for(rg in regulatoryElements){
      arearegs[[rg]]$tfmode <- listOfRegulonsAndModeGmm[[rg]]$gmm
      arearegs[[rg]]$likelihood <- listOfRegulonsAndModeGmm[[rg]]$mi
    }
    if (verbose) {
      cat("Running aREA algorithm...\n")
    }
    nes <- t(aREA(eset=phenotypes, regulon=arearegs, minsize=0, 
                  verbose=FALSE)$nes)
    nes <- nes[samples,regulatoryElements]
    colnames(nes) <- names(regulatoryElements)
    
    #-- for compatibility, wrap-up results into the same format
    regulonActivity <- list(differential=nes)
    regulonActivity <- .tni.stratification.area(regulonActivity)
    regulonActivity$regulatoryElements <- regulatoryElements
    object@results$regulonActivity <- regulonActivity
    
    #---
    if(verbose)cat("-GSEA2 complete! \n\n")
    object@status["Activity"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
## Constructor of TNA Class objects
## Entry point for the TNA pipeline, including pre-processing
setMethod(
  "tni2tna.preprocess",
  "TNI",
  function(object, phenotype=NULL, hits=NULL, phenoIDs=NULL, 
           duplicateRemoverMethod="max", verbose=TRUE) {
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    if(object@status["Preprocess"]!="[x]")
      stop("TNI object requires preprocessing!")
    if(object@status["Permutation"]!="[x]")
      stop("TNI object requires permutation/bootstrap and DPI filter!")  
    if(object@status["DPI.filter"]!="[x]")
      stop("TNI object requires DPI filter!")
    
    ##-----check input arguments
    tnai.checks(name="TNI",para=object)
    tnai.checks(name="phenotype",para=phenotype)
    tnai.checks(name="hits",para=hits)
    phenoIDs <- tnai.checks(name="phenoIDs",para=phenoIDs)
    tnai.checks(name="duplicateRemoverMethod",para=duplicateRemoverMethod)
    tnai.checks(name="verbose",para=verbose)
    ##-----generate a new object of class TNA
    .object <- new("TNA",
                   tni=object,
                   phenotype=phenotype,
                   hits=hits)
    if(!is.null(object@results$conditional) && 
       length(object@results$conditional)>0){
      cdt <- tni.get(object,what="cdt.list")
      lmod <- lapply(cdt,function(reg){
        if(nrow(reg)>0){
          tp <- reg$Mode
          names(tp) <- rownames(reg)
        } else {
          tp=character()
        }
        tp
      })
      .object@listOfModulators <- lmod
    }
    .object <- tna.preprocess(.object, phenoIDs=phenoIDs,
                              duplicateRemoverMethod=duplicateRemoverMethod,
                              verbose=verbose)
    return(.object)
  }
)

##------------------------------------------------------------------------------
##run conditional mutual information analysis
setMethod(
  "tni.conditional",
  "TNI",
  function(object, modulators, tfs=NULL, sampling=35, pValueCutoff=0.01, 
           pAdjustMethod="bonferroni", minRegulonSize=15, minIntersectSize=5, 
           miThreshold="md", prob=0.99, medianEffect=FALSE, 
           iConstraint=TRUE, verbose=TRUE, mdStability=FALSE){
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    ##-----check input arguments
    if(object@status["DPI.filter"]!="[x]")
      stop("input 'object' needs dpi analysis!")
    if(missing(modulators))
      stop("please provide a character vector with 'modulators'!")
    tnai.checks(name="modulators",para=modulators)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="sampling",para=sampling)
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="minIntersectSize",para=minIntersectSize)
    tnai.checks(name="miThreshold",para=miThreshold)
    tnai.checks(name="medianEffect",para=medianEffect)
    tnai.checks(name="prob",para=prob)
    tnai.checks(name="verbose",para=verbose)
    tnai.checks(name="iConstraint",para=iConstraint)
    #check additional (experimental) args
    if(is.logical(mdStability)){
      mrkboot<-NULL
    } else {
      mrkboot<-tnai.checks(name="mdStability.custom",para=mdStability)
      mdStability<-TRUE
    }
    ##-----par info
    object@para$cdt<-list(sampling=sampling, pValueCutoff=pValueCutoff,
                          pAdjustMethod=pAdjustMethod, 
                          minRegulonSize=minRegulonSize, 
                          minIntersectSize=minIntersectSize, 
                          miThreshold=NA, prob=prob,
                          iConstraint=iConstraint)
    ##-----summary info
    cdt<-unlist(object@para$cdt)
    object@summary$para$cdt<-matrix(cdt,nrow=1,ncol=8)
    rownames(object@summary$para$cdt)<-"Parameter"
    colnames(object@summary$para$cdt)<-names(cdt)

    if(verbose)cat("-Preprocessing for input data...\n")
    ##-----make sure all 'tfs' are valid
    if(is.null(tfs)){
      tfs<-object@regulatoryElements
    } else {
      if(verbose)cat("--Checking TFs in the dataset...\n")
      tfs<-as.character(tfs)
      idx<-which(!tfs%in%object@regulatoryElements & 
                   !tfs%in%names(object@regulatoryElements))
      if(length(idx)>0){
        message(paste("Note: input 'tfs' contains", length(idx),
                      " element(s) not listed in the network!\n"))
      }
      idx<-which(object@regulatoryElements%in%tfs | 
                   names(object@regulatoryElements)%in%tfs)
      if(length(idx)<1){
        stop(paste("NOTE: input 'tfs' contains no useful data!\n"))
      }      
      tfs<-object@regulatoryElements[idx]
    }
    
    ##-----make sure all 'modulators' are valid
    if(verbose)cat("--Checking modulators in the dataset...\n")
    modulators <- as.character(modulators)
    col1 <- sapply(1:ncol(object@rowAnnotation),function(i){
      sum(modulators%in%object@rowAnnotation[,i],na.rm=TRUE)
    })
    col1 <- which(col1==max(col1))[1]
    idx<-which(!modulators%in%object@rowAnnotation[[col1]])
    if(length(idx)>0){
      message(paste("Note: input 'modulators' contains", length(idx),
                    " element(s) not listed in the network!\n"))
    }
    idx <- object@rowAnnotation[[col1]] %in% modulators
    modulators <- rownames(object@rowAnnotation)[idx]
    if(length(modulators)==0){
      stop("input 'modulators' contains no useful data!\n")
    }
    ##-----make sure 'modulators' is a named vector
    if(!is.null(object@rowAnnotation$SYMBOL)){
      names(modulators)<-object@rowAnnotation[modulators,"SYMBOL"]
    } else {
      names(modulators)<-modulators
    }
    ##-----remove any empty space or NA from md names
    mdnames<-names(modulators)
    idx<-mdnames==""|mdnames=="NA"
    names(modulators)[idx]<-modulators[idx]
    object@modulators<-modulators
    
    ##-----get TF-targets from tnet
    if(verbose)cat("--Extracting TF-targets...\n")
    tfTargets<-list()
    tfAllTargets<-list()
    for(tf in tfs){
      idx<-object@results$tn.dpi[,tf]!=0
      tfTargets[[tf]]<-rownames(object@results$tn.dpi)[idx]
      idx<-object@results$tn.ref[,tf]!=0
      tfAllTargets[[tf]]<-rownames(object@results$tn.ref)[idx]
    }
    ##-----check regulon size
    gs.size <- unlist(
      lapply(tfTargets, length)
    )
    tfs<-tfs[tfs%in%names(gs.size[gs.size>=minRegulonSize])]
    tfTargets<-tfTargets[tfs]
    tfAllTargets<-tfAllTargets[tfs]
    tnetAllTargets<-rownames(object@results$tn.dpi)[
      rowSums(object@results$tn.dpi!=0)>0]
    ##-----Checking independence of modulators and TFs
    IConstraintList<-list()
    if(iConstraint){
      if(verbose)cat("--Applying modulator independence constraint...\n")
      temp_obj<-tni.dpi.filter(object, eps=0, verbose=FALSE)
      for(tf in tfs){
        idx<-temp_obj@results$tn.ref[,tf]!=0
        IConstraintList[[tf]]<-c(tf,rownames(temp_obj@results$tn.ref)[idx])
      }
    } else {
      for(tf in tfs){
        IConstraintList[[tf]]<-NA
      }
    }
    ##-----set sub-sample idx
    gxtemp<-object@gexp
    spsz<-round(ncol(gxtemp)*sampling/100,0)
    idxLow<-1:spsz
    idxHigh<-(ncol(gxtemp)-spsz+1):ncol(gxtemp)
    ##-----start filtering
    if(verbose)cat("--Applying modulator range constraint...\n")
    if(.isUnloggedData(gxtemp)) gxtemp <- .log2transform(gxtemp)
    gxtemp[is.na(gxtemp)]<-median(gxtemp,na.rm=TRUE)
    gxtemp<-t(apply(
      gxtemp[modulators,,drop=FALSE],1,sort))[,c(idxLow,idxHigh),drop=FALSE]
    if(length(modulators)==1){
      gxtemp<-rbind(gxtemp,gxtemp)
    }
    ##--run limma (for modulator range constraint)
    t <- factor(c(rep("low",spsz),rep("high",spsz)))
    design <- model.matrix(~0+t)
    fit <- lmFit(gxtemp,design)
    thigh=tlow=NULL
    contrasts <- makeContrasts(thigh-tlow, levels=design)
    ct.fit <- eBayes(contrasts.fit(fit, contrasts))
    res.fit<-unclass(decideTests(ct.fit, 
                                 adjust.method=object@para$cdt$pAdjustMethod, 
                                 p.value=object@para$cdt$pValueCutoff))
    RConstraintList<-rownames(res.fit)[res.fit<=0]
    ##--get samples (sorted index) for each pre-selected modulator
    if(verbose)cat("--Selecting subsamples...\n")
    gxtemp<-object@gexp
    gxtemp[is.na(gxtemp)]<-median(gxtemp,na.rm=TRUE)    
    idx<-t(apply(gxtemp[modulators,,drop=FALSE],1,sort.list))
    idxLow<-idx[,idxLow,drop=FALSE]
    idxHigh<-idx[,idxHigh,drop=FALSE]
    
    ##-----estimate mutual information threshold
    if(is.character(miThreshold)){
      if(verbose)cat("\n")
      if(verbose)cat("-Estimating mutual information threshold...\n")
      if(miThreshold=="md.tf"){
        mimark<-miThresholdMdTf(gxtemp[object@targetElements,],tfs=tfs,
                                nsamples=spsz,prob=prob,
                                nPermutations=object@para$perm$nPermutations, 
                                estimator=object@para$perm$estimator,
                                verbose=verbose)
      } else {
        mimark<-miThresholdMd(gxtemp[object@targetElements,],nsamples=spsz,
                              prob=prob,
                              nPermutations=object@para$perm$nPermutations, 
                              estimator=object@para$perm$estimator,
                              verbose=verbose)
      }
    } else {
      mimark<-sort(miThreshold)
      miThreshold<-"md"
      if(length(mimark)==1){
        mimark<-abs(mimark)
        mimark<-c(-mimark,mimark)
      } else {
        if(sum(mimark>0)!=1)
          stop("'miThreshold' upper and lower bounds should have different signals!")
      }
      object@summary$para$cdt[,"prob"]<-"custom"
      object@para$cdt$prob<-NA
    }
    ##-----update miThreshold
    object@summary$para$cdt[,"miThreshold"]<-miThreshold
    object@para$cdt$miThreshold<-mimark
    
    ##-----set data object to save results
    reseffect<-lapply(tfTargets,function(tar){
      data.frame(targets=tar,stringsAsFactors=FALSE)
    })
    rescount<-lapply(tfTargets,function(tar){
      res<-data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                      stringsAsFactors=FALSE)
      colnames(res)<-c("Modulator","irConstraint","nConstraint","TF",
                       "UniverseSize","EffectSize","RegulonSize","Expected",
                       "Observed","Negative","Positive",
                       "Mode","PvFET","AdjPvFET","KS","PvKS","AdjPvKS","R")
      res
    })
    
    ##-----start conditional mutual information analysis
    modregulons<-list()
    glstat<-list()
    if(verbose)cat("\n")
    if(verbose)cat("-Performing conditional mutual information analysis...\n")
    if(verbose)cat("--For", length(tfs), "tfs and" , length(modulators), 
                   "candidate modulator(s) \n")
    if(verbose && !mdStability) pb <- txtProgressBar(style=3)
    for(i in 1:length(modulators)){
      md<-modulators[i]
      #get sample ordering
      lw<-idxLow[md,]
      hg<-idxHigh[md,]
      #compute mi on both tails
      milow<-tni.pmin(gxtemp[object@targetElements,lw],tfs,
                      estimator=object@para$perm$estimator)
      mihigh<-tni.pmin(gxtemp[object@targetElements,hg],tfs,
                       estimator=object@para$perm$estimator)
      milow[is.na(milow)]<-0
      mihigh[is.na(mihigh)]<-0
      #get mi delta
      miDelta<-mihigh-milow
      #identify modulations above mi threshold
      if(miThreshold=="md.tf" && length(tfs)>1){
        sigDelta<-t(apply(miDelta,1,"<",mimark[,1])) | 
          t(apply(miDelta,1,">",mimark[,2]))
      } else {
        sigDelta<- miDelta<mimark[1] | miDelta>mimark[2]
      }
      miDelta<-miDelta/(mihigh+milow)
      
      #check modulator stability (experimental)
      #computational cost forbids default use or advanced customization
      #better only for final verification, and one modulator each time!
      if(mdStability){
        if(is.null(mrkboot)){
          if(verbose)cat("\n-Estimating stability threshold...\n")
          mrkboot<-miThresholdMd(gxtemp[object@targetElements,],nsamples=spsz,
                                 prob=c(0.05, 0.95),
                                 nPermutations=1000,
                                 estimator=object@para$perm$estimator,
                                 verbose=verbose)
        }
        if(verbose)cat("-Checking modulation stability for", names(md), "\n")
        stabt<-cdt.stability(gxtemp,object@para$perm$estimator,
                             mrkboot,miThreshold,
                             spsz,md,tfs,object@targetElements,sigDelta,
                             nboot=100,consensus=75,verbose=verbose)
        sigDelta[!stabt]<-FALSE
      } else {
        if(verbose) setTxtProgressBar(pb, i/length(modulators))
      }
      
      #decision
      miDelta[!sigDelta]<-0
      
      #run main analysis
      for(tf in tfs){
        
        dtvec<-miDelta[,tf]
        tftar<-tfTargets[[tf]]
        tfalltar<-tfAllTargets[[tf]]
        #---get SZs
        #all modulated
        EffectSZ<-sum(dtvec[tnetAllTargets]!=0)
        #all tested
        UniSZ<-length(tnetAllTargets)
        #others
        RegSZ<-length(tftar)
        ExpOV<-(EffectSZ*RegSZ)/UniSZ
        ObsOV<-sum(dtvec[tftar]!=0)
        
        #---run stats
        ObsPos<-sum(dtvec[tftar]>0)
        ObsNeg<-sum(dtvec[tftar]<0)
        #check exclusion list (independence/range constraint)
        irconst<-md%in%IConstraintList[[tf]] || md%in%RConstraintList
        #check minimum number of modulated targets for testing (n constraint)
        nconst<-(ObsOV/RegSZ*100)<minIntersectSize # || ExpOV<1
        if(irconst || nconst ){
          Mode<-NA;pvfet<-NA; dks<-NA;pvks<-NA;dtvec[]<-0; R <- NA
        } else {
          #---get mode of actions (and R, for reference)
          Mode <- if(ObsNeg>ObsPos) -1 else if(ObsNeg<ObsPos) 1 else 0
          R <- cor(gxtemp[tf,],gxtemp[md,], method = "spearman")
          #---set obs to predicted mode
          Obs<-if(Mode==-1) abs(ObsNeg) else if(Mode==1) ObsPos else ObsOV
          #---run fet with phyper (obs-1)
          pvfet <- phyper(Obs-1, RegSZ, UniSZ-RegSZ, EffectSZ, 
                          lower.tail=FALSE)
          #---run ks test
          #pheno
          pheno<-abs(object@results$tn.ref[,tf])
          pheno<-pheno[pheno!=0]
          #hits
          hits<-dtvec[tftar]
          hits<-hits[hits!=0]
          hits<-which(names(pheno)%in%names(hits))
          if(length(hits)>length(pheno)/2){
            dks<-1
            pvks<-0
          } else {
            kst<-suppressWarnings(ks.test(pheno[-hits],pheno[hits],
                                          alternative="greater"))
            dks<-kst$statistic
            pvks<-kst$p.value
          }
          #count (+) and (-) tf-targets in the modulated set
          #expressed by the ratio of (+) or (-) targets, respectively
          if(Mode>=0){
            mdtf.tar<-dtvec[tftar]
            mdtf.tar<-names(mdtf.tar)[mdtf.tar>0]
            tf.tar<-object@results$tn.dpi[tftar,tf]
            mdtf.tar<-tf.tar[mdtf.tar]
            p1<-sum(mdtf.tar>0)/sum(tf.tar>0)
            p2<-sum(mdtf.tar<0)/sum(tf.tar<0)
            bl<-!is.nan(p1) && !is.nan(p2)
          } else {
            mdtf.tar<-dtvec[tftar]
            mdtf.tar<-names(mdtf.tar)[mdtf.tar<0]
            tf.tar<-object@results$tn.dpi[tftar,tf]
            mdtf.tar<-tf.tar[mdtf.tar]
            p1<-sum(mdtf.tar>0)/sum(tf.tar>0)
            p2<-sum(mdtf.tar<0)/sum(tf.tar<0)
            bl<-!is.nan(p1) && !is.nan(p2)
          }
        }
        
        #---add results to a list
        reseffect[[tf]][[md]]<-dtvec[tftar]
        rescount[[tf]][md,]<-c(NA,NA,NA,NA,UniSZ,EffectSZ,RegSZ,ExpOV,ObsOV,
                               ObsNeg,ObsPos,Mode,pvfet,NA,dks,pvks,NA,R)
        rescount[[tf]][md,c(1,4)]<-c(md,tf)
        rescount[[tf]][md,c(2,3)]<-c(irconst,nconst)
        #---retain modulated targets
        mdtftar<-tftar[ dtvec[tftar]!=0]
        if(length(mdtftar)>1){
          modregulons[[md]][[tf]]<-mdtftar
        } else {
          modregulons[[md]][[tf]]<-c(NA,NA)
          modregulons[[md]][[tf]]<-mdtftar
        }
        
      }
      
      #compute mi differential score for each regulon (signal-to-noise ratio)
      #this is a global stats, only used to assess the median effect 
      #for the selected regulons
      if(medianEffect){
        sig2noise<-sapply(names(modregulons[[md]]),function(tf){
          tftar<-modregulons[[md]][[tf]]
          h<-mihigh[tftar,tf]
          l<-milow[tftar,tf]
          (median(h)-median(l))/(sd(h)+sd(l)) 
        })
        sig2noise[is.na(sig2noise)]<-0
        glstat$observed[[md]]$sig2noise<-sig2noise
      }
      
    }
    if(verbose && !mdStability) close(pb)
    
    #set data format
    for(tf in names(rescount)){
      results<-rescount[[tf]][-1,,drop=FALSE]
      if(nrow(results)>0){
        results[,"Expected"]<-round(results[,"Expected"],2)
        results[,"KS"]<-round(results[,"KS"],2)
      }
      rescount[[tf]]<-results
    }
    rescount<-sortblock.cdt.list(cdt=rescount,coln="PvFET")
    
    ##update summary
    object@results$conditional$count<-rescount
    object@results$conditional$effect<-reseffect
    
    #compute null based on each regulon's distribution
    #this is a global stats, only used to assess the median effect on regulons
    #...not use to infer the modulated targets
    if(medianEffect){
      if(verbose)cat("\n")
      if(verbose)cat("-Checking median modulation effect...\n") 
      modulatedTFs<-tni.get(object,what="cdt.list")
      modulatedTFs<-unlist(lapply(modulatedTFs,nrow))
      modulatedTFs<-names(modulatedTFs)[modulatedTFs>0]      
      if(length(modulatedTFs)>0){
        if(verbose)cat("--For", length(modulators), "candidate modulator(s) \n")
        res<-checkModuationEffect(gxtemp,tfs,modregulons,modulatedTFs,glstat,
                                  spsz, minRegulonSize,pValueCutoff,
                                  nPermutations=object@para$perm$nPermutations,
                                  estimator=object@para$perm$estimator,
                                  pAdjustMethod=pAdjustMethod,
                                  count=object@results$conditional$count,
                                  verbose)
        res$md2tf$count<-p.adjust.cdt.list(cdt=res$md2tf$count,
                                           pAdjustMethod=pAdjustMethod,
                                           p.name="PvSNR",
                                           adjp.name="AdjPvSNR",
                                           roundpv=FALSE, global=FALSE)       
        object@results$conditional$mdeffect$md2tf$null<-res$md2tf$null
        object@results$conditional$mdeffect$md2tf$observed<-res$md2tf$observed
        object@results$conditional$mdeffect$tf2md$null<-res$tf2md$null
        object@results$conditional$mdeffect$tf2md$observed<-res$tf2md$observed       
        object@results$conditional$count<-res$md2tf$count
        if(verbose)cat("\n")
      }
  }
  object@status["Conditional"] <- "[x]"
  if(verbose)cat("-Conditional analysis complete! \n\n")
  return(object)
  }
)
#supplementary information: get simple correlation between tfs and modulators
# tni.tfmdcor<-function(x,tfs, mds, estimator="pearson",dg=0, asInteger=FALSE){
#   ids<-unique(c(tfs,setdiff(mds,tfs)))
#   x=x[ids,]
#   x=t(x)
#   #--
#   pcorm=cor(x[,tfs],x[,mds], method=estimator,use="complete.obs")
#   if(asInteger){
#     pcorm[pcorm<0]=-1
#     pcorm[pcorm>0]=1
#   }
#   #--
#   pcorm<-t(pcorm)
#   colnames(pcorm)<-tfs
#   pcorm
# }
#rnet<-tni.tfmdcor(object@gexp,tfs, modulators)

##------------------------------------------------------------------------------
##get graph from TNI
## experimental args:
## mask: a logical value specifying to apply a mask on the 'amapFilter', 
## keeping at least the 
## ......best weighted edge (when verbose=TRUE) or not (when verbose=FALSE).
## hcl: an hclust object with TF's IDs
## overlap: overlapping nodes used for the Jaccard 
## (options: 'all', 'pos', 'neg')
## TODO: revise 'tnai.checks' for new args!
setMethod(
  "tni.graph",
  "TNI",
  function(object, tnet="dpi", gtype="rmap", minRegulonSize=15, regulatoryElements=NULL,
           amapFilter="quantile", amapCutoff=NULL, ntop=NULL, mask=FALSE, 
           hcl=NULL, overlap="all", xlim=c(30,80,5), nquant=5, breaks=NULL, 
           mds=NULL, nbottom=NULL){
    # chech igraph compatibility
    b1<-"package:igraph0" %in% search()
    b2<- "igraph0" %in%  loadedNamespaces()
    if( b1 || b2) {
      stop("\n\n ...conflict with 'igraph0': please use the new 'igraph' package!")
    }
    #---check compatibility
    object <- upgradeTNI(object)
    
    ##-----check input arguments
    if(object@status["DPI.filter"]!="[x]")
      stop("input 'object' needs dpi analysis!")
    tnai.checks(name="tnet",para=tnet)
    tnai.checks(name="tni.gtype",para=gtype)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="regulatoryElements",para=regulatoryElements)
    tnai.checks(name="ntop",para=ntop)
    tnai.checks(name="mds",para=mds)
    tnai.checks(name="amapFilter",para=amapFilter)
    tnai.checks(name="amapCutoff",para=amapCutoff)
    tnai.checks(name="mask",para=mask)
    #---
    if(gtype=="mmap" || gtype=="mmapDetailed")
      tnet="dpi"
    if(tnet=="ref"){
      tnet<-object@results$tn.ref
    } else {
      tnet<-object@results$tn.dpi
    }
    #---
    if(!is.null(hcl)){
      gtype="amapDend"
      tfs<-object@regulatoryElements
      if(!all(hcl$labels%in%tfs | hcl$labels%in%names(tfs)))
        stop("all labels in the 'hclust' object should be listed as 'regulatoryElements'!")
      idx1<-match(hcl$labels,tfs)
      idx2<-match(hcl$labels,names(tfs))
      check<-which(is.na(idx1))
      idx1[check]<-idx2[check]
      tfs<-tfs[idx1]
      hcl$labels<-tfs
    } else if(is.null(regulatoryElements)){
      tfs<-object@regulatoryElements
      minsz<-colnames(tnet)[colSums(tnet!=0)>=minRegulonSize]
      tfs<-tfs[tfs%in%minsz]
    } else {
      tfs<-as.character(regulatoryElements)
      idx<-which(names(object@regulatoryElements)%in%tfs | 
                   object@regulatoryElements%in%tfs)
      if(length(idx)==0)stop("input 'tfs' contains no useful data!\n")
      tfs<-object@regulatoryElements[idx]
    }
    
    #-----------------------------------------
    #-----------------------------------------
    
    if(gtype=="mmap" || gtype=="mmapDetailed"){ #get modulatory maps
      
      ##-----check input arguments
      if(object@status["Conditional"]!="[x]")
        stop("input needs conditional analysis!")
      #get tfs and modulators
      cdt<-tni.get(object,what="cdt.list")
      if(length(cdt)==0)
        stop("input conditional analysis is empty")
      testedtfs<-names(cdt)
      testedtfs<-object@regulatoryElements[
        object@regulatoryElements%in%testedtfs]
      testedtfs<-testedtfs[testedtfs%in%tfs]
      if(length(testedtfs)==0)
        stop("input 'tfs' contains no useful data!\n")
      modulators<-sapply(testedtfs,function(tf){
        rownames(cdt[[tf]])
      })
      modulators<-unlist(modulators)
      modulators<-object@modulators[object@modulators%in%modulators]
      othertfs<-object@regulatoryElements
      othertfs<-othertfs[!othertfs%in%testedtfs]
      othertfs<-othertfs[othertfs%in%modulators]
      #get adjmt
      if(!all(modulators %in% rownames(tnet))){
        stop("for 'mmap', all 'modulators' should be 'targetElements'!")
      }
      tnet<-tnet[unique(c(testedtfs,setdiff(modulators,testedtfs))),
                 testedtfs,drop=FALSE]
      mnet<-tnet;mnet[,]=0
      for(i in colnames(mnet)){
        tp<-cdt[[i]]
        mnet[rownames(tp),i]<-tp$Mode
      }
      pvnet<-tnet;pvnet[,]=1
      for(i in colnames(mnet)){
        tp<-cdt[[i]]
        pvnet[rownames(tp),i]<-tp$PvKS
      }
      #---
      if(gtype=="mmapDetailed"){
        #---experimental!!!
        #return a lista with:
        #1st: a TF
        #2nd: all MDs of a TF
        #3rd: a graph
        if(is.null(mds)){
          mds <- modulators
        } else {
          idx1<-mds%in%modulators
          idx2<-mds%in%names(modulators)
          if(!all(idx1) & !all(idx2)){
            stop("one or more modutors in 'mds' not listed in the TNI object!")
          } else {
            if(sum(idx1)>sum(idx2)){
              mds<-modulators[modulators%in%mds]
            } else {
              mds<-modulators[names(modulators)%in%mds]
            }
          }
        }
        g<-tni.mmap.detailed(object,mnet,testedtfs, mds, ntop=ntop,
                             nbottom=nbottom)
      } else {
        g<-tni.mmap(object,mnet,tnet,pvnet,othertfs,testedtfs,modulators)
      }
      return(g)
      
    } else if(gtype=="rmap"){ 
      
      tnet<-tnet[,tfs,drop=FALSE]
      g<-tni.rmap(tnet)
      #add rowAnnotation
      if(nrow(object@rowAnnotation)>0 & ncol(object@rowAnnotation)>1){
        g<-att.mapv(g=g,dat=object@rowAnnotation,refcol=1)
      }
      #set target names if available
      if(!is.null(V(g)$SYMBOL)){
        g<-att.setv(g=g, from="SYMBOL", to='nodeAlias')
      } else {
        V(g)$nodeAlias<-V(g)$name
      }
      #set TF names
      V(g)$tfs<-as.numeric(V(g)$name%in%tfs)
      idx<-match(tfs,V(g)$name)
      V(g)$nodeAlias[idx]<-names(tfs)
      V(g)$nodeColor<-"black"
      V(g)$nodeLineColor<-"black"
      g<-att.setv(g=g, from="tfs", to='nodeShape',title="")
      g$legNodeShape$shape <- rev(g$legNodeShape$shape)
      g$legNodeShape$legend<-c("Regulator","Target")
      g<-att.setv(g=g, from="tfs", to='nodeSize', xlim=c(20,50,1))
      g<-att.setv(g=g, from="tfs", to='nodeFontSize',xlim=c(10,32,1))
      #remove non-usefull legends
      g<-remove.graph.attribute(g,"legNodeSize")
      g<-remove.graph.attribute(g,"legNodeFontSize")
      if(ecount(g)>0){
        #set edge attr
        g<-att.sete(g=g, from="modeOfAction", to='edgeColor',
                    cols=c("#96D1FF","grey80","#FF8E91"), 
                    title="Interaction",categvec=-1:1)
        g$legEdgeColor$legend<-c("(-)","NA","(+)")
        E(g)$edgeWidth<-1.5
        #map modeOfAction to node attribute (compute average interaction)
        el<-data.frame(get.edgelist(g),E(g)$modeOfAction,stringsAsFactors=FALSE)
        nid<-V(g)$name
        mdmode<-sapply(nid,function(id){
          idx<-el[,2]==id
          median(el[idx,3])
        })
        mdmode[V(g)$tfs==1]=NA
        V(g)$medianModeOfAction<-as.integer(mdmode)
        #assign mode to targets
        g<-att.setv(g=g, from="medianModeOfAction", to='nodeColor',
                    cols=c("#96D1FF","grey80","#FF8E91"), 
                    title="ModeOfAction",categvec=-1:1,pal=1,na.col="grey80")
        V(g)$nodeLineColor<-V(g)$nodeColor
        g<-remove.graph.attribute(g,"legNodeColor")
      }
      return(g)
      
    } else if(gtype=="amap"){
      
      tnet<-tnet[,tfs,drop=FALSE]
      adjmt<-tni.amap(tnet,overlap)
      #-------------------filter J.C.
      if(mask){
        #set a mask to keep at least the best weighted edge
        mask<-sapply(1:ncol(adjmt),function(i){
          tp<-adjmt[,i]
          tp==max(tp)
        })
        nc<-ncol(mask);nr<-nrow(mask)
        mask<-mask+mask[rev(nr:1),rev(nc:1)]>0
      } else {
        mask<-array(0,dim=dim(adjmt))
      }
      if(amapFilter=="phyper"){
        #filter based phyper distribution (remove non-significant overlaps)
        if(is.null(amapCutoff))amapCutoff=0.01
        pvalue<-amapCutoff
        pmat<-tni.phyper(tnet)
        adjmt[pmat>pvalue & mask==0]=0
      } else if(amapFilter=="quantile"){
        #filter based on quantile distribution
        if(is.null(amapCutoff))amapCutoff=0.75
        jc<-as.integer(amapCutoff*100)+1
        tp<-as.numeric(adjmt)
        jc<-quantile(tp[tp>0],probs = seq(0, 1, 0.01), na.rm=TRUE)[jc]
        adjmt[adjmt<jc & mask==0]=0
      } else {
        #custom filter
        if(is.null(amapCutoff))amapCutoff=0
        adjmt[adjmt<amapCutoff & mask==0]=0
      }
      #-------------------
      g<-igraph::graph.adjacency(adjmt, diag=FALSE, mode="undirected", 
                                 weighted=TRUE)
      if(nrow(object@rowAnnotation)>0 & ncol(object@rowAnnotation)>1){
        g<-att.mapv(g=g,dat=object@rowAnnotation,refcol=1)
      }
      sz<-apply(tnet!=0, 2, sum)
      idx<-match(V(g)$name,tfs)
      V(g)$nodeAlias<-names(tfs)[idx]
      V(g)$degree<-sz[idx]
      #---set main attribs
      if(ecount(g)>0)g<-att.sete(g=g, from="weight", to='edgeWidth', 
                                 nquant=nquant, xlim=c(1,15,1),roundleg=2)
      g<-att.setv(g=g, from="degree", to='nodeSize', xlim=xlim, 
                  nquant=nquant, breaks=breaks, roundleg=1,title="Regulon size")
      V(g)$nodeFontSize<-20
      return(g)
      
    } else if(gtype=="amapDend"){
      #---experimental!!!
      #returns a lista with:
      #1st: an igraph
      #2nd: a list with nests
      #3rd: an hclust object
      if(!is.null(hcl)){
        gg<-hclust2igraph(hcl)
      } else {
        x<-tni.amap(tnet[,tfs], overlap)
        diag(x)=1
        hcl <- hclust(as.dist(1-cor(x)), method='complete')
        gg<-hclust2igraph(hcl)
      }
      gg$hcl<-hcl
      #---set alias
      idx<-match(V(gg$g)$name,tfs)
      V(gg$g)$nodeAlias<-names(tfs)[idx]
      V(gg$g)$nodeAlias[is.na(idx)]<-"$hcnode"
      #---set node degree
      V(gg$g)$degree<-2
      sz<-apply(tnet!=0, 2, sum)
      idx<-match(V(gg$g)$name,names(sz))
      V(gg$g)$degree<-sz[idx]
      #---set nest size
      V(gg$g)$nestSize<-V(gg$g)$degree
      nestsz<-sapply(names(gg$nest),function(nest){
        #length(gg$nest[[nest]]) #..count only TFs
        sum(rowSums(tnet[,gg$nest[[nest]]]!=0)>=1)
      })
      idx<-match(names(nestsz),V(gg$g)$name)
      V(gg$g)$nestSize[idx]<-nestsz
      #---set main attribs
      gg$g<-att.setv(g=gg$g, from="degree", to='nodeSize', xlim=xlim, 
                     breaks=breaks, nquant=nquant, roundleg=1, 
                     title="Regulon size")
      V(gg$g)$internalNode <- FALSE
      V(gg$g)$internalNode[is.na(V(gg$g)$degree)] <- TRUE
      V(gg$g)$nodeSize[V(gg$g)$internalNode] <- 10
      E(gg$g)$edgeWidth <- 10
      V(gg$g)$nodeFontSize<-20
      V(gg$g)$nodeFontSize[V(gg$g)$nodeAlias=="$hcnode"]<-1
      V(gg$g)$nodeColor<-"black"
      V(gg$g)$nodeLineColor<-"black"
      E(gg$g)$edgeColor<-"black"
      return(gg)
      
    }
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "tni.replace.samples",
  "TNI",
  function(object, expData, rowAnnotation = NULL, colAnnotation = NULL,
           removeRegNotAnnotated = TRUE, verbose = TRUE){
    
    #--- check compatibility
    object <- upgradeTNI(object)
    
    #--- check arguments
    if(object@status["DPI.filter"]!="[x]")
      stop("input 'object' needs dpi analysis!")
    tnai.checks("verbose",verbose)
    tnai.checks("removeRegNotAnnotated",removeRegNotAnnotated)
    
    #--- check if expData is a summarizedExperiment
    if (is(expData, "SummarizedExperiment") || 
        is(expData, "RangedSummarizedExperiment")) {
      if( !is.null(rowAnnotation) || !is.null(colAnnotation) ){
        tp1 <- "When using a 'SummarizedExperiment' container, "
        tp2 <- "row and col annotations are used from that object!"
        warning(tp1,tp2,call.=FALSE)
      }
      if (length(assays(expData)) > 1) {
        stop("please input a SummarizedExperiment with only one assay")
      }
      rowAnnotation <- as.data.frame(rowData(expData))
      colAnnotation <- as.data.frame(colData(expData))
      expData <- assays(expData)[[1]]
    }
    
    #--- check main arguments
    if(verbose)
      cat("--Mapping 'expData' to 'rowAnnotation' and 'colAnnotation'...\n")
    tmp <- .expDataChecks(expData, rowAnnotation, colAnnotation, verbose=verbose)
    new_expData <- tmp$expData
    new_rowAnnotation <- tmp$rowAnnotation
    new_colAnnotation <- tmp$colAnnotation
    
    #--- match annotation
    if(verbose)cat("--Checking new and old 'rowAnnotation'...")
    old_rowAnnotation <- tni.get(object, "rowAnnotation")
    if(ncol(new_rowAnnotation)==1) 
      new_rowAnnotation <- cbind(new_rowAnnotation, 
                                 'NA'=rep(NA, nrow(new_rowAnnotation)))
    if(ncol(old_rowAnnotation)==1) 
      old_rowAnnotation <- cbind(old_rowAnnotation, 
                                 'NA'=rep(NA, nrow(old_rowAnnotation)))  
    ijcol <- sapply(1:ncol(old_rowAnnotation),function(i){
      sapply(1:ncol(new_rowAnnotation),function(j){
        tp1 <- old_rowAnnotation[,i]; tp1 <- tp1[!is.na(tp1)]
        tp2 <- new_rowAnnotation[,j]; tp2 <- tp2[!is.na(tp2)]
        sum(tp1%in%tp2)
      })
    })
    ijcol <- which(ijcol==max(ijcol), arr.ind = TRUE)[1,]
    checkmatch <- old_rowAnnotation[[ijcol[2]]]%in%new_rowAnnotation[[ijcol[1]]]
    checkmatch <- sum(checkmatch)/nrow(old_rowAnnotation)*100
    if(checkmatch>75){
      if(verbose)cat(paste0(round(checkmatch,2),"% agreement! \n"))
    } else if(checkmatch >= 50 && checkmatch < 75){
      warning("NOTE: ",100-checkmatch,
              "% of the TNI is not represented in new 'rowAnnotation'!")
    } else if(checkmatch < 50){
      stop("\n",round(100-checkmatch,2),
           "% of the TNI is not represented in the new 'rowAnnotation'!")
    }
    
    #--- Simplify new_rowAnnotation
    #--- keep only what will be used to match the datasets
    idx <- new_rowAnnotation[[ijcol[1]]] %in% old_rowAnnotation[[ijcol[2]]]
    new_rowAnnotation <- new_rowAnnotation[idx,,drop=FALSE]
    new_expData <- new_expData[rownames(new_rowAnnotation),,drop=FALSE]
    new_rowAnnotation <- data.frame(ID_NEW=rownames(new_rowAnnotation),
                                    ID_KEY=new_rowAnnotation[[ijcol[1]]], 
                                    row.names = rownames(new_rowAnnotation),
                                    stringsAsFactors = FALSE)
    
    #--- Remove NAs in the key
    idx <- is.na(new_rowAnnotation$ID_KEY) | new_rowAnnotation$ID_KEY=="" | 
      new_rowAnnotation$ID_KEY=="NA"
    new_rowAnnotation <- new_rowAnnotation[!idx,]
    new_expData <- new_expData[rownames(new_rowAnnotation),]
    
    #--- Remove duplications in the key
    if(any(duplicated(new_rowAnnotation$ID_KEY))){
      if(verbose)
        cat("--Removing duplicated genes (keep max coefficient of variation!)...\n")
      #e.g. col1=id, col2=symbol (collapse cv by col2)
      cvres <- cv.filter(new_expData, new_rowAnnotation)
      new_expData <- cvres$gexp
      new_rowAnnotation <- cvres$ids
    }
    
    #--- Final check
    idxFinal <- match(new_rowAnnotation$ID_KEY, old_rowAnnotation[[ijcol[2]]])
    check <- all(new_rowAnnotation$ID_KEY==old_rowAnnotation[[ijcol[2]]][idxFinal])
    if(!check){
      stop("unpredicted exception found in the input data! 
           ...a possible cause is annotation mismatch!")
    }
    
    if(verbose)cat("--Replacing samples and annotation...\n")
    
    #--- Combine keys
    new_rowAnnotation <- cbind(new_rowAnnotation, 
                               old_rowAnnotation[idxFinal,,drop=FALSE])
    rownames(new_rowAnnotation) <- rownames(old_rowAnnotation)[idxFinal]
    new_rowAnnotation <- new_rowAnnotation[,!names(new_rowAnnotation)=="NA", 
                                           drop=FALSE]
    
    #--- Update'targetElements'  listed in the new_rowAnnotation
    object@results$tn.ref <- object@results$tn.ref[
      rownames(new_rowAnnotation),,drop=FALSE]
    object@results$tn.dpi <- object@results$tn.dpi[
      rownames(new_rowAnnotation),,drop=FALSE]
    rownames(object@results$tn.ref) <- new_rowAnnotation$ID_NEW
    rownames(object@results$tn.dpi) <- new_rowAnnotation$ID_NEW
    object@targetElements <- new_rowAnnotation$ID_NEW
    
    #--- Remove 'regulatoryElements' not listed in the new_rowAnnotation
    if(removeRegNotAnnotated){
      idx <- object@regulatoryElements %in% rownames(new_rowAnnotation)
      if(sum(!idx)>0){
        if(verbose) cat("--Removing 'regulatoryElements' not annotated...\n")
        object@regulatoryElements <- object@regulatoryElements[idx]
        object@results$tn.ref <- object@results$tn.ref[
          ,object@regulatoryElements,drop=FALSE]
        object@results$tn.dpi <- object@results$tn.dpi[
          ,object@regulatoryElements,drop=FALSE]
      }
    }
    #--- Update 'regulatoryElements' annotation
    idx <- match(object@regulatoryElements, rownames(new_rowAnnotation))
    object@regulatoryElements[] <- new_rowAnnotation$ID_NEW[idx]
    colnames(object@results$tn.ref) <- as.character(object@regulatoryElements)
    colnames(object@results$tn.dpi) <- as.character(object@regulatoryElements)
    
    #--- Update TNI with new samples and annotation
    object@gexp <- new_expData
    object@rowAnnotation <- new_rowAnnotation
    object@colAnnotation <- new_colAnnotation
    rownames(object@gexp) <- new_rowAnnotation$ID_NEW
    rownames(object@rowAnnotation) <- new_rowAnnotation$ID_NEW
    
    return(object)
  }
)

#-------------------------------------------------------------------------------
## tni.regulon.summary returns a summary of useful information about a 
## particular regulon(s) or the network, to aid in interpretation
setMethod(
  "tni.regulon.summary",
  "TNI",
  function(object, regulatoryElements = NULL, verbose = TRUE) {
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    #-- Basic checks
    if(object@status["DPI.filter"]!="[x]")
      stop("input 'object' needs dpi analysis!")
    if(!is.null(regulatoryElements))
      tnai.checks("regulatoryElements",regulatoryElements)
    tnai.checks("verbose",verbose)
    
    #-- if regulatoryElements = NULL, get a summary of network as a whole
    if (is.null(regulatoryElements)){
      networkSummary <- tni.get(object)$results
      if(verbose){
        nRegulators <- paste(
          "Regulatory network comprised of", 
          networkSummary$tnet["tnet.dpi", "regulatoryElements"],
          "regulons. \n")
        cat(nRegulators)
        message("-- DPI-filtered network: ")
        print(networkSummary$tnet["tnet.dpi",], quote = FALSE)
        print(networkSummary$regulonSize["tnet.dpi",], 
              quote = FALSE, digits = 3)
        message("-- Reference network: ")
        print(networkSummary$tnet["tnet.ref",], quote = FALSE) 
        print(networkSummary$regulonSize["tnet.ref",], 
              quote = FALSE, digits = 3)
        cat("---\n")
      }
      invisible(networkSummary)
    } else { #-- Otherwise, get TF summaries
      regnames <- tni.get(object, "regulatoryElements")
      if(sum(regulatoryElements%in%regnames) > 
         sum(regulatoryElements%in%names(regnames))){
        regulatoryElements <- regnames[regnames%in%regulatoryElements]
      } else {
        regulatoryElements <- regnames[names(regnames)%in%
                                         regulatoryElements]
      }
      if(length(regulatoryElements)==0)
        stop("'regulatoryElements' argument has no valid names!")
      
      #-- Get info of all regulons and refregulons
      tnet <- tni.get(object, "tnet")
      refnet <-  tni.get(object, "refnet")
      
      #-- Get summaries
      allRegulonSummary <- lapply(regulatoryElements, .regulon.summary, 
                                  refnet, tnet, regnames)
      #-- Print
      if(verbose){
        for(tf in names(regulatoryElements)) {
          regSummary <- allRegulonSummary[[tf]]
          
          if (length(regSummary$regulatorsMI) <= 10) {
            textRegsMI <- paste0(paste(names(regSummary$regulatorsMI),
                                       collapse = ", "), "\n", "\n")
          } else {
            textRegsMI <- paste0(paste(names(regSummary$regulatorsMI)[1:10],
                                       collapse = ", "),
                                 "...[", 
                                 length(regSummary$regulatorsMI) - 10, 
                                 " more]", "\n", "\n")
          }
          nTars <- regSummary$targets["DPInet", "Total"]
          
          #-- Size info
          if(nTars < 50) regsize <- "small"
          else if (nTars < 200) regsize <- "medium-sized"
          else regsize <- "large"
          #-- Balance info
          posTars <- regSummary$targets["DPInet","Positive"]
          regbalance <- ifelse(posTars > 0.75*nTars || posTars < 0.25*nTars,
                               "unbalanced", "balanced")
          #-- Print
          cat(paste("The", tf, "regulon", "has", nTars, 
                    "targets, it's a", regsize, 
                    "and", regbalance,
                    "regulon.", "\n"))
          message("-- DPI filtered network targets:")
          print(regSummary$targets["DPInet",], quote = FALSE)
          message("-- Reference network targets:")
          print(regSummary$targets["Refnet",], quote = FALSE)
          message("-- Regulators with mutual information:")
          cat(textRegsMI)
          #-- Warning for < 15 targets in a cloud
          if (regSummary$targets["DPInet","Positive"] < 15) {
            tp <- paste0("This regulon has less than 15 positive targets; ",
                         "regulon activity readings may be unreliable.\n")
            warning(tp)
          } else if (regSummary$targets["DPInet","Negative"] < 15) {
            tp <- paste("This regulon has less than 15 negative targets; ",
                        "regulon activity readings may be unreliable.\n")
            warning(tp)
          }
          cat("---\n")
        }
      }
      invisible(allRegulonSummary)
    }
  }
)

##------------------------------------------------------------------------------
##show summary information on screen
setMethod(
  "show",
  "TNI",
  function(object) {
    cat("A TNI (Transcriptional Network Inference) object:\n")
    message("--status:")
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    print(tni.get(object, what=c("status")), quote=FALSE)
  }
)

##------------------------------------------------------------------------------
##get slots from TNI 
setMethod(
  "tni.get",
  "TNI",
  function(object, what="summary", order=TRUE, ntop=NULL, reportNames=TRUE, 
           idkey=NULL) {
    
    #---check compatibility
    object <- upgradeTNI(object)
    
    ##-----check input arguments
    tnai.checks(name="tni.what",para=what)
    tnai.checks(name="ntop",para=ntop)
    tnai.checks(name="idkey",para=idkey)
    tnai.checks(name="reportNames",para=reportNames)
    ##-----get query
    query <- NULL
    if(what=="regulonActivity"){
      query <- object@results$regulonActivity
    } else if(what=="subgroupEnrichment"){
      query<-object@results$subgroupEnrichment
    } else if(what=="regulonSize"){
      query <- .regulonCounts(tni.get(object, what="regulons.and.mode"))
      if(!is.null(idkey))
        query<-translateQuery(query,idkey,object,"dataframeAndNames",
                              reportNames)
    } else if(what=="refregulonSize"){
      query <- .regulonCounts(tni.get(object, what="refregulons.and.mode"))
      query<-translateQuery(query,idkey,object,"dataframeAndNames",
                            reportNames)
    } else if(what=="gexp"){
      query<-object@gexp
      if(!is.null(idkey))
        query<-translateQuery(query,idkey,object,"gexpAndNames",reportNames)
    } else if(what=="regulatoryElements"){
      query<-object@regulatoryElements
      if(!is.null(idkey))
        query<-translateQuery(query,idkey,object,"vecAndContent",reportNames)
    } else if(what=="modulators"){
      query<-object@modulators
      if(!is.null(idkey))
        query<-translateQuery(query,idkey,object,"vecAndContent",reportNames)
    } else if(what=="targetElements"){
      query<-object@targetElements
      if(!is.null(idkey))
        query<-translateQuery(query,idkey,object,"vecAndContent",reportNames)
    } else if(what=="para"){
      query<-object@para
    } else if(what=="refnet"){
      query<-object@results$tn.ref
      if(is.null(query))stop("empty slot!",call.=FALSE)
      if(!is.null(idkey))query<-translateQuery(query,idkey,object,
                                               "rtnetAndNames",reportNames)
    } else if(what=="tnet"){
      query<-object@results$tn.dpi
      if(is.null(query))stop("empty slot!",call.=FALSE)
      if(!is.null(idkey))query<-translateQuery(query,idkey,object,
                                               "rtnetAndNames",reportNames)
    } else if(what=="refregulons" || what=="refregulons.and.mode"){
      query<-list()
      for(i in object@regulatoryElements){
        idx<-object@results$tn.ref[,i]!=0
        query[[i]]<-rownames(object@results$tn.ref)[idx]
      }
      if(what=="refregulons.and.mode"){
        for(i in names(query)){
          tp<-object@results$tn.ref[query[[i]],i]
          names(tp)<-query[[i]]
          query[[i]]<-tp
        }
        if(!is.null(idkey))query<-translateQuery(query,idkey,object,
                                                 "listAndNames",reportNames)
      } else {
        if(!is.null(idkey))query<-translateQuery(query,idkey,object,
                                                 "listAndContent",reportNames)
      }
    } else if(what=="regulons" || what=="regulons.and.mode"){
      query<-list()
      for(i in object@regulatoryElements){
        idx<-object@results$tn.dpi[,i]!=0
        query[[i]]<-rownames(object@results$tn.dpi)[idx]
      }
      if(what=="regulons.and.mode"){
        for(i in names(query)){
          tp<-object@results$tn.dpi[query[[i]],i]
          names(tp)<-query[[i]]
          query[[i]]<-tp
        }
        if(!is.null(idkey))query<-translateQuery(query,idkey,object,
                                                 "listAndNames",reportNames)
      } else {
        if(!is.null(idkey))query<-translateQuery(query,idkey,object,
                                                 "listAndContent",reportNames)
      }
    } else if(what=="cdt.list" || what=="cdt.table"){
      query <- object@results$conditional$count
      for(nm in names(query)){
        qry <- query[[nm]]
        qry <- qry[ !qry[,2 ] & !qry[,3 ],-c(2,3),drop=FALSE]
        query[[nm]]<-qry
      }
      if(what=="cdt.table"){
        query <- cdt.table(query)
        query <- p.adjust.cdt.table(query,object@para$cdt$pAdjustMethod)
        query <- query[query$AdjPvFET<= object@para$cdt$pValueCutoff,]
        query <- query[query$AdjPvKS <= object@para$cdt$pValueCutoff,]
        if(is.null(query$AdjPvSNR)){
          idx <- order(query$PvKS, query$PvFET, decreasing = F)
        } else {
          query <- query[query$AdjPvSNR <= object@para$cdt$pValueCutoff,]
          idx <- order(query$PvSNR, query$PvKS, query$PvFET, decreasing = F)
        }
        query <- query[idx,]
        rownames(query)<-NULL
        query <- p.format.cdt.table(query)
        if(reportNames){
          idx <- match(query$Modulator,object@modulators)
          query$Modulator <- names(object@modulators)[idx]
          idx <- match(query$TF,object@regulatoryElements)
          query$TF <- names(object@regulatoryElements)[idx]
        }
      } else {
        query <- cdt.list(query,object@para$cdt$pAdjustMethod)
        for(nm in names(query)){
          qry<-query[[nm]]
          qry <- qry[qry$AdjPvFET<= object@para$cdt$pValueCutoff,]
          qry <- qry[qry$AdjPvKS <= object@para$cdt$pValueCutoff,]
          if(!is.null(qry$AdjPvSNR)){
            qry <- qry[qry$AdjPvSNR <= object@para$cdt$pValueCutoff,]
          }
          query[[nm]]<-qry
        }
        query<-query[unlist(lapply(query,nrow))>0]
        if(reportNames){
          for(nm in names(query)){
            if(nrow(query[[nm]])>0){
              idx<-match(query[[nm]][,"Modulator"],object@modulators)
              query[[nm]][,"Modulator"]<-names(object@modulators)[idx]
              idx<-match(query[[nm]][,"TF"],object@regulatoryElements)
              query[[nm]][,"TF"]<-names(object@regulatoryElements)[idx]
            }
          }
        }
        if(length(query)>0)query<-sortblock.cdt.list(query)
      }
      if(!is.null(ntop))
        warning("'ntop' argument has no effect on this query!")
      if(!is.null(idkey))
        warning("'idkey' argument has no effect on this query!")
    } else if(what=="summary"){
      query<-object@summary
      query$results$tnet <- .get.tnet.summary(object)
      query$results$regulonSize <- .get.regulon.summary(object)
    } else if(what=="rowAnnotation"){
      query<-object@rowAnnotation
    } else if(what=="colAnnotation"){
      query<-object@colAnnotation      
    } else if(what=="status"){
      query<-object@status
    } else if(what=="gsea2"){
      getqs<-function(query,order=TRUE,reportNames=TRUE,ntop=NULL){
        if(is.data.frame(query) && nrow(query)>0 ){
          if(is.null(ntop)){
            query<-query[
              query[,"Adjusted.Pvalue"]<=object@para$gsea2$pValueCutoff,,
              drop=FALSE]
          } else {
            if(ntop>nrow(query)|| ntop<0)ntop=nrow(query)
            if(nrow(query)>1){
              idx<-sort.list(query[,"Pvalue"]) 
              query<-query[idx[1:ntop],,drop=FALSE]
            }
          }
          if(order){
            if(nrow(query)>1) query<-query[order(query[,"Observed.Score"]),,
                                           drop=FALSE]
          }
          if(reportNames){
            idx<-match(query[,1],object@regulatoryElements)
            query[,1]<-names(object@regulatoryElements)[idx]
          }
        }
        query
      }
      query<-list()
      if(is.null(ntop)){
        tp<-rownames(getqs(object@results$GSEA2.results$differential))
        tp<-intersect(tp,rownames(getqs(object@results$GSEA2.results$positive)))
        tp<-intersect(tp,rownames(getqs(object@results$GSEA2.results$negative)))
        dft<-getqs(object@results$GSEA2.results$differential,order,reportNames)
        dft<-dft[rownames(dft)%in%tp,,drop=FALSE]
        query$differential<-dft
        query$positive<-object@results$GSEA2.results$positive[rownames(dft),,
                                                              drop=FALSE]
        query$negative<-object@results$GSEA2.results$negative[rownames(dft),,
                                                              drop=FALSE]
      } else {
        query$differential<-getqs(
          object@results$GSEA2.results$differential,order,reportNames,ntop)
        query$positive<-object@results$GSEA2.results$positive[
          rownames(query$differential),,drop=FALSE]
        query$negative<-object@results$GSEA2.results$negative[
          rownames(query$differential),,drop=FALSE]
      }
    } else if(what=="regulons.and.mode.gmm"){
      query <- .mi2gmm.dpi(object, idkey)
    } else if(what=="refregulons.and.mode.gmm"){
      query <- .mi2gmm.ref(object, idkey)
    }
    return(query)
  }
)

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##-------------------------TNI INTERNAL FUNCTIONS-------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
#---check compatibility and upgrade tni objects
upgradeTNI <- function(object){
  if(is(object, "TNI")){
    if(.hasSlot(object, "transcriptionFactors") && 
       !.hasSlot(object, "regulatoryElements")){
      object@regulatoryElements <- object@transcriptionFactors
      object@rowAnnotation <- object@annotation
      ID <- colnames(object@gexp)
      object@colAnnotation <- data.frame(ID, row.names = ID, 
                                         stringsAsFactors = FALSE)
    }
    if(length(object@status)!=6){
      status <- rep("[ ]", 1, 6)
      names(status) <- c("Preprocess", "Permutation", 
                         "Bootstrap", "DPI.filter", 
                         "Conditional","Activity")
      status[names(object@status)] <- object@status
      object@status <- status
    }
    if(length(object@rowAnnotation)==0){
      tp <- rownames(object@gexp)
      object@rowAnnotation <- data.frame(ID=tp, row.names=tp, 
                                         stringsAsFactors = FALSE)
    }
    if(length(object@colAnnotation)==0){
      tp <- colnames(object@gexp)
      object@colAnnotation <- data.frame(ID=tp, row.names=tp, 
                                         stringsAsFactors = FALSE)
    }
    if(!.hasSlot(object, "targetElements")){
      object@targetElements <- rownames(object@rowAnnotation)
    }
    if(is.null(object@summary$targetElements)){
      sum.info.targetElements<-matrix(,1,2)
      rownames(sum.info.targetElements)<-"targetElements"
      colnames(sum.info.targetElements)<-c("input","valid") 
      object@summary$targetElements <- sum.info.targetElements
    }
    if(is.null(object@summary$regulatoryElements)){
      sum.info.regulatoryElements<-matrix(,1,2)
      rownames(sum.info.regulatoryElements)<-"regulatoryElements"
      colnames(sum.info.regulatoryElements)<-c("input","valid") 
      object@summary$regulatoryElements <- sum.info.regulatoryElements
    }
    sum.info.results <- object@summary$results
    colnames(sum.info.results$tnet)<-c("regulatoryElements","Targets","Edges")
    object@summary$results <- sum.info.results
  }
  return(object)
}

##------------------------------------------------------------------------------
##This function returns alternative annotations for get.tni/get.tna methods
translateQuery<-function(query,idkey,object,annottype,reportNames){
  rowAnnotation<-object@rowAnnotation
  if(is.null(query))return(query)
  cnames<-colnames(rowAnnotation)
  if(!idkey%in%cnames){
    tp1<-"'NOTE: <idkey> not available! please use one of: "
    tp2<-paste(cnames,collapse=", ")
    stop(tp1,tp2,call.=FALSE)
  }
  # get TF's lab
  tfs<-object@regulatoryElements
  if(!reportNames){
    idx<-tfs%in%rownames(rowAnnotation)
    names(tfs)[idx]<-rowAnnotation[tfs[idx],idkey]
  }
  if(annottype=="dataframeAndNames"){
    n <- ncol(query)
    idx <- match(rownames(query),rownames(rowAnnotation))
    query[[idkey]] <- rowAnnotation[[idkey]][idx]
    query <- query[,c(n+1,1:n),drop=FALSE]
    query <- query[!is.na(rownames(query)),,drop=FALSE]
  } else if(annottype=="gexpAndNames"){
    idx<-rownames(query)%in%rownames(rowAnnotation)
    rownames(query)[idx]<-rowAnnotation[rownames(query)[idx],idkey]
    query <- query[!is.na(rownames(query)),,drop=FALSE]
  } else if(annottype=="rtnetAndNames"){
    idx<-colnames(query)%in%tfs
    colnames(query)[idx]<-names(tfs)[idx]
    idx<-rownames(query)%in%rownames(rowAnnotation)
    rownames(query)[idx]<-rowAnnotation[rownames(query)[idx],idkey]
    query <- query[!is.na(rownames(query)),,drop=FALSE]
  } else if(annottype=="listAndNames"){
    idx<-names(query)%in%tfs
    names(query)[idx]<-names(tfs)[idx]
    query<-lapply(query,function(qry){
      idx<-names(qry)%in%rownames(rowAnnotation)
      names(qry)[idx]<-rowAnnotation[names(qry)[idx],idkey]
      qry<-qry[!is.na(names(qry))]
      qry
    })
  } else if(annottype=="listAndContent"){
    idx<-names(query)%in%tfs
    names(query)[idx]<-names(tfs)[idx]
    query<-lapply(query,function(qry){
      nms<-names(qry)
      idx<-qry%in%rownames(rowAnnotation)
      qry[idx]<-rowAnnotation[qry[idx],idkey]
      names(qry)<-nms
      qry<-qry[!is.na(qry)]
      unique(qry)
    })
  } else if(annottype=="vecAndContent"){
    nms<-names(query)
    idx<-query%in%rownames(rowAnnotation)
    query[idx]<-rowAnnotation[query[idx],idkey]
    query<-query[!is.na(query)]
  } else if(annottype=="vecAndNames"){
    idx<-names(query)%in%rownames(rowAnnotation)
    names(query)[idx]<-rowAnnotation[names(query)[idx],idkey]
    query<-query[!is.na(names(query))]
  }
  return(query)
}



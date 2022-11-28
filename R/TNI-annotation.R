################################################################################
##################### TNI Functional Annotation Methods ########################
################################################################################

##------------------------------------------------------------------------------
setMethod(
  "tni.annotate.samples",
  "TNI",
  function(object, geneSetList, minSetSize = 15, exponent = 1, 
           samples=NULL, verbose = TRUE){
    
    #--- check compatibility
    object <- upgradeTNI(object)
    
    #-- Basic checks
    if(object@status["Preprocess"]!="[x]")
      stop("input 'object' needs preprocessing!")
    if(missing(geneSetList))
      stop("missing 'geneSetList'.")
    tnai.checks("geneSetList", geneSetList)
    tnai.checks("minSetSize", minSetSize)
    tnai.checks("exponent", exponent)
    tnai.checks(name="samples",para=samples)
    tnai.checks("verbose",verbose)
    
    #--- start preprocessing
    if(verbose)cat("-Preprocessing for input data...\n")
    
    #--- check geneSetList
    if(is.null(names(geneSetList)))
      stop("'geneSetList' should be named (unique names)!")
    if(any(duplicated(names(geneSetList))))
      stop("'geneSetList' should have unique names!")
    
    #--- check geneSetList
    rowAnnotation <- tni.get(object, 'rowAnnotation')
    geneSetList <- .preprocess.genesets(object, geneSetList, 
                                        minSetSize, verbose)
    
    ##----- get gexp and set samples
    gexp <- tni.get(object, "gexp")
    if(!is.null(samples)){
      idx <- samples %in% colnames(gexp)
      if(!all(idx)){
        stop("'samples' should list only valid names!")
      }
      samples <- colnames(gexp)[colnames(gexp) %in% samples]
      gexp <- gexp[,samples, drop=FALSE]
    }
    
    if(verbose) cat("--Checking log space... ")
    if(.isUnloggedData(gexp)){
      if(verbose) cat("applying log2 transformation!\n")
      gexp <- .log2transform(gexp)
    } else {
      if(verbose)cat("OK!\n")
    }
    
    results <- .annotate.samples.gsea(geneSetList, gexp, exponent, verbose)
    
    return(results)
    
  }
)

##------------------------------------------------------------------------------
setMethod(
  "tni.annotate.regulons",
  "TNI",
  function(object, geneSetList, sampleSetList = NULL, regulatoryElements = NULL, 
           minSetSize = 15, sizeFilterMethod="posORneg",
           exponent = 1, verbose = TRUE){
    
    #-- Basic checks
    if(object@status["DPI.filter"]!="[x]")
      stop("input 'object' needs dpi analysis!")
    if(is.null(sampleSetList)){
      if(missing(geneSetList)) stop("missing 'geneSetList'.")
      tnai.checks("geneSetList", geneSetList)
    } else {
      tnai.checks("sampleSetList",sampleSetList)
    }
    tnai.checks("regulatoryElements",regulatoryElements)
    tnai.checks("minSetSize", minSetSize)
    tnai.checks("sizeFilterMethod", sizeFilterMethod)
    tnai.checks("exponent", exponent)
    tnai.checks("verbose",verbose)
  
    #--- check compatibility
    object <- upgradeTNI(object)
    
    #--- start preprocessing
    if(verbose)cat("-Preprocessing for input data...\n")
    
    #--- check geneSetList
    if(is.null(sampleSetList)){
      if(is.null(names(geneSetList)))
        stop("'geneSetList' should be named (unique names)!")
      if(any(duplicated(names(geneSetList))))
        stop("'geneSetList' should have unique names!")
      geneSetList <- .preprocess.genesets(object, geneSetList, 
                                          minSetSize, verbose)
    } else {
      if(is.null(names(sampleSetList)))
        stop("'sampleSetList' should be named (unique names)!")
      if(any(duplicated(names(sampleSetList))))
        stop("'sampleSetList' should have unique names!")
      sampleSetList <- .preprocess.sampsets(object, sampleSetList, 
                                            minSetSize, verbose)
    }
    
    #--- check regulatoryElements
    if(!is.null(regulatoryElements)){
      regnames <- tni.get(object, "regulatoryElements")
      if(sum(regulatoryElements%in%regnames) > 
         sum(regulatoryElements%in%names(regnames))){
        regulatoryElements <- regnames[regnames%in%regulatoryElements]
      } else {
        regulatoryElements <- regnames[names(regnames)%in%regulatoryElements]
      }
      if(length(regulatoryElements)==0)
        stop("'regulatoryElements' argument has no valid names!")
    } else {
      regulatoryElements <- tni.get(object, "regulatoryElements")
    }
    listOfRegulonsAndMode <- tni.get(object, 'regulons.and.mode')
    listOfRegulonsAndMode <- listOfRegulonsAndMode[regulatoryElements]
    
    ##-----check regulon size
    regcounts <- .regulonCounts(listOfRegulonsAndMode)
    if(sizeFilterMethod=="posANDneg"){
      idx <- regcounts$Positive >= minSetSize & 
        regcounts$Negative >= minSetSize
    } else if(sizeFilterMethod=="posORneg"){
      idx <- regcounts$Positive >= minSetSize | 
        regcounts$Negative >= minSetSize
    } else {
      idx <- regcounts$Size >= minSetSize
    }
    regulatoryElements <- regulatoryElements[
      regulatoryElements%in%rownames(regcounts)[idx]]
    listOfRegulonsAndMode <- listOfRegulonsAndMode[regulatoryElements]
    
    ##-----stop when no regulon passes the size requirement
    if(length(listOfRegulonsAndMode)==0){
      stop("no regulon passed the 'minSetSize' requirement!")
    }
    
    if(verbose) cat("--Checking log space... ")
    gexp <- tni.get(object, "gexp")
    if(.isUnloggedData(gexp)){
      if(verbose) cat("applying log2 transformation!\n")
      gexp <- .log2transform(gexp)
    } else {
      if(verbose)cat("OK!\n")
    }
    if(is.null(sampleSetList)){
      results <- .annotate.regulons.gsea2.1(listOfRegulonsAndMode, 
                                            geneSetList, gexp, 
                                            regulatoryElements, 
                                            exponent, verbose)
    } else {
      results <- .annotate.regulons.gsea2.2(listOfRegulonsAndMode, 
                                            sampleSetList, gexp, 
                                            regulatoryElements, 
                                            exponent, verbose)
    }
    results <- t(results$differential)
    
    return(results)
    
  }
)

##------------------------------------------------------------------------------
setMethod(
  "tni.overlap.genesets",
  "TNI",
  function(object, geneSetList, regulatoryElements = NULL, 
           minSetSize = 15, sizeFilterMethod="posORneg",
           method = c("HT","JC"), pValueCutoff = 0.05, 
           pAdjustMethod = "BH", verbose = TRUE){
    
    #-- Basic checks
    if(object@status["DPI.filter"]!="[x]")
      stop("input 'object' needs dpi analysis!")
    if(missing(geneSetList))
      stop("missing 'geneSetList'.")
    tnai.checks("geneSetList", geneSetList)
    tnai.checks("regulatoryElements",regulatoryElements)
    method <- match.arg(method)
    tnai.checks("minSetSize", minSetSize)
    tnai.checks("sizeFilterMethod", sizeFilterMethod)
    tnai.checks("pValueCutoff", pValueCutoff)
    tnai.checks("pAdjustMethod", pAdjustMethod)
    tnai.checks("verbose",verbose)
    
    #--- check compatibility
    object <- upgradeTNI(object)
    
    #--- start preprocessing
    if(verbose)cat("-Preprocessing for input data...\n")
    
    #--- check geneSetList
    if(is.null(names(geneSetList)))
      stop("'geneSetList' should be named (unique names)!")
    if(any(duplicated(names(geneSetList))))
      stop("'geneSetList' should have unique names!")
    
    #--- check regulatoryElements
    if(!is.null(regulatoryElements)){
      regnames <- tni.get(object, "regulatoryElements")
      if(sum(regulatoryElements%in%regnames) > 
         sum(regulatoryElements%in%names(regnames))){
        regulatoryElements <- regnames[regnames%in%regulatoryElements]
      } else {
        regulatoryElements <- regnames[names(regnames)%in%regulatoryElements]
      }
      if(length(regulatoryElements)==0)
        stop("'regulatoryElements' argument has no valid names!")
    } else {
      regulatoryElements <- tni.get(object, "regulatoryElements")
    }
    
    #--- check geneSetList
    rowAnnotation <- tni.get(object, 'rowAnnotation')
    geneSetList <- .preprocess.genesets(object, geneSetList, 
                                        minSetSize, verbose)
    listOfRegulonsAndMode <- tni.get(object, 'regulons.and.mode')
    listOfRegulonsAndMode <- listOfRegulonsAndMode[regulatoryElements]
    
    ##-----check regulon size
    regcounts <- .regulonCounts(listOfRegulonsAndMode)
    if(sizeFilterMethod=="posANDneg"){
      idx <- regcounts$Positive >= minSetSize & 
        regcounts$Negative >= minSetSize
    } else if(sizeFilterMethod=="posORneg"){
      idx <- regcounts$Positive >= minSetSize | 
        regcounts$Negative >= minSetSize
    } else {
      idx <- regcounts$Size >= minSetSize
    }
    regulatoryElements <- regulatoryElements[
      regulatoryElements%in%rownames(regcounts)[idx]]
    listOfRegulonsAndMode <- listOfRegulonsAndMode[regulatoryElements]
    
    ##-----stop when no regulon passes the size requirement
    if(length(listOfRegulonsAndMode)==0){
      stop("no regulon passed the 'minSetSize' requirement!")
    }
    
    # Get regulons' targets
    listOfRegulons <- lapply(names(listOfRegulonsAndMode), function(reg){
      names(listOfRegulonsAndMode[[reg]])
    })
    names(listOfRegulons) <- names(regulatoryElements)
    
    if(method=="JC"){
      results <- .annotate.regulons.jaccard(geneSetList, listOfRegulons, 
                                            verbose=verbose)
      
      jmat <- matrix(0, length(geneSetList), length(listOfRegulons), 
                     dimnames = list(names(geneSetList), names(listOfRegulons)))
      jmat[as.matrix(
        results[c("GeneSet1", "GeneSet2")])] <- results[["Jaccard"]]
      results <- t(jmat)
    } else {
      results <- .annotate.regulons.hypergeo(geneSetList, listOfRegulons, 
                                             universe=rownames(rowAnnotation), 
                                             pAdjustMethod=pAdjustMethod,
                                             verbose=verbose)
      results <- results[results[["Adjusted.Pvalue"]]<pValueCutoff,]
    }
    
    return(results)
    
  }
)

##------------------------------------------------------------------------------
.annotate.samples.gsea <- function(geneSetList, gexp, exponent, verbose){
  
  #-----get phenotypes
  phenotypes <- gexp-apply(gexp,1,mean)
  
  #-----reset names to integer values
  for(i in names(geneSetList)){
    gs <- geneSetList[[i]]
    geneSetList[[i]] <- match(gs,rownames(phenotypes))
  }
  rnames_phenotypes <- rownames(phenotypes)
  rownames(phenotypes)<-1:nrow(phenotypes)
  
  ##-----get ranked phenotypes
  phenoranks <- apply(-phenotypes, 2, rank)
  colnames(phenoranks) <- colnames(phenotypes)
  rownames(phenoranks) <- rownames(phenotypes)
  genesets <- names(geneSetList)
  samples <- colnames(phenotypes)
  
  if(verbose)cat("-Performing single-sample GSEA...\n")
  if(verbose)cat("--For", length(genesets), "gene set(s) and",
                 length(samples),'sample(s)...\n')
  if(verbose)pb <- txtProgressBar(style=3)
  ES <- NULL
  for(i in 1:length(samples)){
    res <- .run.tni.gsea1.alternative(
      geneSetList=geneSetList,
      phenotype=phenotypes[, samples[i]],
      phenorank=phenoranks[, samples[i]],
      exponent=exponent,
      alternative="two.sided"
    )
    ES <- rbind(ES,res)
    if(verbose) setTxtProgressBar(pb, i/length(samples))
  }
  if(verbose) close(pb)
  rownames(ES) <- samples
  colnames(ES) <- genesets
  return(ES)
}

##------------------------------------------------------------------------
##This function applies a hypergeometric test over two lists with sets of
##genes, and returns a data frame with overlaping statistics
.annotate.regulons.hypergeo <- function(geneSetList, listOfRegulons, universe,
                                        pAdjustMethod = "BH", verbose = TRUE) {
  
  if(verbose) cat("-Performing overlap analysis...\n")
  if(verbose) cat("--For", length(listOfRegulons), 
                  "regulon(s) and", 
                  length(geneSetList),'gene set(s)...\n')
  
  #--- Get geneset pairs
  gspairs <- expand.grid(GeneSet1=names(geneSetList), 
                         GeneSet2=names(listOfRegulons), 
                         stringsAsFactors=FALSE)
  gspairs <- gspairs[gspairs$GeneSet1!=gspairs$GeneSet2, ]
  
  if(nrow(gspairs) > 0){
    
    #--- Run .hypergeo
    if(verbose) pb <- txtProgressBar(style=3)
    reshg <- t(sapply(1:nrow(gspairs),function(i){
      if(verbose) setTxtProgressBar(pb, i/nrow(gspairs))
      .hypergeo.genesets(geneSetList[[gspairs$GeneSet1[i]]],
                         listOfRegulons[[gspairs$GeneSet2[i]]],
                         universe)
    }))
    if(verbose) close(pb)
    results <- cbind(gspairs, reshg)
    
    ##Adjust pvalues
    adjPvals <- p.adjust(results[, "Pvalue"], method = pAdjustMethod)
    results <- cbind(results, 'Adjusted.Pvalue'=adjPvals)
  } else {
    results <- matrix(NA, nrow=0, ncol=10)
  }
  colnames(results) <- c("GeneSet","Regulon","Universe.Size", 
                         "GeneSet.Size", "Regulon.Size", 
                         "Expected.Overlap", "Observed.Overlap", 
                         "Pvalue", "Adjusted.Pvalue")
  return(results)
}

##------------------------------------------------------------------------
##This function returns a data frame with Jaccard statistics
.annotate.regulons.jaccard <- function(geneSetList, listOfRegulons, 
                                       verbose = TRUE) {
  
  #--- Get geneset pairs
  gspairs <- expand.grid(GeneSet1=names(geneSetList), 
                         GeneSet2=names(listOfRegulons), 
                         stringsAsFactors=FALSE)
  gspairs <- gspairs[gspairs$GeneSet1!=gspairs$GeneSet2, ]
  
  #--- Run .jaccard
  if(verbose) pb <- txtProgressBar(style=3)
  results <- sapply(1:nrow(gspairs),function(i){
    if(verbose) setTxtProgressBar(pb, i/nrow(gspairs))
    .jaccard.genesets(geneSetList[[gspairs$GeneSet1[i]]],
                      listOfRegulons[[gspairs$GeneSet2[i]]])
  })
  if(verbose) close(pb)
  results <- cbind(gspairs, Jaccard=results)
  return(results)
}


##------------------------------------------------------------------------
##This function takes two gene sets (e.g. regulons), a vector 
##containing the size of the gene universe, and compute the number of 
##genes expected to occur in both gene sets, the actual observed overlap, 
##and the pvalue from a hypergeometric test.
.hypergeo.genesets <- function(geneSet1, geneSet2, universe) {
  ##number of genes in universe
  N <- length(universe)			
  ##remove genes from gene set that are not in universe			
  geneSet1 <- intersect(geneSet1, universe) 
  geneSet2 <- intersect(geneSet2, universe)
  ##size of gene set	
  m <- length(geneSet1) 							
  Nm <- N-m	
  ##geneSet2 in gene set
  overlap <- intersect(geneSet1, geneSet2) 	
  ##number of hits between regulons		
  k <- length(overlap) 							
  n <- length(geneSet2)	
  HGTresults <- phyper(k-1, m, Nm, n, lower.tail = FALSE)
  ex <- (n/N)*m
  if(m == 0 | n == 0) HGTresults <- 1
  hyp.vec <- c(N, m, n, ex, k, HGTresults)
  names(hyp.vec) <- c("Universe.Size", "GeneSet1.Size", "GeneSet2.Size", 
                      "Expected.Overlap", "Observed.Overlap", "Pvalue")
  return(hyp.vec)
}

##------------------------------------------------------------------------
.jaccard.genesets <- function(geneSet1, geneSet2){
  length(intersect(geneSet1,geneSet2))/length(union(geneSet1,geneSet2))
}

##------------------------------------------------------------------------------
.annotate.regulons.gsea2.1 <- function(listOfRegulonsAndMode, geneSetList, gexp,
                                       regulatoryElements, exponent, verbose){
 
  #-----get phenotypes
  genesets <- names(geneSetList)
  phenotypes <- sapply(genesets, function(gs){
    gs <- geneSetList[[gs]]
    samps <- apply(gexp[gs,], 2, mean)
    sp1 <- names(samps[samps>median(samps)])
    sp2 <- names(samps[samps<median(samps)])
    apply(gexp[,sp1], 1, mean) - apply(gexp[,sp2], 1, mean)
  })
  
  #-----reset names to integer values
  listOfRegulons <- lapply(listOfRegulonsAndMode, names)
  for(i in names(listOfRegulonsAndMode)){
    reg <- listOfRegulonsAndMode[[i]]
    names(listOfRegulonsAndMode[[i]]) <- match(names(reg),rownames(phenotypes))
  }
  rownames(phenotypes)<-1:nrow(phenotypes)
  
  ##-----get ranked phenotypes
  phenoranks <- apply(-phenotypes, 2, rank)
  colnames(phenoranks) <- colnames(phenotypes)
  rownames(phenoranks) <- rownames(phenotypes)
  
  if(verbose)cat("-Performing two-tailed GSEA...\n")
  if(verbose)cat("--For", length(listOfRegulonsAndMode), "regulon(s) and",
                 length(genesets),'gene set(s)...\n')
  if(verbose)pb <- txtProgressBar(style=3)
  regulonActivity<-list()
  for(i in 1:length(genesets)){
    res <- .run.tni.gsea2.alternative(
      listOfRegulonsAndMode=listOfRegulonsAndMode,
      phenotype=phenotypes[, genesets[i]],
      phenorank=phenoranks[, genesets[i]],
      exponent=exponent,
      alternative="two.sided"
    )
    regulonActivity$differential<-rbind(regulonActivity$differential,
                                        res$differential[regulatoryElements])
    regulonActivity$positive<-rbind(regulonActivity$positive,
                                    res$positive[regulatoryElements])
    regulonActivity$negative<-rbind(regulonActivity$negative,
                                    res$negative[regulatoryElements])
    if(verbose) setTxtProgressBar(pb, i/length(genesets))
  }
  if(verbose) close(pb)
  rownames(regulonActivity$differential) <- genesets
  rownames(regulonActivity$positive) <- genesets
  rownames(regulonActivity$negative) <- genesets
  colnames(regulonActivity$differential)<-names(regulatoryElements)
  colnames(regulonActivity$positive)<-names(regulatoryElements)
  colnames(regulonActivity$negative)<-names(regulatoryElements)
  regulonActivity$regulatoryElements <- regulatoryElements
  return(regulonActivity)
}

##------------------------------------------------------------------------------
.annotate.regulons.gsea2.2 <- function(listOfRegulonsAndMode, sampleSetList, 
                                       gexp, regulatoryElements, 
                                       exponent, verbose){
  
  #-----get phenotypes
  sgroups <- names(sampleSetList)
  phenotypes <- sapply(sgroups, function(spg){
    samps <- sampleSetList[[spg]]
    sp1 <- names(samps[samps==1])
    sp2 <- names(samps[samps==0])
    apply(gexp[,sp1], 1, mean) - apply(gexp[,sp2], 1, mean)
  })
  
  #-----reset names to integer values
  listOfRegulons <- lapply(listOfRegulonsAndMode, names)
  for(i in names(listOfRegulonsAndMode)){
    reg <- listOfRegulonsAndMode[[i]]
    names(listOfRegulonsAndMode[[i]]) <- match(names(reg),rownames(phenotypes))
  }
  rownames(phenotypes)<-1:nrow(phenotypes)
  
  ##-----get ranked phenotypes
  phenoranks <- apply(-phenotypes, 2, rank)
  colnames(phenoranks) <- colnames(phenotypes)
  rownames(phenoranks) <- rownames(phenotypes)
  
  if(verbose)cat("-Performing two-tailed GSEA...\n")
  if(verbose)cat("--For", length(listOfRegulonsAndMode), "regulon(s) and",
                 length(sgroups),'sample set(s)...\n')
  if(verbose)pb <- txtProgressBar(style=3)
  regulonActivity<-list()
  for(i in 1:length(sgroups)){
    res <- .run.tni.gsea2.alternative(
      listOfRegulonsAndMode=listOfRegulonsAndMode,
      phenotype=phenotypes[, sgroups[i]],
      phenorank=phenoranks[, sgroups[i]],
      exponent=exponent,
      alternative="two.sided"
    )
    regulonActivity$differential<-rbind(regulonActivity$differential,
                                        res$differential[regulatoryElements])
    regulonActivity$positive<-rbind(regulonActivity$positive,
                                    res$positive[regulatoryElements])
    regulonActivity$negative<-rbind(regulonActivity$negative,
                                    res$negative[regulatoryElements])
    if(verbose) setTxtProgressBar(pb, i/length(sgroups))
  }
  if(verbose) close(pb)
  rownames(regulonActivity$differential) <- sgroups
  rownames(regulonActivity$positive) <- sgroups
  rownames(regulonActivity$negative) <- sgroups
  colnames(regulonActivity$differential)<-names(regulatoryElements)
  colnames(regulonActivity$positive)<-names(regulatoryElements)
  colnames(regulonActivity$negative)<-names(regulatoryElements)
  regulonActivity$regulatoryElements <- regulatoryElements
  return(regulonActivity)
}

##------------------------------------------------------------------------------
.preprocess.genesets <- function(object, geneSetList, minSetSize, verbose){
  ##----Checking geneSetList
  ids <- unique(unlist(geneSetList, use.names = FALSE))
  rowAnnotation <- tni.get(object, 'rowAnnotation')
  if(verbose) cat("--Checking agreement between 'geneSetList' and regulons...")
  refcol <- sapply(1:ncol(rowAnnotation),function(i){
    sum(ids%in%rowAnnotation[,i],na.rm=TRUE)
  })
  agreement<-max(refcol)/length(ids)*100
  if(verbose)cat(paste(round(agreement,1),"% ! \n",sep=""))
  if(agreement<30){
    idiff<-round(100-agreement,1)
    tp<-paste("NOTE: ",idiff,
              "% of 'geneSetList' IDs not represented in the 'rowAnnotation'!",
              sep="")
    stop(tp,call.=FALSE)
  } else if(agreement<50){
    idiff<-round(100-agreement,1)
    tp<-paste("NOTE: ",idiff,
              "% of 'geneSetList' IDs not represented in the 'rowAnnotation'!",
              sep="")
    warning(tp,call.=FALSE)
  }
  refcol <- which(refcol==max(refcol))[1]
  geneSetList <- lapply(geneSetList, function(gs){
    idx <- rowAnnotation[[refcol]] %in% gs
    rownames(rowAnnotation)[idx]
  })
  sz <- unlist(lapply(geneSetList, length))
  geneSetList <- geneSetList[sz>=minSetSize]
  if(length(geneSetList)==0)
    stop("input sets in the 'geneSetList' contains no useful data!\n", 
         call.=FALSE)
  return(geneSetList)
}

##------------------------------------------------------------------------------
.preprocess.sampsets <- function(object, sampleSetList, minSetSize, verbose){
  ##----Checking sampleSetList
  ids <- unique(unlist(sampleSetList, use.names = FALSE))
  colAnnotation <- tni.get(object, 'colAnnotation')
  if(verbose) cat("--Checking 'sampleSetList'...\n")
  sampleSetList <- lapply(sampleSetList, function(lt){
    if(!is.vector(lt) || !all.binaryValues(lt)){
      stop("'sampleSetList' should list numerical or integer vectors, with '0s' and '1s'.")
    }
    if(is.null(names(lt))){
      stop("Vectors listed in 'sampleSetList' should be named.")
    }
    lt <- lt[!is.na(lt)]
    lt <- lt[lt%in%c(0,1)]
    lt[names(lt)%in%rownames(colAnnotation)]
  })
  sz0 <- unlist(lapply(sampleSetList, function(lt){
    sum(lt==0)
  }))
  sz1 <- unlist(lapply(sampleSetList, function(lt){
    sum(lt==1)
  }))
  sampleSetList <- sampleSetList[sz0>=minSetSize & sz1>=minSetSize]
  if(length(sampleSetList)==0)
    stop("input sets in the 'sampleSetList' contains no useful data!\n", 
         call.=FALSE)
  
  return(sampleSetList)
}



################################################################################
##########################     TNA Class Methods    ############################
################################################################################

##------------------------------------------------------------------------------
##initialization method
setMethod("initialize",
		"TNA",
		function(.Object, transcriptionalNetwork, referenceNetwork, 
             regulatoryElements, phenotype=NULL, hits=NULL) {
		  
		  #---check compatibility
		  .Object <- upgradeTNA(.Object)
		  
			##-----check arguments
			if(missing(transcriptionalNetwork))
			  stop("NOTE: 'transcriptionalNetwork' is missing!",call.=FALSE)
			if(missing(referenceNetwork))referenceNetwork=transcriptionalNetwork
			if(missing(regulatoryElements))
			  stop("NOTE: 'regulatoryElements' is missing!",call.=FALSE)
			tnai.checks(name="transcriptionalNetwork",transcriptionalNetwork)
			tnai.checks(name="referenceNetwork",referenceNetwork)
			tnai.checks(name="phenotype",phenotype)
			tnai.checks(name="hits",hits)
			tnai.checks(name="regulatoryElements",regulatoryElements)    
      if(is.null(phenotype) && is.null(hits)){
       stop("NOTE: either 'phenotype' or 'hits' should be available!",call.=FALSE)
      }
      b1<-sum(!colnames(referenceNetwork)==colnames(transcriptionalNetwork)) > 0
			b2<-sum(!rownames(referenceNetwork)==rownames(transcriptionalNetwork)) > 0
			if(b1 || b2) 
			  stop("NOTE: col and row names in 'referenceNetwork' should match 'transcriptionalNetwork'!",
			       call.=FALSE)
			if(sum(!regulatoryElements%in%colnames(transcriptionalNetwork))>0)
			  stop("NOTE: one or more 'regulatoryElements' missing in the 'transcriptionalNetwork'!",
			       call.=FALSE)
      if(is.null(names(regulatoryElements)))names(regulatoryElements) <- regulatoryElements      
			##-----initialization
			.Object@transcriptionalNetwork<-transcriptionalNetwork
			.Object@referenceNetwork<-referenceNetwork
			.Object@regulatoryElements<-regulatoryElements
			.Object@phenotype<-phenotype
			.Object@hits<-hits
			.Object@rowAnnotation<-data.frame()
			.Object@listOfRegulons<-list()
			.Object@listOfReferenceRegulons<-list()
			.Object@listOfModulators<-list()
			.Object@para<-list()
			.Object@results<-list()
			#######summary info######
			##-----tnet targets
			sum.info.tar<-matrix(,1,3)
			colnames(sum.info.tar)<-c("input","valid","duplicate.removed")
			rownames(sum.info.tar)<-"Tagets"
			##-----regulon collection
			sum.info.rgc<-matrix(,1,3)
			colnames(sum.info.rgc)<-c("input","valid","above.min.size")
			rownames(sum.info.rgc)<-"Regulons"
			##-----gene hits
			sum.info.hits<-matrix(,1,3)
			colnames(sum.info.hits)<-c("input","valid","duplicate.removed")
			rownames(sum.info.hits)<-"Hits"
			##-----gene list
			sum.info.gl<-matrix(,1,3)
			colnames(sum.info.gl)<-c("input","valid","duplicate.removed")
			rownames(sum.info.gl)<-"Genes"
			##-----parameters
			sum.info.para <- list()
			sum.info.para$mra <- matrix(,1,4)
			colnames(sum.info.para$mra) <- c("pValueCutoff", "pAdjustMethod",
			                                 "minRegulonSize", "tnet")
			rownames(sum.info.para$mra) <- "Parameters"
			sum.info.para$overlap <- matrix(,1,4)
			colnames(sum.info.para$overlap) <- c("pValueCutoff", "pAdjustMethod", 
			                                     "minRegulonSize", "tnet")
			rownames(sum.info.para$overlap) <- "Parameters"
			sum.info.para$gsea1 <- matrix(,1,7)
			colnames(sum.info.para$gsea1) <- c("pValueCutoff", "pAdjustMethod", 
			                                   "minRegulonSize", "nPermutations", 
                                        "exponent", "tnet","orderAbsValue") 
			rownames(sum.info.para$gsea1) <- "Parameters"
			sum.info.para$gsea2 <- matrix(,1,6)
			colnames(sum.info.para$gsea2) <- c("pValueCutoff", "pAdjustMethod", 
			                                   "minRegulonSize", "nPermutations", 
			                                  "exponent", "tnet")
			rownames(sum.info.para$gsea2) <- "Parameters"
			##-----results
			sum.info.results<-matrix(,3,1)
			colnames(sum.info.results)<-"TNA"
			rownames(sum.info.results)<-c("MRA", "GSEA1", "GSEA2")
			##-----summary info
			.Object@summary<-list(tar=sum.info.tar, rgc=sum.info.rgc, hts=sum.info.hits, 
			                      gl=sum.info.gl, para=sum.info.para, 
			                      results=sum.info.results)
			##-----status info
			.Object@status<-list()
			.Object@status$preprocess <- rep("[ ]", 1, 3)
			names(.Object@status$preprocess) <- c("integration","phenotype", "hits")
			.Object@status$analysis <- rep("[ ]", 1, 3)
			names(.Object@status$analysis) <- c("MRA", "GSEA1", "GSEA2")
			return(.Object)
		}
)

##------------------------------------------------------------------------------
##get slots from TNA
setMethod(
  "tna.get",
  "TNA",
  function(object, what="summary", order=TRUE, ntop=NULL, 
           reportNames=TRUE, idkey=NULL) {
    
    #---check compatibility
    object <- upgradeTNA(object)
    
    ##-----check input arguments
    tnai.checks(name="tna.what",para=what)
    tnai.checks(name="order",para=order)
    tnai.checks(name="ntop",para=ntop)
    tnai.checks(name="reportNames",para=reportNames)
    tnai.checks(name="idkey",para=idkey)
    ##-----get query
    query<-NULL
    if(what=="regulonSize"){
      query <- .regulonCounts(tna.get(object, what="regulons.and.mode"))  
    } else if(what=="refregulonSize"){
      query <- .regulonCounts(tna.get(object, what="refregulons.and.mode"))
    } else if(what=="tnet"){
      query<-object@transcriptionalNetwork
      if(!is.null(idkey))
        query<-translateQuery(query,idkey,object,"matrixAndNames",reportNames)
    } else if(what=="refnet"){
      query<-object@referenceNetwork
      if(!is.null(idkey))query<-translateQuery(query,idkey,object,"matrixAndNames")
    } else if(what=="regulatoryElements"){
      query<-object@regulatoryElements
      if(!is.null(idkey))query<-translateQuery(query,idkey,object,"vecAndContent")
    } else if(what=="pheno"){
      query<-object@phenotype
      if(!is.null(idkey))warning("'idkey' argument has no effect on phenotype data!")
    } else if(what=="regulons" || what=="regulons.and.pheno"){
      query<-object@listOfRegulons
      if( what=="regulons.and.pheno" && !is.null(object@phenotype) ){
        pheno<-object@phenotype
        for(rg in names(query)){
          tp<-query[[rg]]
          idx<-match(tp,names(pheno))
          idx<-idx[!is.na(idx)]
          tpp<-pheno[idx]
          query[[rg]]<-tpp
        }
        if(!is.null(idkey))
          query<-translateQuery(query,idkey,object,"listAndNames",reportNames)
      } else {
        if(!is.null(idkey))
          query<-translateQuery(query,idkey,object,"listAndContent",reportNames)
      }
    } else if(what=="refregulons" || what=="refregulons.and.pheno"){
      query<-object@listOfReferenceRegulons
      if( what=="refregulons.and.pheno" && !is.null(object@phenotype) ){
        pheno<-object@phenotype
        for(rg in names(query)){
          tp<-query[[rg]]
          idx<-match(tp,names(pheno))
          idx<-idx[!is.na(idx)]
          tpp<-pheno[idx]
          query[[rg]]<-tpp
        }
        if(!is.null(idkey))
          query<-translateQuery(query,idkey,object,"listAndNames",reportNames)
      } else {
        if(!is.null(idkey))
          query<-translateQuery(query,idkey,object,"listAndContent",reportNames)
      }
    } else if(what=="para"){
      query<-object@para
    } else if(what=="mra"){
      query<-object@results$MRA.results
      if(is.data.frame(query) && nrow(query)>0){
        if(is.null(ntop)){
          query<-query[query[,"Adjusted.Pvalue"] <= object@para$mra$pValueCutoff,,
                       drop=FALSE]
          if(order && nrow(query)>1) 
            query<-query[order(query[,"Pvalue"]),,drop=FALSE]
        } else {
          if(ntop>nrow(query) || ntop<0)ntop=nrow(query)
          if(nrow(query)>1){
            idx<-sort.list(query[,"Pvalue"]) 
            query<-query[idx[1:ntop],,drop=FALSE]
          }
        }
        if(reportNames){
          idx<-match(query[,1],object@regulatoryElements)
          query[,1]<-names(object@regulatoryElements)[idx]
        }
      }
    }
    else if(what=="gsea1"){
      query <- object@results$GSEA1.results
      if(is.data.frame(query) && nrow(query)>0 ){
        if(is.null(ntop)){
          query<-query[query[,"Adjusted.Pvalue"] <= object@para$gsea1$pValueCutoff,,
                       drop=FALSE]
          if(order && nrow(query)>1)
            query <- query[with(query, order(Pvalue, -abs(Observed.Score))), ,drop=FALSE]
        } else {
          if(ntop>nrow(query)|| ntop<0)ntop=nrow(query)
          if(order && nrow(query)>1){
            query <-query[with(query, order(Pvalue, -abs(Observed.Score))), ,drop=FALSE]
            query <- query[1:ntop,,drop=FALSE]
          }
        }
        if(reportNames){
          idx<-match(query[,1],object@regulatoryElements)
          query[,1]<-names(object@regulatoryElements)[idx]
        }
      }
      if(!is.null(idkey))warning("'idkey' argument has no effect on this table!")
    } else if(what%in%c("gsea2","gsea2summary")){
      getqs<-function(query,order=TRUE,reportNames=TRUE,ntop=NULL){
        if(is.data.frame(query) && nrow(query)>0 ){
          if(is.null(ntop)){
            query<-query[query[,"Adjusted.Pvalue"] <= object@para$gsea2$pValueCutoff,,
                         drop=FALSE]
            if(order && nrow(query)>1)query<-query[order(query[,"Pvalue"]),,drop=FALSE]
          } else {
            if(ntop>nrow(query) || ntop<0)ntop=nrow(query)
            if(order && nrow(query)>1){
              idx<-sort.list(query[,"Pvalue"]) 
              query<-query[idx[1:ntop],,drop=FALSE]
            }
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
        dft<-getqs(object@results$GSEA2.results$differential,order,reportNames)
        dft<-dft[rownames(dft)%in%tp,,drop=FALSE]
        query$differential<-dft
        query$positive<-object@results$GSEA2.results$positive[
          rownames(dft),,drop=FALSE]
        query$negative<-object@results$GSEA2.results$negative[
          rownames(dft),,drop=FALSE]
      } else {
        query$differential<-getqs(object@results$GSEA2.results$differential,
                                  order,reportNames,ntop)
        query$positive<-object@results$GSEA2.results$positive[
          rownames(query$differential),,drop=FALSE]
        query$negative<-object@results$GSEA2.results$negative[
          rownames(query$differential),,drop=FALSE]
      }
      #simplify results for gsea2
      if(what=="gsea2summary"){
        query$differential$Size.Positive <- query$positive$Regulon.Size
        query$differential$Size.Negative <- query$negative$Regulon.Size
        query <- query$differential
        tp <- c("Regulon","Regulon.Size","Size.Positive","Size.Negative",
                "Observed.Score","Pvalue","Adjusted.Pvalue")
        query <- query[,tp]
      }
      if(!is.null(idkey))
        warning("'idkey' argument has no effect on consolidated tables!")
    } else if(what=="regulons.and.mode"){
      query<-list()
      for(i in names(object@listOfRegulons)){
        tp<-object@transcriptionalNetwork[object@listOfRegulons[[i]],i]
        names(tp)<-object@listOfRegulons[[i]]
        query[[i]]<-tp
      }
      if(!is.null(idkey))
        query<-translateQuery(query,idkey,object,"listAndNames",reportNames)
    } else if(what=="refregulons.and.mode"){
      query<-list()
      for(i in names(object@listOfReferenceRegulons)){
        tp<-object@referenceNetwork[object@listOfReferenceRegulons[[i]],i]
        names(tp)<-object@listOfReferenceRegulons[[i]]
        query[[i]]<-tp
      }
      if(!is.null(idkey))
        query<-translateQuery(query,idkey,object,"listAndNames",reportNames)
    } else if(what=="summary"){
      query<-object@summary
    } else if(what=="rowAnnotation"){
      query<-object@rowAnnotation
    } else if(what=="colAnnotation"){
      query<-object@colAnnotation      
    } else if(what=="status"){
      query<-object@status
    } else if(what=="hits"){
      query<-object@hits
    }
    return(query)
  }
)
##------------------------------------------------------------------------------
##show summary information on screen
setMethod(
  "show",
  "TNA",
  function(object) {
    
    #---check compatibility
    object <- upgradeTNA(object)
    
    status<-tna.get(object, what=c("status"))
    cat("A TNA (transcriptional network analysis) object:\n")
    message("--preprocessing status:")
    print(status$preprocess, quote=FALSE)
    message("--analysis status:")
    print(status$analysis, quote=FALSE)
  }
)

##------------------------------------------------------------------------------
##Master regulator analysis
setMethod(
  "tna.mra",
  "TNA",
  function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15,
           tnet="dpi", tfs=NULL, verbose=TRUE) {
    
    #---check compatibility
    object <- upgradeTNA(object)
    
    if(object@status$preprocess["integration"]!="[x]")
      stop("NOTE: input 'object' needs preprocessing!")
    if(object@status$preprocess["hits"]!="[x]")
      stop("NOTE: input 'hits' is empty and/or needs preprocessing!")
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="tnet",para=tnet)   
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="verbose",para=verbose)
    object@para$mra<-list(pValueCutoff=pValueCutoff, 
                          pAdjustMethod=pAdjustMethod, 
                          minRegulonSize=minRegulonSize, tnet=tnet)
    object@summary$para$mra[1,]<-c(pValueCutoff, pAdjustMethod, 
                                   minRegulonSize, tnet)
    
    ##-----get tnet
    if(tnet=="dpi"){
      rgcs<-object@listOfRegulons
    } else {
      rgcs<-object@listOfReferenceRegulons
    }
    if(!is.null(tfs)){
      if(sum(tfs%in%object@regulatoryElements) > 
         sum(tfs%in%names(object@regulatoryElements) ) ){
        tfs<-object@regulatoryElements[object@regulatoryElements%in%tfs]
      } else {
        tfs<-object@regulatoryElements[names(object@regulatoryElements)%in%tfs]
      }
      if(length(tfs)==0)stop("NOTE: 'tfs' argument has no valid name!")
      rgcs<-rgcs[tfs]
    }
    tnet.universe<-rownames(object@referenceNetwork)
    rgcs<-lapply(rgcs, intersect, y=tnet.universe)
    gs.size <- unlist(lapply(rgcs, length))
    object@summary$rgc[,"above.min.size"]<-sum(gs.size>=minRegulonSize)
    
    ##-----run mra analysis
    if(verbose)cat("-Performing master regulatory analysis...\n")
    if(verbose)cat("--For", object@summary$rgc[,"above.min.size"], "regulons...\n")
    if(object@summary$rgc[,"above.min.size"] > 0){
      MRA.results <- hyperGeoTest4RTN(rgcs, universe=tnet.universe, 
                                      hits=object@hits, 
                                      minGeneSetSize=object@para$mra$minRegulonSize, 
                                      pAdjustMethod=object@para$mra$pAdjustMethod, 
                                      verbose=verbose)
    } else {
      MRA.results <- matrix(, nrow=0, ncol=7)
    }
    colnames(MRA.results) <- c("Universe.Size", "Regulon.Size", "Total.Hits", 
                               "Expected.Hits", "Observed.Hits", "Pvalue", 
                               "Adjusted.Pvalue")
    MRA.results<-data.frame(Regulon=rownames(MRA.results),
                            MRA.results,stringsAsFactors=FALSE)
    #-----format results 
    MRA.results$Expected.Hits<-round(MRA.results$Expected.Hits,2)
    MRA.results$Pvalue<-signif(MRA.results$Pvalue, digits=2)
    MRA.results$Adjusted.Pvalue<-signif(MRA.results$Adjusted.Pvalue, digits=2)
    ##-----add results
    object@results$MRA.results<-MRA.results
    MRA.results<-tna.get(object,what="mra", reportNames=FALSE)
    object@summary$results["MRA",]<-ifelse(is.data.frame(MRA.results),
                                           nrow(MRA.results),0)      
    if(verbose)cat("-Master regulatory analysis complete\n\n")    
    ##-----update status and return results
    object@status$analysis["MRA"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##GSEA1
setMethod(
  "tna.gsea1",
  "TNA",
  function(object, pValueCutoff=0.05, pAdjustMethod="BH",  minRegulonSize=15, 
           sizeFilterMethod="posORneg", nPermutations=1000, exponent=1, tnet="dpi", 
           orderAbsValue=TRUE, tfs=NULL, verbose=TRUE){
    
    #---check compatibility
    object <- upgradeTNA(object)
    
    if(object@status$preprocess["integration"]!="[x]")
      stop("NOTE: input 'object' needs preprocessing!")
    if(object@status$preprocess["phenotype"]!="[x]")
      stop("NOTE: input 'phenotype' is empty and/or needs preprocessing!")
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="sizeFilterMethod",para=sizeFilterMethod)
    tnai.checks(name="nPermutations",para=nPermutations)
    tnai.checks(name="exponent",para=exponent)
    tnai.checks(name="gsea.tnet",para=tnet)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="orderAbsValue",para=orderAbsValue)
    tnai.checks(name="verbose",para=verbose)    
    object@para$gsea1<-list(pValueCutoff=pValueCutoff, 
                            pAdjustMethod=pAdjustMethod, 
                            minRegulonSize=minRegulonSize, 
                            nPermutations=nPermutations, 
                            exponent=exponent,tnet=tnet,
                            orderAbsValue=orderAbsValue)
    object@summary$para$gsea1[1,]<-c(pValueCutoff, pAdjustMethod, 
                                     minRegulonSize, 
                                     nPermutations,exponent,
                                     tnet,orderAbsValue)
    ##-----get regulons
    if(tnet=="ref"){
      rgcs<-object@listOfReferenceRegulons
    } else {
      rgcs<-object@listOfRegulons
    }
    
    ##-----
    if(!is.null(tfs)){
      if(sum(tfs%in%object@regulatoryElements) > 
         sum(tfs%in%names(object@regulatoryElements) ) ){
        tfs<-object@regulatoryElements[object@regulatoryElements%in%tfs]
      } else {
        tfs<-object@regulatoryElements[names(object@regulatoryElements)%in%tfs]
      }
      if(length(tfs)==0)stop("NOTE: 'tfs' argument has no valid name!")
      rgcs<-rgcs[tfs]
    } else {
      tfs <- tna.get(object, what = "regulatoryElements")
    }
    
    ##-----check regulon size
    regulonSize <- tna.get(object, what = "regulonSize")
    if(sizeFilterMethod=="posANDneg"){
      idx <- regulonSize$Positive>=minRegulonSize & regulonSize$Negative>=minRegulonSize 
    } else if(sizeFilterMethod=="posORneg"){
      idx <- regulonSize$Positive>=minRegulonSize | regulonSize$Negative>=minRegulonSize 
    } else {
      idx <- regulonSize$Size >= minRegulonSize
    }
    tfs <- tfs[tfs%in%rownames(regulonSize)[idx]]
    rgcs<-rgcs[tfs]
    
    ##----count above min size (input)
    object@summary$rgc[,"above.min.size"]<-length(tfs)
    
    ##-----filter 'rgcs' by 'minRegulonSize'
    if(length(tfs)<1){ 
      stop("NOTE: no regulon passed the 'minRegulonSize' requirement!")
    }
    
    ##-----remove genes not listed in the phenotype
    for(i in names(rgcs)){
      regs <- rgcs[[i]]
      rgcs[[i]] <- regs[regs %in% names(object@phenotype)]
    }
    
    ##-----run gsea1
    GSEA1.results<-run.gsea1(
      listOfRegulons=rgcs,
      phenotype=object@phenotype,
      pAdjustMethod=object@para$gsea1$pAdjustMethod,
      nPermutations=object@para$gsea1$nPermutations,
      exponent=object@para$gsea1$exponent,
      orderAbsValue=object@para$gsea1$orderAbsValue,
      verbose=verbose
    )
    
    ##-----if orderAbsValue, check rounding problems
    if(orderAbsValue){
      idx<-GSEA1.results$Observed.Score<=0
      GSEA1.results$Pvalue[idx]<-1
      GSEA1.results$Adjusted.Pvalue[idx]<-1
    }
    ##-----add results
    object@results$GSEA1.results<-GSEA1.results     
    GSEA1.results<-tna.get(object,what="gsea1", reportNames=FALSE)
    object@summary$results["GSEA1",]<-ifelse(is.data.frame(GSEA1.results),
                                             nrow(GSEA1.results),0)       
    ##-----update status and return results
    object@status$analysis["GSEA1"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##GSEA2
setMethod(
  "tna.gsea2",
  "TNA",
  function(object, pValueCutoff=0.05, pAdjustMethod="BH", minRegulonSize=15, 
           sizeFilterMethod="posORneg", nPermutations=1000, exponent=1, tnet="dpi", 
           tfs=NULL, verbose=TRUE, doSizeFilter=NULL){
    
    #---check compatibility
    object <- upgradeTNA(object)
    
    if(object@status$preprocess["integration"]!="[x]")
      stop("NOTE: input 'object' needs preprocessing!")
    if(object@status$preprocess["phenotype"]!="[x]")
      stop("NOTE: input 'phenotype' is empty and/or needs preprocessing!")
    ##-----check and assign parameters
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="minRegulonSize",para=minRegulonSize)
    tnai.checks(name="sizeFilterMethod",para=sizeFilterMethod)
    tnai.checks(name="nPermutations",para=nPermutations)
    tnai.checks(name="exponent",para=exponent)
    tnai.checks(name="gsea.tnet",para=tnet)
    tnai.checks(name="tfs",para=tfs)
    tnai.checks(name="verbose",para=verbose)
    object@para$gsea2<-list(pValueCutoff=pValueCutoff, 
                            pAdjustMethod=pAdjustMethod, 
                           minRegulonSize=minRegulonSize, 
                           sizeFilterMethod=sizeFilterMethod,
                           nPermutations=nPermutations, 
                           exponent=exponent,tnet=tnet)
    object@summary$para$gsea2[1,]<-c(pValueCutoff, 
                                     pAdjustMethod, minRegulonSize, 
                                     nPermutations, exponent, tnet)
    
    if(!is.null(doSizeFilter)){
      warning("'doSizeFilter' is deprecated, please use the 'sizeFilterMethod' parameter.")
      tnai.checks(name="doSizeFilter",para=doSizeFilter)
      if(doSizeFilter){
        sizeFilterMethod="posANDneg"
      } else {
        sizeFilterMethod="posORneg"
      }
    }
    
    ##------check phenotype for gsea2
    if(!min(object@phenotype)<0 || !max(object@phenotype)>0){
      warning("NOTE: it is expected 'phenotype' data as differential expression values (e.g. logFC)!")
    }
    ##-----get tnet and regulons
    if(tnet=="ref"){
      tnet <- object@referenceNetwork
      listOfRegulonsAndMode <- tna.get(object,what="refregulons.and.mode")
    } else {
      tnet<-object@transcriptionalNetwork
      listOfRegulonsAndMode <- tna.get(object,what="regulons.and.mode")
    }
    ##-----
    if(!is.null(tfs)){
      if(sum(tfs%in%object@regulatoryElements) > 
         sum(tfs%in%names(object@regulatoryElements) ) ){
        tfs<-object@regulatoryElements[object@regulatoryElements%in%tfs]
      } else {
        tfs<-object@regulatoryElements[names(object@regulatoryElements)%in%tfs]
      }
      if(length(tfs)==0)stop("NOTE: 'tfs' argument has no valid names!")
    } else {
      tfs<-object@regulatoryElements
    }
    listOfRegulonsAndMode<-listOfRegulonsAndMode[tfs]

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
    tfs <- tfs[tfs%in%rownames(regcounts)[idx]]
    listOfRegulonsAndMode <- listOfRegulonsAndMode[tfs]
    
    ##-----stop when no regulon passes the size requirement
    if(length(listOfRegulonsAndMode)==0){
      stop("NOTE: no regulon passed the 'minRegulonSize' requirement!")
    }
    object@summary$rgc[,"above.min.size"] <- length(listOfRegulonsAndMode)
    
    ##-----remove genes not listed in phenotype
    for(i in names(listOfRegulonsAndMode)){
      regs <- listOfRegulonsAndMode[[i]] 
      listOfRegulonsAndMode[[i]] <- regs[names(regs)%in%names(object@phenotype)]
    }
    
    ##-----run gsea2
    GSEA2.results<-run.gsea2(
      listOfRegulonsAndMode=listOfRegulonsAndMode,
      phenotype=object@phenotype,
      pAdjustMethod=object@para$gsea2$pAdjustMethod,
      nPermutations=object@para$gsea2$nPermutations, 
      exponent=object@para$gsea2$exponent,
      verbose=verbose
    )
    
    ##-----add results
    object@results$GSEA2.results<-GSEA2.results
    GSEA2.results<-tna.get(object,what="gsea2", reportNames=FALSE)
    object@summary$results["GSEA2",]<-ifelse(is.data.frame(GSEA2.results),
                                             nrow(GSEA2.results),0)       
    ##-----update status and return results
    object@status$analysis["GSEA2"] <- "[x]"
    if(verbose)cat("-GSEA2 analysis complete \n\n")
    return(object)
  }
)


##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##-------------------------TNA INTERNAL FUNCTIONS-------------------------------
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
#---check compatibility and upgrade tna objects
upgradeTNA <- function(object){
  if(is(object,"TNA")){
    if(.hasSlot(object, "transcriptionFactors")){
      object@regulatoryElements <- object@transcriptionFactors
      object@rowAnnotation <- object@annotation
      IDs <- colnames(object@gexp)
      object@colAnnotation <- data.frame(IDs, row.names = IDs, 
                                         stringsAsFactors = FALSE)
    }
  }
  return(object)
}

##------------------------------------------------------------------------------
##internal pre-processing (input via tni2tna.preprocess)
tna.preprocess<-function(object, phenoIDs=NULL, duplicateRemoverMethod="max", 
                         verbose=TRUE) {
  ##-----data preprocessing
  if(verbose)cat("-Preprocessing for input data...\n")
  ##-----data preprocessing: phenotype
  if(!is.null(object@phenotype))
    object<-pheno.preprocess(object, phenoIDs, duplicateRemoverMethod, verbose)
  ##-----data preprocessing: hits
  if(!is.null(object@hits))
    object<-hits.preprocess(object, phenoIDs, verbose)
  ##-----data preprocessing: integration
  object<-data.integration(object, verbose)
  ##-----update and return preprocessing
  if(verbose)cat("-Preprocessing complete!\n\n")
  object@status$analysis[c("MRA", "Overlap", "GSEA1", "GSEA2")] <- "[ ]"
  return(object)
}
##This function returns a preprocessed tna object
pheno.preprocess<-function(object, phenoIDs, duplicateRemoverMethod, verbose){
  ##-----input summary         
  object@summary$gl[,"input"]<-length(object@phenotype)    
  ##-----check phenoIDs if available
  if(!is.null(phenoIDs)){
    if(verbose)cat("--Mapping 'phenotype' to 'phenoIDs'...\n")
    ids<-phenoIDs[,2]
    names(ids)<-phenoIDs[,1]
    if( !all(names(object@phenotype) %in% names(ids)) ){
      stop("NOTE: all names in 'phenotype' should be available in col1 of 'phenoIDs'!",
           call.=FALSE)
    }
    names(object@phenotype)<-ids[names(object@phenotype)]
  }
  ##-----remove NAs in phenotype
  pheno<-object@phenotype
  idx<-!is.na(pheno) & names(pheno)!="" & !is.na(names(pheno))
  if(any(!idx)){
    if(verbose) cat("--Removing genes without names or values in 'phenotype'...\n")
    pheno<-pheno[idx]
  }
  object@summary$gl[,"valid"]<-length(pheno) #genes with valid values
  if(length(pheno)==0)stop("NOTE: input 'phenotype' contains no useful data!\n",
                           call.=FALSE)
  ##-----duplicate remover in phenotype
  uninames<-unique(names(pheno))
  if(length(names(pheno))>length(uninames)){
    if(verbose) cat("--Removing duplicated genes...\n")
    pheno<-tna.duplicate.remover(phenotype=pheno,method=duplicateRemoverMethod)
  }
  object@summary$gl[,"duplicate.removed"]<-length(pheno)	#after removing duplicates
  if(length(pheno)==0)stop("NOTE: input 'phenotype' contains no useful data!\n",
                           call.=FALSE)
  object@phenotype<-pheno

  ##-----check phenotype names in transcriptionalNetwork
  #if(verbose)cat("--Checking 'transcriptionalNetwork' targets in 'phenotype'...")
  #idx<-rownames(object@transcriptionalNetwork) %in% names(object@phenotype)
  #checkmatch<-sum(idx)/length(idx)
  #if(checkmatch==0){
  #  stop("NOTE: 'no agreement between names in 'transcriptionalNetwork' targets and 'phenotype'!")
  #} else if(checkmatch<0.9){
   # warning(paste("Only",round(checkmatch*100,1),
   #               "% of 'transcriptionalNetwork' targets can be mapped to 'phenotype'!"))
  #} else {
  #  if(verbose)cat(paste(round(checkmatch*100,1),"% agreement! \n"))
  #}
  
  ##-----update and return
  object@status$preprocess["phenotype"] <- "[x]"
  return(object)
}

##------------------------------------------------------------------------------
##This function returns a preprocessed tna object
hits.preprocess<-function(object, phenoIDs, verbose){
  ##-----make sure vector 'hits' is set to character
  object@hits<-as.character(object@hits)
  ##-----input summary         
  object@summary$hts[,"input"]<-length(object@hits)
  ##-----check phenoIDs if available
  if(!is.null(phenoIDs)){
    if(verbose)cat("--Mapping 'hits' to 'phenoIDs'...\n")
    ids<-phenoIDs[,2]
    names(ids)<-phenoIDs[,1]
    if(sum( !(names(object@hits) %in% names(ids)) )>0 ){
      stop("NOTE: all names in 'hits' should be available in 'phenoIDs'!",
           call.=FALSE)
    }
    object@hits<-ids[object@hits]
  }
  ##-----remove duplicated hits
  if(verbose) cat("--Removing duplicated hits...\n")
  object@hits<-unique(object@hits)
  object@hits<-object@hits[!is.na(object@hits)]
  object@summary$hts[,"duplicate.removed"]<-length(object@hits)
  if(length(object@hits)==0)
    stop("NOTE: input 'hits' contains no useful data!\n",call.=FALSE)
  
  ##-----update and return
  object@status$preprocess["hits"] <- "[x]"
  return(object)
}

##------------------------------------------------------------------------------
##This function returns the final preprocessed tna object
data.integration<-function(object, verbose){
  
  ##-----input summary
  object@summary$tar[,"input"]<-nrow(object@transcriptionalNetwork)
  object@summary$rgc[,"input"]<-length(object@regulatoryElements)
  
  ##-----check rowAnnotation if available  
  if(nrow(object@rowAnnotation)>0){
    if(verbose)
      cat("--Mapping 'transcriptionalNetwork' annotation to 'phenotype'...\n")
    annot<-object@rowAnnotation
    #col with possible current refids
    col0<-sapply(1:ncol(annot),function(i){
      sum(rownames(annot)%in%annot[,i],na.rm=TRUE)
    })
    if(max(col0)==nrow(annot)){
      col0 <- which(col0==max(col0))[1]
    } else {
      col0 <- 0
    }
    #col with possible phenoIDs (phenotype and hits)
    phenoIDs<-unique(c(names(object@phenotype),object@hits))
    col1<-sapply(1:ncol(annot),function(i){
      sum(phenoIDs%in%annot[,i],na.rm=TRUE)
    })
    col1<-which(col1==max(col1))[1]
    #col with possible gene symbols
    #..obs. posso chamar SYMBOL aqui sem problema! a entrada e' unica, via objetos TNI,
    #..ja verificados exaustivamente em pipeline anterior.
    col2<-which(colnames(annot)=="SYMBOL")[1]
    col2<-ifelse(!is.na(col2),col2,0)
    col2<-ifelse(col1!=col2,col2,0)
    #other cols
    othercols<-1:ncol(annot)
    othercols<-othercols[!othercols%in%c(col0,col1,col2)]
    if(col0!=col1){
      #set rowAnnotation to correct order
      object@rowAnnotation<-annot[,c(col1,col2,col0,othercols),drop=FALSE]
      #update rownames 
      #(na anotacao somente mais adiante por nao aceitar nomes duplicados)
      ids<-object@rowAnnotation[,1]
      names(ids)<-rownames(object@rowAnnotation)
      rownames(object@referenceNetwork)<-ids[
        rownames(object@referenceNetwork)]
      rownames(object@transcriptionalNetwork)<-ids[
        rownames(object@transcriptionalNetwork)]
      #update TFs and colnames
      tfs<-object@regulatoryElements
      coltf<-sapply(1:ncol(object@rowAnnotation),function(i){
        sum(tfs%in%object@rowAnnotation[,i],na.rm=TRUE)
      })
      coltf<-which(coltf==max(coltf))[1]
      idx<-match(tfs,object@rowAnnotation[,coltf])
      tnames<-object@rowAnnotation[idx,1]
      names(tnames)<-names(tfs)
      object@regulatoryElements<-tnames
      #update transcriptionalNetwork colnames
      tfs<-colnames(object@transcriptionalNetwork)
      coltf<-sapply(1:ncol(object@rowAnnotation),function(i){
        sum(tfs%in%object@rowAnnotation[,i],na.rm=TRUE)
      })
      coltf<-which(coltf==max(coltf))[1]
      idx<-match(tfs,object@rowAnnotation[,coltf])
      tnames<-object@rowAnnotation[idx,1]
      names(tnames)<-names(tfs)  
      colnames(object@transcriptionalNetwork)<-tnames
      #update referenceNetwork colnames
      tfs<-colnames(object@referenceNetwork)
      coltf<-sapply(1:ncol(object@rowAnnotation),function(i){
        sum(tfs%in%object@rowAnnotation[,i],na.rm=TRUE)
      })
      coltf<-which(coltf==max(coltf))[1]      
      idx<-match(tfs,object@rowAnnotation[,coltf])
      tnames<-object@rowAnnotation[idx,1]
      names(tnames)<-names(tfs)  
      colnames(object@referenceNetwork)<-tnames
      ##-----remove unnamed nodes in the tnets
      ##tnet
      tnames<-rownames(object@transcriptionalNetwork)
      tnames<-!is.na(tnames) & tnames!=""
      object@transcriptionalNetwork<-object@transcriptionalNetwork[
        tnames,,drop=FALSE]
      ##refnet
      tnames<-rownames(object@referenceNetwork)
      tnames<-!is.na(tnames) & tnames!=""
      object@referenceNetwork<-object@referenceNetwork[tnames,,drop=FALSE]
      ##rowAnnotation
      tnames<-object@rowAnnotation[,1]
      tnames<-!is.na(tnames) & tnames!=""
      object@rowAnnotation<-object@rowAnnotation[tnames,,drop=FALSE]
      ##-----remove duplicate nodes in tnet
      uninames<-unique(rownames(object@transcriptionalNetwork))
      if(length(rownames(object@transcriptionalNetwork))>length(uninames)){
        if(verbose) cat("--Removing duplicated targets...\n")
        abmax<-function(x){
          imax<-max(x);imin<-min(x)
          ifelse(imax>abs(imin),imax,imin)
        }
        tnet<-object@transcriptionalNetwork
        idx<-rownames(tnet)[duplicated(rownames(tnet))]
        idx<-which(rownames(tnet)%in%idx)
        tnetdp<-tnet[idx,,drop=FALSE];tnet<-tnet[-idx,,drop=FALSE]
        #---
        tnetdp<-aggregate(tnetdp,by=list(rownames(tnetdp)),abmax, simplify=TRUE)
        rownames(tnetdp)<-tnetdp[,1]
        tnetdp<-as.matrix(tnetdp[,-1,drop=FALSE])
        #---
        tnetdp<-tnetdp[,object@regulatoryElements,drop=FALSE]
        tnet<-tnet[,object@regulatoryElements,drop=FALSE]
        tnet<-rbind(tnet,tnetdp)
        #---
        object@transcriptionalNetwork<-tnet    
      }
      object@summary$tar[,"duplicate.removed"]<-nrow(object@transcriptionalNetwork)
      if(prod(dim(object@transcriptionalNetwork))==0)
        stop("NOTE: input 'transcriptionalNetwork' contains no useful data!\n",
             call.=FALSE)
      ##-----duplicate remover in refnet
      uninames<-unique(rownames(object@referenceNetwork))
      if(length(rownames(object@referenceNetwork))>length(uninames)){
        abmax<-function(x){
          imax<-max(x);imin<-min(x)
          ifelse(imax>abs(imin),imax,imin)
        }
        tnet<-object@referenceNetwork
        idx<-rownames(tnet)[duplicated(rownames(tnet))]
        idx<-which(rownames(tnet)%in%idx)
        tnetdp<-tnet[idx,,drop=FALSE];tnet<-tnet[-idx,,drop=FALSE]
        #---
        tnetdp<-aggregate(tnetdp,by=list(rownames(tnetdp)),abmax, simplify=TRUE)
        rownames(tnetdp)<-tnetdp[,1]
        tnetdp<-as.matrix(tnetdp[,-1,drop=FALSE])
        #---        
        tnetdp<-tnetdp[,object@regulatoryElements,drop=FALSE]
        tnet<-tnet[,object@regulatoryElements,drop=FALSE]
        tnet<-rbind(tnet,tnetdp)
        #---
        object@referenceNetwork<-tnet    
      }
      ##-----update modulator list if available
      if(length(object@listOfModulators)>0){
        tfs<-names(object@listOfModulators)
        coltf<-sapply(1:ncol(object@rowAnnotation),function(i){
          sum(tfs%in%object@rowAnnotation[,i],na.rm=TRUE)
        })
        coltf<-which(coltf==max(coltf))[1]    
        idx<-match(tfs,object@rowAnnotation[,coltf])
        tnames<-object@rowAnnotation[idx,1] 
        names(object@listOfModulators)<-tnames
        lmod<-sapply(object@listOfModulators,function(reg){
          if(length(reg)>0){
            idx<-match(names(reg),object@rowAnnotation[,coltf])
            mnames<-object@rowAnnotation[idx,1]
            mnames<-aggregate(reg,by=list(mnames),max, simplify=TRUE)
            reg<-mnames[,2]
            names(reg)<-mnames[,1]
          }
          reg
        })
        object@listOfModulators<-lmod
      }
      ##-----duplicate remover in rowAnnotation
      #..get current refids col
      col0<-sapply(1:ncol(object@rowAnnotation),function(i){
        sum(rownames(object@rowAnnotation)%in%object@rowAnnotation[,i],
            na.rm=TRUE)
      })
      if(max(col0)==nrow(object@rowAnnotation)){
        col0 <- which(col0==max(col0))[1]
      } else {
        col0<-0
      }
      uninames<-unique(object@rowAnnotation[,1])
      idx<-match(uninames,object@rowAnnotation[,1])
      object@rowAnnotation<-object@rowAnnotation[idx,,drop=FALSE]
      rownames(object@rowAnnotation)<-object@rowAnnotation[,1]
      object@rowAnnotation<-object@rowAnnotation[,-col0,drop=FALSE] #agora ok p/ retirar!
      ##-----check ordering
      tp1<-rownames(object@rowAnnotation)
      tp2<-rownames(object@transcriptionalNetwork)
      tp3<-rownames(object@referenceNetwork)
      b1<-all(tp1%in%tp2) & all(tp1%in%tp3)
      b2<-length(tp1)==length(tp2) && length(tp1)==length(tp3)
      if(b1 && b2){
        object@transcriptionalNetwork<-object@transcriptionalNetwork[
          rownames(object@rowAnnotation),,drop=FALSE]
        object@referenceNetwork<-object@referenceNetwork[
          rownames(object@rowAnnotation),,drop=FALSE]
      } else {
        warning("NOTE: possible mismatched names between 'transcriptionalNetwork' and 'rowAnnotation'!",
                call.=FALSE)
      }
    }
  }
  ##update
  object@summary$tar[,"valid"]<-nrow(object@transcriptionalNetwork)
  if(nrow(object@transcriptionalNetwork)==0)
    stop("NOTE: input 'transcriptionalNetwork' contains no useful data!\n",
         call.=FALSE)
  
  #-----check phenotype names in transcriptionalNetwork
  if(!is.null(object@phenotype)){
    if(verbose)cat("--Checking agreement between 'transcriptionalNetwork' and 'phenotype'... ")
    phenoIDs<-unique(c(names(object@phenotype),object@hits))
    idx<-rownames(object@transcriptionalNetwork) %in% phenoIDs
    agreement<-sum(idx)/length(idx)*100
    if(verbose)cat(paste(round(agreement,1),"% ! \n",sep=""))
    if(agreement<30){
      idiff<-round(100-agreement,1)
      tp<-paste("NOTE: ",idiff,
                "% of 'transcriptionalNetwork' targets not represented in the 'phenotype'!",
                sep="")
      stop(tp,call.=FALSE)
    } else if(agreement<90){
      idiff<-round(100-agreement,1)
      tp<-paste("NOTE: ",idiff,
                "% of 'transcriptionalNetwork' targets not represented in the 'phenotype'!",
                sep="")
      warning(tp,call.=FALSE)
    }
  }
  
  #-----check and remove hits not listed in the universe
  if(!is.null(object@hits)){
    hits.int<-intersect(object@hits,rownames(object@referenceNetwork))
    if(length(hits.int)<length(object@hits)){
      if(verbose) 
        cat("--Removing 'hits' not listed in 'transcriptionalNetwork' universe...\n")
      object@hits<-hits.int
    }
    object@summary$hts[,"valid"]<-length(object@hits)
    if(length(object@hits)==0)
      stop("NOTE: input 'hits' contains no useful data!\n", call.=FALSE)
  }
  
  ##-----extracting regulons from 'transcriptionalNetwork' 
  if(verbose) cat("--Extracting regulons...\n")
  #Regulons from tnet
  listOfRegulons<-list()
  for(i in object@regulatoryElements){
    idx<-object@transcriptionalNetwork[,i]!=0
    listOfRegulons[[i]]<-rownames(object@transcriptionalNetwork)[idx]
  }
  object@listOfRegulons<-listOfRegulons
  
  #Regulons refnet
  listOfReferenceRegulons<-list()
  for(i in object@regulatoryElements){
    idx<-object@referenceNetwork[,i]!=0
    listOfReferenceRegulons[[i]]<-rownames(object@referenceNetwork)[idx]
  }
  object@listOfReferenceRegulons<-listOfReferenceRegulons
  if(length(object@listOfRegulons)==0)
    stop("NOTE: derived 'listOfRegulons' contains no useful data!\n",
         call.=FALSE)
  
  ##-----update and return
  object@status$preprocess["integration"] <- "[x]"
  return(object)
}

##------------------------------------------------------------------------------
##This function takes a list of gene sets (regulons), a named phenotype 
##vector, and returns the results of gene set enrichment analysis for all 
##regulons (with multiple hypothesis testing correction).
run.gsea1 <- function(listOfRegulons, phenotype, pAdjustMethod="BH", 
                      nPermutations=1000, exponent=1, 
                      orderAbsValue=TRUE, verbose=TRUE) {
  
  ##-----get ordered phenotype
  if(orderAbsValue) phenotype <- abs(phenotype)
  phenotype <- phenotype[order(phenotype,decreasing=TRUE)]
  
  ##-----reset names to integer values
  listOfRegulons<-lapply(listOfRegulons, match, table=names(phenotype))
  names(phenotype)<-1:length(phenotype)
  
  
  ##-----calculate enrichment scores for all regulons
  test.collection<-list()
  if(length(listOfRegulons) > 0){
    test.collection <- gsea1tna(listOfRegulons=listOfRegulons,phenotype=phenotype,
                                exponent=exponent,nPermutations=nPermutations,
                                verbose=verbose)
  } else {
    test.collection<-list(Observed.scores=NULL, Permutation.scores=NULL)
  }
  
  ##-----compute pvals
  if(length(test.collection$Observed.scores) > 0) {
    gs.size <- unlist(lapply(listOfRegulons, length))
    test.pvalues.collection <- tna.perm.pvalues(
      permScores = test.collection$Permutation.scores,
      observedScores = test.collection$Observed.scores)  	
    gsea.adjust.pval <- p.adjust(test.pvalues.collection,
                                 method = pAdjustMethod)
    GSEA1.results <- cbind(gs.size, test.collection$Observed.scores,
                          test.pvalues.collection, gsea.adjust.pval)			
    colnames(GSEA1.results) <- c("Regulon.Size", "Observed.Score", 
                                 "Pvalue", "Adjusted.Pvalue")
  } else {
    GSEA1.results <- matrix(, nrow=0, ncol=4)
    colnames(GSEA1.results) <- c("Regulon.Size", "Observed.Score", 
                                 "Pvalue", "Adjusted.Pvalue")
  }
  
  #-----format results 
  GSEA1.results[,"Observed.Score"]<-round(GSEA1.results[,"Observed.Score"],2)
  GSEA1.results[,"Pvalue"]<-signif(GSEA1.results[,"Pvalue"], digits=5)
  GSEA1.results[,"Adjusted.Pvalue"]<-signif(GSEA1.results[,"Adjusted.Pvalue"], 
                                            digits=5)
  GSEA1.results<-data.frame(Regulon=rownames(GSEA1.results),GSEA1.results,
                            stringsAsFactors=FALSE)
  if(verbose)cat("-Gene set enrichment analysis complete \n\n")
  return( GSEA1.results ) 
}

##------------------------------------------------------------------------------
##This function takes a list of gene sets (regulons), a named phenotype 
##vector, and returns the results of gene set enrichment analysis for all 
##regulons (with multiple hypothesis testing correction).
# fgseaRes <- fgsea(pathways=listOfRegulonsAndMode, nperm=10000, stats=object@phenotype, gseaParam=1)
run.gsea2 <- function(listOfRegulonsAndMode, phenotype, pAdjustMethod="BH", 
                      nPermutations=1000, exponent=1, verbose=TRUE) {
  
  ##-----get ordered phenotype
  phenotype<-phenotype[order(phenotype,decreasing=TRUE)]
  
  ##---reset names to integer values
  for(i in names(listOfRegulonsAndMode)){
    reg<-listOfRegulonsAndMode[[i]]
    names(listOfRegulonsAndMode[[i]]) <- match(names(reg),names(phenotype))
  }
  names(phenotype) <- 1:length(phenotype)
  
  ##-----calculate enrichment scores for all regulons
  test.collection.up <- list()
  test.collection.down <- list()
  if(length(listOfRegulonsAndMode) > 0){
    listOfRegulonsUp <- lapply(names(listOfRegulonsAndMode), function(reg){
      tp<-listOfRegulonsAndMode[[reg]]
      names(tp[tp>0])
    })
    names(listOfRegulonsUp) <- names(listOfRegulonsAndMode)
    listOfRegulonsDown <- lapply(names(listOfRegulonsAndMode), function(reg){
      tp<-listOfRegulonsAndMode[[reg]]
      names(tp[tp<0])
    })
    names(listOfRegulonsDown) <- names(listOfRegulonsAndMode)
    test.collection.up <- gsea2tna(listOfRegulons=listOfRegulonsUp, 
                                   phenotype=phenotype,
                                   exponent=exponent, 
                                   nPermutations=nPermutations, 
                                   verbose1=verbose, verbose2=TRUE)
    test.collection.down <- gsea2tna(listOfRegulons=listOfRegulonsDown, 
                                     phenotype=phenotype,
                                     exponent=exponent, 
                                     nPermutations=nPermutations, 
                                     verbose1=verbose, verbose2=FALSE)
  } else {
    test.collection.up<-list(Observed.scores=NULL, Permutation.scores=NULL)
    test.collection.down<-list(Observed.scores=NULL, Permutation.scores=NULL)
  }
  ##-----compute pvals
  b1<-length(test.collection.up$Observed.scores)>0 && 
    length(test.collection.down$Observed.scores)>0
  b2<-length(test.collection.up$Observed.scores) == 
    length(test.collection.down$Observed.scores)
  if(b1 && b2) {
    gs.size.up <- unlist(lapply(listOfRegulonsUp, length))
    gs.size.down <- unlist(lapply(listOfRegulonsDown, length))
    pvalues.up <- tna.perm.pvalues(permScores = test.collection.up$Permutation.scores,
                                   observedScores = test.collection.up$Observed.scores)    
    adjust.pval.up <- p.adjust(pvalues.up, method = pAdjustMethod, 
                               n=length(pvalues.up)*2)
    GSEA2.results.up <- cbind(gs.size.up, test.collection.up$Observed.scores, 
                              pvalues.up, adjust.pval.up)
    pvalues.down <- tna.perm.pvalues(permScores = test.collection.down$Permutation.scores,
                                     observedScores = test.collection.down$Observed.scores)    
    adjust.pval.down <- p.adjust(pvalues.down, method = pAdjustMethod, 
                                 n=length(pvalues.down)*2)
    GSEA2.results.down <- cbind(gs.size.down,test.collection.down$Observed.scores,
                                pvalues.down, adjust.pval.down)
    test.collection.both<-list(Permutation.scores=test.collection.up$Permutation.scores-
                               test.collection.down$Permutation.scores,
                               Observed.scores=test.collection.up$Observed.scores-
                               test.collection.down$Observed.scores)
    pvalues.both <- tna.perm.pvalues(permScores = test.collection.both$Permutation.scores,
                                     observedScores = test.collection.both$Observed.scores) 
    adjust.pval.both <- p.adjust(pvalues.both, method = pAdjustMethod)
    GSEA2.results.both <- cbind(gs.size.up+gs.size.down,
                                test.collection.both$Observed.scores,
                                pvalues.both, adjust.pval.both)
  } else {
    GSEA2.results.up <- matrix(nrow=0, ncol=4)
    GSEA2.results.down <- matrix(nrow=0, ncol=4)
    GSEA2.results.both <- matrix(nrow=0, ncol=4)
  }
  colnames(GSEA2.results.up) <- c("Regulon.Size", "Observed.Score", 
                                  "Pvalue", "Adjusted.Pvalue")
  colnames(GSEA2.results.down) <- c("Regulon.Size", "Observed.Score", 
                                    "Pvalue", "Adjusted.Pvalue")
  colnames(GSEA2.results.both) <- c("Regulon.Size", "Observed.Score", 
                                    "Pvalue", "Adjusted.Pvalue")
  #-----format results
  GSEA2.results.up[,"Observed.Score"]<-round(GSEA2.results.up[,"Observed.Score"],2)
  GSEA2.results.up[,"Pvalue"]<-signif(GSEA2.results.up[,"Pvalue"], digits=5)
  GSEA2.results.up[,"Adjusted.Pvalue"]<-signif(GSEA2.results.up[,"Adjusted.Pvalue"], 
                                               digits=5)
  GSEA2.results.down[,"Observed.Score"]<-round(GSEA2.results.down[,"Observed.Score"],2)
  GSEA2.results.down[,"Pvalue"]<-signif(GSEA2.results.down[,"Pvalue"], digits=5)
  GSEA2.results.down[,"Adjusted.Pvalue"]<-signif(GSEA2.results.down[,"Adjusted.Pvalue"], 
                                                 digits=5)
  GSEA2.results.both[,"Observed.Score"]<-round(GSEA2.results.both[,"Observed.Score"],2)
  GSEA2.results.both[,"Pvalue"]<-signif(GSEA2.results.both[,"Pvalue"], digits=5)
  GSEA2.results.both[,"Adjusted.Pvalue"]<-signif(GSEA2.results.both[,"Adjusted.Pvalue"], 
                                                 digits=5)
  #-----pack and return
  GSEA2.results.up<-data.frame(Regulon=rownames(GSEA2.results.up),
                               GSEA2.results.up, stringsAsFactors=FALSE)
  GSEA2.results.down<-data.frame(Regulon=rownames(GSEA2.results.down),
                                 GSEA2.results.down, stringsAsFactors=FALSE)
  GSEA2.results.both<-data.frame(Regulon=rownames(GSEA2.results.both),
                                 GSEA2.results.both, stringsAsFactors=FALSE)
  GSEA2.results<-list(positive=GSEA2.results.up,negative=GSEA2.results.down,
                      differential=GSEA2.results.both)
  return( GSEA2.results ) 
}


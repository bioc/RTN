################################################################################
##########################     AVS Class Methods    ############################
################################################################################

##------------------------------------------------------------------------------
##initialization method
setMethod("initialize",
          "AVS",
          function(.Object, markers) {
            ##-----check arguments
            if(missing(markers))stop("NOTE: 'markers' is missing!",call.=FALSE) 
            markers <- .avs.checks(name="markers",markers)
            ##-----initialization
            .Object@validatedMarkers<-markers
            .Object@markers<-markers$rsid
            .Object@variantSet<-list()
            .Object@randomSet<-list()
            .Object@results<-list()
            ##-----status matrix
            .Object@status <- rep("[ ]", 1, 5)
            names(.Object@status) <- c("Preprocess", "VSE", "EVSE", "PEVSE","RVSE")
            ##-----summary info
            ##-----markers
            sum.info.markers<-matrix(,1,4)
            rownames(sum.info.markers)<-"Marker"
            colnames(sum.info.markers)<-c("input","valid","universe.removed","colinked.removed")
            ##-----parameters
            sum.info.para <- list()
            sum.info.para$avs<-matrix(,1,3)
            colnames(sum.info.para$avs)<-c("nrand","reldata","snpop")
            rownames(sum.info.para$avs)<-"Parameter"
            sum.info.para$vse<-matrix(,1,3)
            colnames(sum.info.para$vse)<-c("maxgap","pValueCutoff","pAdjustMethod")
            rownames(sum.info.para$vse)<-"Parameter"
            sum.info.para$evse<-matrix(,1,3)
            colnames(sum.info.para$evse)<-c("maxgap","pValueCutoff","pAdjustMethod")
            rownames(sum.info.para$evse)<-"Parameter"      
            sum.info.para$pevse<-matrix(,1,3)
            colnames(sum.info.para$pevse)<-c("maxgap","pValueCutoff","pAdjustMethod")
            rownames(sum.info.para$pevse)<-"Parameter" 
            sum.info.para$rvse<-matrix(,1,3)
            colnames(sum.info.para$rvse)<-c("maxgap","pValueCutoff","pAdjustMethod")
            rownames(sum.info.para$rvse)<-"Parameter" 
            ##-----results
            sum.info.results<-matrix(,4,1)
            colnames(sum.info.results)<-"Annotation"
            rownames(sum.info.results)<-c("VSE","EVSE","PEVSE","RVSE")
            .Object@summary<-list(markers=sum.info.markers,para=sum.info.para,results=sum.info.results)			
            .Object
          }
)

##------------------------------------------------------------------------------
setMethod(
  "avs.vse",
  "AVS",
  function(object, annotation, glist=NULL, maxgap=0, minSize=100, 
           pValueCutoff=0.05, pAdjustMethod="bonferroni", 
           boxcox=TRUE, verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")stop("NOTE: input data need preprocessing!",call.=FALSE)
    
    #---initial checks
    if(ncol(annotation)<3 && !is.null(glist)){
      stop("'annotation' input should also provide the IDs available the 'glist'! ",call.=FALSE)
    }
    annotation<-tnai.checks(name="annotation.vse",para=annotation)
    tnai.checks(name="maxgap",para=maxgap)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="boxcox",para=boxcox)
    glist<-tnai.checks(name="glist",para=glist)
    tnai.checks(name="minSize",para=minSize)
    tnai.checks(name="verbose",para=verbose)
    object@summary$para$vse[1,]<-c(maxgap,pValueCutoff,pAdjustMethod)
    object@para$vse<-list(maxgap=maxgap,pValueCutoff=pValueCutoff,
                          pAdjustMethod=pAdjustMethod)
    maxgap <- maxgap*1000 #set to bp
    
    #---check glist agreement with annotation
    if(!is.null(glist)){
      gnames<-unique(unlist(glist))
      if(verbose)cat("-Checking agreement between 'glist' and the 'annotation' dataset... ")
      agreement<-sum(gnames%in%annotation$ID)/length(gnames)*100
      if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
      if(agreement<90){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        warning(tp,call.=FALSE)
      } else if(agreement<50){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",sep="")
        stop(tp,call.=FALSE)
      }
      glist<-lapply(glist,intersect,y=annotation$ID)
      gsz<-unlist(lapply(glist,length))
      glist<-glist[gsz>minSize]
      if(length(glist)==0){
        stop("NOTE: no gene set > 'minSize' in the 'glist'!",call.=FALSE)
      }
      #map names to integer values
      annot<-data.table(aid=annotation$ID,ord=1:nrow(annotation))
      setkey(annot,'aid')
      glist<-lapply(glist,function(gl){
        annot[gl,nomatch=0]$ord
      })
    } else {
      glist<-list(annotation$ID)
      names(glist)<-"annotation"
    }
    vSet<-object@variantSet
    rSet<-object@randomSet
    
    #---start vse analysis
    if(isParallel()){
      if(verbose)cat("-Running VSE analysis (parallel version - ProgressBar disabled)...\n")
      getTree=FALSE
    } else {
      if(verbose)cat("-Running VSE analysis...\n")
      getTree=TRUE
    }
    n <- length(glist); labs <- names(glist)
    for(lab in labs){
      if(verbose){
        tp <- paste0("... (",which(labs==lab),"/",n,")")
        cat("--For ",lab, tp,"\n",sep="")
      }
      annot<-getAnnotRanges(annotation[glist[[lab]],],
                            maxgap=maxgap,getTree=getTree,getReduced=TRUE)
      vse<-vsea(vSet,rSet,annot=annot) 
      object@results$vse[[lab]]<-vse
    }
    
    #---map avs to annotation
    annot <- getAnnotRanges(annotation,maxgap=maxgap,getTree=FALSE,
                            getReduced=FALSE)
    annotdist <- getAnnotOverlap1(vSet,annot)
    annotation$OverlapAVS <- FALSE
    annotation[names(annotdist),"OverlapAVS"] <- annotdist
    annotdist <- getAnnotOverlap2(vSet,annot)
    annotation$OverlapCluster <- NA
    annotation[names(annotdist),"OverlapCluster"] <- annotdist
    object@results$annotation$vse <- annotation
    
    #---compute enrichment stats
    object@results$stats$vse<-vseformat(object@results$vse,
                                        pValueCutoff=pValueCutoff,
                                        pAdjustMethod=pAdjustMethod, 
                                        boxcox=boxcox)
    
    #get universe counts (marker and annotation counts)
    # REVISAR: contagem de anotacao nao relevante p/ VSE, talvez seja 
    # desnecessaria quando nao entrar com glist... 
    #revisar correspondente no EVSE!!!
    universeCounts<-getUniverseCounts1(vSet,annotation,maxgap)
    object@results$counts$vse<-universeCounts
    
    ##-----update status and return results
    object@status["VSE"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
setMethod(
  "avs.evse",
  "AVS",
  function(object, annotation, gxdata, snpdata, glist=NULL, maxgap=250, 
           minSize=100, pValueCutoff=0.05, pAdjustMethod="bonferroni", 
           boxcox=TRUE, fineMapping=TRUE, verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")
      stop("NOTE: input data need preprocessing!",call.=FALSE)
    
    #---initial checks
    annotation<-tnai.checks(name="annotation.evse",para=annotation)
    tnai.checks(name="gxdata",para=gxdata)
    tnai.checks(name="maxgap",para=maxgap)
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="boxcox",para=boxcox)
    glist<-tnai.checks(name="glist",para=glist)
    minSize=tnai.checks(name="evse.minSize",para=minSize)
    tnai.checks(name="fineMapping",para=fineMapping)
    tnai.checks(name="verbose",para=verbose)
    object@summary$para$evse[1,]<-c(maxgap,pValueCutoff,pAdjustMethod)
    object@para$evse<-list(maxgap=maxgap,pValueCutoff=pValueCutoff,
                           pAdjustMethod=pAdjustMethod)
    maxgap <- maxgap*1000 #set to bp
    
    #---check gxdata agreement with annotation
    if(verbose)
      cat("-Checking agreement between 'gxdata' and the 'annotation' dataset... ")
    agreement<-sum(rownames(gxdata)%in%annotation$ID)
    agreement<-agreement/nrow(gxdata)*100
    if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
    if(agreement<90){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,
                "% of the ids in 'gxdata' are not represented in the 'annotation' dataset!",
                sep="")
      warning(tp,call.=FALSE)
    } else if(agreement<50){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,
                "% of the ids in 'gxdata' are not represented in the 'annotation' dataset!",
                sep="")
      stop(tp,call.=FALSE)
    }
    annotation<-annotation[annotation$ID%in%rownames(gxdata),,drop=FALSE]
    
    #---check glist agreement with annotation
    if(!is.null(glist)){
      gnames<-unique(unlist(glist))
      if(verbose)
        cat("-Checking agreement between 'glist' and the 'annotation' dataset...  ")
      agreement<-sum(gnames%in%annotation$ID)/length(gnames)*100
      if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
      if(agreement<90){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,
                  "% of the ids in 'glist' are not represented in the 'annotation' dataset!",
                  sep="")
        warning(tp,call.=FALSE)
      } else if(agreement<50){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,
                  "% of the ids in 'glist' are not represented in the 'annotation' dataset!",
                  sep="")
        stop(tp,call.=FALSE)
      }
      glist<-lapply(glist,intersect,y=annotation$ID)
      gsz<-unlist(lapply(glist,length))
      glist<-glist[gsz>minSize[1]]
      if(length(glist)==0){
        stop("NOTE: no gene set > 'minSize' in the 'glist'!",call.=FALSE)
      }
      ##if not fine mapping, get a proxy for the nulls with pre-predefined sizes
      if(!fineMapping){
        gsz<-unlist(lapply(glist,length))
        maxSize<-round(max(gsz)/minSize[2])*minSize[2]
        nproxy<-rev(round(seq(minSize[2],maxSize,by=minSize[2])))
        check<-sapply(gsz,function(gz){
          tp<-abs(nproxy-gz)
          nproxy[which(tp==min(tp))[1]]
        })
        nproxy<-nproxy[nproxy%in%unique(check)]
        names(nproxy)<-1:length(nproxy)
        nproxyids<-sapply(gsz,function(gz){
          tp<-abs(nproxy-gz)
          which(tp==min(tp))[1]
        })
        names(nproxyids)<-names(gsz)
      }
    } else {
      glist<-list(annotation$ID)
      names(glist)<-"annotation"
    }
    
    #---check snpdata matrix
    b1<-!is.matrix(snpdata) && !inherits(snpdata, "ff")
    b2<-!is.integer(snpdata[1,])
    if( b1 && b2){
      stop("'snpdata' should be a matrix (or ff matrix) of integer values!",
           call.=FALSE)
    }
    b1<-is.null(colnames(snpdata)) || is.null(rownames(snpdata))
    b2<-length(unique(rownames(snpdata))) < length(rownames(snpdata))
    b3<-length(unique(colnames(snpdata))) < length(colnames(snpdata))   
    if(  b1 || b2 || b3 ){
      stop("'snpdata' matrix should be named with unique names on rows and cols!",
           call.=FALSE)
    }
    
    #---check gxdata/snpdata matching
    b1<-!all(colnames(gxdata)%in%colnames(snpdata))
    b2<-ncol(gxdata)!=ncol(snpdata)
    if(b1 || b2){
      stop("inconsistent 'gxdata' and 'snpdata' colnames!",call.=FALSE)
    }
    gxdata<-gxdata[,colnames(snpdata)]
    
    #---check avs agreement with snpdata
    if(verbose)
      cat("-Checking agreement between 'AVS' and 'snpdata' datasets... ")
    rMarkers <- avs.get(object,what="randomMarkers")
    lMarkers <- avs.get(object,what="linkedMarkers")
    allMarkers <- unique(c(lMarkers,rMarkers))
    agreement<-sum(allMarkers%in%rownames(snpdata))/length(allMarkers)*100
    if(verbose)cat(paste(round(agreement,digits=1),"% !\n\n",sep=""))
    if(agreement<50){
      idiff<-round(100-agreement,digits=1)
      tp1 <- paste("NOTE: ",idiff,
                   "% of the SNPs in the 'AVS' are not represented in the 'snpdata'!\n",
                   sep="")
      tp2 <- "Although the ideal case would be a perfect matching, it is common\n"
      tp3 <- "to see large GWAS studies interrogating a fraction of the annotated\n"
      tp4 <- "variation. So, given that the annotated variation in the 'AVS' object\n"
      tp5 <- "represents the target population (universe), it is expected\n"
      tp6 <- "a certain level of underepresation of the 'snpdata' in the 'AVS'.\n"
      tp7 <- "Please evaluate whether this number is acceptable for your study.\n"
      warning(tp1,tp2,tp3,tp4,tp5,tp6,tp7,call.=FALSE)
    }
    snpdata <- snpdata[rownames(snpdata)%in%allMarkers,,drop=FALSE]
    if(nrow(snpdata)==0)
      stop("rownames in the 'snpdata' are not compatible with SNPs listed in the 'AVS'")
    
    #---set marker ids to integer in order to improve computational performance 
    if(verbose)cat("-Mapping marker ids to 'snpdata'...\n")
    vSet<-object@variantSet
    rSet<-object@randomSet
    vSet<-mapvset(vSet,snpnames=rownames(snpdata))
    rSet<-maprset(rSet,snpnames=rownames(snpdata),verbose=verbose)
    cat("\n")
    
    #--- start evse analysis
    if(isParallel()){
      getTree=FALSE
      if(fineMapping){
        if(verbose)
          cat("-Running EVSE analysis (parallel version - ProgressBar disabled)...\n")
      } else {
        if(verbose)
          cat("-Running EVSE analysis - pooled null (parallel version)...\n")
      }
    } else {
      getTree=TRUE
      if(fineMapping){
        if(verbose)cat("-Running EVSE analysis...\n")
      } else {
        if(verbose)cat("-Running EVSE analysis - pooled null...\n")
      }
    }
    if(fineMapping){
      n <- length(glist); labs <- names(glist)
      for(lab in labs){
        if(verbose){
          tp <- paste0("... (",which(labs==lab),"/",n,")")
          cat("--For ",lab, tp,"\n",sep="")
        }
        #---run default evsea
        annot <- getAnnotRanges(annotation[glist[[lab]],], maxgap=maxgap,
                              getTree=getTree, getReduced=FALSE)
        evse <- evsea(vSet,rSet,annot,gxdata,snpdata,pValueCutoff,verbose)
        #---get individual eqtls
        annot <- getAnnotRanges(annotation[glist[[lab]],], maxgap=maxgap,
                              getTree=FALSE, getReduced=FALSE)
        eqtls <- eqtlExtract(vSet,annot,gxdata,snpdata,pValueCutoff)
        #---check
        mtally <- names(evse$mtally[evse$mtally])
        bl <- all(unique(eqtls$RiskSNP)%in%mtally)
        if(!bl){warning("...mismatched 'mtally' counts for ", lab, call.=FALSE)}
        evse$eqtls <- eqtls
        object@results$evse[[lab]] <- evse
      }
      # #...for testing only: to update 'eqtls' summary
      # for(lab in labs){
      #   cat("--For ",lab, "\n",sep="")
      #   annot <- getAnnotRanges(annotation[glist[[lab]],], maxgap=maxgap,
      #                           getTree=getTree, getReduced=FALSE)
      #   evse <- object@results$evse[[lab]]
      #   annot <- getAnnotRanges(annotation[glist[[lab]],], maxgap=maxgap,
      #                           getTree=FALSE, getReduced=FALSE)
      #   eqtls <- eqtlExtract(vSet,annot,gxdata,snpdata,pValueCutoff)
      #   mtally <- names(evse$mtally[evse$mtally])
      #   bl <- all(unique(eqtls$RiskSNP)%in%mtally)
      #   if(!bl){warning("...mismatched 'mtally' counts for ", lab, call.=FALSE)}
      #   evse$eqtls <- eqtls
      #   object@results$evse[[lab]] <- evse
      # }
      #---
      object@results$evse.matrix$probs <- getEvseMatrix(object,"P.Multi")
      object@results$evse.matrix$fstat <- getEvseMatrix(object,"F.Multi")
      object@results$evse.matrix$rstat <- getEvseMatrix(object,"R")
    } else {
      #---run evsea null
      nullproxy<-sapply(1:length(nproxy),function(i){
        if(verbose)cat("-- ",i,"/",length(nproxy),"...\n",sep="")
        annot<-sample(1:nrow(annotation),nproxy[i])
        annot<-getAnnotRanges(annotation[annot,],maxgap=maxgap,getTree=getTree,
                              getReduced=FALSE)
        nullproxy<-evseaproxy(rSet,annot,gxdata,snpdata,pValueCutoff,
                              verbose=verbose)
        return(nullproxy)
      })
      #---run evsea
      if(verbose)cat("-- concluding batch processing...\n")
      if(verbose) pb <- txtProgressBar(style=3)
      for (i in 1:length(glist)){
        if(verbose) setTxtProgressBar(pb, i/length(glist))
        lab<-names(glist)[i]
        annot<-getAnnotRanges(annotation[glist[[lab]],],maxgap=maxgap,
                              getTree=getTree,getReduced=FALSE)
        evse<-list()
        evse$mtally<-get.eqtldist.evsea(vSet, annot, gxdata, snpdata, pValueCutoff)
        evse$nulldist<-nullproxy[,nproxyids[lab]]
        evse$nclusters<-length(evse$mtally)
        object@results$evse[[lab]]<-evse
      }
      if(verbose) close(pb)
    }
    
    #---map avs to annotation
    annot<-getAnnotRanges(annotation,maxgap=maxgap, getTree=FALSE, 
                          getReduced=FALSE)
    annotdist<-getAnnotOverlap1(vSet,annot)
    annotation$OverlapAVS <- FALSE
    annotation[names(annotdist),"OverlapAVS"] <- annotdist
    annotdist <- getAnnotOverlap2(vSet,annot)
    annotation$OverlapCluster <- NA
    annotation[names(annotdist),"OverlapCluster"] <- annotdist
    object@results$annotation$evse <- annotation
    
    #---compute enrichment stats
    object@results$stats$evse<-vseformat(object@results$evse,
                                         pValueCutoff=pValueCutoff,
                                         pAdjustMethod=pAdjustMethod, 
                                         boxcox=boxcox)
    
    #get universe counts (marker and gene counts)
    universeCounts<-getUniverseCounts2(vSet,annotation,maxgap)
    object@results$counts$evse<-universeCounts
    
    ##-----update status and return results
    object@status["EVSE"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
setMethod(
  "avs.pevse",
  "AVS",
  function(object, annotation, eqtls, glist, maxgap=250, minSize=100, 
           pValueCutoff=0.05, pAdjustMethod="bonferroni", 
           boxcox=TRUE, verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")
      stop("NOTE: input data need preprocessing!",call.=FALSE)
    
    #---initial checks
    annotation<-tnai.checks(name="annotation.evse",para=annotation)
    eqtls<-tnai.checks(name="eqtls",para=eqtls)
    tnai.checks(name="maxgap",para=maxgap)
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="boxcox",para=boxcox)
    glist<-tnai.checks(name="glist",para=glist)
    tnai.checks(name="minSize",para=minSize)
    tnai.checks(name="verbose",para=verbose)
    object@summary$para$pevse<-matrix(,1,3)
    colnames(object@summary$para$pevse)<-c("maxgap","pValueCutoff","pAdjustMethod")
    rownames(object@summary$para$pevse)<-"Parameter" 
    object@summary$para$pevse[1,]<-c(maxgap,pValueCutoff,pAdjustMethod)
    object@para$pevse<-list(maxgap=maxgap,pValueCutoff=pValueCutoff,
                            pAdjustMethod=pAdjustMethod)
    maxgap <- maxgap*1000 #set to bp
    
    #---check glist agreement with annotation
    if(!is.null(glist)){
      gnames<-unique(unlist(glist))
      if(verbose)
        cat("-Checking agreement between 'glist' and the 'annotation' dataset...  ")
      agreement<-sum(gnames%in%annotation$ID)/length(gnames)*100
      if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
      if(agreement<90){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,"% of the ids in 'glist' are not represented in the 'annotation' dataset!",
                  sep="")
        warning(tp,call.=FALSE)
      } else if(agreement<50){
        idiff<-round(100-agreement,digits=1)
        tp<-paste("NOTE: ",idiff,
                  "% of the ids in 'glist' are not represented in the 'annotation' dataset!",
                  sep="")
        stop(tp,call.=FALSE)
      }
      glist<-lapply(glist,intersect,y=annotation$ID)
      gsz<-unlist(lapply(glist,length))
      glist<-glist[gsz>minSize]
      if(length(glist)==0){
        stop("NOTE: no gene set > 'minSize' in the 'glist'!",call.=FALSE)
      }
    } else {
      glist<-list(annotation$ID)
      names(glist)<-"annotation"
    }
    
    #---check avs agreement with eqtls
    if(verbose)cat("-Checking agreement between 'eqtls' and the 'annotation' dataset...  ")
    gnames <- unique(eqtls$GENEID)
    agreement<-sum(gnames%in%annotation$ID)/length(gnames)*100
    if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
    if(agreement<90){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,
                "% of the ids in 'eqtls' are not represented in the 'annotation' dataset!",
                sep="")
      warning(tp,call.=FALSE)
    } else if(agreement<50){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,
                "% of the ids in 'eqtls' are not represented in the 'annotation' dataset!",
                sep="")
      stop(tp,call.=FALSE)
    }
    eqtls <- paste(eqtls$RSID, eqtls$GENEID, sep="~")
    eqtls <- sort(unique(eqtls))
    vSet<-object@variantSet
    rSet<-object@randomSet
    
    #--- start pevse analysis
    if(isParallel()){
      getTree=FALSE
      if(verbose)
        cat("-Running pEVSE analysis (parallel version - ProgressBar disabled)...\n")
    } else {
      getTree=TRUE
      if(verbose)cat("-Running pEVSE analysis...\n")
    }
    n <- length(glist); labs <- names(glist)
    for(lab in labs){
      if(verbose){
        tp <- paste0("... (",which(labs==lab),"/",n,")")
        cat("--For ",lab, tp,"\n",sep="")
      }
      #---run evsea with pre-defined eQTLs
      annot<-getAnnotRanges(annotation[glist[[lab]],],maxgap=maxgap,
                            getTree=getTree, getReduced=FALSE)
      pevse<-pre_evsea(vSet,rSet,annot,eqtls,verbose=verbose)
      object@results$pevse[[lab]]<-pevse
    }
    
    #---map avs to annotation
    annot<-getAnnotRanges(annotation,maxgap=maxgap, getTree=FALSE, 
                          getReduced=FALSE)
    annotdist<-getAnnotOverlap1(vSet,annot)
    annotation$OverlapAVS <- FALSE
    annotation[names(annotdist),"OverlapAVS"] <- annotdist
    annotdist <- getAnnotOverlap2(vSet,annot)
    annotation$OverlapCluster <- NA
    annotation[names(annotdist),"OverlapCluster"] <- annotdist
    object@results$annotation$pevse <- annotation
    
    #---compute enrichment stats
    object@results$stats$pevse<-vseformat(object@results$pevse,
                                          pValueCutoff=pValueCutoff,
                                          pAdjustMethod=pAdjustMethod, 
                                          boxcox=boxcox)
    
    #get universe counts (marker and gene counts)
    universeCounts<-getUniverseCounts1(vSet,annotation,maxgap)
    object@results$counts$pevse<-universeCounts
    
    ##-----update status and return results
    object@status["PEVSE"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
setMethod(
  "avs.rvse",
  "AVS",
  function(object, annotation, regdata, snpdata, glist, maxgap=250,
           minSize=100, pValueCutoff=0.05, pAdjustMethod="bonferroni", 
           boxcox=TRUE, verbose=TRUE){
    if(object@status["Preprocess"]!="[x]")
      stop("NOTE: input data need preprocessing!",call.=FALSE)
    
    #---initial checks rvse
    annotation<-tnai.checks(name="annotation.evse",para=annotation)
    tnai.checks(name="regdata",para=regdata)
    tnai.checks(name="maxgap",para=maxgap)
    tnai.checks(name="pValueCutoff",para=pValueCutoff)
    tnai.checks(name="pAdjustMethod",para=pAdjustMethod)
    tnai.checks(name="boxcox",para=boxcox)
    glist<-tnai.checks(name="glist",para=glist)
    tnai.checks(name="minSize",para=minSize)
    tnai.checks(name="verbose",para=verbose)
    object@summary$para$rvse<-matrix(,1,3)
    colnames(object@summary$para$rvse)<-c("maxgap","pValueCutoff","pAdjustMethod")
    rownames(object@summary$para$rvse)<-"Parameter" 
    object@summary$para$rvse[1,]<-c(maxgap,pValueCutoff,pAdjustMethod)
    object@para$rvse<-list(maxgap=maxgap,pValueCutoff=pValueCutoff,
                           pAdjustMethod=pAdjustMethod)
    maxgap <- maxgap*1000 #set to bp
    
    #---check glist agreement with annotation
    gnames<-unique(unlist(glist))
    if(verbose)
      cat("-Checking agreement between 'glist' and the 'annotation' dataset...  ")
    agreement<-sum(gnames%in%annotation$ID)/length(gnames)*100
    if(verbose)cat(paste(round(agreement,digits=1),"% !\n",sep=""))
    if(agreement<90){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,
                "% of the ids in 'glist' are not represented in the 'annotation' dataset!",
                sep="")
      warning(tp,call.=FALSE)
    } else if(agreement<50){
      idiff<-round(100-agreement,digits=1)
      tp<-paste("NOTE: ",idiff,
                "% of the ids in 'glist' are not represented in the 'annotation' dataset!",
                sep="")
      stop(tp,call.=FALSE)
    }
    glist<-lapply(glist,intersect,y=annotation$ID)
    gsz<-unlist(lapply(glist,length))
    glist<-glist[gsz>minSize]
    if(length(glist)==0){
      stop("NOTE: no gene set > 'minSize' in the 'glist'!",call.=FALSE)
    }
    
    #---check snpdata matrix
    b1<-!is.matrix(snpdata) && !inherits(snpdata, "ff")
    b2<-!is.integer(snpdata[1,])
    if( b1 && b2){
      stop("'snpdata' should be a matrix (or ff matrix) of integer values!",
           call.=FALSE)
    }
    b1<-is.null(colnames(snpdata)) || is.null(rownames(snpdata))
    b2<-length(unique(rownames(snpdata))) < length(rownames(snpdata))
    b3<-length(unique(colnames(snpdata))) < length(colnames(snpdata))   
    if(  b1 || b2 || b3 ){
      stop("'snpdata' matrix should be named with unique names on rows and cols!",
           call.=FALSE)
    }
    
    #---check regdata/snpdata matching
    b1<-!all(colnames(regdata)%in%colnames(snpdata))
    b2<-ncol(regdata)!=ncol(snpdata)
    if(b1 || b2){
      stop("inconsistent 'regdata' colnames with 'snpdata' colnames!",call.=FALSE)
    }
    regdata<-regdata[,colnames(snpdata)]
    
    #---check regdata/glist matching
    b1<-!all(names(glist)%in%rownames(regdata))
    if(b1){
      stop("inconsistent 'glist' names with 'regdata' rownames!",call.=FALSE)
    }
    regdata<-regdata[names(glist),]
    
    #---check avs agreement with snpdata
    if(verbose)
      cat("-Checking agreement between 'AVS' and 'snpdata' datasets... ")
    rMarkers <- avs.get(object,what="randomMarkers")
    lMarkers <- avs.get(object,what="linkedMarkers")
    allMarkers <- unique(c(lMarkers,rMarkers))
    agreement<-sum(allMarkers%in%rownames(snpdata))/length(allMarkers)*100
    if(verbose)cat(paste(round(agreement,digits=1),"% !\n\n",sep=""))
    if(agreement<50){
      idiff<-round(100-agreement,digits=1)
      tp1 <- paste("NOTE: ",idiff,
                   "% of the SNPs in the 'AVS' are not represented in the 'snpdata'!\n",
                   sep="")
      tp2 <- "Although the ideal case would be a perfect matching, it is common\n"
      tp3 <- "to see large GWAS studies interrogating a fraction of the annotated\n"
      tp4 <- "variation. So, given that the annotated variation in the 'AVS' object\n"
      tp5 <- "represents the target population (universe), it is expected\n"
      tp6 <- "a certain level of underepresation of the 'snpdata' in the 'AVS'.\n"
      tp7 <- "Please evaluate whether this number is acceptable for your study.\n"
      warning(tp1,tp2,tp3,tp4,tp5,tp6,tp7,call.=FALSE)
    }
    snpdata <- snpdata[rownames(snpdata)%in%allMarkers,,drop=FALSE]
    if(nrow(snpdata)==0)
      stop("rownames in the 'snpdata' are not compatible with SNPs listed in the 'AVS'")
    
    #---set marker ids to integer in order to improve computational performance 
    if(verbose)cat("-Mapping marker ids to 'snpdata'...\n")
    vSet<-object@variantSet
    rSet<-object@randomSet
    vSet<-mapvset(vSet,snpnames=rownames(snpdata))
    rSet<-maprset(rSet,snpnames=rownames(snpdata),verbose=verbose)
    cat("\n")
    
    #--- start rvse analysis
    if(isParallel()){
      getTree=FALSE
      if(verbose)
        cat("-Running RVSE analysis (parallel version - ProgressBar disabled)...\n")
    } else {
      getTree=TRUE
      if(verbose)cat("-Running RVSE analysis...\n")
    }
    n <- nrow(regdata); labs <- rownames(regdata)
    for(lab in labs){
      if(verbose){
        tp <- paste0("... (",which(labs==lab),"/",n,")")
        cat("--For ",lab, tp,"\n",sep="")
      }
      annot <- getAnnotRanges(annotation[glist[[lab]],], maxgap=maxgap,
                              getTree=getTree, getReduced=FALSE)
      regact <- regdata[lab,]
      rvse <- rvsea(vSet, rSet, annot, regact, snpdata, pValueCutoff, verbose)
      object@results$rvse[[lab]]<-rvse
    }
    
    #---compute enrichment stats
    object@results$stats$rvse<-vseformat(object@results$rvse,
                                         pValueCutoff=pValueCutoff,
                                         pAdjustMethod=pAdjustMethod, 
                                         boxcox=boxcox)
    
    #get universe counts (marker and gene counts)
    universeCounts<-getUniverseCounts3(vSet)
    object@results$counts$rvse<-universeCounts
    
    ##-----update status and return results
    object@status["RVSE"] <- "[x]"
    return(object)
  }
)

##------------------------------------------------------------------------------
##get slots from AVS 
setMethod(
  "avs.get",
  "AVS",
  function(object, what="summary", report=FALSE, pValueCutoff=NULL){
    ##-----check input arguments
    tnai.checks(name="avs.what",para=what)
    if(!is.null(pValueCutoff))tnai.checks(name="pValueCutoff",para=pValueCutoff)
    ##-----get query
    query<-NULL
    if(what=="markers"){
      query<-object@markers
    } else if(what=="validatedMarkers"){
      query<-object@validatedMarkers      
    } else if(what=="variantSet"){
      if(report){
        query<-report.vset(object@variantSet)
      } else {
        query<-object@variantSet        
      }
    } else if(what=="randomSet"){
      query<-object@randomSet
    } else if(what=="randomMarkers"){
      query <- getMarkers.rset(object@randomSet)
    } else if(what=="linkedMarkers"){
      query <- getMarkers.vset(object@variantSet)
    } else if(what=="evse"){
      if(is.null(pValueCutoff))pValueCutoff<-object@para$evse$pValueCutoff
      if(!is.null(object@results$evse)){
        query<-vseformat(object@results$evse, pValueCutoff=pValueCutoff, 
                         pAdjustMethod=object@para$evse$pAdjustMethod, 
                         boxcox=TRUE)
        if(report)query<-vsereport(query)
      }
    } else if(what=="rvse"){
      if(is.null(pValueCutoff))pValueCutoff<-object@para$rvse$pValueCutoff
      if(!is.null(object@results$rvse)){
        query<-vseformat(object@results$rvse, pValueCutoff=pValueCutoff, 
                         pAdjustMethod=object@para$rvse$pAdjustMethod, 
                         boxcox=TRUE)
        if(report)query<-vsereport(query)
      }
    } else if(what=="pevse"){
      if(is.null(pValueCutoff))pValueCutoff<-object@para$pevse$pValueCutoff
      if(!is.null(object@results$pevse)){
        query<-vseformat(object@results$pevse, pValueCutoff=pValueCutoff, 
                         pAdjustMethod=object@para$pevse$pAdjustMethod, 
                         boxcox=TRUE)
        if(report)query<-vsereport(query)
      }
    } else if(what=="vse"){
      if(is.null(pValueCutoff))pValueCutoff<-object@para$vse$pValueCutoff
      if(!is.null(object@results$vse)){
        query<-vseformat(object@results$vse, pValueCutoff=pValueCutoff, 
                         pAdjustMethod=object@para$vse$pAdjustMethod, 
                         boxcox=TRUE)
        if(report)query<-vsereport(query)
      }
    } else if(what=="summary"){
      query<-object@summary
    } else if(what=="status"){
      query<-object@status
    } else if(what=="annotation.vse"){
      query<-object@results$annotation$vse
    } else if(what=="annotation.evse"){
      query<-object@results$annotation$evse
    } else if(what=="annotation.pevse"){
      query<-object@results$annotation$pevse
    }
    return(query)
  }
)

##------------------------------------------------------------------------------
##show summary information on screen
setMethod(
  "show",
  "AVS",
  function(object) {
    cat("An AVS (Associated Variant Set) object:\n")
    message("--status:")
    print(avs.get(object, what=c("status")), quote=FALSE)
  }
)

##------------------------------------------------------------------------------
##This function is used for argument checking
.avs.checks <- function(name, para){
  if(name=="markers") {
    if( !is.data.frame(para) || ncol(para)<4 ){
      stop("'markers' should be a dataframe, a 'BED file' format with ncol >= 4 !",
           call.=FALSE)
    }
    para<-para[,1:4]
    b1 <- !is.numeric(para[,2]) && !is.integer(para[,2])
    b2 <- !is.numeric(para[,3]) && !is.integer(para[,3])
    if(b1 || b2){
      stop("'markers' should have a 'BED file' format, with chromosomal positions as integer values!",
           call.=FALSE)
    }
    para$start<-as.integer(para$start)
    para$end<-as.integer(para$end)
    colnames(para)<-c("chrom","start","end","rsid")
    if(is.numeric(para$chrom) || is.integer(para$chrom)){
      para$chrom <- paste("chr",para$chrom,sep="")
    }
    para$chrom<-as.character(para$chrom)
    para$rsid<-as.character(para$rsid)
    return(para)
  }
}

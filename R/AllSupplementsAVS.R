

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##                      -- AVS supplements --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##-------------------------------------------------------------------------
## sort annotation
sortAnnotation<-function(annotation){
  #sort annotation
  chrs<-c(paste("chr",1:22,sep=""),"chrX","chrY")
  sortannot<-data.frame(stringsAsFactors=FALSE)
  for(chr in chrs){
    annot<-annotation[annotation$chrom==chr,,drop=FALSE]
    if(nrow(annot)>0){
      idx<-sort.list(annot$start)
      sortannot<-rbind(sortannot,annot[idx,,drop=FALSE])
    }
  }
  rownames(sortannot)<-NULL
  sortannot
}

##------------------------------------------------------------------------
##get markers from computed AVS
getMarkers.vset<-function(variantSet,getlinked=TRUE){
  markers<-NULL
  lkmarkers<-unique(unlist(
    lapply(variantSet,function(vset){
      if(class(vset)=="IRanges"){
        markers<<-c(markers,names(vset@metadata$blocks))
        tp<-lapply(vset@metadata$blocks,names)
      } else {
        markers<<-c(markers,names(vset))
        tp<-lapply(vset,names)
      }
      tp
    })
  ))
  if(getlinked){
    return(lkmarkers)
  } else {
    return(markers)
  }
}

##------------------------------------------------------------------------
##get markers from computed random AVS
getMarkers.rset<-function(randomSet,getlinked=TRUE){
  lkmarkers<-lapply(1:length(randomSet),function(i){
    res<-getMarkers.vset(randomSet[[i]],getlinked)
  })
  unique(unlist(lkmarkers))
}

##-------------------------------------------------------------------------
##map AVS to snpdate (speed-up the permutation step)
mapvset<-function(vSet,snpnames){
  snpnames <- data.table(x=snpnames, y=1:length(snpnames),key="x")
  y=NULL
  for(i in 1:length(vSet)){
    tp<-.mtdata(vSet[[i]])
    mappedMarkers<-list()
    for(j in 1:length(tp$blocks)){
      tpp<-tp$blocks[[j]]
      mapm<-snpnames[data.table(names(tpp))][,y]
      mapm<-mapm[!is.na(mapm)]
      mappedMarkers[[names(tp$blocks[j])]]<-mapm
    }
    vSet[[i]]@metadata$mappedMarkers<-mappedMarkers
  }
  vSet
}
.mtdata<-function(x) {
  if (is.null(x@metadata) || is.character(x@metadata)){
    mdt<-list(metadata = x@metadata)
  } else {
    mdt<-x@metadata
  }
  mdt
}

##-------------------------------------------------------------------------
##map ramdom AVS to snpdate (speed-up the permutation step)
maprset<-function(rSet,snpnames,verbose=TRUE){
  snpnames <- data.table(x=snpnames, y=1:length(snpnames),key="x")
  nr<-length(rSet)
  if(verbose) pb <- txtProgressBar(style=3)
  y=NULL
  resrset<-lapply(1:nr,function(i){
    vSet<-rSet[[i]]
    if(verbose) setTxtProgressBar(pb, i/nr)
    for(i in 1:length(vSet)){
      tp<-.mtdata(vSet[[i]])
      mappedMarkers<-list()
      for(j in 1:length(tp$blocks)){
        tpp<-tp$blocks[[j]]
        mapm<-snpnames[data.table(names(tpp))][,y]
        mapm<-mapm[!is.na(mapm)]
        mappedMarkers[[names(tp$blocks[j])]]<-mapm
      }
      vSet[[i]]@metadata$mappedMarkers<-mappedMarkers
    }
    vSet
  })
  resrset
}

##-------------------------------------------------------------------------
##get IRanges for the AVS
getAvsRanges<-function(vSet){
  clustersRanges<-sapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    blocks<-vSet[[i]]
    clRanges<-IRanges()
    index<-NULL
    for(j in 1:length(blocks)){
      pos<-as.integer(blocks[[j]])
      query<-IRanges(start=pos, end=pos,names=names(blocks[[j]]))
      index<-c(index,rep(j,length(query)))
      clRanges<-c(clRanges,query)
    }
    clRanges@metadata$chr<-chr
    clRanges@metadata$markers<-names(blocks)
    clRanges@metadata$blocks<-blocks
    clRanges@metadata$index<-index
    clRanges
  })
  names(clustersRanges)<-names(vSet)
  return(clustersRanges)
}
##get IRanges for the random AVS
getRandomAvsRanges<-function(rSet, verbose=TRUE){
  if(verbose) pb <- txtProgressBar(style=3)
  clustersRanges<-lapply(1:length(rSet),function(i){
    if(verbose) setTxtProgressBar(pb, i/length(rSet)) 
    getAvsRanges(rSet[[i]])
  })
  if(verbose)close(pb)
  clustersRanges
}
##get IRanges for annotation
getAnnotRanges<-function(annotation,maxgap=0, getTree=TRUE, getReduced=FALSE){
  annot<-annotation
  idx<-annot$START>annot$END
  annot$START[idx]<-annotation$END
  annot$END[idx]<-annotation$START
  annotation<-annot
  chrs<-c(paste("chr",1:22,sep=""),"chrX")
  chrs<-chrs[chrs%in%unique(annotation$CHROM)]
  annotRange<-sapply(chrs,function(chr){
    annot<-annotation[annotation$CHROM==chr,c("ID","START","END")]
    start<-as.integer(annot$START-maxgap)
    start[start<0]<-0
    end<-as.integer(annot$END+maxgap)
    subject<-IRanges(start, end, names=annot$ID)
    if(getReduced)subject<-reduce(subject)
    if(getTree)subject<-NCList(subject)
    subject@metadata<-list(mappedAnnotations=annot$ID)
    subject
  })
  annotRange
}

##------------------------------------------------------------------------
##get markers and genes from computed evse
##return a named character vector with the RiskAssociatedSNPs mapped in the VSE analysis
##..names indicate the LD cluster (i.e. the RiskSNP assignment) 
getMappedClusters<-function(object){
  MarkerIDs<-NULL
  for(vset in object@variantSet){
    tp<-lapply(vset@metadata$blocks,names)
    tp<-unlist(tp)
    names(tp)<-vset@metadata$markers[vset@metadata$index]
    MarkerIDs<-c(MarkerIDs,tp)
  }
  mappedIds<-getMappedMarkers(object)
  MarkerIDs<-MarkerIDs[names(MarkerIDs)%in%mappedIds$RiskSNP]
  MarkerIDs<-MarkerIDs[MarkerIDs%in%mappedIds$RiskAssociatedSNP]
  return(MarkerIDs)
}
##return a list with RiskSNPs and RiskAssociatedSNPs mapped in the VSE analysis
getMappedMarkers<-function(object){
  RiskSNP<-NULL
  RiskAssociatedSNP<-NULL
  for(reg in names(object@results$evse)){
    eqtls<-object@results$evse[[reg]]$eqtls
    RiskSNP<-c(RiskSNP,eqtls$RiskSNP)
    RiskAssociatedSNP<-c(RiskAssociatedSNP,eqtls$RiskAssociatedSNP)
  }
  rsid<-object@validatedMarkers$rsid
  RiskSNP<-rsid[rsid%in%RiskSNP]
  RiskAssociatedSNP<-unique(RiskAssociatedSNP)
  return(list(RiskSNP=RiskSNP,RiskAssociatedSNP=RiskAssociatedSNP))
}

###########################################################################
## VSE analysis
###########################################################################
##-------------------------------------------------------------------------
vsea<-function(vSet,rSet,annot,verbose=TRUE){
  #compute vse
  resvset<-get.avsdist(vSet=vSet,annot=annot)
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("get.avsdist","IRanges","overlapsAny",".mtdata"),
                        envir=environment())
    resrset <- snow::parSapply(cl, 1:length(rSet), function(i) {
      sum(get.avsdist(rSet[[i]],annot))
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    resrset<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet)) 
      sum(get.avsdist(rSet[[i]],annot))
    })
    if(verbose)close(pb)
  }
  return(list(mtally=resvset,nulldist=resrset,nclusters=length(resvset)))
}

##-------------------------------------------------------------------------
##get avs overlap for a variant set
get.avsdist<-function(vSet,annot){
  clusterMapping<-lapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]   
    if(!is.null(subject)){
      res<-rep(FALSE,length(query@metadata$markers))
      ov<-query@metadata$index[overlapsAny(query,subject)]
      res[ov]<-TRUE
    } else {
      res<-rep(FALSE,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    res
  })
  resvset<-unlist(clusterMapping)
  return(resvset)
}

##-------------------------------------------------------------------------
##get annotation overlap
##annot should be named!
getAnnotOverlap1<-function(vSet,annot){
  annotMapping<-lapply(1:length(annot),function(i){
    chr<-names(annot[i])
    query<-annot[[i]]
    subject<-vSet[[chr]]  
    if(!is.null(subject)){
      res<-overlapsAny(query,subject)
    } else {
      res<-rep(FALSE,length(query))
    }
    names(res)<-names(query)
    res
  })
  annotMapping<-unlist(annotMapping)
  return(annotMapping)
}
getAnnotOverlap2<-function(vSet,annot){
  annotMapping <- NULL
  for(i in 1:length(annot)){
    chr<-names(annot[i])
    query<-annot[[i]]
    subject<-vSet[[chr]] 
    if(!is.null(subject)){
      tp <- findOverlaps(query,subject)
      sb <- subject@metadata$index[to(tp)]
      sb <- subject@metadata$markers[sb]
      qr <- names(query)[from(tp)]
      idx <- !duplicated(paste(qr,sb))
      res <- cbind(qr[idx],sb[idx])
      annotMapping <- rbind(annotMapping,res)
    }
  }
  if(!is.null(annotMapping) && nrow(annotMapping)>0){
    ids <- unique(annotMapping[,1])
    annotMapping <- sapply(ids,function(id){
      idx <- annotMapping[,1]==id
      paste(annotMapping[idx,2], collapse = ", ")
    })
  } else {
    annotMapping <- NULL
  }
  return(annotMapping)
}

###########################################################################
## pre-EVSE analysis
###########################################################################
##-------------------------------------------------------------------------
##run evsea for observed and random variant sets
pre_evsea<-function(vSet, rSet, annot, eqtls, verbose=TRUE){
  #compute evse
  mtally<-get.pre_eqtldist(vSet, annot, eqtls)
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("get.pre_eqtldist","IRanges","overlapsAny",
                                 "findOverlaps","%chin%",".mtdata"),
                        envir=environment())
    resrset <- snow::parSapply(cl, 1:length(rSet), function(i) {
      sum(get.pre_eqtldist(rSet[[i]], annot, eqtls))
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    resrset<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet))
      sum(get.pre_eqtldist(rSet[[i]], annot, eqtls))
    })
    if(verbose)close(pb)
  }
  return(list(mtally=mtally,nulldist=resrset,nclusters=length(mtally)))
}

##-------------------------------------------------------------------------
##get avs/eqtl dist for a variant set
get.pre_eqtldist<-function(vSet,annot,eqtls){
  # mapping tally
  clusterMapping<-sapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      overlaps<-findOverlaps(query,subject)
      geneList<-.mtdata(subject)$mappedAnnotations
      res<-sapply(1:length(query@metadata$markers),function(j){
        snpList<-names(.mtdata(query)$blocks[[j]])
        ov<-S4Vectors::from(overlaps)%in%which(query@metadata$index==j)
        if(any(ov) && length(snpList)>0){
          ov<-unique(S4Vectors::to(overlaps)[ov])
          gList<-geneList[ov]
          if(length(gList)>0){
            sg <- paste(snpList, gList, sep="~")
            bl <- any(sg %chin% eqtls)
            res <- ifelse(bl,TRUE,FALSE)
          } else {
            res <- FALSE
          }
        } else {
          res <- FALSE
        }
        return(res)
      })
    } else {
      res<-rep(FALSE,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    return(res)
  })
  clusterMapping<-unlist(clusterMapping)
  clusterMapping
}


###########################################################################
## EVSE analysis
###########################################################################

##-------------------------------------------------------------------------
##run evsea for observed and random variant sets
evsea<-function(vSet, rSet, annot, gxdata, snpdata, pValueCutoff=0.01,verbose=TRUE){
  #compute evse
  mtally<-get.eqtldist(vSet, annot, gxdata, snpdata, pValueCutoff)
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("get.eqtldist","eqtlTest","IRanges","overlapsAny",
                                 "findOverlaps",".mtdata"),
                        envir=environment())
    resrset <- snow::parSapply(cl, 1:length(rSet), function(i) {
      sum(get.eqtldist(rSet[[i]], annot, gxdata, snpdata, pValueCutoff))
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    resrset<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet))
      sum(get.eqtldist(rSet[[i]], annot, gxdata, snpdata, pValueCutoff))
    })
    if(verbose)close(pb)
  }
  return(list(mtally=mtally,nulldist=resrset,nclusters=length(mtally)))
}

##-------------------------------------------------------------------------
##run evsea for observed and random variant sets
evseaproxy<-function(rSet, annot, gxdata, snpdata, pValueCutoff=0.01,verbose=TRUE){
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("get.eqtldist","eqtlTest","IRanges",
                                 "overlapsAny","findOverlaps",".mtdata"),
                        envir=environment())
    nulldist <- snow::parSapply(cl, 1:length(rSet), function(i) {
      res<-get.eqtldist(rSet[[i]], annot, gxdata, snpdata, pValueCutoff)
      sum(res)
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    nulldist<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet))
      res<-get.eqtldist(rSet[[i]], annot, gxdata, snpdata, pValueCutoff)
      sum(res)
    })
    if(verbose)close(pb)
  }
  return(nulldist)
}

##-------------------------------------------------------------------------
##get avs/eqtl dist for a variant set
get.eqtldist<-function(vSet,annot,gxdata,snpdata,pValueCutoff=0.01){
  # mapping tally
  clusterMapping<-sapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      overlaps<-findOverlaps(query,subject)
      geneList<-.mtdata(subject)$mappedAnnotations
      res<-sapply(1:length(query@metadata$markers),function(j){
        snpList<-.mtdata(query)$mappedMarkers[[j]]
        ov<-S4Vectors::from(overlaps)%in%which(query@metadata$index==j)
        if(any(ov) && length(snpList)>0){
          ov<-unique(S4Vectors::to(overlaps)[ov])
          gList<-geneList[ov]
          if(length(gList)>0){
            res<-eqtlTest(geneList=as.character(gList),snpList=as.integer(snpList),gxdata,snpdata)
          } else {
            res<-1.0
          }
        } else {
          res<-1.0
        }
        return(res)
      })
    } else {
      res<-rep(1.0,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    return(res)
  })
  clusterMapping<-unlist(clusterMapping)
  clusterMapping<pValueCutoff
}

##-------------------------------------------------------------------------
##eqtl test
##run two-way manova with multiple additive factors (gx as response varible)
eqtlTest<-function(geneList,snpList,gxdata,snpdata){
  gxdata<-t(gxdata[geneList,,drop=FALSE])
  snpdata<-t(snpdata[snpList,,drop=FALSE])
  # set names for formulae
  colnames(gxdata)<-paste(rep("G",ncol(gxdata)),1:ncol(gxdata),sep="")
  colnames(snpdata)<-paste(rep("S",ncol(snpdata)),1:ncol(snpdata),sep="")
  # run lm
  fm1<-paste(colnames(snpdata),collapse="+")
  fmla <- formula( paste("gxdata ~",fm1, collapse=" ") )
  resfit<-lm(fmla, data=as.data.frame(snpdata,stringsAsFactors=FALSE))
  if(ncol(gxdata)>1){
    resf<-summary(manova(resfit))
    resf<-resf$stats[,"Pr(>F)"]
  } else {
    resf<-anova(resfit) 
    resf<-resf[["Pr(>F)"]]
  }
  if(sum(is.na(resf))==length(resf)){
    res=1
  } else {
    res<-min(resf,na.rm=TRUE)
  }
  return(res)
}

#-------------------------------------------------------------------------
getUniverseCounts1<-function(vSet,annotation,maxgap){
  #count markers in hapmap
  clusterCounts<-unlist(lapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$blocks,length))
  }))
  #count tested genes, overlap
  annot<-getAnnotRanges(annotation,maxgap=maxgap,getTree=FALSE, getReduced=FALSE)
  # mapping tally
  geneCounts<-lapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      geneList<-.mtdata(subject)$mappedAnnotations
      res<-sapply(1:length(query@metadata$markers),function(j){
        gList<-geneList[overlapsAny(subject,query[query@metadata$index==j])]
        length(gList)
      })
    } else {
      res<-rep(0,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    return(res)
  })
  geneCounts<-unlist(geneCounts)
  #merge counts
  counts<-t(rbind(clusterCounts,geneCounts))
  colnames(counts)<-c("markers","annotation")
  counts
}

#-------------------------------------------------------------------------
getUniverseCounts2<-function(vSet,annotation,maxgap){
  #count markers in hapmap
  clusterCounts<-unlist(sapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$blocks,length))
  }))
  #count markers mapped to the genotype data
  clusterCounts<-rbind(clusterCounts,unlist(sapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$mappedMarkers,length))
  })))
  #count tested genes, overlap
  annot<-getAnnotRanges(annotation=annotation,maxgap=maxgap,getTree=FALSE, getReduced=FALSE)
  # mapping tally
  geneCounts<-sapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      geneList<-.mtdata(subject)$mappedAnnotations
      res<-sapply(1:length(query@metadata$markers),function(j){
        gList<-geneList[overlapsAny(subject,query[query@metadata$index==j])]
        length(gList)
      })
    } else {
      res<-rep(0,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    return(res)
  })
  geneCounts<-unlist(geneCounts)
  #merge counts
  counts<-t(rbind(clusterCounts,geneCounts))
  colnames(counts)<-c("totalMarkers","markers","annotation")
  counts
}

#-------------------------------------------------------------------------
vseformat<-function(resavs, pValueCutoff=0.01, 
                    pAdjustMethod="bonferroni", boxcox=TRUE){
  groups<-rep(1,length(resavs)) #'groups' nao ativo!
  ntests<-length(resavs)
  #get mtally
  mtally<-sapply(names(resavs),function(i){
    resavs[[i]]$mtally
  })
  #get mtally (transformed/normalized if possible)
  isnormlzd<-rep(FALSE,length(resavs))
  names(isnormlzd)<-names(resavs)
  nulldist<-sapply(names(resavs),function(i){
    null<-resavs[[i]]$nulldist
    obs<-sum(resavs[[i]]$mtally)
    tp<-c(null,obs)
    if(sd(null)>0 && median(null)>0){
      isnormlzd[i]<<-TRUE
      tp<-(tp-median(tp))/sd(tp)
    }
    tp
  })
  score<-nulldist[nrow(nulldist),]
  nulldist<-nulldist[-nrow(nulldist),,drop=FALSE]
  #----get ci
  pvals <- pnorm(score, lower.tail=FALSE)
  if(pAdjustMethod=="bonferroni"){
    ci <- qnorm(1-(pValueCutoff/length(pvals)))
  } else {
    ci <- qnorm(1-p.threshold(pvals, pValueCutoff, pAdjustMethod))
  }
  
  #----powerTransform
  if(boxcox){
    ptdist <- sapply(1:ncol(nulldist),function(i){
      null <- nulldist[,i]
      obs <- score[i] 
      if(isnormlzd[i] && shtest(null)){
        minval <- min(c(nulldist[,i],score[i]))
        minval <- ifelse(minval<=0,abs(minval)+1,minval)
        nullm <- null+minval
        obsm <- obs+minval
        obsm <- round(obsm,digits=5)
        # l <- coef(powerTransform(c(nullm,obsm)), round=TRUE)
        # ptdat <- bcPower(c(nullm,obsm),l)
        l <- round(.estimate.bcpower(c(nullm,obsm)),digits=5)
        ptdat <- .bcpower(c(nullm,obsm),l)
        ptdat <- (ptdat-median(ptdat))/sd(ptdat)
        return(ptdat)
      } else {
        return(c(null,obs))
      }
    })
    colnames(ptdist) <- names(score)
    score <- ptdist[nrow(ptdist),]
    nulldist <- ptdist[-nrow(ptdist),,drop=FALSE]
    pvals <- pnorm(score, lower.tail=FALSE)
    # (NEW) it corrects distributions not able of transformation and
    # obvious non-significant cases introduced by distortions of 
    # very sparse null distributions or absence of observations
    # in the mapping tally
    for (i in names(isnormlzd)){
      if(!isnormlzd[i]){
        p <- (1 + sum(score[i]<=nulldist[,i]) ) / (1 + nrow(nulldist))
        p <- min(0.5,p)
        pvals[i]<-p
        score[i]<-qnorm(p,lower.tail=FALSE)
      }
    }
  }
  #----reorder
  if(length(groups)>1){
    gs<-unique(groups)
    ord<-lapply(1:length(gs),function(i){
      idx<-sort.list(score[groups==i],decreasing=TRUE)
      labs<-names(score)[groups==i]
      labs[idx]
    })
    ord<-unlist(ord)
    score<-score[ord]
    pvals<-pvals[ord]
    mtally<-mtally[,ord]
    nulldist<-nulldist[,ord]
  }
  return(list(mtally=mtally,nulldist=nulldist,score=score,pvalue=pvals,ci=ci))
}

#-------------------------------------------------------------------------
# this internal function was required to fix installation of
# some dependencies (see "car" package for original implamation)
.bcpower <- function(U, lambda, jacobian.adjusted=FALSE){
  bc1 <- function(U, lambda){
    if(any(U[!is.na(U)] <= 0)) 
      stop("First argument must be strictly positive.")
    if (abs(lambda) <= 1e-06){
      z <- log(U)
    } else {
      z <- ((U^lambda) - 1)/lambda
    }
    if (jacobian.adjusted == TRUE){
      z <- z * (exp(mean(log(U), na.rm = TRUE)))^(1 - lambda)
    }
    z
  }
  out <- U
  if(is.matrix(out) | is.data.frame(out)){
    if (is.null(colnames(out))) 
      colnames(out) <- paste("Z", 1:dim(out)[2], sep = "")
    for (j in 1:ncol(out)){
      out[, j] <- bc1(out[, j], lambda[j])
    }
    colnames(out) <- paste(colnames(out), round(lambda, 2), sep = "^")
  } else {
    out <- bc1(out, lambda)
  }
  out
}
.estimate.bcpower <- function(obj){
  fam <- .bcpower
  Y <- as.matrix(obj)
  X <- matrix(rep(1, dim(Y)[1]), ncol=1) 
  w <- 1
  nc <- dim(Y)[2]
  nr <- nrow(Y)
  xqr <- qr(w * X)
  llik <- function(lambda){
    (nr/2)*log(((nr - 1)/nr) * det(var(qr.resid(xqr, w*fam(Y, lambda, j=TRUE)))))
  }
  llik1d <- function(lambda,Y){
    (nr/2)*log(((nr - 1)/nr) * var(qr.resid(xqr, w*fam(Y, lambda, j=TRUE))))
  }
  start <- rep(1, nc)
  for(j in 1:nc){
    res<- suppressWarnings(
      optimize(
        f = function(lambda) llik1d(lambda,Y[ , j, drop=FALSE]),
        lower=-3, upper=+3)
    )
    start[j] <- res$minimum
  }
  start
}

#-------------------------------------------------------------------------
vsereport<-function(obj){
  if(any(names(obj)=="escore")){
    obj$score<-obj$escore #just to correct a label!
  }
  mtally<-t(obj$mtally)
  mtally[,]<-as.integer(mtally)
  score<-obj$score
  pvalue<-obj$pvalue
  null<-t(boxplot(obj$nulldist, plot = FALSE)$stats)
  colnames(null)<-paste("q",1:5,sep="")
  report<-data.frame(Annotation=names(score),Pvalue=format(round(-log10(pvalue),3)), Score=format(round(score,3)), format(round(null,3)), mtally, stringsAsFactors = FALSE)
  rownames(report)<-NULL
  tp<-rowSums(mtally);tp<-tp/sum(tp)
  idx<-sort.list(score+tp,decreasing=TRUE)
  report<-report[idx,]
  return(report)
}

#-------------------------------------------------------------------------
shtest<-function(null){
  nnull<-length(null)
  nd<-as.integer(nnull*0.05)
  nd<-max(nd,10)
  qt<-quantile(null,probs=seq(0,1,length.out=nd),names=FALSE)
  shapiro.test(qt)$p.value<0.05
}

##-------------------------------------------------------------------------
##check if snow cluster is loaded
isParallel<-function(){
  b1<-"package:snow" %in% search()
  b2<-tryCatch({
    cl<-getOption("cluster")
    cl.check<-FALSE
    if(is(cl, "cluster")){
      cl.check <- all( sapply(1:length(cl),function(i)isOpen(cl[[i]]$con) ) == TRUE )
    }
    cl.check
  }, error=function(e){ FALSE 
  })
  all(c(b1,b2))
}

###########################################################################
## Methods (under development) to extract consolidated results
###########################################################################

#-------------------------------------------------------------------------
# extract eqtls, return consolidated results
eqtlExtract<-function(vSet,annot,gxdata,snpdata,pValueCutoff){
  eqtls<-eqtlExtractFull(vSet,annot,gxdata,snpdata)
  eqtls<-eqtls[eqtls$'Pr(>F)'<pValueCutoff,]
  rownames(eqtls)<-NULL
  eqtls
}
eqtlExtractFull<-function(vSet,annot,gxdata,snpdata){ 
  # mapping tally
  clusterMapping<-lapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      overlaps<-findOverlaps(query,subject)
      geneList<-.mtdata(subject)$mappedAnnotations
      res<-lapply(1:length(query@metadata$markers),function(j){
        snpList<-.mtdata(query)$mappedMarkers[[j]]
        ov<-S4Vectors::from(overlaps)%in%which(query@metadata$index==j)
        ov<-unique(S4Vectors::to(overlaps)[ov])
        gList<-geneList[ov]
        if(length(gList)>0 && length(snpList)>0){
          res<-eqtlTestDetailed(geneList=as.character(gList),snpList=as.integer(snpList),gxdata,snpdata)
          rownames(res)<-rownames(snpdata)[snpList]
        } else {
          res<-matrix(NA,nrow=length(snpList),ncol=2)
          rownames(res)<-snpList;colnames(res)<-c("F","Pr(>F)")
        }
        return(res)
      })
      names(res)<-query@metadata$markers
    } else {
      res<-NA
    }
    return(res)
  })
  #---simplify list
  res<-list()
  for(i in 1:length(clusterMapping)){
    tp<-clusterMapping[[i]]
    if(is.list(tp)){
      for(j in names(tp)){
        tpp<-tp[[j]]
        tpp<-tpp[!is.na(tpp[,2]),]
        if(ncol(tpp)>2 && nrow(tpp)>0){
          res[[j]]<-tpp
        }
      }
    }
  }
  clusterMapping<-res
  #---get summary
  summ<-data.frame(NULL,stringsAsFactors=FALSE)
  for(riskSNP in names(clusterMapping)){
    tp<-clusterMapping[[riskSNP]]
    Fstat<-tp[,1,drop=FALSE]
    Pstat<-tp[,2,drop=FALSE]
    geneid<-colnames(tp)[-c(1,2)]
    for(associatedSNP in rownames(Fstat)){
      tpp<-data.frame(riskSNP,associatedSNP,geneid,Fstat[associatedSNP,],Pstat[associatedSNP,],stringsAsFactors=FALSE)
      summ<-rbind(summ,tpp)
    }
  }
  if(nrow(summ)>0){
    colnames(summ)<-c("RiskSNP","RiskAssociatedSNP","GeneID","F","Pr(>F)")
  }
  summ
}

#-------------------------------------------------------------------------
##eqtl test for the 'eqtlExtract' function
##run two-way manova with multiple additive factors (gx as response varible)
eqtlTestDetailed<-function(geneList,snpList,gxdata,snpdata){
  gxdata<-t(gxdata[geneList,,drop=FALSE])
  snpdata<-t(snpdata[snpList,,drop=FALSE])
  #set names for formulae
  names(snpList)<-paste(rep("S",ncol(snpdata)),1:ncol(snpdata),sep="")
  names(geneList)<-paste(rep("G",ncol(gxdata)),1:ncol(gxdata),sep="")
  colnames(gxdata)<-names(geneList)
  colnames(snpdata)<-names(snpList)
  #run lm
  fm1<-paste(colnames(snpdata),collapse="+")
  fmla <- formula( paste("gxdata ~",fm1, collapse=" ") )
  resfit<-lm(fmla, data=as.data.frame(snpdata,stringsAsFactors=FALSE))
  if(ncol(gxdata)>1){
    resf<-summary(manova(resfit))
    resf<-as.data.frame(resf$stats)
    resf<-resf[1:(nrow(resf)-1),]
    resf<-resf[,c("approx F","Pr(>F)")]
  } else {
    resf<-anova(resfit)
    resf<-resf[1:(nrow(resf)-1),]
    resf<-as.data.frame(resf)
    resf<-resf[,c("F value","Pr(>F)")]
  }
  #library('gplots')
  #plotmeans(G2 ~ S6, data=cbind(gxdata,snpdata), col="red", barcol="blue",connect=FALSE,pch=15, las=1)
  #get SNP stats
  colnames(resf)<-c("F","Pr(>F)")
  #get gene stats
  resaov<-summary(aov(resfit))
  resp<-sapply(1:length(resaov),function(i){
    tp<-resaov[[i]][["Pr(>F)"]]
    tp[-length(tp)]
  })
  resp<-matrix(resp,ncol=length(geneList))
  colnames(resp)<-geneList
  #combine results
  res<-cbind(resf,resp)
  res<-res[names(snpList),]
  rownames(res)<-snpList
  return(res)
}

#-------------------------------------------------------------------------
eqtlExtractAnova<-function(vSet,annot,gxdata,snpdata){ 
  # mapping tally
  clusterMapping<-lapply(1:length(vSet),function(i){
    chr<-names(vSet[i])
    query<-vSet[[i]]
    subject<-annot[[chr]]
    if(!is.null(subject)){
      overlaps<-findOverlaps(query,subject)
      geneList<-.mtdata(subject)$mappedAnnotations
      resfit<-NULL
      for(j in 1:length(query@metadata$markers)){
        snpList<-.mtdata(query)$mappedMarkers[[j]]
        ov<-S4Vectors::from(overlaps)%in%which(query@metadata$index==j)
        ov<-unique(S4Vectors::to(overlaps)[ov])
        gList<-geneList[ov]
        if(length(gList)>0 && length(snpList)>0){
          res<-eqtlTestDetailedAnova(geneList=as.character(gList),snpList=as.integer(snpList),gxdata,snpdata)
          res$RiskAssociatedSNP<-rownames(snpdata)[res$RiskAssociatedSNP]
          res<-data.frame(RiskSNP=query@metadata$markers[j],res,stringsAsFactors=FALSE)
          resfit<-rbind(resfit,res)
        }
      }
    } else {
      resfit<-NA
    }
    return(resfit)
  })
  #---simplify list
  summ<-NULL
  for(i in 1:length(clusterMapping)){
    summ<-rbind(summ,clusterMapping[[i]])
  }
  if(nrow(summ)>0){
    colnames(summ)<-c(".RiskSNP",".RiskAssociatedSNP",".GeneID",".F",".Pr(>F)",".Coef",".R2")
  }
  summ
}

#-------------------------------------------------------------------------
eqtlTestDetailedAnova<-function(geneList,snpList,gxdata,snpdata){
  gxdata<-t(gxdata[geneList,,drop=FALSE])
  snpdata<-t(snpdata[snpList,,drop=FALSE])
  #set names for formulae
  names(snpList)<-paste(rep("S",ncol(snpdata)),1:ncol(snpdata),sep="")
  names(geneList)<-paste(rep("G",ncol(gxdata)),1:ncol(gxdata),sep="")
  colnames(gxdata)<-names(geneList)
  colnames(snpdata)<-names(snpList)
  #run lm
  resfit<-NULL
  if(ncol(snpdata)>1 && ncol(gxdata)>1){
    for(S in colnames(snpdata)){
      fmla <- formula( paste("gxdata ~",S, collapse=" ") )
      raov <- aov(fmla, data=as.data.frame(snpdata,stringsAsFactors=FALSE))
      cf<-raov$coefficients[S,]
      sf<-sapply(summary(raov),function(rv){
        tp<-as.data.frame(rv)
        tp<-tp[1:(nrow(tp)-1),4:5]
        as.numeric(tp)
      })
      rownames(sf)<-c("F","Prob")
      colnames(sf)<-names(cf)
      sf<-t(sf)
      raov<-data.frame(RiskAssociatedSNP=S,GeneID=rownames(sf),sf,Coef=cf,stringsAsFactors=FALSE)
      resfit<-rbind(resfit,raov)
    }
  } else {
    for(S in colnames(snpdata)){
      fmla <- formula( paste("gxdata ~",S, collapse=" ") )
      raov <- aov(fmla, data=as.data.frame(snpdata,stringsAsFactors=FALSE))
      cf<-raov$coefficients[2]
      sf<-as.data.frame(summary(raov)[[1]])
      sf<-sf[1,4:5]
      raov<-data.frame(RiskAssociatedSNP=S,GeneID=colnames(gxdata),F=sf[,1],Prob=sf[,2],Coef=cf,stringsAsFactors=FALSE)
      resfit<-rbind(resfit,raov)
    }
  }
  rownames(resfit)<-NULL
  resfit$RiskAssociatedSNP<-snpList[resfit$RiskAssociatedSNP]
  resfit$GeneID<-geneList[resfit$GeneID]
  #--
  #run cor
  R2<-cor(gxdata,snpdata)
  colnames(R2)<-snpList
  rownames(R2)<-geneList
  summ<-NULL
  for(i in colnames(R2)){
    tp<-data.frame(RiskAssociatedSNP=i,GeneID=rownames(R2),R2=R2[,i],stringsAsFactors=FALSE)
    summ<-rbind(summ,tp)
  }
  rownames(summ)<-NULL
  #---
  resfit<-cbind(resfit,R2=summ$R2)
  return(resfit)
}

##-------------------------------------------------------------------------
#return consolidated results in a matrix
#ps."object" should be an "avs" already evaluated by the "avs.evse" method
getEvseMatrix<-function(object,what="probs"){
  mappedIds<-getMappedMarkers(object)
  RiskSNP<-mappedIds$RiskSNP
  evsemtx<-NULL
  if(what=="probs"){
    vl=1;cl="Pr(>F)";efun=min
  } else if(what=="fstat"){
    vl=0;cl="F";efun=max
  }
  for(reg in names(object@results$evse)){
    eqtls<-object@results$evse[[reg]]$eqtls
    rvec<-rep(vl,length(RiskSNP))
    names(rvec)<-RiskSNP
    for(rs in RiskSNP){
      idx<-which(eqtls$RiskSNP==rs)
      if(length(idx)>0){
        rvec[rs]<-efun(eqtls[idx,cl])
      }
    }
    evsemtx<-rbind(evsemtx,rvec)
  }
  rownames(evsemtx)<-names(object@results$evse)
  evsemtx
}

##-------------------------------------------------------------------------
#another table to extract eqtls, return consolidated results
#ps."object" should be an "avs" already evaluated by the "avs.evse" method
getEvseEqtls<-function(object,tfs=NULL){
  if(is.null(tfs))tfs<-colnames(object@results$stats$evse$mtally)
  tp<-object@results$evse[tfs]
  res<-NULL
  for(tf in tfs){
    tpp<-tp[[tf]]$eqtls
    tpp<-data.frame(Regulon=tf,tpp,check.names = FALSE,stringsAsFactors = FALSE)
    res<-rbind(res,tpp)
  }
  res
}

##------------------------------------------------------------------------
##report markers and linked markers from computed variantSet
report.vset<-function(variantSet){
  lkmarkers<-lapply(variantSet,function(vset){
    if(class(vset)=="IRanges"){
      res<-lapply(names(vset@metadata$blocks),function(rs){
        linked_rs<-names(vset@metadata$blocks[[rs]])
        cbind(rs,rev(linked_rs))
      })
    } else {
      stop("Please, check 'vset' class! Method implemented for 'IRanges' objects only!")
    }
    res
  })
  summ<-NULL
  for(lt in lkmarkers){
    for(ltt in lt){
      summ<-rbind(summ,ltt)
    }
  }
  idx<-which(summ[,1]!=summ[,2])
  summ<-summ[idx,]
  summ<-data.frame(summ,stringsAsFactors = FALSE)
  colnames(summ) <- c("rs","linked_rs")
  return(summ)
}

##------------------------------------------------------------------------
##returns rejection threshold for methods in 'p.adjust'
p.threshold <- function (pvals, alpha=0.05, method="BH") {
  pvals <- sort(pvals)
  padj <- p.adjust(pvals, method = method)
  thset <- which(padj <= alpha)
  if(length(thset)>0){
    mx1 <- mx2 <- which.max(thset)
    if(mx2<length(padj)) mx2 <- mx2 + 1
    th <- (pvals[mx1] + min(pvals[mx2],alpha) ) / 2
  } else {
    th <- min(c(alpha,pvals))
  }
  return(th)
}
# p.threshold <- function (pvals, alpha=0.05, method="BH") {
#   pvals <- sort(pvals)
#   padj <- p.adjust(pvals, method = method)
#   thset <- which(padj <= alpha)
#   ifelse(length(thset) == 0, 0, pvals[thset[which.max(thset)]])
# }




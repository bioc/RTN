

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
getMarkers.vset<-function(variantSet){
  lkblocks <- lapply(variantSet,function(vset){
    if(is(vset,"IRanges")){
      tp <- lapply(vset@metadata$blocks,names)
    } else {
      tp <- lapply(vset,names)
    }
    tp
  })
  lkmarkers <- NULL
  for(i in 1:length(lkblocks)){
    tp <- lkblocks[[i]]
    for(j in 1:length(tp)){
      tpp <- tp[[j]]
      names(tpp) <- rep(names(tp)[j],length(tpp))
      lkmarkers <- c(lkmarkers, tpp)
    }
  }
  lkmarkers <- lkmarkers[!duplicated(lkmarkers)]
  return(lkmarkers)
}

##------------------------------------------------------------------------
##get markers from computed random AVS
getMarkers.rset<-function(randomSet){
  lkmarkers<-lapply(1:length(randomSet),function(i){
    res <- getMarkers.vset(randomSet[[i]])
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
##return a named character vector with the RiskAssociatedSNPs mapped in the 
##VSE analysis names indicate the LD cluster (i.e. the RiskSNP assignment) 
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
    snow::clusterExport(cl, list("get.avsdist","IRanges",
                                 "overlapsAny",".mtdata"),
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
evsea<-function(vSet, rSet, annot, gxdata, snpdata, pValueCutoff, 
                verbose=TRUE){
  
  #compute evse
  mtally<-get.eqtldist.evsea(vSet, annot, gxdata, snpdata, pValueCutoff)
  
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("get.eqtldist.evsea","eqtlTest","IRanges",
                                 "overlapsAny","findOverlaps",".mtdata"),
                        envir=environment())
    resrset <- snow::parSapply(cl, 1:length(rSet), function(i) {
      sum(get.eqtldist.evsea(rSet[[i]], annot, gxdata, snpdata, pValueCutoff))
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    resrset<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet))
      sum(get.eqtldist.evsea(rSet[[i]], annot, gxdata, snpdata, pValueCutoff))
    })
    if(verbose)close(pb)
  }
  return(list(mtally=mtally,nulldist=resrset,nclusters=length(mtally)))
}

##-------------------------------------------------------------------------
##run evsea for observed and random variant sets
evseaproxy<-function(rSet, annot, gxdata, snpdata, pValueCutoff, 
                     verbose=TRUE){
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("get.eqtldist.evsea","eqtlTest","IRanges",
                                 "overlapsAny","findOverlaps",".mtdata"),
                        envir=environment())
    nulldist <- snow::parSapply(cl, 1:length(rSet), function(i) {
      res<-get.eqtldist.evsea(rSet[[i]], annot, gxdata, snpdata, pValueCutoff)
      sum(res)
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    nulldist<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet))
      res<-get.eqtldist.evsea(rSet[[i]], annot, gxdata, snpdata, pValueCutoff)
      sum(res)
    })
    if(verbose)close(pb)
  }
  return(nulldist)
}

##-------------------------------------------------------------------------
##get avs/eqtl dist for a variant set
get.eqtldist.evsea<-function(vSet,annot,gxdata,snpdata,pValueCutoff){
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
            res<-eqtlTest(geneList=gList, snpList, gxdata, snpdata)
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
##run two-way manova with multiple additive factors (gx as response variable)
eqtlTest<-function(geneList,snpList,gxdata,snpdata){
  gxdata<-t(gxdata[geneList,,drop=FALSE])
  snpdata<-t(snpdata[snpList,,drop=FALSE])
  # set names for formulae
  colnames(gxdata)<-paste0(rep("G",ncol(gxdata)),1:ncol(gxdata))
  colnames(snpdata)<-paste0(rep("S",ncol(snpdata)),1:ncol(snpdata))
  # run lm
  fm1<-paste(colnames(snpdata),collapse="+")
  fmla <- formula(paste("gxdata ~",fm1, collapse=" "))
  resfit<-lm(fmla, data=as.data.frame(snpdata))
  if(ncol(gxdata)>1){
    resf<-summary(manova(resfit))
    resf<-resf$stats[,"Pr(>F)"]
  } else {
    resf<-anova(resfit) 
    resf<-resf[["Pr(>F)"]]
  }
  if(all(is.na(resf))){
    resf=1
  } else {
    resf<-min(resf,na.rm=TRUE)
  }
  return(resf)
}

###########################################################################
## RVSE analysis
###########################################################################

##-------------------------------------------------------------------------
##run rvsea for observed and random variant sets
rvsea<-function(vSet, rSet, annot, regact, snpdata, 
                pValueCutoff, verbose=TRUE){
  
  #compute evse
  mtally <- get.eqtldist.rvsea(vSet, annot, regact, 
                               snpdata, pValueCutoff)
  #compute null
  if(isParallel()){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("get.eqtldist.rvsea","eqtlTest.rvsea",".mtdata"),
                        envir=environment())
    resrset <- snow::parSapply(cl, 1:length(rSet), function(i) {
      sum(get.eqtldist.rvsea(rSet[[i]], annot, regact, 
                             snpdata, pValueCutoff))
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    resrset<-sapply(1:length(rSet),function(i){
      if(verbose) setTxtProgressBar(pb, i/length(rSet))
      sum(get.eqtldist.rvsea(rSet[[i]], annot, regact, 
                             snpdata, pValueCutoff))
    })
    if(verbose)close(pb)
  }
  return(list(mtally=mtally,nulldist=resrset,nclusters=length(mtally)))
}

##-------------------------------------------------------------------------
##get avs/rqtl dist for a variant set
get.eqtldist.rvsea<-function(vSet,annot,regact,snpdata,pValueCutoff){
  clusterMapping<-lapply(1:length(vSet),function(i){
    chr <- names(vSet[i])
    query <- vSet[[i]]
    subject <- annot[[chr]]
    if(!is.null(subject)){
      ov <- rep(FALSE,length(query@metadata$markers))
      ov[query@metadata$index[overlapsAny(query,subject)]] <- TRUE
      res <- sapply(1:length(query@metadata$markers),function(j){
        snpList <- .mtdata(query)$mappedMarkers[[j]]
        if(ov[j] && length(snpList)>0){
          res <- eqtlTest.rvsea(regact, snpdata, snpList)
        } else {
          res <- 1.0
        }
        return(res)
      })
    } else {
      res<-rep(1,length(query@metadata$markers))
    }
    names(res)<-query@metadata$markers
    return(res)
  })
  clusterMapping<-unlist(clusterMapping)
  clusterMapping<pValueCutoff
}

##-------------------------------------------------------------------------
##run lm with multiple additive factors (regact as response variable)
eqtlTest.rvsea<-function(regact, snpdata, snpList){
  snpdata <- t(snpdata[snpList,,drop=FALSE])
  colnames(snpdata) <- paste0(rep("S",ncol(snpdata)),1:ncol(snpdata))
  fm1 <- paste(colnames(snpdata),collapse="+")
  fmla <- formula(paste("regact ~",fm1, collapse=" "))
  resf <- summary(aov(fmla, data=as.data.frame(snpdata)))
  resf <- resf[[1]][["Pr(>F)"]]
  if(all(is.na(resf))){
    resf=1
  } else {
    resf<-min(resf,na.rm=TRUE)
  }
  return(resf)
}

##-------------------------------------------------------------------------
# eqtlTest.rvsea<-function(regact, snpdata, snpList){
#   snpdata <- t(snpdata[snpList,,drop=FALSE])
#   N <- nrow(snpdata)
#   R <- as.numeric(cor(regact,snpdata, method = "spearman", use="pairwise.complete.obs"))
#   resr <- 2 * pt( -abs( R*sqrt((N-2)/(1-R^2)) ), N-2)
#   if(all(is.na(resr))){
#     resr=1
#   } else {
#     resr<-min(resr,na.rm=TRUE)
#   }
#   return(resr)
# }

##-------------------------------------------------------------------------
##get avs/rqtl dist for a variant set
# get.eqtldist.rvsea<-function(vSet, annot, regact, snpdata, pValueCutoff){
#   clusterMapping<-sapply(1:length(vSet),function(i){
#     mappedMarkers <- .mtdata(vSet[[i]])$mappedMarkers
#     res<-sapply(1:length(mappedMarkers),function(j){
#       snpList<-mappedMarkers[[j]]
#       if(length(snpList)>0){
#         res<-eqtlTest.rvsea(regact, snpdata, snpList)
#       } else {
#         res<-1.0
#       }
#       return(res)
#     })
#     names(res)<-names(mappedMarkers)
#     return(res)
#   })
#   clusterMapping<-unlist(clusterMapping)
#   clusterMapping<pValueCutoff
# }
##-------------------------------------------------------------------------
# eqtlTest.rvsea <- function(regact, snpdata, snpList){
#   snpdata <- t(snpdata[snpList,,drop=FALSE])
#   fit <- aov(regact~snpdata)
#   asgn <- fit$assign[fit$qr$pivot[1L:fit$rank]]
#   uasgn <- unique(asgn)
#   effects <- fit$effects
#   effects <- as.matrix(effects)[seq_along(asgn), , drop = FALSE]
#   rdf <- fit$df.residual
#   resid <- as.matrix(fit$residuals)
#   ai <- (asgn == uasgn[2])
#   df <- sum(ai)
#   ss <- sum(effects[ai, 1]^2)
#   df <- c(df, rdf)
#   ss <- c(ss, sum(resid[, 1]^2))
#   ms <- ss/df
#   TT <- ms/ms[length(df)]
#   p <- pf(TT[1], df[1], rdf, lower.tail = FALSE)
#   p <- ifelse(is.na(p),1,p)
#   return(p)
# }

###########################################################################
## Others
###########################################################################

#-------------------------------------------------------------------------
getUniverseCounts1<-function(vSet,annotation,maxgap){
  #count markers in hapmap
  clusterCounts<-unlist(lapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$blocks,length))
  }))
  #count tested genes, overlap
  annot<-getAnnotRanges(annotation,maxgap=maxgap,getTree=FALSE, 
                        getReduced=FALSE)
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
  geneCounts <- unlist(geneCounts)
  #merge counts
  counts <- cbind(clusterCounts,geneCounts)
  colnames(counts)<-c("markers","annotation")
  counts
}

#-------------------------------------------------------------------------
getUniverseCounts2<-function(vSet,annotation,maxgap){
  #count markers
  totalClusterCounts<-unlist(sapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$blocks,length))
  }))
  #count markers mapped to the genotype data
  mappedClusterCounts<-unlist(sapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$mappedMarkers,length))
  }))
  #count tested genes, overlap
  annot<-getAnnotRanges(annotation=annotation,maxgap=maxgap,getTree=FALSE, 
                        getReduced=FALSE)
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
  counts <- cbind(totalClusterCounts,mappedClusterCounts,geneCounts)
  colnames(counts)<-c("totalMarkers","markers","annotation")
  return(counts)
}

#-------------------------------------------------------------------------
getUniverseCounts3<-function(vSet){
  #count total markers
  totalClusterCounts<-unlist(lapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$blocks,length))
  }))
  #count markers mapped to the genotype data
  mappedClusterCounts<-unlist(sapply(1:length(vSet),function(i){
    unlist(lapply(vSet[[i]]@metadata$mappedMarkers,length))
  }))
  counts<-cbind(totalClusterCounts,mappedClusterCounts)
  colnames(counts)<-c("totalMarkers","markers")
  return(counts)
}

#-------------------------------------------------------------------------
vseformat<-function(resavs, pValueCutoff, pAdjustMethod="bonferroni", 
                    boxcox=TRUE){
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
    (nr/2)*log(((nr - 1)/nr) * det(var(
      qr.resid(xqr, w*fam(Y, lambda, j=TRUE)))))
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
  report<-data.frame(Annotation=names(score),
                     Pvalue=format(round(-log10(pvalue),3)), 
                     Score=format(round(score,3)), 
                     format(round(null,3)), mtally, stringsAsFactors = FALSE)
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
      cl.check <- all( sapply(
        1:length(cl),function(i)isOpen(cl[[i]]$con) ) == TRUE )
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
  eqtls<-eqtls[eqtls$P.Multi<pValueCutoff,]
  rownames(eqtls)<-NULL
  eqtls
}
eqtlExtractFull<-function(vSet,annot,gxdata,snpdata){ 
  # mapping
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
          res <- eqtlTestDetailed(geneList=as.character(gList),
                                  snpList=snpList, gxdata, snpdata)
        } else {
          res <- NA
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
      tpp <- tp[!is.na(tp)]
      res <- c(res,tpp)
    }
  }
  clusterMapping<-res
  #---get summary
  summ<-data.frame(NULL,stringsAsFactors=FALSE)
  for(riskSNP in names(clusterMapping)){
    tp <- clusterMapping[[riskSNP]]
    nr <- nrow(tp$f.uni.stat)
    nc <- ncol(tp$f.uni.stat)
    asnp <- rep(rownames(tp$f.uni.stat),nc)
    genes <- rep(colnames(tp$f.uni.stat),each=nr)
    #---
    fmulti <- rep(tp$f.multi.stat[,"F"],nc)
    pmulti <- rep(tp$f.multi.stat[, "Pr(>F)"],nc) 
    funi <- as.numeric(tp$f.uni.stat)
    puni <- as.numeric(tp$p.uni.stat)
    r <- as.numeric(tp$r.stat)
    tpp <- data.frame(riskSNP, asnp, genes, fmulti, pmulti, 
                      funi, puni, r, stringsAsFactors=FALSE)
    summ<-rbind(summ,tpp)
  }
  if(nrow(summ)>0){
    colnames(summ)<-c("RiskSNP","RiskAssociatedSNP","GeneID","F.Multi","P.Multi","F.Uni","P.Uni","R")
  }
  summ
}

#-------------------------------------------------------------------------
##for 'eqtlExtract' function
##run two-way manova with multiple additive factors (gx as response variable)
eqtlTestDetailed<-function(geneList,snpList,gxdata,snpdata){
  snpListLabs <- snpList; names(snpListLabs) <- rownames(snpdata)[snpList]
  gxdt <- t(gxdata[geneList,,drop=FALSE])
  snpdt <- t(snpdata[snpList,,drop=FALSE])
  #set names for formulae
  names(snpList) <- paste(rep("S",ncol(snpdt)),1:ncol(snpdt),sep="")
  names(geneList) <- paste(rep("G",ncol(gxdt)),1:ncol(gxdt),sep="")
  colnames(gxdt) <- names(geneList)
  colnames(snpdt) <- names(snpList)
  #---run lm multivar
  fm1 <- paste(colnames(snpdt),collapse="+")
  fmla <- formula( paste("gxdt ~",fm1, collapse=" ") )
  resfit <- lm(fmla, data=as.data.frame(snpdt,stringsAsFactors=FALSE))
  if(ncol(gxdt)>1){
    resf.multi <- summary(manova(resfit))
    resf.multi <- as.data.frame(resf.multi$stats)
    resf.multi <- resf.multi[1:(nrow(resf.multi)-1),]
    resf.multi <- resf.multi[,c("approx F","Pr(>F)")]
  } else {
    resf.multi <- anova(resfit)
    resf.multi <- resf.multi[1:(nrow(resf.multi)-1),]
    resf.multi <- as.data.frame(resf.multi)
    resf.multi <- resf.multi[,c("F value","Pr(>F)")]
  }
  colnames(resf.multi) <- c("F","Pr(>F)")
  resf.multi <- resf.multi[names(snpList),,drop=FALSE]
  resf.multi <- as.matrix(resf.multi)
  # library('gplots')
  # plotmeans(G1 ~ S3, data=cbind(gxdt,snpdt), col="red",
  # barcol="blue",connect=FALSE,pch=15, las=1)
  #---run lm univar
  resf.uni <- resp.uni <- matrix(NA,ncol=length(geneList), nrow=length(snpList))
  for(i in 1:ncol(snpdt)){
    for(j in 1:ncol(gxdt)){
      tp <- summary(aov(snpdt[,i]~gxdt[,j]))[[1]]
      resf.uni[i,j] <- tp[["F value"]][1]
      resp.uni[i,j] <- tp[["Pr(>F)"]][1]
    }
  }
  colnames(resf.uni) <- colnames(resp.uni) <- as.character(geneList)
  rownames(resf.uni) <- rownames(resp.uni) <- names(snpList)
  #---run cor
  rescor <- t(cor(gxdt,snpdt, method = "spearman"))
  colnames(rescor) <- geneList
  rescor <- rescor[names(snpList),,drop=FALSE]
  #---update names and remove NAs
  rownames(resf.multi) <- names(snpListLabs)
  rownames(resf.uni) <- names(snpListLabs)
  rownames(resp.uni) <- names(snpListLabs)
  rownames(rescor) <- names(snpListLabs)
  idx <- rowSums(is.na(resf.multi))==0
  resf.multi <- resf.multi[idx,,drop=FALSE]
  resf.uni <- resf.uni[idx,,drop=FALSE]
  resp.uni <- resp.uni[idx,,drop=FALSE]
  rescor <- rescor[idx,,drop=FALSE]
  res <- list(f.multi.stat=resf.multi, f.uni.stat=resf.uni,
              p.uni.stat=resp.uni, r.stat=rescor)
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
          res < -eqtlTestDetailedAnova(geneList=as.character(gList),
                                     snpList=as.integer(snpList),gxdata,snpdata)
          res$RiskAssociatedSNP<-rownames(snpdata)[res$RiskAssociatedSNP]
          res<-data.frame(RiskSNP=query@metadata$markers[j],res, stringsAsFactors=FALSE)
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
    colnames(summ)<-c(".RiskSNP",".RiskAssociatedSNP",".GeneID",
                      ".F",".Pr(>F)",".Coef",".R2")
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
      raov<-data.frame(RiskAssociatedSNP=S,GeneID=rownames(sf),sf,
                       Coef=cf,stringsAsFactors=FALSE)
      resfit<-rbind(resfit,raov)
    }
  } else {
    for(S in colnames(snpdata)){
      fmla <- formula( paste("gxdata ~",S, collapse=" ") )
      raov <- aov(fmla, data=as.data.frame(snpdata,stringsAsFactors=FALSE))
      cf<-raov$coefficients[2]
      sf<-as.data.frame(summary(raov)[[1]])
      sf<-sf[1,4:5]
      raov<-data.frame(RiskAssociatedSNP=S,GeneID=colnames(gxdata),
                       F=sf[,1],Prob=sf[,2],Coef=cf,stringsAsFactors=FALSE)
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
    tp<-data.frame(RiskAssociatedSNP=i,GeneID=rownames(R2),R2=R2[,i],
                   stringsAsFactors=FALSE)
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
getEvseMatrix<-function(object, what="P.Multi"){
  mappedIds<-getMappedMarkers(object)
  RiskSNP<-mappedIds$RiskSNP
  evsemtx<-NULL
  if(what=="P.Multi"){
    vl=1;cl="P.Multi";efun=min
  } else if(what=="F.Multi"){
    vl=0;cl="F.Multi";efun=max
  } else if(what=="F.Uni"){
    vl=0;cl=c("F.Uni","P.Uni");
    efun=function(x){ x[which.min(x[,2]),1] }
  } else if(what=="R"){
    vl=0;cl=c("R","P.Uni");
    efun=function(x){ x[which.min(x[,2]),1] }
  }
  for(reg in names(object@results$evse)){
    eqtls<-object@results$evse[[reg]]$eqtls
    rvec<-rep(vl,length(RiskSNP))
    names(rvec)<-RiskSNP
    for(rs in RiskSNP){
      idx <- which(eqtls$RiskSNP==rs)
      if(length(idx)>0){
        rvec[rs] <- efun(eqtls[idx,cl])
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
    if(is(vset,"IRanges")){
      res<-lapply(names(vset@metadata$blocks),function(rs){
        linked_rs<-names(vset@metadata$blocks[[rs]])
        cbind(rs,rev(linked_rs))
      })
    } else {
      stop(
        "Please, check 'vset' class! Method implemented for 'IRanges' 
        objects only!"
        )
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




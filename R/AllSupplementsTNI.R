

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##       OPTIMIZED MI FUNCTIONS FOR LARGE-SCALE TNETS
##      -- IMPLEMENTS AN INTERFACE FOR ARACNE/MINET --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##------------------------------------------------------------------------
##This function takes a gene expression matrix (x), a list of TFs
##and computes a partial mutal information matrix, i.e., between each TF 
##and all potential targets. Sig. mi values are inferred by permutation 
##analysis.
tni.pmin<-function(x, tfs, estimator="spearman", 
                   boxcox=FALSE, isboot=FALSE, 
                   setdg=0, simplified=FALSE, getadj=FALSE){
  x <- t(x)
  pmim <- suppressWarnings(cor(x[,tfs],x, method=estimator, 
                               use="complete.obs"))
  if(boxcox) pmim <- .boxcox.adj(pmim, isboot)
  pmim <- pmim^2
  pmim[is.na(pmim)] <- 0
  if(length(tfs)>1){
    diag(pmim[,tfs])=setdg
  } else {
    pmim[,tfs]=setdg
  }
  maxi=0.999999
  pmim[which(pmim > maxi)]=maxi
  pmim = -0.5 * log(1 - pmim)
  pmim[pmim<0]=0
  if(simplified){
    return(t(pmim))
  } else if(getadj){
    mim=matrix(setdg,ncol=ncol(x),nrow=ncol(x))
    colnames(mim)=colnames(x)
    rownames(mim)=colnames(x)  
    mim[tfs,]=pmim
    mim[,tfs]=t(pmim)
    return(mim)
  } else {
    pmim<-t(pmim)
    colnames(pmim)<-tfs
    return(pmim)
  }
}
#We have observed that transcription factor (TF) regulons reconstructed from 
#the RTN package exhibit different proportions of positive and negative targets. 
#While the proportion can vary between different regulons, we have observed a 
#consistent higher proportion of positive targets, specially when using RNA-seq 
#data. RTN uses mutual information (MI) to assess TF-target associations, 
#assigning the direction of the inferred associations by Spearman's correlation. 
#Dam et al. (2017) have acknowledged that different normalization methods 
#introduce different biases in co-expression analysis, usually towards positive 
#correlation, possibly affected by read-depth differences between samples and 
#the large abundance of 0 values present in RNA-seq-derived expression matrices. 
#In order to correct this positive correlation bias we suggest using this 
#box-cox correction strategy.
.boxcox.adj <- function(xmat, isboot=FALSE){
  xvec <- as.numeric(xmat)
  if(isboot){
    #only used to speed the bootstrep rounds
    nn <- max( min(1e+5,length(xvec)), length(xvec)*0.1)
    xvec <- sample(xvec, round(nn))
  }
  l <- coef(powerTransform(xvec,  family="yjPower"))
  xmat[,] <- yjPower(as.numeric(xmat),l)
  return(xmat)
}
##local pmin function for tni.perm.separate
.perm.pmin.separate<-function(x,tf,estimator="spearman",nPerm=1000){
  x=t(x)
  x=cor(replicate(nPerm,sample(x[,tf])),x[,-tf], 
        method=estimator,use="complete.obs")^2
  maxi=0.999999
  x[which(x > maxi)]=maxi
  x = -0.5 * log(1 - x)
  x[x<0]=0
  t(x)
}
##local pmin function for pooled tni.perm.pooled
.perm.pmin.pooled<-function(x, n, estimator="spearman"){
  x<-matrix(sample(x),ncol=nrow(x),nrow=ncol(x))
  x<-cor(x[,1:n],x, method=estimator, use="complete.obs")^2
  if(n>1){
    diag(x[,1:n])=NA
  } else {
    x[,1:n]=NA
  }
  maxi=0.999999
  x[which(x>maxi)]=maxi
  x = -0.5*log(1-x)
  x[x<0]=0
  t(x)
}

##------------------------------------------------------------------------
##This function takes a partial mutal information matrix and apply the dpi 
##filter (aracne algorithm). TFs should be in cols and potential targets 
##in rows.
tni.dpi<-function(pmim,eps=0){
  tfs<-match(colnames(pmim),rownames(pmim))
  tar<-setdiff(1:nrow(pmim),tfs)
  x=pmim
  for(i in tar){
    idx<-pmim[i,]>0
    if(sum(idx)>1){
      mi<-pmim[c(i,tfs[idx]),idx]
      mi<-cbind(c(0,mi[1,]),mi)
      mi<-aracne(mi,eps=eps)
      tarvec<-mi[-1,1]
      x[i,idx]<-tarvec
    }
  }
  if(ncol(pmim)>1)x[tfs,]<-aracne(pmim[tfs,],eps=eps)
  return(x)
}

##------------------------------------------------------------------------
## permutation analysis with separated null distribution
tni.perm.separate<-function(object,verbose=TRUE){
  ##check if package snow has been loaded and a cluster object has been created
  if(isParallel() && length(object@regulatoryElements)>1){
    ispar <- TRUE
    if(verbose)
      cat("-Performing permutation analysis (parallel version - ProgressBar disabled)...\n")
    if(verbose)
      cat("--For", length(object@regulatoryElements), "regulons...\n")
  } else {
    ispar <- FALSE
    if(verbose)cat("-Performing permutation analysis...\n")
    if(verbose)cat("--For", length(object@regulatoryElements), "regulons...\n")
  }
  ##compute partial mi matrix
  gexp <- object@gexp[object@targetElements,]
  pmim<-tni.pmin(gexp,object@regulatoryElements,object@para$perm$estimator, 
                 object@para$perm$boxcox)
  ##compute null distributions
  if(ispar){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list(".perm.pmin.separate"),envir=environment())
    mipval <- snow::parSapply(cl, object@regulatoryElements, function(tf){
      pi<-which(tf==rownames(gexp))
      midist <- .perm.pmin.separate(gexp, pi, object@para$perm$estimator, 
                                    object@para$perm$nPermutations)
      midist<-sort(midist, na.last = NA)
      np<-length(midist)
      midist<-np-findInterval(pmim[,tf],midist)
      midist <- (1 + midist)/(1 + np)
      midist
    })
  } else {
    if(verbose) pb <- txtProgressBar(style=3)
    mipval<-sapply(1:length(object@regulatoryElements), function(i){
      tf<-object@regulatoryElements[i]
      pi<-which(tf==rownames(gexp))
      midist <- .perm.pmin.separate(gexp, pi, object@para$perm$estimator, 
                                    object@para$perm$nPermutations)
      midist<-sort(midist, na.last = NA)
      np<-length(midist)
      midist<-np-findInterval(pmim[,tf],midist)
      midist <- (1 + midist)/(1 + np)
      if(verbose) setTxtProgressBar(pb, i/length(object@regulatoryElements))      
      midist   
    })
    if(verbose)close(pb)
  }
  ##if globalAdjustment...
  if(object@para$perm$globalAdjustment){
    miadjpv<-p.adjust(mipval,method=object@para$perm$pAdjustMethod)
    miadjpv<-matrix(miadjpv,nrow=nrow(pmim),ncol=ncol(pmim))    
  } else {
    miadjpv<-apply(mipval,2,p.adjust,method=object@para$perm$pAdjustMethod)
  }
  dimnames(mipval)<-dimnames(miadjpv)<-dimnames(pmim)
  ##decide on the significance and return results
  pmim[miadjpv>object@para$perm$pValueCutoff]<-0
  return(list(tn.ref=pmim,mipval=mipval, miadjpv=miadjpv))
}

##------------------------------------------------------------------------
##permutation with pooled null distribution
tni.perm.pooled<-function(object, parChunks=10, verbose=TRUE){
  ##check if package snow has been loaded and a cluster object has been created
  if(isParallel() && length(object@regulatoryElements)>1){
    ispar <- TRUE
    if(verbose)cat("-Performing permutation analysis (parallel version)...\n")
    if(verbose)cat("--For", length(object@regulatoryElements), "regulons...\n")
  } else {
    ispar <- FALSE
    if(verbose)cat("-Performing permutation analysis...\n")
    if(verbose)cat("--For", length(object@regulatoryElements), "regulons...\n")
  }
  ##compute partial mi matrix and get unique values
  gexp <- object@gexp[object@targetElements,]
  uniqueVec<-tni.pmin(gexp,object@regulatoryElements,
                      object@para$perm$estimator,
                      object@para$perm$boxcox) 
  uniqueVec<-sort(unique(as.numeric(uniqueVec)))
  ##initialize permutation count
  ctsum<-numeric(length(uniqueVec))
  ##build the null distribution via permutation
  if(ispar){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list(".perm.pmin.pooled"),envir=environment())
    nper <- object@para$perm$nPermutations
    if(is.null(parChunks)){
      parChunks <- seq.int(length(cl), nper, by=length(cl)*3)
      if(!nper%in%parChunks)parChunks <- c(parChunks, nper)
      parChunks <- length(parChunks)
    }
    if(parChunks>(nper/2)){
      parChunks<-as.integer(nper/2)
    }
    nperChunks<-1:nper
    nperChunks<-split(nperChunks, cut(nperChunks, parChunks))
    nperChunks<-unlist(lapply(nperChunks,length))
    if(verbose)pb<-txtProgressBar(style=3)
    for(i in 1:parChunks){
      permdist <- snow::parLapply(cl,1:nperChunks[i],function(j){
        permt <- .perm.pmin.pooled(gexp,length(object@regulatoryElements),
                                   object@para$perm$estimator)
        permt <- sort(permt, na.last = NA)
        length(permt)-findInterval(uniqueVec,permt)
      })
      for(j in 1:length(permdist) ){
        ctsum <- ctsum + permdist[[j]]
      }
      rm(permdist);gc()
      if(verbose)setTxtProgressBar(pb, i/parChunks)
    }
  } else {
    if(verbose)pb<-txtProgressBar(style=3)
    for(i in 1:object@para$perm$nPermutations){
      permt <- .perm.pmin.pooled(gexp,length(object@regulatoryElements),
                                 object@para$perm$estimator)
      permt <- sort(permt, na.last = NA)
      permt <- length(permt)-findInterval(uniqueVec,permt)
      ctsum <- ctsum + permt
      if(verbose)setTxtProgressBar(pb, i/object@para$perm$nPermutations)
    }
  }
  if(verbose)close(pb)
  ##compute pvals
  pmim<-tni.pmin(gexp,object@regulatoryElements,object@para$perm$estimator,
                 object@para$perm$boxcox)
  np <- object@para$perm$nPermutations * 
    ( prod(dim(pmim)) - length(object@regulatoryElements) )
  mipval <- (1 + ctsum)/(1 + np)
  mipval<-mipval[match(as.numeric(pmim),uniqueVec)]
  mipval<-matrix(mipval,nrow=nrow(pmim),ncol=ncol(pmim))  
  ##adjust pvals
  if(object@para$perm$globalAdjustment){
    miadjpv<-p.adjust(mipval,method=object@para$perm$pAdjustMethod)
    miadjpv<-matrix(miadjpv,nrow=nrow(pmim),ncol=ncol(pmim))  
  } else {
    miadjpv<-apply(mipval,2,p.adjust,method=object@para$perm$pAdjustMethod)
  }
  dimnames(mipval) <- dimnames(miadjpv) <- dimnames(pmim)
  ##decide on the significance
  if(object@para$perm$pValueCutoff <= 1/np){
    tp1 <- "The input 'pValueCutoff' is below the permutation resolution (>"
    tp2 <- ")\n...a larger number of permutations is required to accurately 
    estimate this threshold!"
    warning(tp1,signif(1/np,1),tp2,call.=FALSE)
  }
  pmim[miadjpv>object@para$perm$pValueCutoff]=0.0
  
  return(list(tn.ref=pmim, mipval=mipval, miadjpv=miadjpv))
}

##------------------------------------------------------------------------
#bootstrap analysis
tni.boot<-function(object, parChunks=10, verbose=TRUE){
  ##check if package snow has been loaded and a cluster object has been created
  if(isParallel() && length(object@regulatoryElements)>1){
    ispar <- TRUE
    if(verbose)cat("-Performing bootstrap analysis (parallel version)...\n")
    if(verbose)cat("--For", length(object@regulatoryElements), "regulons...\n")
  } else {
    ispar <- FALSE
    if(verbose)cat("-Performing bootstrap analysis...\n")
    if(verbose)cat("--For", length(object@regulatoryElements), "regulons...\n")
  }
  ##identify empirical boundaries from tn.ref
  ##..this strategy turns the bootstrap more tractable for large 
  ##..scale datasets.
  mimark<-sapply(1:ncol(object@results$tn.ref), function(i){
    tp<-abs(object@results$tn.ref[object@results$tn.ref[,i]!=0,i])
    if(length(tp)>0) min(tp) else 0
  })
  mimark<-rep(median(mimark),length(mimark))
  ##initialize bootstrap matrix
  bcount<-matrix(0,ncol=ncol(object@results$tn.ref),
                 nrow=nrow(object@results$tn.ref))
  gexp <- object@gexp[object@targetElements,]
  ##run bootstrap
  if(ispar){
    cl<-getOption("cluster")
    snow::clusterExport(cl, list("tni.pmin",".boxcox.adj",
                                 "powerTransform","yjPower"),
                        envir=environment())
    nboot <- object@para$boot$nBootstraps
    if(is.null(parChunks)){
      parChunks <- seq.int(length(cl), nboot, by=length(cl))
      if(!nboot%in%parChunks)parChunks <- c(parChunks, nboot)
      parChunks <- length(parChunks)
    }
    if(parChunks>(nboot/2)){
      parChunks<-as.integer(nboot/2)
    }
    nbootChunks<-1:nboot
    nbootChunks<-split(nbootChunks, cut(nbootChunks, parChunks))
    nbootChunks<-unlist(lapply(nbootChunks,length))
    if(verbose)pb<-txtProgressBar(style=3)
    for(i in 1:parChunks){
      bootdist <- snow::parLapply(cl,1:nbootChunks[i],function(j){
        #--- get bootstrap sample
        x <- gexp[,sample(ncol(gexp),replace=TRUE)]
        #--- add a null dist to deal with eventual sd==0
        idx <- which(apply(x,1,sd)==0)
        if(length(idx)>0) x[idx,] <- sample(x,length(idx)*ncol(x))
        #--- compute bootdist
        boott<-tni.pmin(x, tfs=object@regulatoryElements, 
                        estimator=object@para$boot$estimator,
                        boxcox=object@para$boot$boxcox, 
                        isboot=TRUE, simplified=TRUE)
        sapply(1:ncol(boott),function(k){as.numeric(boott[,k]>mimark[k])})
      })
      for(j in 1:length(bootdist)){
        bcount<-bcount+bootdist[[j]]
      }
      rm(bootdist);gc()
      if(verbose)setTxtProgressBar(pb, i/parChunks)
    }
  } else {
    if(verbose) pb <- txtProgressBar(style=3)  
    for(i in 1:object@para$boot$nBootstraps){
      if(verbose) setTxtProgressBar(pb, i/object@para$boot$nBootstraps)
      #--- get bootstrap sample
      x <- gexp[,sample(ncol(gexp),replace=TRUE)]
      #--- add a null dist to deal with eventual sd==0
      idx <- which(apply(x,1,sd)==0)
      if(length(idx)>0) x[idx,] <- sample(x,length(idx)*ncol(x))
      #--- compute bootdist
      bootdist<-tni.pmin(x, tfs=object@regulatoryElements,
                         estimator=object@para$boot$estimator,
                         boxcox=object@para$boot$boxcox, 
                         isboot=TRUE, simplified=TRUE)
      bootdist<-sapply(1:ncol(bootdist),function(j){
        as.numeric(bootdist[,j]>mimark[j])
        })
      bcount <- bcount + bootdist
    }
  }
  if(verbose)close(pb)
  ##decide on the stability of the inferred associations
  cons<-as.integer(object@para$boot$nBootstraps * 
                     (object@para$boot$consensus/100) )
  object@results$tn.ref[bcount<=cons]<-0
  return(object@results$tn.ref)
}


##------------------------------------------------------------------------
##------------------------------------------------------------------------
##   -- COMPLEMENTARY FUNCTIONS TO DERIVE DIRECTIONALITY FOR TNETs --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##------------------------------------------------------------------------
##This function adds (-/+) correlation signals to inferred associations.
##It takes a gene expression matrix (x) and a tnet, and computes a simple 
##correlation between each TF and all targets; i.e. it adds direction 
##to all pre-computed associations.
tni.cor<-function(x,tnet,estimator="spearman",dg=0, asInteger=TRUE, 
                  mapAssignedAssociation=TRUE){
  tnet <- abs(tnet)
  tfs <- colnames(tnet)
  tar <- rownames(tnet)
  ids <- unique(c(tfs,setdiff(tar,tfs)))
  x=x[ids,]
  x=t(x)
  #--
  pcorm=cor(x[,tfs],x[,tar], method=estimator,use="complete.obs")
  if(asInteger){
    pcorm[pcorm<0]=-1
    pcorm[pcorm>0]=1
  }
  if(length(tfs)>1)diag(pcorm[,tfs])=dg
  #--
  pcorm<-t(pcorm)
  colnames(pcorm)<-tfs
  if(mapAssignedAssociation)pcorm[tnet==0]=0
  return(pcorm)
}

##------------------------------------------------------------------------
# .alpha.adjust(nB=90, nA=900, alphaA = 0.001, betaA = 0.5)
.alpha.adjust <- function(nB, nA, alphaA, betaA){
  #--- find the effect size given nA, alphaA, and betaA
  pw <- 1-betaA
  rA <- tryCatch({
    pwr.r.test(n=nA, r=NULL, sig.level=alphaA, power=pw)$r
  }, error=function(e){ NULL })
  if(is.null(rA)){
    tp1 <- "NOTE: 'uniroot' could not find a root for the input arguments; "
    tp2 <- "please see 'tni.alpha.adjust' function documentation."
    stop(tp1, tp2, call. = FALSE)
  }
  #--- get a sequence of alpha values
  alpha <- c(seq(0.01, 0.1, by = 0.001), seq(0.11, 0.99, by = 0.01) )
  if(alphaA < 0.01){
    alpha <- c(alpha, seq(0.1, 0.5, by = 0.1) %o% 10^(-2:log10(alphaA)))
  }
  alpha <- sort(unique(alpha), decreasing = FALSE)
  #--- find size for each alpha value
  size <- round(sapply(alpha, function(a){
    tryCatch({
      pwr.r.test(n=NULL, r=rA, sig.level=a, power=pw)$n
    }, error=function(e){ NA })
  }))
  alpha <- alpha[!is.na(size)]
  size <- size[!is.na(size)]
  if(length(alpha)==0){
    tp1 <- "NOTE: 'uniroot' could not find a root for the input arguments; "
    tp2 <- "please see 'tni.alpha.adjust' function documentation."
    stop(tp1, tp2, call. = FALSE)
  }
  #--- find alpha for nB
  alphaB <- alpha[which.min(abs(size-nB))]
  return(alphaB)
}

# .alpha.adjust <- function(nB, nA, alphaA, betaA){
#   esizeA <- pwr.r.test(n=nA, r = NULL, sig.level=alphaA, power=1-betaA)$r
#   alphaB <- pwr.r.test(n=nB, r = esizeA, sig.level=NULL, 
#             power=1-betaA)$sig.level
#   return(alphaB)
# }

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##        -- COMPLEMENTARY FUNCTIONS FOR CDT MI ANALYSIS --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##------------------------------------------------------------------------
##compute delta mi from pooled null distributions, returns global mi markers
miThresholdMd <- function(gexp, nsamples, prob=c(0.05,0.95), 
                          nPermutations=1000, estimator="spearman", 
                          verbose=TRUE){
  if(length(prob)==1)prob=sort(c(1-prob,prob))
  if(nsamples>ncol(gexp))nsamples=ncol(gexp)
  ##build null distribution
  nsc<-10 #scale-up the rounds of this global permutation
  if(verbose)pb<-txtProgressBar(style=3)
  nullmark<-sapply(1:nPermutations, function(i){
    if(verbose)setTxtProgressBar(pb, i/nPermutations)
    rmil<-tni.pmin(gexp[,sample(ncol(gexp),nsamples)],
                   tfs=sample(rownames(gexp),nsc),
                   estimator=estimator, simplified=TRUE)
    rmih<-tni.pmin(gexp[,sample(ncol(gexp),nsamples)],
                   tfs=sample(rownames(gexp),nsc),
                   estimator=estimator, simplified=TRUE)
    quantile(rmih-rmil,probs=prob,na.rm=TRUE,names=FALSE)
  })
  if(verbose)close(pb)
  mimark<-c(median(nullmark[1,],na.rm=TRUE),
            median(nullmark[2,],na.rm=TRUE))
  return(mimark)
}

##------------------------------------------------------------------------
##compute delta mi, returns mi markers for each tf individually
miThresholdMdTf<-function(gexp, tfs, nsamples, prob=c(0.05,0.95), 
                          nPermutations=1000, 
                          estimator="spearman", verbose=TRUE){
  if(length(prob)==1)prob=sort(c(1-prob,prob))
  #---get low/high sample idxs
  nc<-ncol(gexp)
  idxl<-1:nsamples
  idxh<-(nc-nsamples-1):nc
  #---run permutation
  if(verbose)pb<-txtProgressBar(style=3)
  rmid<-sapply(1:nPermutations,function(i){
    if(verbose)setTxtProgressBar(pb, i/nPermutations)
    ridx<-sample(nc)
    rmil<-tni.pmin(gexp[,ridx[idxl]],tfs=tfs,estimator=estimator,
                   simplified=TRUE)
    rmih<-tni.pmin(gexp[,ridx[idxh]],tfs=tfs,estimator=estimator,
                   simplified=TRUE)
    apply(rmih-rmil,2,quantile,probs=prob,na.rm=TRUE,names=FALSE)
  })
  if(verbose)close(pb)
  idx<-c(1:(length(tfs)*2))%%2
  rmil<-rmid[idx==1,]
  rmih<-rmid[idx==0,]
  if(length(tfs)>1){
    dmimark<-cbind(apply(rmil,1,median,na.rm=TRUE),
                   apply(rmih,1,median,na.rm=TRUE))
    rownames(dmimark)<-tfs
  } else {
    dmimark<-c(median(rmid,na.rm=TRUE),median(rmid,na.rm=TRUE))
  }
  dmimark
}

##------------------------------------------------------------------------
##compute delta mi, returns mi null for each tf 
##randomization applyed on modulators (very conservative)
miMdTfStats<-function(gexp, regulons, nsamples, nPermutations=1000, 
                      estimator="spearman", minRegulonSize,
                      verbose=TRUE){
  tfs<-names(regulons)
  #---set minimun size to compute the null
  #---modulated targets < minRegulonSize should score zero
  sz<-unlist(lapply(regulons,length))
  if(any(sz<minRegulonSize)){
    for(tf in names(sz[sz==0])){
      regulons[[tf]]<-sample(rownames(gexp),minRegulonSize)
    }
  }
  #---get low/high sample idxs
  nc<-ncol(gexp)
  idxl<-1:nsamples
  idxh<-(nc-nsamples-1):nc
  idxp<-unique(c(tfs,as.character(unlist(regulons))))
  gexp<-gexp[idxp,,drop=FALSE]
  #---run permutation
  if(verbose)pb<-txtProgressBar(style=3)
  rmid<-sapply(1:nPermutations,function(i){
    if(verbose)setTxtProgressBar(pb, i/nPermutations)
    ridx<-sample(nc)
    rmil<-tni.pmin(gexp[,ridx[idxl]],tfs=tfs,estimator=estimator,
                   simplified=TRUE)
    rmih<-tni.pmin(gexp[,ridx[idxh]],tfs=tfs,estimator=estimator,
                   simplified=TRUE)
    sapply(tfs,function(tf){
      l<-rmil[regulons[[tf]],tf]
      h<-rmih[regulons[[tf]],tf]
      #median(h-l)/(sd(h)+sd(l)) 
      (median(h)-median(l))/(sd(h)+sd(l)) 
    })
  })
  if(verbose)close(pb)
  if(length(tfs)==1){
    rmid<-array(rmid,dim=c(1,length(rmid)))
    rownames(rmid)<-tfs
    med<-apply(rmid,1,median)
    std<-apply(rmid,1,sd)
    zscore<-rmid-med
    zscore<-zscore/std
  } else {
    rownames(rmid)<-tfs
    med<-apply(rmid,1,median)
    std<-apply(rmid,1,sd)
    zscore<-apply(rmid,2,'-',med)
    zscore<-apply( zscore, 2 ,'/',std)
  }
  list(dist=rmid,median=med,sd=std,zscore=zscore)
}

##------------------------------------------------------------------------
###check modulator stability
#1) computational cost forbids default use in the current implementation
#2) only one modulator each time
#3) intends to remove evident unstable modulations of selected modulators
cdt.stability<-function(gexp,estimator,mrkboot,miThreshold,spsz,md,tfs,
                        targetElements, sigDelta, nboot=100, consensus=95, 
                        verbose=TRUE){
  bcount<-sigDelta;bcount[,]<-0
  if(verbose)pb<-txtProgressBar(style=3)
  for(i in 1:nboot){
    if(verbose)setTxtProgressBar(pb, i/nboot)
    gxtemp<-gexp[,sample(ncol(gexp),replace=TRUE)]
    idxl<-1:spsz
    idxh<-(ncol(gxtemp)-spsz+1):ncol(gxtemp)
    idx<-sort.list(gxtemp[md,])
    idxl<-idx[idxl]
    idxh<-idx[idxh]
    mil<-tni.pmin(gxtemp[targetElements,idxl],tfs,estimator=estimator,
                  simplified=TRUE)
    mih<-tni.pmin(gxtemp[targetElements,idxh],tfs,estimator=estimator,
                  simplified=TRUE)
    mid<-mih-mil
    if(miThreshold=="md.tf" && length(tfs)>1){
      stabt<-t(apply(mid,1,"<",mrkboot[,1])) | 
        t(apply(mid,1,">",mrkboot[,2]))
    } else {
      stabt<- mid<mrkboot[1] | mid>mrkboot[2]
    }
    bcount <- bcount + stabt
  }
  if(verbose)close(pb)
  bcount>=as.integer( nboot * (consensus/100) )
}

##------------------------------------------------------------------------
checkModuationEffect<-function(gexp,tfs,modregulons,modulatedTFs,glstat,spsz,
                               minRegulonSize, pValueCutoff,nPermutations,
                               estimator,pAdjustMethod, count,verbose){
  for(md in names(modregulons)){
    mdreg<-modregulons[[md]][modulatedTFs]
    glstat$observed[[md]]$sig2noise<-
      glstat$observed[[md]]$sig2noise[modulatedTFs]
    #null dist
    glstat$null[[md]]<-miMdTfStats(gexp,regulons=mdreg, nsamples=spsz,
                                   nPermutations=nPermutations, 
                                   estimator=estimator, 
                                   minRegulonSize=minRegulonSize,
                                   verbose=verbose)
    glstat$null[[md]]$ci<-qnorm(1-(pValueCutoff/length(modulatedTFs)))
    #observed scores
    zscore<-glstat$observed[[md]]$sig2noise
    zscore<-(zscore-glstat$null[[md]]$median)/glstat$null[[md]]$sd
    glstat$observed[[md]]$zscore<-zscore
    pvals<-pnorm(abs(zscore), lower.tail=FALSE)
    glstat$observed[[md]]$pvals<-pvals
    glstat$observed[[md]]$adjpvals<-p.adjust(pvals, method=pAdjustMethod)
    glstat$observed[[md]]<-data.frame(glstat$observed[[md]])
  }
  #get glstatrev
  glstatrev<-list()
  mds<-names(glstat$observed)
  for(md in mds){
    observed<-glstat$observed[[md]]
    null<-glstat$null[[md]]
    for(tf in modulatedTFs){
      glstatrev$observed[[tf]]<-rbind(glstatrev$observed[[tf]],observed[tf,])
      glstatrev$null[[tf]]$dist<-rbind(glstatrev$null[[tf]]$dist,null$dist[tf,])
      glstatrev$null[[tf]]$median<-c(
        glstatrev$null[[tf]]$median,null$median[tf])
      glstatrev$null[[tf]]$sd<-c(glstatrev$null[[tf]]$sd,null$sd[tf])
      glstatrev$null[[tf]]$zscore<-rbind(glstatrev$null[[tf]]$zscore,
                                         null$zscore[tf,])
      glstatrev$null[[tf]]$ci<-qnorm(1-(pValueCutoff/length(mds)))
    }
  }
  for(tf in modulatedTFs){
    rownames(glstatrev$observed[[tf]])<-mds
    rownames(glstatrev$null[[tf]]$dist)<-mds
    names(glstatrev$null[[tf]]$median)<-mds
    names(glstatrev$null[[tf]]$sd)<-mds
    rownames(glstatrev$null[[tf]]$zscore)<-mds
    glstatrev$observed[[tf]]$adjpvals<-p.adjust(glstatrev$observed[[tf]]$pvals, 
                                                method=pAdjustMethod)
  }
  #update summary
  tp<-glstat$observed
  addstats<-lapply(tfs,function(tf){
    tpp<-count[[tf]]
    if(tf%in%modulatedTFs && nrow(tpp)>0){
      res<-t(sapply(rownames(tpp), function(md){
        c(tp[[md]][tf,"zscore"],tp[[md]][tf,"pvals"],NA)
      }))
      colnames(res)<-c("SNR","PvSNR","AdjPvSNR")
      res[,"SNR"]<-signif(res[,"SNR"],digits=3)
      tpp<-data.frame(tpp,res)
    } else {
      tpp<-data.frame(tpp,SNR=NA,PvSNR=NA,AdjPvSNR=NA)
    }
    tpp
  })
  names(addstats)<-tfs
  glstat$count<-addstats
  #--
  list(md2tf=glstat,tf2md=glstatrev)
}

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##          -- COMPLEMENTARY FUNCTIONS FOR PREPROCESSING --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##------------------------------------------------------------------------
##This function takes a row gene expression matrix, a named vector with
##matched ids and remove duplicated genes using the coefficient
##of variation (CV) to decide for the most informative probes.
cv.filter<-function(gexp, ids){ 
  # compute CV
  meangx<-apply(gexp,1,mean,na.rm=TRUE)
  sdgx<-apply(gexp,1,sd,na.rm=TRUE)
  cvgx<-abs(sdgx/meangx)
  # align ids
  ids <- data.frame(ids[rownames(gexp),,drop=FALSE], CV=cvgx, 
                    stringsAsFactors=FALSE)
  # remove probes without rowAnnotation
  ids <- ids[!is.na(ids[,2]),,drop=FALSE]
  ids <- ids[ids[,2]!="",,drop=FALSE]
  # remove duplicated genes by CV
  ids <- ids[sort.list(ids$CV,decreasing=TRUE),,drop=FALSE]
  unids<-unique(ids[,2])
  idx<-match(unids,ids[,2])
  ids<-ids[idx,,drop=FALSE]
  ids <- ids[sort.list(ids$CV,decreasing=FALSE),,drop=FALSE]
  # update and return gexp matrix
  gexp <- gexp[rownames(ids),,drop=FALSE]
  ids <- ids[,-ncol(ids),drop=FALSE]
  return(list(gexp=gexp,ids=ids))
}

##------------------------------------------------------------------------
##Simplified version of cv.filter
##This function takes a gene expression matrix and remove duplicated genes 
##using the coeff. variation (CV) to decide for the most informative probes.
##Input format -> rownames: probeID; col[1]: geneID, Col[2...n]:numeric data
# cv.filter.simp<-function(gexp){
#   x<-as.matrix(gexp[,-1])
#   ids<-data.frame(ProbeID=rownames(gexp),geneID=gexp[,1],
#                   stringsAsFactors=FALSE)
#   rownames(ids)<-ids$ID1
#   # compute CV
#   meangx<-apply(x,1,mean,na.rm=TRUE)
#   sdgx<-apply(x,1,sd,na.rm=TRUE)
#   cvgx<-sdgx/meangx
#   # aline ids
#   ids <- data.frame(ids, CV=cvgx, stringsAsFactors=FALSE)
#   # remove proves without rowAnnotation
#   ids <- ids[!is.na(ids[,2]),]
#   ids <- ids[ids[,2]!="",]
#   # remove duplicated genes by CV
#   ids <- ids[sort.list(ids$CV,decreasing=TRUE),]
#   unids<-unique(ids[,2])
#   idx<-match(unids,ids[,2])
#   ids<-ids[idx,]
#   # update and return gexp matrix
#   x<-x[rownames(ids),]
#   ids<-ids[,-ncol(ids)]
#   rownames(x)<-ids$geneID
#   x
# }


##------------------------------------------------------------------------
##------------------------------------------------------------------------
##              FUNCTIONS FOR REGULONS AND TNETS
##          -- IMPLEMENTS METHODS FOR REDER/IGRAPH --
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##------------------------------------------------------------------------
##return an adjacence matrix with regulons (tfs and targets)
tni.rmap<-function(tnet){
  tnet[tnet>0]=1;tnet[tnet<0]=-1
  tfs<-colnames(tnet)
  if(sum(!tfs%in%rownames(tnet))>0){
    idx<-!tfs%in%rownames(tnet)
    addtfs<-matrix(0,nrow=sum(idx),ncol=ncol(tnet))
    rownames(addtfs)<-tfs[idx]
    tnet<-rbind(addtfs,tnet)
  }
  idx<-apply(tnet!=0,1,sum)>0
  onlytar<-rownames(tnet)[idx]
  onlytar<-onlytar[!onlytar%in%tfs]
  tnet<-tnet[c(tfs,onlytar),,drop=FALSE]
  tmat<-t(tnet[onlytar,,drop=FALSE])
  zmat<-matrix(0,length(onlytar),length(onlytar))
  bmat<-rbind(tmat,zmat)
  rmap<-cbind(tnet,bmat)
  rmap[onlytar,tfs]<-0
  g<-igraph::graph.adjacency(rmap, diag=FALSE, mode="directed",
                             weighted=TRUE)
  edl<-igraph::get.edgelist(g)
  if(nrow(edl)>0){
    etype<-edl[,1]%in%tfs+edl[,2]%in%tfs
    etype[etype==2]<-0
    E(g)$modeOfAction<-E(g)$weight
    E(g)$modeOfAction[etype==0]=0
    E(g)$arrowType<-E(g)$modeOfAction
  }
  return(g)
}

##-----------------------------------------------------------------------------
##return an adjacence matrix for regulons
##using jaccard coefficient
tni.amap<-function(tnet, overlap="all"){
  if(overlap=="agreement"){
    tnet[tnet>0]<-1;tnet[tnet<0]<--1
    jc<-function(x,xmat){
      c<-x*xmat
      a<-colSums(c==1)
      c<-abs(x)+abs(xmat)
      b<-colSums(c!=0)
      b[b==0]=1
      a/b
    }
  } else {
    tnet[tnet!=0]=1
    jc<-function(x,xmat){
      c<-x+xmat
      a<-colSums(c==2)
      b<-colSums(c>0)
      b[b==0]=1
      a/b
    }
  }
  amap<-apply(tnet,2,jc,xmat=tnet)
  colnames(amap)<-colnames(tnet)
  diag(amap)=0
  return(amap)
}

##-----------------------------------------------------------------------------
##return an adjacence matrix
tni.mmap<-function(object,mnet,tnet,pvnet,othertfs,testedtfs,modulators){
  tnet[tnet<0]=-1;tnet[tnet>0]=1
  mnet<-mnet[,colnames(tnet),drop=FALSE]
  mnet<-mnet[rownames(tnet),,drop=FALSE]
  pvnet<-pvnet[,colnames(tnet),drop=FALSE]
  pvnet<-pvnet[rownames(tnet),,drop=FALSE]
  #---merge/simplify tnet and mnet
  mcnet<-mnet*10
  mcnet[mcnet==0]<-tnet[mcnet==0]
  #---get ids
  tfs<-colnames(mcnet)
  mod<-rownames(mcnet)
  onlytar<-setdiff(mod,tfs)
  #---check tfs
  if(sum(!tfs%in%rownames(mcnet))>0){
    idx<-!tfs%in%rownames(mcnet)
    addtfs<-matrix(0,nrow=sum(idx),ncol=ncol(mcnet))
    rownames(addtfs)<-tfs[idx]
    mcnet<-rbind(addtfs,mcnet)
    pvnet<-rbind(addtfs,pvnet)
  }
  #---set ordering
  mcnet<-mcnet[c(tfs,onlytar),,drop=FALSE]
  pvnet<-pvnet[c(tfs,onlytar),,drop=FALSE]
  #---
  elist<-NULL
  emode<-NULL
  etype<-NULL
  pvalue<-NULL
  tpnet<-mcnet[tfs,tfs,drop=FALSE]
  tpnet[tpnet!=0]=1
  for(i in 1:nrow(tpnet)){
    for(j in 1:ncol(tpnet)){
      tp<-tpnet[i,j] + tpnet[j,i]
      if(j>i){
        if(tp!=0){
          emode<-c(emode,0)
          etype<-c(etype,"TF.TF")
          elist<-rbind(elist,c(tfs[i],tfs[j]))
          pvalue<-c(pvalue,NA)
          emode<-c(emode,0)
          etype<-c(etype,"TF.TF")
          elist<-rbind(elist,c(tfs[j],tfs[i]))
          pvalue<-c(pvalue,NA)
        }
      }
    }
  }
  #others
  if(length(onlytar)>0){
    tpnet<-mcnet[onlytar,tfs,drop=FALSE]
    tpvnet<-pvnet[onlytar,tfs,drop=FALSE]
    for(i in 1:nrow(tpnet)){
      for(j in 1:ncol(tpnet)){
        tp<-tpnet[i,j]
        tpv<-tpvnet[i,j]
        if(abs(tp)==10){
          emode<-c(emode,tp)
          tpp<-ifelse(tp<0,"-Md.TF","+Md.TF")
          etype<-c(etype,tpp)
          elist<-rbind(elist,c(onlytar[i],tfs[j]))
          pvalue<-c(pvalue,tpv)
        } else if(abs(tp)==1){
          emode<-c(emode,tp)
          tpp<-ifelse(tp<0,"-TF.Target","+TF.Target")
          etype<-c(etype,tpp)
          elist<-rbind(elist,c(tfs[j],onlytar[i]))
          pvalue<-c(pvalue,tpv)
        }
      }
    }
  }
  #Correct TF-TF assigments derived from the MI analysis
  if(length(othertfs)>0){
    alltfs<-unique(c(tfs,othertfs))
    idx<-elist[,1]%in%alltfs+elist[,2]%in%alltfs == 2
    idx<-idx & abs(emode)==1
    emode[idx]<-0
    etype[idx]<-"TF.TF"
  }
  if(!is.null(emode)){
    idx<-abs(emode)==10
    emode[idx]<-emode[idx]/10
  }
  #---get graph and set mode
  if(is.null(elist)){
    g<-graph.empty(n=0, directed=TRUE)
  } else {
    g<-igraph::graph.edgelist(elist, directed=TRUE)
    E(g)$arrowType<-emode
    E(g)$pvalueModulation<-pvalue
    E(g)$interactionType<-etype
    E(g)$modeOfAction<-NA
    E(g)$modeOfAction[E(g)$interactionType=="-Md.TF"]<- -1
    E(g)$modeOfAction[E(g)$interactionType=="+Md.TF"]<- 1
    E(g)$modeOfAction[E(g)$interactionType=="TF.TF"]<- 0
    E(g)$modeOfAction[E(g)$interactionType=="-TF.Target"]<- -1
    E(g)$modeOfAction[E(g)$interactionType=="+TF.Target"]<- 1
  }
  #---add disconnected tfs 
  idx<-!tfs%in%V(g)$name
  if(any(idx))g<-g+vertices(tfs[idx])
  #---
  g<-setregs1(object,g,testedtfs,modulators)
  return(g)
}
setregs1<-function(object,g,testedtfs,modulators){
  #add rowAnnotation
  if(.hasSlot(object, "rowAnnotation")){
    tp <- dim(object@rowAnnotation)
    if(tp[1]>0 && tp[2]>1)g<-att.mapv(g=g,dat=object@rowAnnotation,refcol=1)
  }
  #set names if available
  if(!is.null(V(g)$SYMBOL)){
    g<-att.setv(g=g, from="SYMBOL", to='nodeAlias')
  } else {
    V(g)$nodeAlias<-V(g)$name
  }
  #map tested tfs and modulators
  V(g)$tfs<-as.numeric(V(g)$name%in%testedtfs)
  V(g)$modulators<-as.numeric(V(g)$name%in%modulators)
  #get regulon size (downstream targets)
  rsize<-tni.get(object,what="regulons")
  rsize<-sapply(rsize,length)
  rsize<-rsize[names(rsize)%in%V(g)$name]
  rsize<-rsize[names(rsize)%in%testedtfs]
  idx<-match(names(rsize),V(g)$name)
  V(g)$regulonSize<-NA
  V(g)$regulonSize[idx]<-rsize
  #---interno, so p/ customizar -- linearSize!!
  customSz=TRUE
  if(customSz){
    nquant=NULL
    tp1<-min(V(g)$regulonSize,na.rm=TRUE)
    tp3<-max(V(g)$regulonSize,na.rm=TRUE)
    tp1<-pretty(tp1,n=1,shrink=0.5)[1]
    tp3<-pretty(tp3,n=1,shrink=0.5)[2]
    tp2<-round((tp1+tp3)/2,0)
    breaks<-c(tp1,tp2,tp3)
  } else {
    breaks=NULL
    nq<-unique(V(g)$regulonSize)
    nq<-sum(!is.na(nq))
    nquant<-ifelse(nq>5,5,nq)
    nq<-length(unique(quantile(V(g)$regulonSize,na.rm=TRUE)))
    if(nq<3)nquant=NULL
  }
  xlim = c(25, 70, 10)
  #---
  g<-att.setv(g=g, from="regulonSize", to='nodeSize',xlim=xlim,roundleg=0,
              nquant=nquant,title="Downstream targets (n)",breaks=breaks)
  #set node attr
  V(g)$nodeFontSize=15
  V(g)$nodeLineWidth<-1.5
  V(g)$nodeColor<-"grey"
  V(g)$nodeLineColor<-"grey"
  V(g)$nodeFontSize[V(g)$tfs==1]=25
  #set node shape and legend
  g<-att.setv(g=g, from="tfs", to='nodeShape',shapes=c("DIAMOND","ELLIPSE"),
              title="", 
              categvec=c(0,1))
  g$legNodeShape$legend<-c("Modulator","TF")
  if(ecount(g)>0){
    #set edge attr
    E(g)$edgeWidth<-1.5
    if(any(E(g)$modeOfAction==0)){
      g<-att.sete(g=g, from="modeOfAction", to='edgeColor',
                  cols=c("#96D1FF","grey80","#FF8E91"), 
                  title="Mode",categvec=c(-1,0,1))
      g$legEdgeColor$legend<-c("Down","NA","Up")
    } else {
      g<-att.sete(g=g, from="modeOfAction", to='edgeColor',
                  cols=c("#96D1FF","grey80","#FF8E91"), 
                  title="Mode",categvec=c(-1,1))
      g$legEdgeColor$legend<-c("Down","Up")
    }
    #map upstring edges
    dg<-degree(g)
    el<-get.edgelist(g,names=FALSE)
    ew<-sapply(1:nrow(el),function(i){
      min(dg[el[i,1]],dg[el[i,2]])
    })
    E(g)$terminalEdge<-0
    E(g)$terminalEdge[ew==1]<-1
    #map modeOfAction to node attribute (compute the median mode of action)
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
    V(g)$nodeColor[V(g)$tfs==1]<-"white"
    V(g)$nodeLineColor[V(g)$tfs==1]<-"black"
  }
  g
}

##-----------------------------------------------------------------------------
##experimental code!!! (opem mmaps)
tni.mmap.detailed<-function(object,mnet,testedtfs,mds,ntop,nbottom){
  mnet.targets<-list()
  for(tf in testedtfs){
    tp1<-object@results$conditional$effect[[tf]]
    tp2<-mnet[,tf];tp2<-tp2[tp2!=0]
    tp3<-tp1[,names(tp2),drop=FALSE]
    if(ncol(tp3)>0){
      res<-sapply(names(tp2),function(md){
        tp<-tp3[,md]
        if(tp2[md]>0) tp[tp<0]<-0 else tp[tp>0]<-0
        tp
      })
      if("targets"%in%colnames(tp1)){
        rownames(res)<-tp1$targets
      } else {
        rownames(res)<-tp1$tagets
      }
      mnet.targets[[tf]]<-res
    }
  }
  #---
  gLists<-lapply(names(mnet.targets),function(tf){
    tar<-object@results$tn.dpi[,tf]
    mtar<-mnet.targets[[tf]]
    tpmds<-mds[mds%in%colnames(mtar)]
    glists<-lapply(tpmds,function(md){
      mtartp<-mtar[,md]
      mode<-mnet[md,tf]
      tartp<-tar[names(mtartp)]
      elist<-cbind(A=c(md,tf,rep("TF|Md",length(tartp))),
                   B=c("TF|Md","TF|Md",names(tartp)))
      rownames(elist)<-NULL
      #---
      modeOfAction<-c(NA,NA,tartp)
      modulationEffect<-c(mode,1,mtartp)
      modulationType<-modulationEffect
      modulationType[modulationType!=0]<-1
      modulationType[1:2]<-2
      #---
      interactionType<-rep("TF.Target",length(modeOfAction))
      interactionType[1:2]<-"TF|Md"
      elist<-data.frame(elist,interactionType,modeOfAction,modulationEffect,
                        modulationType,stringsAsFactors=FALSE)
      idx<-sort.list(abs(elist$modeOfAction),na.last=FALSE,decreasing=FALSE)
      elist<-elist[idx,,drop=FALSE]
      if(!is.null(ntop)){
        tp<-abs(elist$modulationEffect)
        tp[1:2]<-NA
        idx<-sort.list(tp,na.last=FALSE,decreasing=FALSE)
        elist<-elist[idx,,drop=FALSE]
        idx<-c(1,2,rev(3:nrow(elist)))[1:(ntop+2)]
        if(!is.null(nbottom)){
          idx<-c(idx, 3:(nbottom+2) )
        }
        elist<-elist[sort(idx),,drop=FALSE]
      }
      #---
      g<-igraph::graph.edgelist(as.matrix(elist[,1:2]), directed=TRUE)
      #---
      xy<-drcircle(nv=vcount(g)-3)
      hd<-data.frame(x=c(-0.14,0,0.2),y=c(-0.9,-0.7,-0.8))
      xy<-rbind(hd,xy)
      xy$name<-V(g)$name
      #---resort
      idx<-abs(elist$modeOfAction)
      idx[elist$modulationType==0]<-NA
      idx<-sort.list(idx,na.last=FALSE,decreasing=FALSE)
      elist<-elist[idx,,drop=FALSE]
      #---
      g<-igraph::graph.edgelist(as.matrix(elist[,1:2]), directed=TRUE)
      #---set node xy   
      V(g)[xy$name]$coordX<-xy$x
      V(g)[xy$name]$coordY<-xy$y
      #---set node type
      V(g)$nodeType<-rep("Target",length(V(g)$name))
      idx<-match(c(md,"TF|Md",tf),V(g)$name)
      V(g)$nodeType[idx]<-c("Md","TF|Md","TF")
      #---
      E(g)$interactionType<-elist$interactionType
      E(g)$modeOfAction<-elist$modeOfAction
      E(g)$modulationEffect<-elist$modulationEffect
      E(g)$modulationType<-elist$modulationType
      #---
      g
    })
    names(glists)<-names(tpmds)
    glists
  })
  # names(gLists)<-names(mnet.targets)
  names(gLists)<-names(testedtfs)
  #---
  gLists<-lapply(gLists,function(gg){
    gg<-lapply(gg,function(g){
      if(.hasSlot(object, "rowAnnotation") && nrow(object@rowAnnotation)>0){
        g<-att.mapv(g=g,dat=object@rowAnnotation,refcol=1)
      } else {
        V(g)$SYMBOL<-V(g)$name
      }
      if(ecount(g)>0)g<-setregs2(g)
    })
    gg
  })
  gLists
}
setregs2<-function(g){
  g<-att.setv(g=g, from="SYMBOL", to='nodeAlias')
  V(g)$nodeAlias[V(g)$name=="TF|Md"]<-"TF | Md"
  V(g)$nodeColor<-"white"
  g<-att.setv(g=g, from="nodeType", to='nodeLineColor',
              categvec=c("Md","TF","TF|Md","Target"),
              cols=c("#FF6666","#66CCFF","#FF9900","#00CC99"),title="")
  g$legNodeLineColor$legend<-c("Md","TF","TF | Md", "Mod. target")
  #--set legNodeColor
  g$legNodeColor<-g$legNodeLineColor
  g<-remove.graph.attribute(g,"legNodeLineColor")
  #--
  V(g)$nodeSize<-30
  V(g)$nodeSize[V(g)$nodeType=="Target"]<-20
  V(g)$nodeFontSize<-12
  V(g)$nodeFontSize[V(g)$nodeType!="Target"]<-16
  V(g)$nodeLineWidth<-2.8
  #set edge attr
  E(g)$edgeColor<-"grey65"
  E(g)$edgeWidth<-2
  E(g)$edgeWidth[E(g)$modulationType==0]<-2
  #--
  g<-att.sete(g=g, from="modulationType", to='edgeColor', title="", 
              categvec=c(2,1,0),cols=c("black","grey50","grey70"))
  #--
  #E(g)$absModulationEffect<-sign(E(g)$modulationEffect)
  # g<-att.sete(g=g, from="absModulationEffect", to='edgeColor', title="",  
  #             categvec=-1:1, cols=c("grey65","grey95","grey65"))
  #--
  E(g)$edgeColor[E(g)$modulationType==2]<-"black"
  g$legEdgeColor$legend<-c("conditioned","modulated","non-modulated")
  #--
  el<-data.frame(get.edgelist(g),E(g)$modulationType,stringsAsFactors=FALSE)
  ndm<-el[el[,3]==0,2]
  V(g)$nodeLineColor[V(g)$name%in%ndm]<-"grey60"
  V(g)$nodeFontColor<-"black"
  V(g)$nodeFontColor[V(g)$name%in%ndm]<-"grey40"
  g
}
drcircle<-function (nv, ybend=1, xbend=1, ang=1, rotate=0, plot=FALSE) {
  angle.inc <- ang * pi/(nv-1)
  angles <- seq(0, ang*pi, by=angle.inc)+rotate*pi
  radius<-0.5
  xv <- cos(angles) * radius * xbend
  yv <- sin(angles) * radius * ybend
  if(plot){
    plot(-1:1,-1:1,type="n",xlab="",ylab="",main="")
    polygon(xv, yv)
  }
  data.frame(x=xv,y=yv)
}
#drcircle(nv=10,ang=1,plot=TRUE)

#---------------------------------------------------------------
#test overlap among regulons in a tnet via phyper function
tni.phyper<-function(tnet){
  tnet[tnet!=0]<-1
  tnet<-tnet[rowSums(tnet)>0,,drop=FALSE]
  phtest<-function(x,xmat){
    c<-x+xmat
    pop<-length(x)
    sample1<-sum(x>0)
    sample2<-colSums(xmat>0)
    overlap<-colSums(c==2)
    resph<-phyper(overlap-1, sample1, pop - sample1, sample2, 
                  lower.tail=FALSE)
    resph
  }
  pmat<-apply(tnet,2,phtest,xmat=tnet)
  colnames(pmat)<-colnames(tnet)
  diag(pmat)=1
  pmat<-matrix(p.adjust(pmat),ncol=ncol(pmat),nrow=nrow(pmat),
               dimnames=list(rownames(pmat),colnames(pmat)))
  return(pmat)
}

#---------------------------------------------------------------
cdt.table <- function(cdt){
  cdt.tb <- NULL
  for(i in names(cdt)){
    tp <- cdt[[i]]
    cdt.tb <- rbind(cdt.tb,tp)
  }
  rownames(cdt.tb)<-NULL
  return(cdt.tb)
}
p.adjust.cdt.table <- function(cdt.tb, pAdjustMethod="bonferroni"){
  cdt.tb$AdjPvFET <- p.adjust(cdt.tb$PvFET, method=pAdjustMethod)
  cdt.tb$AdjPvKS <- p.adjust(cdt.tb$PvKS, method=pAdjustMethod)
  if(!is.null(cdt.tb$PvSNR)){
    cdt.tb$AdjPvSNR <- p.adjust(cdt.tb$PvSNR, method=pAdjustMethod)
  }
  return(cdt.tb)
}
p.format.cdt.table <- function(cdt.tb, pAdjustMethod="bonferroni"){
  cdt.tb$PvFET <- signif(cdt.tb$PvFET, digits = 3)
  cdt.tb$AdjPvFET <- signif(cdt.tb$AdjPvFET, digits = 3)
  cdt.tb$PvKS <- signif(cdt.tb$PvKS, digits = 3)
  cdt.tb$AdjPvKS <- signif(cdt.tb$AdjPvKS, digits = 3)
  if(!is.null(cdt.tb$R))cdt.tb$R <- round(cdt.tb$R, digits = 4)
  if(!is.null(cdt.tb$PvSNR)){
    cdt.tb$PvSNR <- signif(cdt.tb$PvSNR, digits = 3)
    cdt.tb$AdjPvSNR <- signif(cdt.tb$AdjPvSNR, digits = 3)
  }
  return(cdt.tb)
}
#---------------------------------------------------------------
cdt.list<-function(cdt,pAdjustMethod="bonferroni"){
  cdt<-p.adjust.cdt.list(cdt=cdt,pAdjustMethod=pAdjustMethod,
                         p.name="PvFET",adjp.name="AdjPvFET")
  cdt<-p.adjust.cdt.list(cdt=cdt,pAdjustMethod=pAdjustMethod,
                         p.name="PvKS",adjp.name="AdjPvKS",sort.name="PvKS")
  if(!is.null(cdt$PvSNR)){
    cdt<-p.adjust.cdt.list(cdt=cdt,pAdjustMethod=pAdjustMethod,
                           p.name="PvSNR",adjp.name="AdjPvSNR",
                           sort.name="PvSNR",
                           global=FALSE)
  }
  cdt
}
p.adjust.cdt.list<-function(cdt,pAdjustMethod="bonferroni",p.name="Pvalue",
                       adjp.name="AdjustedPvalue",sort.name=NULL,
                       decreasing=FALSE, roundpv=TRUE, global=TRUE){
  if(global){
    #get pvalues
    pvals<-array(,dim=c(1,2))
    for(i in 1:length(cdt)){
      tp<-cdt[[i]]
      if(nrow(tp)>0 && !is.null(tp[[p.name]]))
        pvals<-rbind(pvals,cbind(i,tp[[p.name]]))
    }
    pvals<-pvals[-1,,drop=FALSE]
    colnames(pvals)<-c("idx","p")
    if(nrow(pvals)==0)return(cdt)
    #adjust pvals
    adjp<-pvals
    adjp[,"p"]<-p.adjust(adjp[,"p"],method=pAdjustMethod)
    adjp[,"p"]<-signif(adjp[,"p"], digits=3)
    pvals[,"p"]<-signif(pvals[,"p"], digits=3)
    #update adjusted pvalues
    idx<-unique(adjp[,1])
    for(i in idx){
      tp<-cdt[[i]]
      tp[[adjp.name]]<-adjp[adjp[,1]==i,2]
      if(roundpv)tp[[p.name]]<-pvals[pvals[,1]==i,2]
      cdt[[i]]<-tp
    }
  } else {
    for(i in 1:length(cdt)){
      tp<-cdt[[i]]
      if(nrow(tp)>0 && !is.null(tp[[p.name]])){
        tp[[adjp.name]]<-p.adjust(tp[[p.name]],method=pAdjustMethod)
        tp[[adjp.name]]<-signif(tp[[adjp.name]], digits=3)
        if(roundpv)tp[[p.name]]<-signif(tp[[p.name]], digits=3)
      }
      cdt[[i]]<-tp
    }
  }
  #sort
  if(is.character(sort.name)){
    for(i in 1:length(cdt)){
      tp<-cdt[[i]]
      if(nrow(tp)>0 && !is.null(tp[[p.name]])){
        tp<-tp[sort.list(tp[[sort.name]],decreasing=decreasing),]
        cdt[[i]]<-tp
      }
    }
  }
  cdt
}
sortblock.cdt.list<-function(cdt,coln="PvFET"){
  #sort blocks
  idx1<-unlist(lapply(cdt,nrow))
  idx1<-idx1/max(idx1)+1
  idx2<-lapply(1:length(cdt),function(i){
    tp<-cdt[[i]][[coln]]
    if(length(tp)>0 && sum(tp,na.rm=TRUE)>0){
      tp<--log10(mean(tp,na.rm=TRUE))
    } else {
      tp<-0
    }
    tp
  })
  idx2<-unlist(idx2)
  idx2<-idx2/sum(idx2)
  idx<-sort.list(idx1+idx2,decreasing=TRUE)
  cdt<-cdt[idx]
  cdt
}

##------------------------------------------------------------------------------
# GEO2R auto-detect checks and log2 transformation
.isUnloggedData <- function(gexp){
  qx <- as.numeric(quantile(gexp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), 
                            na.rm=TRUE))
  LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || 
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  LogC
}
.log2transform<-function(gexp){
  gexp[which(gexp <= 0)] <- NaN
  gexp <- log2(gexp) 
  gexp[is.nan(gexp)] <- 0
  gexp
}

#---------------------------------------------------------------
#reverse results from conditional analyses, from TF-MD to MD-TF
# cdt.list.rev<-function(cdt,pAdjustMethod="bonferroni"){
#   cdtrev<-list()
#   for(i in 1:length(cdt)){
#     tf<-names(cdt)[i]
#     tp<-cdt[[i]]
#     if(nrow(tp)>0){
#       mds<-rownames(tp)
#       for(md in mds){
#         tpp<-tp[md,c(2,1,3:ncol(tp))]
#         rownames(tpp)<-tf
#         cdtrev[[md]]<-rbind(cdtrev[[md]],tpp)
#       }
#     }
#   }
#   if(length(cdtrev)>0){
#     cdtrev<-p.adjust.cdt.list(cdt=cdtrev,pAdjustMethod=pAdjustMethod, 
#                               p.name="PvFET",adjp.name="AdjPvFET")
#     cdtrev<-p.adjust.cdt.list(cdt=cdtrev,pAdjustMethod=pAdjustMethod, 
#                               p.name="PvKS",adjp.name="AdjPvKS",
#                               sort.name="PvKS")
#     cdtrev<-p.adjust.cdt.list(cdt=cdtrev,pAdjustMethod=pAdjustMethod, 
#                               p.name="PvSNR",adjp.name="AdjPvSNR",
#                               sort.name="PvSNR",global=FALSE)
#   }
#   cdtrev
# }

##-----------------------------------------------------------------------------
##build an igraph object from hclust
hclust2igraph<-function(hc){
  if(!is(hc,"hclust"))stop("'hc' should be an 'hclust' object!")
  if(is.null(hc$labels))hc$labels=as.character(sort(hc$order))
  #options to remove nests
  rmnodes<-NULL
  if(!is.null(hc$cutoff)){
    #remove nests based on length cutoff
    tmap<-treemap(hc)
    hcNodes<-tmap$hcNodes
    hcNests<-hcNodes[hcNodes$type=="nest",,drop=FALSE]
    hcEdges<-tmap$hcEdges
    nestList<-tmap$nest
    rmnodes<-hcEdges[hcEdges$edgeLength<hc$cutoff,,drop=FALSE]
    rmnodes<-unique(rmnodes$parentNode)
    rmnodes<-hcNodes$hcId[hcNodes$node%in%rmnodes]
    mergeBlocks=FALSE
  } else if(!is.null(hc$keep)){
    #option to assign parent nodes to be preserved
    ids<-as.numeric(hc$merge)
    ids<-sort(ids[ids>0])
    rmnodes<-ids[!ids%in%hc$keep]
    mergeBlocks=TRUE
  } else if(!is.null(hc$remove)){
    #option to assign parent nodes to be removed
    ids<-as.numeric(hc$merge)
    ids<-sort(ids[ids>0])
    rmnodes<-ids[ids%in%hc$remove]
    mergeBlocks=FALSE
  }
  #build igraph and return assigments
  .hclust2igraph(hc,rmnodes,mergeBlocks)
}
##-----------------------------------------------------------------------------
##build an igraph object from pvclust
# pvclust2igraph<-function(hc,alpha=0.95,max.only=FALSE){
#   if(!is(hc,"pvclust"))stop("'hc' should be an 'pvclust' object!")
#   if(is.null(hc$hclust$labels))hc$hclust$labels=as.character(
#        sort(hc$hclust$order))
#   keep<-pvpick(hc,alpha=alpha,max.only=max.only)$edges
#   ids<-as.numeric(hc$hclust$merge)
#   ids<-sort(ids[ids>0])
#   rmnodes<-ids[!ids%in%keep]
#   .hclust2igraph(hc$hclust,rmnodes,mergeBlocks=TRUE)
# }
##-----------------------------------------------------------------------------
##wrap-up  for 'hclust2igraph' function
.hclust2igraph<-function(hc,rmnodes=NULL, mergeBlocks=TRUE){
  #get treemap
  tmap<-treemap(hc)
  hcNodes<-tmap$hcNodes
  hcNests<-hcNodes[hcNodes$type=="nest",,drop=FALSE]
  hcEdges<-tmap$hcEdges
  nestList<-tmap$nest
  nestmap<-tmap$nestmap
  lineage<-tmap$lineage
  rootid<-tmap$rootid
  #set nests do be removed
  if(!is.null(rmnodes) && length(rmnodes)>0){
    rmnodes<-hcNodes[hcNodes$mergeId%in%rmnodes,,drop=FALSE]
    rmnodes<-rmnodes$node[!rmnodes$node==rootid]
    #update hcEdges and nestList
    if(length(rmnodes)>0){
      hcEdges<-hcEdges.filter(hcEdges,hcNodes,rmnodes,lineage,
                              rootid,mergeBlocks)
      nestList<-nestList[names(nestList)%in%hcEdges$parentNode]
    }
  }
  #build igraph and return assigments
  g<-graph.edgelist(as.matrix(hcEdges[,1:2]), directed=TRUE)
  #E(g)$edgeWeight<-hcEdges$edgeLength
  #tp<-hcEdges$parentHeight-min(hcEdges$parentHeight)
  #E(g)$edgeWeight<-(max(tp)-tp)/max(tp) * 100
  E(g)$arrowType<-0
  list(g=g,nest=nestList)
}
hcEdges.filter<-function(hcEdges,hcNodes,rmnodes,lineage,rootid, 
                         mergeBlocks=TRUE){
  getpar<-function(hcEdges,childNode){
    idx<-which(hcEdges$childNode==childNode)
    hcEdges$parentNode[idx]
  }
  #update rmnodes with any continous bottom-up struture not in rmnodes
  if(mergeBlocks){
    nids<-hcNodes$node[hcNodes$type=="nest"]
    nids<-nids[!nids%in%rmnodes]
    for(i in 1:length(nids)){
      pn<-getpar(hcEdges,nids[i])
      if(length(pn)>0 && !pn%in%rmnodes && !pn%in%rootid){
        rmnodes<-c(rmnodes,nids[i])
      }
    } 
  }
  #update lineage
  lineage<-lapply(lineage,setdiff,rmnodes)
  #get the correct order
  rmNests<-hcNodes[hcNodes$node%in%rmnodes,,drop=FALSE]
  rmEdges<-hcEdges[hcEdges$childNode%in%rmnodes,,drop=FALSE]
  #get filtered hcEdges
  hcEdges<-hcEdges[!hcEdges$childNode%in%rmNests$node,,drop=FALSE]
  #update parents'ids
  for(i in 1:nrow(hcEdges)){
    newpn<-lineage[[hcEdges$childNode[i]]][1]
    hcEdges$parentNode[i]<-newpn
  }
  rownames(hcEdges)<-NULL
  hcEdges
}
##-----------------------------------------------------------------------------
##build nested maps from hclust objects
treemap<-function(hc){
  A=hc$merge
  B=list()
  C=list()
  D=list()
  E=list()
  nest=list()
  if(is.null(hc$labels))hc$labels=as.character(sort(hc$order))
  for(i in 1:nrow(A)){
    ai=A[i,1]
    if(ai < 0){
      B[[i]]= -ai
      C[[i]]=1
    } else {
      B[[i]]=B[[ai]]      
      C[[i]]=C[[ai]]+1 
    }
    ai=A[i,2]
    if(ai < 0){
      B[[i]]=sort(c(B[[i]],-ai))
    } else {
      B[[i]]=sort(c(B[[i]],B[[ai]]))
      C[[i]]=max(C[[i]],C[[ai]]+1)
    }
    p=match(i,A)
    D[[i]]=ifelse(p>nrow(A),p-nrow(A),p)
    nest[[i]]=hc$labels[B[[i]]]
  }
  D[[nrow(A)]]=nrow(A)+1
  for(i in 1:nrow(A)){
    step=1
    find=D[[i]]  
    while(find<D[[nrow(A)]]){
      find=D[[find]]
      step=step+1
    }
    E[[i]]=step
  }
  # get dendogram xy position
  nn=nrow(A) + 1
  xaxis=c()
  yaxis=hc$height
  tp=rep(0,2)
  mm=match(1:length(hc$order),hc$order)
  for(i in 1:(nn-1)) {
    ai=A[i,1]
    if(ai < 0){
      tp[1]=mm[-ai]
    } else {
      tp[1]=xaxis[ai]
    }
    ai=A[i,2]
    if(ai < 0){
      tp[2]=mm[-ai]
    } else {
      tp[2]=xaxis[ai]
    }
    xaxis[i]=mean(tp)
  }
  xyaxis=data.frame(xaxis=xaxis,yaxis=yaxis,stringsAsFactors=FALSE)
  # return res
  C=as.numeric(C)
  D=as.numeric(D)
  E=as.numeric(E)
  N=hc$merge>0
  N=N[,1]+N[,2]
  obj<-list(nest=nest,compids=B,labels=hc$labels,parent=D,leafdist=C,
            rootdist=E,height=hc$height,nnest=N, xyaxis=xyaxis)
  #---get unified edges
  N<-nrow(hc$merge);nn<-N+1
  hcEdges<-NULL
  eLength<-NULL
  for(i in 1:N){
    y1<-hc$merge[i,1]
    y2<-hc$merge[i,2]
    if(y1>0){
      l1<-hc$height[i] - hc$height[y1]
    } else {
      l1<-hc$height[i]
    }
    if(y2>0){
      l2<-hc$height[i] - hc$height[y2]
    } else {
      l2<-hc$height[i]
    }    
    tp<-cbind(rbind(c(i,y1),c(i,y2)),c(l1,l2))
    hcEdges<-rbind(hcEdges,tp)
    NULL
  }
  colnames(hcEdges)<-c("parentNode","childNode","edgeLength")
  hcEdges<-data.frame(hcEdges,stringsAsFactors=FALSE)
  hcEdges$parentHeight<-obj$height[hcEdges$parentNode]
  #---get unified nodes
  hcl<-data.frame(node=hc$labels,mergeId=-c(1:nn),hcId=c(1:nn), type="leaf",
                  stringsAsFactors=FALSE)
  hcn<-data.frame(node=paste("N",c(1:N),sep=""),mergeId=c(1:N),hcId=c(1:N),
                  type="nest",stringsAsFactors=FALSE)
  hcNodes<-rbind(hcl,hcn)
  hcEdges$parentNode<-hcNodes$node[match(hcEdges$parentNode,hcNodes$mergeId)]
  hcEdges$childNode<-hcNodes$node[match(hcEdges$childNode,hcNodes$mergeId)]
  obj$hcNodes<-hcNodes;obj$hcEdges<-hcEdges
  names(obj$nest)<-paste("N",c(1:N),sep="")
  #---get mest maps (optional)
  nestmap<-obj$nest
  nestall<-obj$nest
  nts<-names(obj$nest)
  for(i in 1:length(nts)){
    tp<-unlist(lapply(lapply(obj$nest,"%in%",obj$nest[[i]]),all))
    tp<-tp[-i]
    tp<-nts[which(tp)]
    nestmap[[i]]<-tp
    nestall[[i]]<-c(nestall[[i]],tp)
  }
  obj$nestmap<-nestmap
  obj$nestall<-nestall
  #---get node lineage
  getpar<-function(hcEdges,childNode){
    idx<-which(hcEdges$childNode==childNode)
    hcEdges$parentNode[idx]
  }
  root<-hcEdges$parentNode[length(hcEdges$parentNode)]
  lineage<-list()
  for(i in 1:nrow(obj$hcEdges)){
    ch<-obj$hcEdges$childNode[i]
    pn<-obj$hcEdges$parentNode[i]
    lineage[[ch]]<-pn
    while(pn!=root){
      pn<-getpar(obj$hcEdges,pn)
      lineage[[ch]]<-c(lineage[[ch]],pn)
    }
  }
  obj$lineage<-lineage
  obj$rootid<-root
  return(obj)
}

##-----------------------------------------------------------------------------
.run.tni.gsea2.alternative <- function(listOfRegulonsAndMode, 
                                       phenotype, phenorank, exponent, 
                                       alternative){
  ##-----get regulons
  listOfRegulonsUp <- lapply(listOfRegulonsAndMode, function(reg){
    as.numeric(names(reg[reg>0]))
  })
  listOfRegulonsDown <- lapply(listOfRegulonsAndMode, function(reg){
    as.numeric(names(reg[reg<0]))
  })
  ##-----set alternative hypothesis
  ## ...observe that greater/less refer to regulon activity
  if(alternative=="two.sided"){
    alternative_up <- alternative_down <- alternative
  } else {
    alternative_up <- ifelse(alternative=="greater","greater", "less")
    alternative_down <- ifelse(alternative=="greater","less", "greater")
  }
  GSEA2.results.up <- sapply(listOfRegulonsUp, .fgseaScores4TNI, 
                             phenotype=phenotype, phenorank=phenorank, 
                             exponent=exponent, alternative=alternative_up)
  GSEA2.results.down <- sapply(listOfRegulonsDown, .fgseaScores4TNI, 
                               phenotype=phenotype, phenorank=phenorank, 
                               exponent=exponent, alternative=alternative_down)
  b1<-length(GSEA2.results.up)>0 && length(GSEA2.results.down)>0
  b2<-length(GSEA2.results.up)==length(GSEA2.results.down)
  if(b1 && b2) {
    GSEA2.results.both <- GSEA2.results.up-GSEA2.results.down
  } else {
    GSEA2.results.up <- numeric()
    GSEA2.results.down <- numeric()
    GSEA2.results.both <- numeric()
  }
  #-----format, pack and return results
  GSEA2.results.up<-round(GSEA2.results.up,4)
  GSEA2.results.down<-round(GSEA2.results.down,4)
  GSEA2.results.both<-round(GSEA2.results.both,4)
  GSEA2.results<-list(positive=GSEA2.results.up,
                      negative=GSEA2.results.down,
                      differential=GSEA2.results.both)
  return( GSEA2.results )
}

##-----------------------------------------------------------------------------
## Used only for .annotate.samples.gsea
.run.tni.gsea1.alternative <- function(geneSetList, phenotype, phenorank, 
                                       exponent, alternative){
  GSEA2.results <- sapply(geneSetList, .fgseaScores4TNI, 
                          phenotype=phenotype, phenorank=phenorank, 
                          exponent=exponent, alternative=alternative)
  if(length(GSEA2.results)==0) {
    GSEA2.results <- numeric()
  }
  #-----format, pack and return results
  GSEA2.results <- round(GSEA2.results,4)
  return( GSEA2.results )
}

##---------------------------------------------------------
.fgseaScores4TNI <- function(geneSet, phenotype, phenorank, 
                             exponent=1, alternative="two.sided"){
  if(length(geneSet)>0){
    nh <- length(geneSet)
    N <- length(phenotype)
    Phit <- rep(0, N)
    Pmiss <- rep(0, N)
    hits <- phenorank[geneSet]
    Phit[hits] <- abs(phenotype[geneSet])^exponent
    NR <- sum(Phit)
    Pmiss[-hits] <- 1/(N-nh)
    Phit <- cumsum(Phit/NR)
    Pmiss <- cumsum(Pmiss)
    runningES <- Phit-Pmiss
    ESmax <- max(runningES)
    ESmin <- min(runningES)
    if(alternative=="greater"){
      ES <- ESmax
    } else if(alternative=="less"){
      ES <- ESmin
    } else {
      ES <- ifelse(abs(ESmin)>abs(ESmax), ESmin, ESmax)
    }
  } else {
    ES <- 0
  }
  return(ES)
}

#---------------------------------------------------------------
.get.regulon.summary <- function(object){
  ref <- tni.get(object, what = "refregulons")
  if(length(ref)>0){
    ref <- summary(unlist(lapply(ref, length)))
  } else {
    ref <- summary(0)
    ref[] <- NA
  }
  dpi <- tni.get(object, what = "regulons")
  if(length(dpi)>0){
    dpi <- summary(unlist(lapply(dpi, length)))
  } else {
    dpi <- summary(0)
    dpi[] <- NA
  }
  rbind(tnet.ref=ref,tnet.dpi=dpi)
}
#---------------------------------------------------------------
.get.tnet.summary <- function(object){
  tnetsumm <- object@summary$results$tnet
  if(object@status["Permutation"]=="[x]"){
    bin<-object@results$tn.ref
    bin[bin!=0]<-1
    tnetsumm[1,]<-c(ncol(bin),sum(rowSums(bin)>0),sum(bin))
  }
  if(object@status["DPI.filter"]=="[x]"){
    bin<-object@results$tn.dpi
    bin[bin!=0]<-1
    tnetsumm[2,]<-c(ncol(bin),sum(rowSums(bin)>0),sum(bin))
  }
  return(tnetsumm)
}

##------------------------------------------------------------------------------
.tni.stratification.gsea2 <- function(regulonActivity, sections=1, 
                                      center=TRUE){
  regstatus <- sign(regulonActivity$dif)
  for (reg in colnames(regstatus)){
    sq <- c(seq_len(sections))
    pos <- regulonActivity$pos[, reg]
    neg <- regulonActivity$neg[, reg]
    dif <- regulonActivity$dif[, reg]
    #---
    regstatus[sign(pos) == sign(neg), reg] <- 0
    tp <- regstatus[, reg]
    #---
    tp1 <- sort(dif[tp > 0], decreasing = TRUE)
    tp1[] <- rep(sq, each = ceiling(length(tp1)/sections), 
                 length.out = length(tp1))
    regstatus[names(tp1), reg] <- tp1
    #---
    tp2 <- sort(dif[tp < 0], decreasing = TRUE)
    tp2[] <- rep(sq + sections + 1, each = ceiling(length(tp2)/sections), 
                 length.out = length(tp2))
    regstatus[names(tp2), reg] <- tp2
  }
  mid <- sections + 1
  regstatus[regstatus == 0] <- mid
  #---
  if(center){
    regstatus <- -1 * (regstatus - mid)
    mid <- 0
  }
  #---
  regulonActivity$status <- regstatus
  regulonActivity$sections <- sections
  regulonActivity$center <- mid
  return(regulonActivity)
}

##------------------------------------------------------------------------------
.tni.stratification.area <- function(regulonActivity, sections=1, 
                                     center=FALSE){
  regstatus <- sign(regulonActivity$dif)
  for (reg in colnames(regstatus)){
    sq <- c(seq_len(sections))
    dif <- regulonActivity$dif[, reg]
    tp <- regstatus[, reg]
    #---
    tp1 <- sort(dif[tp > 0], decreasing = TRUE)
    tp1[] <- rep(sq, each = ceiling(length(tp1)/sections), 
                 length.out = length(tp1))
    regstatus[names(tp1), reg] <- tp1
    #---
    tp2 <- sort(dif[tp < 0], decreasing = TRUE)
    tp2[] <- rep(sq + sections + 1, each = ceiling(length(tp2)/sections), 
                 length.out = length(tp2))
    regstatus[names(tp2), reg] <- tp2
  }
  #--- obs. this stratification generates a 'midle' group 
  #--- only when there are samples with regulon activity assigned with 0 or NA
  mid <- sections + 1
  regstatus[regstatus == 0 | is.na(regstatus)] <- mid
  #---
  if(center){
    regstatus <- -1 * (regstatus - mid)
    mid <- 0
  }
  #---
  regulonActivity$status <- regstatus
  regulonActivity$sections <- sections
  regulonActivity$center <- mid
  return(regulonActivity)
}

#---------------------------------------------------------------
#--- simple function to generate gmm scores
.mi2gmm.dpi <- function(object, idkey=NULL){
  #--- set input
  regs <- tni.get(object, what="regulatoryElements", idkey=idkey)
  #--- get tf-tar mi
  tftar_mi <- lapply(tni.get(object, what="regulons.and.mode", 
                             idkey=idkey), abs)
  #--- get tf-tar correlation
  object@results$tn.dpi <- tni.cor(object@gexp, object@results$tn.dpi, 
                                   asInteger=FALSE,
                                   estimator=object@para$perm$estimator)
  tftar_cor <- tni.get(object, what="regulons.and.mode", idkey=idkey)
  #--- get tf-tar mixed scores
  tftar_gmm <- .gmm.scores(tftar_cor)
  #--- get a list
  .gmmlist(tftar_gmm, tftar_mi, regs)
}
.mi2gmm.ref <- function(object, idkey=NULL){
  #--- set input
  regs <- tni.get(object, what="regulatoryElements", idkey=idkey)
  #--- get tf-tar mi
  tftar_mi <- lapply(tni.get(object, what="refregulons.and.mode", idkey=idkey), 
                     abs)
  #--- get tf-tar correlation
  object@results$tn.ref <- tni.cor(object@gexp, object@results$tn.ref, 
                                   asInteger=FALSE,
                                   estimator=object@para$perm$estimator)
  tftar_cor <- tni.get(object, what="refregulons.and.mode", idkey=idkey)
  #--- get tf-tar mixed scores
  tftar_gmm <- .gmm.scores(tftar_cor)
  #--- get a list
  .gmmlist(tftar_gmm, tftar_mi, regs)
}
.gmmlist <- function(tftar_gmm, tftar_mi, regs){
  res <- list()
  for(rg in regs){
    res[[rg]]$gmm <- tftar_gmm[[rg]]
    res[[rg]]$mi <- tftar_mi[[rg]]
  }
  return(res)
}

#--- fit a normal mixture model that maximizes the conditional 
#--- expected log-likelihood at each M-step of the algorithm
.gmm.scores <- function(tftar_cor, lambda=c(0.3,0.3,0.3), 
                        mu=c(-0.5,0,0.5), sigma=c(0.15,0.2,0.15)){
  x <- unlist(tftar_cor, use.names=FALSE)
  emf <- mixtools::normalmixEM(x, lambda=lambda, mu=mu, sigma=sigma, 
                               mean.constr=c(NA, 0, NA), verb=FALSE)
  tftar_scores <- lapply(tftar_cor, function(x){
    df1 <- pnorm(x, emf$mu[1], emf$sigma[1], lower.tail = FALSE)
    df2 <- pnorm(x, emf$mu[2], emf$sigma[2], lower.tail = FALSE)
    df3 <- pnorm(x, emf$mu[2], emf$sigma[2], lower.tail = TRUE)
    df4 <- pnorm(x, emf$mu[3], emf$sigma[3], lower.tail = TRUE)
    r1 <- df1/(df1 + df3 + df4) * (x < 0)
    r2 <- df4/(df1 + df2 + df4) * (x >= 0)
    r2 - r1
  })
  return(tftar_scores)
}

#---------------------------------------------------------------
#--- helper function for tni.regulon.summary
.regulon.summary <- function(tf, refnet, tnet, regnames) {
    #-- Initialize
    regulonSummary <- list()
    
    #-- Add target information
    regulonSummary$targets <- as.matrix(data.frame(
        Total = c(sum(tnet[,tf] != 0), sum(refnet[,tf] != 0)),
        `Positive` = c(sum(tnet[,tf] > 0), 
                       `Positive targets` = sum(refnet[,tf] > 0)),
        `Negative` = c(sum(tnet[,tf] < 0),
                       `Negative targets` = sum(refnet[,tf] < 0)),
        row.names = c("DPInet", "Refnet")))
    #-- Add information about MI with other regulators
    tfMI <- refnet[tf,]
    other_regs <- tfMI[order(abs(tfMI), decreasing = TRUE)]
    other_regs <- names(other_regs[other_regs != 0])
    other_regs <- regnames[regnames%in%other_regs]
    regulonSummary$regulatorsMI <- other_regs
    return(regulonSummary)
}

##------------------------------------------------------------------------------
## target contribution for regulon activity
.target.contribution <- function(listOfRegulonsAndMode, regulonActivity, 
                                 phenoranks, phenotypes, exponent, 
                                 alternative, verbose){
  regs <- names(listOfRegulonsAndMode)
  nsamp <- ncol(phenoranks)
  if(isParallel() && length(regs)>1){
    if(verbose)
      cat("-Assessing target contribution (parallel version - ProgressBar disabled)...\n")
    if(verbose)
      cat("--For", length(listOfRegulonsAndMode), "regulon(s) and", 
          nsamp,'sample(s)...\n')
    cl<-getOption("cluster")
    snow::clusterExport(cl, list(".gsea2.target.contribution",
                                 ".run.tni.gsea2.alternative",
                                 ".fgseaScores4TNI"), 
                        envir=environment())
    tcontrib <- snow::parLapply(cl,regs, function(reg){
      regulonAndMode <- listOfRegulonsAndMode[[reg]]
      tc <- sapply(1:length(regulonAndMode), function(i){
        react <- .gsea2.target.contribution(regulonAndMode[-i], phenotypes, 
                                            phenoranks, exponent, alternative)
        sqrt(mean((react - regulonActivity$dif[,reg])^2))
      })
      names(tc) <- names(regulonAndMode)
      return(tc)
    })
  } else {
    if(verbose) 
      cat("-Assessing target contribution...\n")
    if(verbose)
      cat("--For", length(listOfRegulonsAndMode), "regulon(s) and", 
          nsamp,'sample(s)...\n')
    if(verbose) pb <- txtProgressBar(style=3)
    ntotal <- sum(unlist(lapply(listOfRegulonsAndMode, length)))
    ncount <- 0
    tcontrib <- lapply(regs, function(reg){
      regulonAndMode <- listOfRegulonsAndMode[[reg]]
      tc <- sapply(1:length(regulonAndMode), function(i){
        ncount <<- ncount + 1
        if(verbose) setTxtProgressBar(pb, ncount/ntotal)
        react <- .gsea2.target.contribution(regulonAndMode[-i], phenotypes, 
                                            phenoranks, exponent, alternative)
        sqrt(mean((react - regulonActivity$dif[,reg])^2))
      })
      names(tc) <- names(regulonAndMode)
      return(tc)
    })
    if(verbose) close(pb)
  }
  names(tcontrib) <- regs
  return(tcontrib)
}
## gsea2 for target contribution
.gsea2.target.contribution <- function(regulonAndMode, phenotypes, 
                                       phenoranks, exponent, alternative){
  samples <- colnames(phenotypes)
  regActivity <- NULL
  for(samp in samples){
    res <- .run.tni.gsea2.alternative(
      listOfRegulonsAndMode=list(regulonAndMode),
      phenotype=phenotypes[, samp],
      phenorank=phenoranks[, samp],
      exponent=exponent,
      alternative=alternative
    )
    regActivity <- c(regActivity,res$differential)
  }
  names(regActivity) <- samples
  return(regActivity)
}
## Count regulon's targets
.regulonCounts <- function(regulonsAndMode){
  counts <- lapply(regulonsAndMode, function(reg){
    c(length(reg),sum(reg>0),sum(reg<0))
  })
  counts <- data.frame(t(data.frame(counts, check.names = FALSE)), 
                       check.names = FALSE)
  colnames(counts) <- c("Size","Positive", "Negative")
  idx <- counts$Positive > 0.75  *counts$Size | counts$Positive < 0.25*counts$Size
  counts$TargetDistribution <- c("balanced","unbalanced")[1+idx]
  return(counts)
}

#---------------------------------------------------------------
#test Cohen's Kappa agreement between modulons and regulons
# pkappa<-function (Modulon,Regulon){
#   ttab<-data.frame(A=Modulon,B=Regulon)
#   ttab<-table(ttab)
#   #---get tabs
#   nc <- ncol(ttab)
#   ns <- sum(ttab)
#   weighttab <- matrix(0, nrow = nc, ncol = nc)
#   diag(weighttab)<-c(1,1)
#   #---kvalue
#   agreeP <- sum(ttab * weighttab)/ns
#   tm1 <- apply(ttab, 1, sum)
#   tm2 <- apply(ttab, 2, sum)
#   eij <- outer(tm1, tm2)/ns
#   chanceP <- sum(eij * weighttab)/ns
#   kvalue <- (agreeP - chanceP)/(1 - chanceP)
#   return(kvalue)
# }
#test overlap of among regulons in a tnet via phyper function
# tni.phyper<-function(tnet){
#   tnet[tnet!=0]=1
#   tnet<-tnet[rowSums(tnet)>0,]
#   labs<-colnames(tnet)
#   pmat<-matrix(NA,ncol=ncol(tnet),nrow=ncol(tnet), dimnames=list(labs, labs))
#   phtest<-function(c1,c2){
#     c<-c1+c2
#     pop<-length(c)
#     sample1<-sum(c1>0)
#     sample2<-sum(c2>0)
#     overlap<-sum(c==2)    
#     resph<-phyper(overlap, sample1, pop - sample1, sample2, lower.tail=FALSE)
#     resph
#   }
#   for(i in 1:ncol(tnet)){
#     v1<-tnet[,i]
#     for(j in 1:ncol(tnet)){
#       if(j>i){
#         v2<-tnet[,j]
#         ph<-phtest(as.integer(v1),as.integer(v2))
#         pmat[i,j]<-ph
#       }
#     }
#   }
#   lt<-lower.tri(pmat)
#   pmat[lt]=t(pmat)[lt]
#   diag(pmat)=1
#   pmat<-matrix(p.adjust(pmat),ncol=ncol(pmat),nrow=nrow(pmat),
#                dimnames=list(rownames(pmat),colnames(pmat)))           
#   pmat
# }

# ##----------------------------------------------------------------------------
# ##build an igraph object from hclust
# hclust2igraph<-function(hc,length.cutoff=NULL){
#   if(!is(hc,"hclust"))stop("'hc' should be an 'hclust' object!")
#   if(is.null(hc$labels))hc$labels=as.character(sort(hc$order))
#   #get treemap
#   tmap<-treemap(hc)
#   hcNodes<-tmap$hcNodes
#   hcNests<-hcNodes[hcNodes$type=="nest",]
#   hcEdges<-tmap$hcEdges
#   nestList<-tmap$nest
#   #remove nests based on length cutoff
#   if(!is.null(length.cutoff)){
#     #identify rmnodes
#     rmnodes<-hcEdges[hcEdges$edgeLength<length.cutoff,]
#     rmnodes<-rmnodes$childNode[rmnodes$childNode%in%hcNests$node]
#     #update hcEdges and nestList
#     hcEdges<-hcEdges.filter1(hcEdges,hcNodes,rmnodes)
#     nestList<-nestList[!names(nestList)%in%rmnodes]
#   }
#   #build igraph and return assigments
#   g<-graph.edgelist(as.matrix(hcEdges[,1:2]), directed=TRUE)
#   #E(g)$edgeWeight<-hcEdges$edgeLength
#   tp<-hcEdges$parentHeight-min(hcEdges$parentHeight)
#   E(g)$edgeWeight<-(max(tp)-tp)/max(tp) * 100
#   E(g)$arrowType<-0
#   list(g=g,nest=nestList)
# }
# hcEdges.filter1<-function(hcEdges,hcNodes,rmnodes){
#   #get the correct order
#   rmNests<-hcNodes[hcNodes$node%in%rmnodes,]
#   rmEdges<-hcEdges[hcEdges$childNode%in%rmnodes,]
#   #get filtered hcEdges
#   hcEdges<-hcEdges[!hcEdges$childNode%in%rmNests$node,]
#   #update parents'ids
#   for(oldid in rmNests$node){
#     newid<-rmEdges$parentNode[which(rmEdges$childNode==oldid)]
#     idx<-which(hcEdges$parentNode==oldid)
#     hcEdges$parentNode[idx]<-newid
#   }
#   rownames(hcEdges)<-NULL
#   hcEdges
# }

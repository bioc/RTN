

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##      OPTIMIZED GSEA FUNCTIONS FOR REGULONS AND TNETS
##------------------------------------------------------------------------
##------------------------------------------------------------------------

##------------------------------------------------------------------------
##This function computes observed and permutation-based scores from 
##a gene set enrichment analysis for a collection of regulons
gsea1tna <- function(listOfRegulons, phenotype, exponent=1, 
                     nPermutations=1000,verbose=TRUE) {	 
  #OBS1:names provided as integer values!
  #OBS2:max/min sizes checked in the previous functions!
  pheno.names <- as.integer(names(phenotype))
  nRegulons <- length(listOfRegulons)
  if(nRegulons > 0) {
    ##Generate a matrix to store the permutation-based scores, with 
    ##one row for each gene set (that has been tagged) and one column 
    ##for each permutation	
    permScores <- matrix(rep(0, (nPermutations * nRegulons)), nrow=nRegulons)
    rownames(permScores) <- names(listOfRegulons)
    ##Generate a vector to store the experimental scores
    ##one entry for each gene set (that has been tagged)
    observedScores <- rep(0, nRegulons)
    names(observedScores) <- names(listOfRegulons)
    ##Compute the scores	
    ##create permutation gene list
    perm.gL <- sapply(1:nPermutations, function(n) pheno.names[
      sample(1:length(phenotype), length(phenotype),replace=FALSE)])
    perm.gL<-cbind(pheno.names,perm.gL)
    ##check if package snow has been loaded and a cluster object 
    ##has been created for HTSanalyzeR
    if(isParallel() && nRegulons>1) {
      if(verbose)
        cat("-Performing GSEA (parallel version - ProgressBar disabled)...\n")
      if(verbose)
        cat("--For", length(listOfRegulons), "regulons...\n")  
      cl<-getOption("cluster")
      snow::clusterExport(cl, list("gseaScoresBatch4RTN"), envir=environment())
      scores <- snow::parSapply(cl, 1:length(listOfRegulons), function(i){
        gseaScoresBatch4RTN(phenotype=phenotype, geneNames.perm=perm.gL, 
                            geneset=listOfRegulons[[i]], exponent=exponent,
                            nPermutations=nPermutations)
      })
      for(i in 1:nRegulons){
        permScores[i,]<-unlist(scores["permScores",i])
        observedScores[i]<-unlist(scores["observedScores",i])
      }
    } else {
      if(verbose) cat("-Performing gene set enrichment analysis...\n")
      if(verbose) cat("--For", length(listOfRegulons), "regulons...\n")
      if(verbose) pb <- txtProgressBar(style=3)
      for(i in 1:nRegulons) {
        scores <- gseaScoresBatch4RTN(phenotype=phenotype, 
                                      geneNames.perm=perm.gL, 
                                      geneset=listOfRegulons[[i]],
                                      exponent=exponent,
                                      nPermutations=nPermutations)
        observedScores[i] <- scores$observedScores
        permScores[i,] <- scores$permScores
        if(verbose) setTxtProgressBar(pb, i/nRegulons)
      }	
      if(verbose) close(pb)
    }
  } else {
    observedScores <- NULL
    permScores <- NULL
  }
  return(list("Observed.scores" = observedScores , 
              "Permutation.scores" = permScores))	
}

##------------------------------------------------------------------------
##This function computes observed and permutation-based scores 
gsea2tna <- function(listOfRegulons, phenotype, exponent=1, 
                     nPermutations=1000, verbose1=TRUE, 
                     verbose2=TRUE) {   
  pheno.names <- as.integer(names(phenotype))
  nRegulons <- length(listOfRegulons)
  if(nRegulons > 0){
    permScores <- matrix(rep(0, (nPermutations * nRegulons)), nrow=nRegulons)
    rownames(permScores) <- names(listOfRegulons)
    observedScores <- rep(0, nRegulons)
    names(observedScores) <- names(listOfRegulons)
    perm.gL <- sapply(1:nPermutations, function(n) pheno.names[
      sample(1:length(phenotype), length(phenotype),replace=FALSE)])
    perm.gL<-cbind(pheno.names,perm.gL)
    if(isParallel() && nRegulons>1){
      if(verbose1 && verbose2)
        cat("-Performing two-tailed GSEA (parallel version - ProgressBar disabled)...\n")
      if(verbose1 && verbose2)
        cat("--For", length(listOfRegulons), "regulons...\n") 
      cl<-getOption("cluster")
      snow::clusterExport(cl, list("gseaScoresBatch4RTN"), envir=environment())
      scores <- snow::parSapply(cl, 1:length(listOfRegulons), function(i){
        gseaScoresBatch4RTN(phenotype=phenotype, geneNames.perm=perm.gL, 
                            geneset=listOfRegulons[[i]], exponent=exponent,
                            nPermutations=nPermutations)
      })
      for(i in 1:nRegulons){
        permScores[i,]<-unlist(scores["permScores",i])
        observedScores[i]<-unlist(scores["observedScores",i])
      }
    } else {
      if(verbose1 && verbose2) 
        cat("-Performing two-tailed GSEA analysis...\n")
      if(verbose1 && verbose2) 
        cat("--For", length(listOfRegulons), "regulons...\n")
      if(verbose1) pb <- txtProgressBar(style=3)
      for(i in 1:nRegulons) {
        scores <- gseaScoresBatch4RTN(phenotype=phenotype, 
                                      geneNames.perm=perm.gL, 
                                      geneset=listOfRegulons[[i]],
                                      exponent=exponent, 
                                      nPermutations=nPermutations)
        observedScores[i] <- scores$observedScores
        permScores[i,] <- scores$permScores
        if(verbose1) setTxtProgressBar(pb, i/nRegulons)
      }	
      if(verbose1) close(pb)
    }
  } else {
    observedScores <- NULL
    permScores <- NULL
  }
  return(list("Observed.scores" = observedScores , 
              "Permutation.scores" = permScores))	
}

##------------------------------------------------------------------------
##This function gets rid of the duplicates in a gene list.
tna.duplicate.remover <- function(phenotype, method = "max") {
  ##Get the unique names and create a vector that will store the 
  ##processed values corresponding to those names
  pheno.names <- names(phenotype)
  datanames <- unique(pheno.names)
  data.processed <- rep(0, length(datanames))
  names(data.processed) <- datanames
  data.processed2 <- data.processed
  l.data.processed <- length(data.processed)	
  ##If the absolute value of the min is bigger than the absolute 
  ##value of the max, then it is the min that is kept
  if(method=="max") {
    abmax<-function(x){
      imax<-max(x);imin<-min(x)
      ifelse(imax>abs(imin),imax,imin)
    }
    tp<-aggregate(phenotype,by=list(names(phenotype)),abmax)
    data.processed<-tp[,2]
    names(data.processed)<-tp[,1] 
  } else if(method=="min") {
    abmin<-function(x){
      imax<-max(x);imin<-min(x)
      ifelse(imax>abs(imin),imin,imax)
    }
    tp<-aggregate(phenotype,by=list(names(phenotype)),abmin)
    data.processed<-tp[,2]
    names(data.processed)<-tp[,1] 
  } else if(method=="average") {
    tp<-aggregate(phenotype,by=list(names(phenotype)),mean)
    data.processed<-tp[,2]
    names(data.processed)<-tp[,1]
  }	
  data.processed
}

##------------------------------------------------------------------------
## This function returns the nominal p-value for the GSEA results
tna.perm.pvalues <- function(permScores, observedScores, exact=TRUE){
	nPerm <- ncol(permScores)
	pvals <- rep(1, length(observedScores))
	##check how many permScores are higher (in magnitude) than
	##the observedScores; done separately for negative and positive
	##scores; pseudocounts are added to avoid P-values of zero.
	valid.id <- which(!is.na(observedScores) & observedScores!=0)
	for(i in valid.id) {
	  pvals[i] <- ifelse(
			observedScores[i] > 0,
			(sum(permScores[i, ] > observedScores[i])+1)/(nPerm+1),
			(sum(permScores[i, ] < observedScores[i])+1)/(nPerm+1)
		)
	}
	names(pvals) <- names(observedScores)
	if(exact){
	  p.resolution <- (1+1)/(nPerm+1)
	  mu <- apply(permScores, 1, mean)
	  std <- apply(permScores, 1, sd)
	  idx1 <- pvals <= p.resolution & std!=0
	  if(any(idx1)){
	    ep <- .exact_pvalue(mu[idx1], std[idx1], val=observedScores[idx1])
	    idx2 <- ep < pvals[idx1]
	    pvals[idx1][idx2] <- ep[idx2]
	  }
	}
	return(pvals)
}
## exact pvalues using permutational central limit theorem
.exact_pvalue <- function(mu, std, val){
  z <- (val - mu)/std
  p_vals <- pnorm(-abs(z), lower.tail=TRUE)
  return(p_vals)
}

#--------------------------------------------------------------------
##This function performs hypergeometric tests for over-representation 
##of hits, on a list of gene sets.
hyperGeoTest4RTN <- function(collectionOfGeneSets, universe, hits, 
                              minGeneSetSize = 15, pAdjustMethod = "BH", 
                              verbose = TRUE) {
  geneset.size <- lapply(lapply(collectionOfGeneSets, intersect, y=universe), 
                         length)
  geneset.size <- unlist(geneset.size)
  geneset.filtered <- which(geneset.size >= minGeneSetSize)
  l.GeneSet <- length(geneset.filtered)
  if(verbose) pb <- txtProgressBar(style=3)
  results <- t(
    sapply(geneset.filtered, 
           function(i) {
             if(verbose) setTxtProgressBar(pb, i/l.GeneSet)		
             .hyperGeoTest4RTN(collectionOfGeneSets[i], universe, hits)
           }
    )
  )
  if(verbose) close(pb)
  if(length(results) > 0) {
    adjPvals <- p.adjust(results[, "Pvalue"], method = pAdjustMethod)
    results <- cbind(results, adjPvals)
    colnames(results)[ncol(results)] <- "Adjusted.Pvalue"
    results <- results[order(results[, "Adjusted.Pvalue"]), , drop=FALSE]		
  } else {
    reuslts <- matrix(, nrow=0, ncol=7)
    colnames(results) <- c("Universe Size", "Gene Set Size", 
                           "Total Hits", "Expected Hits", "Observed Hits", 
                           "Pvalue", "Adjusted.Pvalue")
  }
  return(results)
}
.hyperGeoTest4RTN <- function(geneset, universe, hits) {
  #number of genes in universe
  N <- length(universe) 			
  #remove genes from gene set that are not in universe			
  geneset <- intersect(geneset[[1]], universe) 
  #size of gene set	
  m <- length(geneset) 							
  Nm <- N-m	
  #hits in gene set
  overlap <- intersect(geneset, hits) 	
  #number of hits in gene set		
  k <- length(overlap) 							
  n <- length(hits)	
  HGTresults <- phyper(k-1, m, Nm, n, lower.tail = FALSE)
  ex <- (n/N)*m
  if(m == 0) HGTresults <- NA
  hyp.vec <- c(N, m, n, ex, k, HGTresults)
  names(hyp.vec) <- c("Universe Size", "Gene Set Size", "Total Hits", 
                      "Expected Hits", "Observed Hits", "Pvalue")
  return(hyp.vec)
}

#--------------------------------------------------------------------
##This function computes observed and permutation enrichment scores 
gseaScoresBatch4RTN <- function(phenotype, geneNames.perm, geneset, 
                                exponent=1, nPermutations=1000){
  nh <- length(geneset)
  N <- length(phenotype)
  if(N>0 && nh>0){
    geneset <- as.numeric(geneset)
    ES <- sapply(1:(nPermutations+1), function(i){
      setidx <- geneNames.perm[geneset,i]
      .fgseaScores4TNA(phenotype,setidx,exponent)
    })
  } else {
    ES <- rep(0,nPermutations+1)
  }
  ES <- list(observedScores=ES[1], permScores=ES[2:(nPermutations+1)])
  return(ES)	
}
.fgseaScores4TNA <- function (phenotype, setidx, exponent=1){
  if(length(setidx)>0){
    setidx <- sort(setidx)
    nh <- length(setidx)
    N <- length(phenotype)
    hits <- abs(phenotype[setidx])^exponent
    NR <- sum(hits)
    hcumsum <- cumsum(hits)/NR
    topES <- hcumsum - (setidx - seq_along(setidx))/(N - nh)
    botES <- topES - hits/NR
    ESmax <- max(topES)
    ESmin <- min(botES)
    ES <- ifelse(abs(ESmin)>abs(ESmax), ESmin, ESmax)
  } else {
    ES <- 0
  }
  return(ES)
}

#--------------------------------------------------------------------
##This function returns enrichment scores, running scores, and 
##position of hits for a gene set; used only in the plotting functions
gseaScores4RTN <- function(phenotype, geneset, exponent=1, 
                           mode=c("score","runningscores")) {
  if( is.character(geneset) || is.null(names(geneset)) ){
    geneSetType<-rep(1,length(geneset))
    names(geneSetType)<-geneset
  } else {
    geneSetType<-geneset
    geneset<-names(geneset)
  }
  phenotype <- phenotype[!is.na(phenotype)]
  geneset<-intersect(names(phenotype), geneset)
  nh <- length(geneset)
  N <- length(phenotype)
  ES <- 0
  Phit <- rep(0, N)
  Pmiss <- rep(0, N)
  runningES <- rep(0, N)
  hits <- rep(FALSE, N)
  hitsType <- rep(0, N)
  hits[which(!is.na(match(names(phenotype), geneset)))] <- TRUE
  hitsType[which(!is.na(
    match(names(phenotype), names(geneSetType[geneSetType>0]))))] <- 1
  hitsType[which(!is.na(
    match(names(phenotype), names(geneSetType[geneSetType<0]))))] <- -1
  if(sum(hits)!=0){
    Phit[which(hits)]<-abs(phenotype[which(hits)])^exponent
    NR=sum(Phit)
    Pmiss[which(!hits)]<-1/(N-nh)
    Phit=cumsum(Phit/NR)
    Pmiss=cumsum(Pmiss)
    runningES<-Phit-Pmiss
    runningES[is.nan(runningES)] <- 0
    ESmax<-max(runningES)
    ESmin<-min(runningES)
    ES<-ifelse(abs(ESmin)>abs(ESmax), ESmin, ESmax)
  }
  if(mode[1]=="score"){
    return(ES)
  } else if(mode[1]=="runningscores"){
    #for a hit at extreme positions, set zero to avoid a misleading plot
    runningES[1] <- 0; runningES[length(runningES)] <- 0
    return(list("enrichmentScore"=ES, "runningScore"=runningES, 
                "positions"=hitsType))
  }
}

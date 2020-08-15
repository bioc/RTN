
#-- tni prune: 'pruned' selected regulons
setMethod(
  "tni.prune",
  "TNI", function(object, regulatoryElements = NULL, minRegCor = 0.95, 
                  tarPriorityMethod = "EC", minPrunedSize = 30, 
                  verbose = TRUE, ...){
    
    object <- upgradeTNI(object)
    
    if(object@status["Preprocess"]!="[x]")
      stop("NOTE: TNI object is not compleate: requires preprocessing!")
    if(object@status["Permutation"]!="[x]")
      stop("NOTE: TNI object is not compleate: requires permutation/bootstrap and DPI filter!")  
    if(object@status["DPI.filter"]!="[x]")
      stop("NOTE: TNI object is not compleate: requires DPI filter!")
    
    #-- checks
    tnai.checks(name = "minRegCor", para = minRegCor)
    tnai.checks(name = "tarPriorityMethod", para = tarPriorityMethod)
    tnai.checks(name = "minPrunedSize", para = minPrunedSize)
    tnai.checks(name = "verbose", para = verbose)
    
    #-- get regulatoryElements
    regnames <- tni.get(object, "regulatoryElements")
    if(!is.null(regulatoryElements)){
      tnai.checks(name="regulatoryElements", para=regulatoryElements)
      if(sum(regulatoryElements%in%regnames) > 
         sum(regulatoryElements%in%names(regnames))){
        regulatoryElements <- regnames[regnames%in%regulatoryElements]
      } else {
        regulatoryElements <- regnames[names(regnames)%in%regulatoryElements]
      }
      if(length(regulatoryElements)==0)
        stop("NOTE: 'regulatoryElements' argument has no valid names!")
    } else {
      regulatoryElements <- regnames
    }
    regnames <- regnames[names(regulatoryElements)]
    names(regulatoryElements) <- regulatoryElements
    
    #-- creating effect list
    if (tarPriorityMethod == "TC"){
      object <- tni.gsea2(object, regulatoryElements = regulatoryElements, 
                          targetContribution=TRUE, verbose=verbose, ...=...)
      dataActivity <- tni.get(object, what = "regulonActivity")
      listOfTargetEffect <- dataActivity$data$listOfTargetContribution
    } else if (tarPriorityMethod == "EC"){
      object <- tni.gsea2(object, regulatoryElements = regulatoryElements, 
                          additionalData=TRUE, verbose=verbose, ...=...)
      dataActivity <- tni.get(object, what = "regulonActivity")
      listOfTargetEffect <- lapply(regulatoryElements, expressionCorrelation, 
                                   object, dataActivity)
    } else if (tarPriorityMethod == "MI"){
      object <- tni.gsea2(object, regulatoryElements = regulatoryElements, 
                          additionalData=TRUE, verbose=verbose, ...=...)
      dataActivity <- tni.get(object, what = "regulonActivity")
      listOfTargetEffect <- lapply(regulatoryElements, mutualInformation,
                             object, dataActivity)
    }
    regulatoryElements <- regulatoryElements[names(listOfTargetEffect)]
    
    #-- pruning algorithm
    nsamp <- ncol(dataActivity$data$phenoranks)
    if(length(regulatoryElements)>0){
      if (isParallel() && length(regulatoryElements) > 1){
        if(verbose)
          cat("-Running pruning algorithm (parallel version - ProgressBar disabled)...\n")
        if(verbose)
          cat("--For", length(regulatoryElements), "regulon(s) and", 
              nsamp,'sample(s)...\n')
        cl<-getOption("cluster")
        snow::clusterExport(cl, list("minRegulonComputation",
                                     ".gsea2.target.contribution",
                               ".run.tni.gsea2.alternative",
                               ".fgseaScores4TNI","checkTars"), 
                            envir=environment())
        prunedTarlist <- snow::parLapply(cl, regulatoryElements, 
                                         minRegulonComputation, 
                                   dataActivity, listOfTargetEffect, 
                                   minPrunedSize, 
                                   minRegCor, verbose=FALSE, ...=...)
      } else {
        if(verbose)cat("-Running pruning algorithm...\n")
        if(verbose)cat("--For", length(regulatoryElements), "regulon(s) and", 
                       nsamp,'sample(s)...\n')
        prunedTarlist <- lapply(regulatoryElements, minRegulonComputation, 
                                dataActivity, listOfTargetEffect, minPrunedSize, 
                                minRegCor, verbose, ...=...)
      }
      names(prunedTarlist) <- regulatoryElements
      
      #-- prune the tnet
      prunedTnet <- sapply(regulatoryElements, function(reg){
        tnet <- tni.get(object, "tnet")
        notTars <- setdiff(rownames(tnet), prunedTarlist[[reg]])
        tnet[notTars,reg] <- 0
        return(tnet[,reg])
      })
      tnet <- object@results$tn.dpi
      tnet[,regulatoryElements] <- prunedTnet
      object@results$tn.dpi <- tnet
      
      #-- add processing info
      object@para$pruning <- list(minRegCor = minRegCor, 
                                  tarPriorityMethod = tarPriorityMethod,
                                  minPrunedSize = minPrunedSize,
                                  prunedRegs = regulatoryElements)
    }
    
    return(object)
    
  })

##------------------------------------------------------------------------------
#-- Auxiliary functions for pruning
##------------------------------------------------------------------------------

#-- pruning algorithm
minRegulonComputation <- function (reg, dataActivity, listOfTargetEffect,
                                   minPrunedSize, minRegCor, verbose){
  #-- get regulon and target contribution metric
  regulonAndMode <- dataActivity$data$listOfRegulonsAndMode[[reg]]
  targetEffect <- sort(listOfTargetEffect[[reg]], decreasing = TRUE)
  #-- start pruning
  e_cor <- 0
  best_ntar <- length(targetEffect)
  new_best_ntar <- trunc(best_ntar/2)
  mintail <- min(checkTars(regulonAndMode))
  while(!(new_best_ntar %in% best_ntar)){
    ntar <- new_best_ntar
    current_tars <- names(targetEffect[1:ntar])
    mintail <- min(checkTars(regulonAndMode[current_tars]))
    #-- compute new activity
    react <- .gsea2.target.contribution(regulonAndMode[current_tars], 
                                        dataActivity$data$phenotypes, 
                                        dataActivity$data$phenoranks, 
                                        dataActivity$data$exponent, 
                                        dataActivity$data$alternative)
    #-- compare by correlation
    e_cor <- cor(react, dataActivity$dif[,reg], method = "spearman")
    #-- check minRegCor and minPrunedSize
    if(e_cor > minRegCor && mintail == minPrunedSize){
      break
    } else if (e_cor > minRegCor && mintail > minPrunedSize){
      best_ntar <- c(best_ntar, ntar)
      new_best_ntar <- trunc(ntar/2)
    } else if (e_cor < minRegCor || mintail < minPrunedSize){
      best_ntar <- c(best_ntar, ntar)
      new_best_ntar <- ntar + trunc(ntar/2)
    }
    if(verbose) cat(".")
  }
  if(verbose) cat("\n")
  current_tars <- dataActivity$data$rnames_phenotypes[as.numeric(current_tars)]
  return(current_tars)
}

#-- old pruning algorithm
# minRegulonComputation_Old <- function (reg, dataActivity, listOfTargetEffect,
#                                    minPrunedSize, minRegCor, verbose){
#   #-- get regulon and target contribution metric
#   regulonAndMode <- dataActivity$data$listOfRegulonsAndMode[[reg]]
#   targetEffect <- sort(listOfTargetEffect[[reg]], decreasing = TRUE)
#   #-- start pruning
#   e_cor <- 0
#   best_ntar <- length(targetEffect)
#   new_best_ntar <- trunc(best_ntar/2)
#   while(!(new_best_ntar %in% best_ntar)){
#     ntar <- new_best_ntar
#     current_tars <- names(targetEffect[1:ntar])
#     #-- compute new activity
#     react <- .gsea2.target.contribution(regulonAndMode[current_tars], 
#                                         dataActivity$data$phenotypes, 
#                                         dataActivity$data$phenoranks, 
#                                         dataActivity$data$exponent, 
#                                         dataActivity$data$alternative)
#     #-- compare by correlation
#     e_cor <- cor(react, dataActivity$dif[,reg], method = "spearman")
#     #-- check minRegCor and minPrunedSize
#     if (e_cor > minRegCor && ntar > minPrunedSize){
#       best_ntar <- c(best_ntar, ntar)
#       new_best_ntar <- trunc(ntar/2)
#     } else if (e_cor < minRegCor){
#       best_ntar <- c(best_ntar, ntar)
#       new_best_ntar <- ntar + trunc(ntar/2)
#     } else if(e_cor > minRegCor && ntar <= minPrunedSize) {
#       current_tars <- names(targetEffect[1:minPrunedSize])
#       break
#     }
#     if(verbose) cat(".")
#   }
#   if(verbose) cat("\n")
#   current_tars <- dataActivity$data$rnames_phenotypes[as.numeric(current_tars)]
#   return(current_tars)
# }

#-- check size of regulon tails
checkTars <- function(regulonAndMode){
  c(pos=sum(regulonAndMode>0),neg=sum(regulonAndMode<0))
}

#-- computation of target effect by correlation expression/activity
expressionCorrelation <- function(reg, tni, dataActivity) {
  regActivity <- dataActivity$dif[,reg]
  tar_ids <- names(dataActivity$data$listOfRegulonsAndMode[[reg]])
  tar_names <- dataActivity$data$listOfRegulons[[reg]]
  gexp <- t(tni@gexp[tar_names,])
  tars_contrib <- abs(as.vector(cor(regActivity, gexp)))
  names(tars_contrib) <- tar_ids
  return(tars_contrib)
}

#-- computation of target effect by mutual information
mutualInformation <- function(reg, tni, dataActivity){
  tar_ids <- names(dataActivity$data$listOfRegulonsAndMode[[reg]])
  tar_names <- dataActivity$data$listOfRegulons[[reg]]
  tars_contrib <- abs(tni@results$tn.dpi[tar_names,reg])
  names(tars_contrib) <- tar_ids
  return(tars_contrib)
}



#---------------------------------------------------------------------
#Subgroup Regulon Enrichment for TNI-class objects
# rtni_tmpLIHCCHOL <- tni.sre(rtni_tmpLIHCCHOL, sampleGroups, pValueCutoff=0.0001)

setMethod(
  "tni.sre",
  "TNI",
  function(object, sampleGroups, regulatoryElements = NULL, pValueCutoff = 0.05, 
           pAdjustMethod = "BH") {
    
    #--- check compatibility
    object <- upgradeTNI(object)
    
    #-- Checks
    if(object@status["Preprocess"]!="[x]")
      stop("TNI object requires preprocessing!")
    if(object@status["Permutation"]!="[x]")
      stop("TNI object requires permutation/bootstrap and DPI filter!")  
    if(object@status["DPI.filter"]!="[x]")
      stop("TNI object requires DPI filter!")
    if(object@status["Activity"]!="[x]")
      stop("TNI object requires regulon activity results! 
           Please see 'tni.gsea2' function.")
    tnai.checks("sampleGroups", sampleGroups)
    tnai.checks("regulatoryElements",regulatoryElements)
    tnai.checks("pValueCutoff",pValueCutoff)
    tnai.checks("pAdjustMethod",pAdjustMethod)
    
    #-- Get data
    colAnnotation <- tni.get(object, "colAnnotation")
    
    #-- check sampleGroups
    if(is.list(sampleGroups)){
      if(!all(unlist(sampleGroups)%in%rownames(colAnnotation)))
        stop("All IDs in 'sampleGroups' should be available the 'colAnnotation' slot of the TNI object.")
      gnames <- names(sampleGroups)
      if(is.null(gnames) || any(duplicated(gnames)))
         stop("'sampleGroups' should be named with unique names.")
    } else {
      if(!sampleGroups%in%colnames(colAnnotation)){
        stop("Varible 'sampleGroups' should be listed in the 'colAnnotation' slot of the TNI object.")
      }
      sampleGroups <- colAnnotation[[sampleGroups]]
      if(!any(duplicated(sampleGroups))) {
        stop("'sampleGroups' column doesn't contain information to divide samples into subgroups")
      }
      sampleGroups <- split(rownames(colAnnotation),sampleGroups)
    }
    
    #-- check regulatoryElements
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
    regulonActivity <- tni.get(object, what = "regulonActivity")
    regulonActivity <- regulonActivity$dif
    
    #-- Create results table
    res_list <- lapply(names(sampleGroups), function(gp) {
      group_res <- rbindlist(
        lapply(colnames(regulonActivity), .regulonGroupFET, regulonActivity, 
               gp, sampleGroups, pValueCutoff))
    })
    res_tb <- rbindlist(res_list)
    res_tb$FET_pAdjusted <- p.adjust(res_tb$FET_pValue, method = pAdjustMethod)
    
    #-- Use pAdjusted for enrichment rest:
    res_tb[res_tb$FET_pAdjusted > pValueCutoff, "Enrichment_mode"] <- 0
    
    #-- Reorder
    res_tb <- res_tb[order(res_tb$FET_pAdjusted),]
    
    #-- Add to tns
    object@results$subgroupEnrichment <- res_tb
    return(object)
  })

#---------------------------------------------------------------------
tni.plot.sre <- function(object, nGroupsEnriched = NULL, nTopEnriched = NULL,
                         colors = c("blue","white","red"), 
                         breaks = seq(-1.5, 1.5, 0.1), 
                         markEnriched = TRUE, ...) {
  
  #-- Checks
  if(!is(object,"TNI"))stop("not a 'TNI' object!")
  tnai.checks("nGroupsEnriched",nGroupsEnriched)
  tnai.checks("nTopEnriched",nTopEnriched)
  tnai.checks("colorvec",colors)
  tnai.checks("breaks",breaks)
  tnai.checks("markEnriched",markEnriched)
  if(!is.null(nTopEnriched) && !is.null(nGroupsEnriched))
    stop("it must be use either the 'nGroupsEnriched' or 'nTopEnriched' argument.")
  
  #-- subgroupEnrichment checks
  fet_res <- tni.get(object, "subgroupEnrichment")
  if(is.null(fet_res)){
    stop("'tni.plot.sre' requires results from the 'tni.sre' analysis!")
  }
  fet_ls <- split(fet_res, fet_res$Regulon)
  breaks <- sort(breaks)
  
  #-- Filter checks
  if(!is.null(nGroupsEnriched)){
    #-- Filter to keep regulons enriched in nGroupsEnriched
    idx <- sapply(fet_ls, function(regtb){ 
      sum(regtb$Enrichment_mode != 0) >= nGroupsEnriched
    })
    fet_ls <- fet_ls[idx]
  } else if(!is.null(nTopEnriched)){
    #-- Filter to keep nTopEnriched regulons in each group
    fet_Gls <- split(fet_res, fet_res$Group)
    idx <- sapply(fet_Gls, function(grouptb) {
      unlist((grouptb[order(grouptb$FET_pValue), "Regulon"][1:nTopEnriched]))
    })
    idx <- unique(as.vector(idx))
    fet_ls <- fet_ls[idx]
  }
  filt_fet_res <- rbindlist(fet_ls)
  if(nrow(filt_fet_res)==0){
    stop(paste0("no regulon passed the '",by,"' threshold"))
  }
  
  #-- Get dESaverage matrix
  dESaverage <- as.data.frame(
    dcast(filt_fet_res, Regulon ~ Group, value.var = "dESaverage"))
  rownames(dESaverage) <- dESaverage$Regulon
  dESaverage <- dESaverage[,-1]
  colors <- colorRampPalette(colors)(length(breaks) - 1)
  
  #-- Get Enrichment_mode matrix
  Enrichment_mode <- as.data.frame(
    dcast(filt_fet_res, Regulon ~ Group, value.var = "Enrichment_mode")
  )
  rownames(Enrichment_mode) <- Enrichment_mode$Regulon
  Enrichment_mode <- Enrichment_mode[,-1]
  
  if(markEnriched){
    
    #-- Spread "*"s
    idx <- Enrichment_mode
    Enrichment_mode[idx==0] <- ""
    Enrichment_mode[idx!=0] <- "*"
    
    #-- Plot
    pheatmap(dESaverage,
             color = colors,
             breaks = breaks,
             border_color = NA, 
             display_numbers = Enrichment_mode,
             ...=...)
  } else {
    pheatmap(dESaverage,
             color = colors,
             breaks = breaks,
             border_color = NA,
             ...=...)
  }
  summary <- list(dESaverage=dESaverage, 
                  Enrichment_mode=Enrichment_mode)
  invisible(summary)
}

#---------------------------------------------------------------------
.regulonGroupFET <- function(reg, regact, group, grouping, pValueCutoff) {
  
    #-- Getting group position in grouping list
    grn <- which(names(grouping) == group)
    
    #-- Initializing result
    res <- list(
        Regulon = reg,
        Group = group,
        Enrichment_mode = 0,
        dESaverage = 0,
        FET_pValue = NA
        )
    
    #-- For one regulon
    reg_dES <- regact[,reg]
    idx <- reg_dES >= 0
    act <- names(reg_dES[idx])
    rep <- names(reg_dES[!idx])
    
    #-- ct table
    ct <- matrix(c(
        length(intersect(act, grouping[[grn]])),
        length(intersect(rep, grouping[[grn]])),
        length(intersect(act, unlist(grouping[-grn]))),
        length(intersect(rep, unlist(grouping[-grn])))
    ), nrow = 2, dimnames = list(Regulon = c("Active", "Repressed"),
                                 Group = c("Belong", "Don't belong")))
    
    #-- Positively enriched test
    ft_act <- fisher.test(ct, alternative = "greater")
    
    #-- Negatively enriched test
    ft_rep <- fisher.test(ct[2:1,], alternative = "greater")
    
    #-- Add results
    ft <- list(ft_act, ft_rep)
    idx <- which.min(c(ft_act$p.value, ft_rep$p.value))
    res$FET_pValue <- ft[[idx]]$p.value
    res$Enrichment_mode <- ifelse(idx == 1, 1, -1)
    res$Enrichment_mode <- ifelse(res$FET_pValue > pValueCutoff, 0, 
                                  res$Enrichment_mode)
    res$dESaverage <- mean(reg_dES[grouping[[grn]]])
    return(res)
}

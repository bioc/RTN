# Unit tests fot TNI-class methods
test_tni <- function(){
  data(tniData)
  tfs4test<-c("PTTG1","FOXM1")
  rtni <- tni.constructor(expData=tniData$expData, 
                          regulatoryElements=tfs4test, 
                          rowAnnotation=tniData$rowAnnotation)
  #tni.permutation
  rtni<-tni.permutation(rtni,nPermutations=10)
  res<-tni.get(rtni,what="refnet")
  checkTrue(is.matrix(res) && ncol(res)==2)
  #tni.bootstrap
  rtni<-tni.bootstrap(rtni,nBootstraps=10)
  res<-tni.get(rtni,what="refnet")
  checkTrue(is.matrix(res) && ncol(res)==2)
  #tni.dpi.filter
  rtni<-tni.dpi.filter(rtni)
  res<-tni.get(rtni,what="tnet")
  checkTrue(is.matrix(res) && ncol(res)==2)
  #tni.conditional
  annot <- res<-tni.get(rtni,what="rowAnnotation")
  idx <- annot$SYMBOL%in%c("FGF2","ERBB2")
  mod4test<-annot$PROBEID[idx]
  rtni<-tni.conditional(rtni, modulators=mod4test, minRegulonSize=1)
  res<-tni.get(rtni,what="cdt.list")
  checkTrue(is.list(res))
  #tni.graph
  #res<-tni.graph(rtni, gtype="rmap", tfs=tfs4test)
  #checkTrue(is.igraph(res))
}
# Unit tests fot TNA-class methods
test_tna <- function(){
  data(tniData)
  data(tnaData)
  tfs4test<-c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")
  rtni <- tni.constructor(expData=tniData$expData, 
                          regulatoryElements=tfs4test, 
                          rowAnnotation=tniData$rowAnnotation)
  rtni<-tni.permutation(rtni,nPermutations=10)
  rtni<-tni.bootstrap(rtni,nBootstraps=10)
  rtni<-tni.dpi.filter(rtni)
  rtna<-tni2tna.preprocess(rtni, phenotype=tnaData$phenotype, hits=tnaData$hits, 
                           phenoIDs=tnaData$phenoIDs)
  #tna.mra
  rtna <- tna.mra(rtna, minRegulonSize=1)
  res<-tna.get(rtna,what="mra")
  checkTrue(is.data.frame(res) && ncol(res)==8)
  #tna.gsea1
  rtna <- tna.gsea1(rtna, nPermutations=10, minRegulonSize=1)
  res <- tna.get(rtna,what="gsea1")
  checkTrue(is.data.frame(res) && ncol(res)==5)
  #tna.gsea2
  rtna <- tna.gsea2(rtna, nPermutations=10, minRegulonSize=1)
  res<-tna.get(rtna,what="gsea2")
  checkTrue(is.list(res) && length(res)==3)
}

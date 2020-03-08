
##------------------------------------------------------------------------------
## Check target distribution
tni.plot.checks <- function(object, minRegulonSize = 15, 
                            option = c("barplot","edf","points")){
  object <- upgradeTNI(object)
  status <- tni.get(object, what="status")
  if(status["DPI.filter"]!="[x]") 
    stop("NOTE: input 'object' needs dpi analysis!")
  rgcounts <- tni.get(object, what="regulonSize")
  option <- match.arg(option)
  if(option=="barplot"){
    .barplot.pos.vs.neg(rgcounts, minRegulonSize)
  } else if(option=="edf"){
    .edf.pos.vs.neg(rgcounts)
  } else if(option=="points"){
    .points.pos.vs.neg(rgcounts)
  }
}

##------------------------------------------------------------------------------
.barplot.pos.vs.neg <- function(rgcounts, minRegulonSize){
  rgcounts$pve_minus_nve <- rgcounts$Positive-rgcounts$Negative
  rgcounts <- rgcounts[sort.list(rgcounts$pve_minus_nve),]
  op <- par(mar=c(5, 5, 4, 4), mgp=c(2.2,0.4,0), tcl=-0.2)
  cols <- adjustcolor( c("#FF8E91","#96D1FF","grey50"), alpha.f = 0.5)
  legcols <- adjustcolor( c("#FF8E91","#96D1FF","grey50"), alpha.f = 1)
  plot.new()
  n <- nrow(rgcounts)
  plot.window(xlim=c(1,n*1.2), 
              ylim=range(c(-rgcounts$Negative,rgcounts$Positive)))
  box(col="grey"); axis(2, las=2, tcl=-0.2)
  barplot(-rgcounts$Negative, col = cols[2], border=cols[2], space=0.2,
          ylab="", xlab="", las=2, axes=FALSE, add = TRUE)
  barplot(rgcounts$Positive, col = cols[1], border=cols[1], space=0.2,
          ylab="", xlab="", las=2, axes=FALSE, add = TRUE)
  segments(y0 = c(minRegulonSize,-minRegulonSize), x0 = c(0,0), cex=0.8,
           x1= c(n,n,n,n)*1.2, col = 'grey50', lty=2)
  points(y = rgcounts$pve_minus_nve, x=0:(n-1)*1.2, pch=16, cex=0.25, 
         col=legcols[3])
  title(ylab="Number of targets", mgp=c(2.5,0.4,0))
  title(xlab=paste0("n = ", length(rgcounts$Size), " regulons"), 
        mgp=c(0.5,0.4,0))
  size_th1 <- sum(rgcounts$Size>=minRegulonSize)
  size_th2 <- sum(rgcounts$Positive>=minRegulonSize | 
                    rgcounts$Negative>=minRegulonSize)
  size_th3 <- sum(rgcounts$Positive>=minRegulonSize & 
                    rgcounts$Negative>=minRegulonSize)
  size_th1 <- paste0("- ", size_th1, " regulons have >= ", 
                     minRegulonSize," total targets")
  size_th2 <- paste0("- ", size_th2, " regulons have >= ", 
                     minRegulonSize," '+' or '-' targets")
  size_th3 <- paste0("- ", size_th3, " regulons have >= ", 
                     minRegulonSize," '+' and '-' targets")
  legend("topleft", x.intersp = 0,
         legend = c("Summary:",size_th1,size_th2,size_th3), 
         cex=0.7, bty = "n", pch=NA)
  legend("bottomright", 
         legend = c("Positive targets","Negative targets", 
                    "Positive - Negative"), 
         inset = c(0.05,0), cex=0.7, col = legcols, 
         bty = "n", pch=c(15,15,16))
  par(op)
}

##------------------------------------------------------------------------------
.edf.pos.vs.neg <- function(rgcounts){
  op <- par(mar=c(5, 5, 4, 4), mgp=c(2.2,0.4,0), tcl=-0.2)
  legcols <- adjustcolor(c("#FF8E91","#96D1FF","grey50"), alpha.f = 0.8)
  cols <- adjustcolor(legcols, alpha.f = 0.8)
  plot.new()
  plot.window(xlim=c(1,max(rgcounts$Size)), ylim=c(0,1))
  segments( x0 = c(15,30,45,60), y0 = c(0,0,0,0), y1= c(1,1,1,1), 
            col = 'grey70', lty=2)
  plot(ecdf(rgcounts$Positive), verticals=FALSE, pch=20, cex=0.8, col.hor=NA, 
       xlim=c(0,max(rgcounts$Size)), 
       col=cols[1], ylim=c(0,1), main="", xlab="", ylab="", las=1, 
       col.01line = 'grey70',
       panel.first=grid(lty=1, lwd=0.25, col="lightgray"), add=TRUE)
  plot(ecdf(rgcounts$Negative), verticals=FALSE, pch=20, cex=0.8, 
       col.hor=NA, xlim=c(0,max(rgcounts$Size)),
       col=cols[2], ylim=c(0,1), main="", xlab="", ylab="", las=1, 
       add=TRUE, col.01line = NA)
  plot(ecdf(rgcounts$Size), verticals=FALSE, pch=20, cex=0.8, col.hor=NA, 
       xlim=c(0,max(rgcounts$Size)), 
       col=cols[3], ylim=c(0,1), main="", xlab="", ylab="", las=1, 
       add=TRUE, col.01line = NA)
  text(x = c(15,30,45,60), y=c(0.02,0.06,0.1,0.14), 
       labels = c(15,30,45,60), cex=0.8)
  box(col="grey"); axis(2, las=2, tcl=-0.2); axis(1); axis(2, las=2)
  title(ylab = paste0("Fraction of ", nrow(rgcounts), " regulons"),
        mgp=c(2.5,0.4,0))
  title(xlab = "Number of regulon targets", mgp=c(1.8,0.4,0))
  legend("bottomright", legend = c("Positive targets",
                                   "Negative targets", "All targets"), 
         col = legcols, pch=19, bty="n", cex=0.7, inset = c(0.05,0.05))
  par(op)
}

##------------------------------------------------------------------------------
.points.pos.vs.neg <- function(rgcounts, minRegulonSize){
  rgcounts$pve_minus_nve <- rgcounts$Positive-rgcounts$Negative
  rgcounts <- rgcounts[sort.list(rgcounts$Size),]
  op <- par(mar=c(5, 5, 4, 4), mgp=c(2.2,0.5,0), tcl=-0.2)
  legcols <- adjustcolor(c("red","royalblue","grey20"), alpha.f = 0.6)
  cols <- adjustcolor(legcols, alpha.f = 0.5)
  plot.new()
  plot.window(xlim=c(1,length(rgcounts$Size)), ylim=c(0,max(rgcounts$Size)))
  x <- 1:nrow(rgcounts)
  points(x=x, y=rgcounts$Positive, pch=20, cex=0.8, 
         col=cols[1], main="", xlab="", 
         ylab="", las=1, panel.first=grid(lty=1, lwd=0.25, col="lightgray"))
  points(x=x, y=rgcounts$Negative, pch=20, cex=0.8, col=cols[2], 
         main="", xlab="", 
         ylab="", las=1, panel.first=grid(lty=1, lwd=0.25, col="lightgray"))
  box(col="grey"); axis(2, las=2, tcl=-0.2); axis(1); axis(2, las=2)
  title(ylab = "Number of regulon targets", mgp=c(3,0.4,0))
  title(mgp=c(2,0.4,0), xlab=paste0("Regulons sorted by total targets ", 
                                    "(n = ",length(rgcounts$Size),
                                    " regulons)"),) 
  legend("topleft", legend = c("(+) targets","(-) targets"), 
         col = legcols, pch=19, bty="n", cex=0.8)
  par(op)
}


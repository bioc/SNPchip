setMethod("chromosomeAnnotation", "AnnotatedSnpCallSet", function(object) object@chromosomeAnnotation)
setReplaceMethod("chromosomeAnnotation", c("AnnotatedSnpCallSet", "data.frame"), function(object, value){
  object@chromosomeAnnotation <- value
  object
})

setMethod("initialize", "AnnotatedSnpCallSet",
          function(.Object,
                   assayData = assayDataNew(
                     calls=calls,
                     callsConfidence = callsConfidence),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=character(),
                   calls=new("matrix"),
                   callsConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                 chromosomeAnnotation = data.frame()){
            .Object@assayData <- assayDataNew(calls = calls, callsConfidence = callsConfidence)
            .Object@phenoData <- phenoData
            .Object@annotation <- annotation
            .Object@featureData <- featureData
            .Object@experimentData <- experimentData
            .Object@chromosomeAnnotation <- chromosomeAnnotation
            .Object
          })

setMethod("show", "AnnotatedSnpCallSet", function(object) {
  cat(class( object ), " (storageMode: ", storageMode(object), ")\n", sep="")
  adim <- dim(object)
  if (length(adim)>1)
    cat("assayData:",
        if (length(adim)>1)
        paste(adim[[1]], "features,",
              adim[[2]], "samples") else NULL,
        "\n")
  cat("  element names:",
      paste(assayDataElementNames(object), collapse=", "), "\n")
  cat("experimentData: use 'experimentData(object)'\n")
  pmids <- pubMedIds(object)
  if (length(pmids) > 0 && all(pmids != ""))
    cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
  cat("Annotation:", annotation(object), "\n")
  cat("phenoData\n")
  if(length(adim) > 1) show(phenoData(object))
  cat("featureData\n")      
  if(length(adim) > 1) show(featureData(object))
  cat("Annotation ")
  show(annotation(object))
  cat("\nchromosomeAnnotation\n")
  adim <- nrow(chromosomeAnnotation(object))
  if(adim > 1){
    idx <- selectSomeIndex(chromosomeAnnotation(object), maxToShow=4)
    pData <- chromosomeAnnotation(object)[c(idx[[1]], idx[[3]]), , drop=FALSE]
    rnms <- rownames(pData)
    nms <- c(rnms[idx[[1]]], idx[[2]],
             if (!is.null(idx[[1]])) rnms[-idx[[1]]] else NULL)
    pData <- pData[nms,]
    pData[nms=="...", ] <- rep("...", ncol(pData))
    rownames(pData)[nms=="..."] <- "..."
    print(pData)
  }
})

##setMethod("plotSnp", "AnnotatedSnpCallSet",
##          function(object, chromosomes, samples, ylim=NULL, xlim=NULL, col="black",
##                   colAA="blue", colAB="red", colNC="green3", colCentromere="bisque",
##                   col.axis="brown", pch=".", cex=1, cexAA=1, cexAB=1, cexNC=1,
##                   cex.axis=0.8, cex.legend=1, cex.chr=0.8, oma=c(5, 3, 4, 0.5),
##                   mar=c(0, 0, 0, 0.2), width.right=NULL, nsummaryPanel=FALSE,
##                   showLayout=TRUE, plotIt=TRUE, digits=3, legend=TRUE, legend.stats="left",
##                   legend.pch="topleft", legend.bty="o", legend.col="white",
##                   alternate.xaxis=TRUE, xaxis=TRUE, jitter=TRUE, factor=0.1, yTicks=5,
##                   xTicks=2, bty="n", bw=FALSE, ...){
##            ##If black and white, change the default color scheme
##            if(bw){
##              colAA <- gray(0.7)
##              colAB <- gray(0)
##              colNC <- gray(0.9)
##              col.axis <- gray(0)
##              colCentromere <- gray(0.7)
##            }
##            chromosomes <- paste("chr", chromosomes, sep="")
##            chromosomes[chromosomes == "chr23"] <- "chrX"
##            chromosomes[chromosomes == "chr24"] <- "chrY"
##            object <- object[chromosome(object) %in% chromosomes, samples]
##            obj <- object
##            
##            samplenames <- sampleNames(object)
##            chrAnn <- chromosomeAnnotation(object)[chromosomes, ]
##            
##            object <- as(object, "AnnotatedSnpCallSetList")            
##            object <- SNPchip::snpSetList(object)
##            names(object) <- chromosomes
##            k <- 1
##            S <- dim(object[[1]])[2]
##            N <- length(object)
##            
##            widths <- chrAnn$chromosomeSize
##            widths <- widths/min(widths)
##            if(summaryPanel){
##              if(is.null(width.right)){
##                width.right <- length(chromosomes)/1.5
##              }
##              N <- N+1
##              widths <- c(widths, width.right)
##            }
##            nf <- layout(matrix(1:(S*N),
##                                nc=N,
##                                byrow=FALSE), widths=widths)
##
##            ##Option to return just the layout
##            if(!plotIt) {
##              layout.show(nf)
##              return()
##            }
##            par(mar=mar, oma=oma)
##            for(chrom in chromosomes){
##              for(i in 1:S){
##                if(k == 1){ yaxis=TRUE; side.last <- 3} else yaxis=FALSE
##                chromosomeSize <- chrAnn[chrom, 3]                
##                xlim <- c(0, chromosomeSize)
##                plotChromosome(object[[chrom]][, i], colAA=colAA, colAB=colAB,
##                               colNC=colNC, cexAA=cexAA, cexAB=cexAB, cexNC=cexNC,
##                               ylim=ylim, xlim=xlim, yaxt="n", xaxt="n",  yaxs="r",
##                               legendStats=FALSE, legendPch=FALSE,
##                               centromereBorder=NA, xlab="", ylab="", bty=bty,
##                               panel.xaxis=FALSE, panel.yaxis=yaxis, yTicks=yTicks,
##                               mar=mar)
##                if(i == 1){
##                  if(length(chromosomes) > 1){
##                    if(chrom == chromosomes[2]){
##                      label <- chrom
##                    } else {
##                      label <- strsplit(chrom, "chr")[[1]][2]
##                    }
##                  }
##                  if(side.last == 1) mtext(label, 3, line=2.5, cex=cex.chr)
##                }
##                if(i == S) {
##                  if(length(chromosomes) <= 6){
##                    probs <- seq(0, 1, by=1/(xTicks+2))
##                    probs <- probs[-c(1, length(probs))]
##                    quants <- quantile(xlim, probs)                                                                                
##                    tcl <- NULL
##                    las <- 1
##                  } else {
##                    probs <- c(0, 0.5, 1)
##                    quants <- quantile(xlim, probs)                                                            
##                    labels <- c("", round(chromosomeSize/1e6, 0), "")
##                    tcl <- 0
##                    las <- 3
##                  }
##                  if(length(chromosomes) <= 6)
##                    labels <- as.character(round(quants/1e6, 0))
##                  if(side.last == 1) side <- 3 else side <- 1
##                  axis(side, at=quants, outer=TRUE, labels=labels, tcl=tcl,
##                       cex.axis=cex.axis, col=col.axis, col.axis=col.axis,
##                       las=las, line=0, lwd=1, mgp=c(2, 0.5, 0))
##                  if(side == 1){
##                    if(chrom == chromosomes[1]){
##                      label <- chrom
##                    } else{
##                      strsplit(chrom, "chr")[[1]][2]
##                    }
##                    mtext(label, 1, line=3.5, cex=cex.chr)
##                  }
##                  side.last <- side
##                  if(!alternate.xaxis) side.last <- 3
##                }
##              }
##              k <- k+1
##            }
##            mtext("Mb ", 1, at=0, line=0, outer=TRUE, cex=cex.axis,
##                  col=col.axis, adj=1, las=las)
##            if(alternate.xaxis & length(chromosomes) >= 2)
##              mtext("Mb ", 3, at=0, line=0, outer=TRUE, cex=cex.axis,
##                    col=col.axis, adj=0, las=las)            
##
##            ###########################################################################
##            ##Plot summary statistics
##            object <- obj
##            if(summaryPanel){
##              if(sum(chromosome(object) != "chrX") > 0){
##                obj <- object[chromosome(object) != "chrX", ]
##              } else { obj <- object}
##              ht <- colMeans(ifelse(calls(obj) == 2, 1, 0))
##              ho <- colMeans(ifelse(calls(obj) == 1 | calls(obj) == 3, 1, 0))
##              stats <- cbind(ht, ho)
##              stats <- round(stats, digits)
##            }            
##            if(summaryPanel){
##              showSummary <- function(x){
##                par(mar=rep(0,4))
##                plot(0:1, 0:1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
##                legend(legend.stats, legend=c(
##                                    paste(x["ho"], " %AA/BB", sep=""),
##                                    paste(x["ht"], " %AB", sep="")), bty="n",
##                       title=substr(x["samplenames"], 1, 10),
##                       y.intersp=1.5,
##                       cex=cex.legend,
##                       text.col=c(colAA, colAB))
##              }
##              stats <- data.frame(stats); stats$samplenames <- samplenames
##              apply(stats, 1, showSummary)
##            }
##            ##plot legend
##            if(legend){
##              op <- par(bg=legend.col)
##              if(summaryPanel){
##                legend.bty <- "n"
##              }
##              legend(legend.pch, pch=20, 
##                     col=c(colAA, colAB),
##                     legend=c("AA/BB", "AB"), cex=cex.legend,
##                     y.intersp=1.5,
##                     bty=legend.bty, pt.cex=cex.legend*1.5)
##            }
##          })
##
##setMethod("plotChromosome", "AnnotatedSnpCallSet",
##          function(object,  col="black", colAA="blue", colAB="red", colNC="green3",
##                   colCentromere="bisque", centromereBorder=NA, pch=".",
##                   cex=1, cexAA=1, cexAB=1, cexNC=1, bty="o", legend.bty="n",
##                   xlim=NULL, ylim=NULL, cex.axis=0.8, cex.legend=1, cex.main=1,
##                   xlab="Mb", ylab="copy number", ps=16, digits=3, lwdCn=2,
##                   legendStats=TRUE, legendPch=TRUE, legend.loc="topleft",
##                   xaxt="s", yaxt="s", xaxs="i", yaxs="r", main="",
##                   mar=c(4, 4, 1, 0.2), panel.xaxis=FALSE, panel.yaxis=FALSE,
##                   yTicks=5, xTicks=5, log="", bw=FALSE, ...){
##            browser()
##            if(bw){
##              colAA <- gray(0.6)
##              colAB <- "black"
##              colNC <- gray(0.8)
##              col.axis <- "black"
##              colCentromere <- gray(0.7)
##            }
##            chrom <- unique(SNPchip::chromosome(object))
##            if(length(chrom) > 1) stop("object must contain only 1 chromosome")
##            if(is.null(ylim)) ylim <- c(floor(min(copyNumber(object), na.rm=TRUE)), ceiling(max(copyNumber(object), na.rm=TRUE)))
##            ##If more than 1 sample is present, it only plots the first
##            if(dim(object)[2] > 1){
##              warning("Only plotting the first sample")
##              object <- object[,1]
##            }
##
##            y <- copyNumber(object)
##            y[y < ylim[1]] <- ylim[1]
##            y[y > ylim[2]] <- ylim[2]
##            x <- position(object)
##
##            chrAnn <- SNPchip::chromosomeAnnotation(object)[chrom, ]
##            chromosomeSize <- chrAnn$chromosomeSize
##            if(is.null(xlim)){
##              xlim <- c(0, chromosomeSize)
##              xlim[2] <- xlim[2]+2e6
##            }
##            if(log == "y") ylimit <- NULL else ylimit <- ylim
##            par(las=1, ps=ps, mar=mar, "ylog")
##            plot(x, y,
##                 pch=pch,
##                 cex.axis=cex.axis,
##                 cex.main=cex.main,
##                 xlab=xlab,
##                 ylab=ylab,
##                 main=main,
##                 xaxs=xaxs,                 
##                 yaxs=yaxs,
##                 xlim=xlim,                 
##                 ylim=ylimit,
##                 xaxt="n",
##                 yaxt=yaxt,
##                 type="n",
##                 bty=bty,
##                 log=log)
##            if(panel.xaxis) axis(side=1, at=pretty(xlim, n=xTicks),
##                                 labels=as.character(round(pretty(xlim, n=xTicks)/1e6,0)),
##                                 cex.axis=cex.axis)
##            if(panel.yaxis) axis(side=2, at=pretty(ylim, n=yTicks), outer=TRUE,
##                                 cex.axis=cex.axis)
##            calls <- calls(object)
##            hom <- calls == 1 | calls == 3
##            het <- calls == 2
##            points(jitter(x[hom], 1), jitter(y[hom], 1),
##                   col = colAA,
##                   cex = cexAA,
##                   pch = pch)
##            points(jitter(x[het],1), jitter(y[het], 1),
##                   col = colAB,
##                   cex = cexAB,
##                   pch = pch)
##            if(any(calls==4)){
##              points(x[calls == 4], y[calls == 4, 1],
##                     col = colNC,
##                     cex = cexNC,
##                     pch = pch)
##            }
##            ######################################################################
##            ##Draw centromere
##            ######################################################################
##            centromere <- chrAnn[chrom, 1:2]
##            rect(xleft=centromere[[1]], ybottom=ylim[1],
##                 xright=centromere[[2]], ytop=ylim[2], col=colCentromere,
##                 border=centromereBorder)
##            cn <- mean(copyNumber(object))
##            cn.sd <- sd(copyNumber(object))
##            ht <- mean(ifelse(calls(object) == 2, 1, 0))
##            ho <- mean(ifelse(calls(object) == 1 | calls(object) == 3, 1, 0))
##            stats <- c(cn, cn.sd, ht, ho)
##            stats <- round(stats, digits)
##            names(stats) <- c("cn", "cn.sd", "ht", "ho")
##            if(legendStats){
##              par(bg="antiquewhite1")
##              legend(legend.loc, legend = c(substr(x["samplenames"], 1, min(nchar(x["samplenames"]),10)),
##                                  paste(stats["ho"], " %AA/BB", sep = ""),
##                                  paste(stats["ht"], " %AB", sep = ""),
##                                  paste(stats["cn"], " avg CN"),
##                                  paste(stats["cn.sd"], " sd")), bty = legend.bty, cex = cex.legend,
##                     text.col = c("black", colAA, colAB, "black", "black"))
##            }
##            if(legendPch){
##              legend("topright",
##                     pch=20,
##                     col=c(colAA, colAB),
##                     legend=c("AA/BB", "AB"),
##                     cex=cex.legend,
##                     bty="n",
##                     pt.cex=cex.legend*1.5)
##              par(bg="white")
##            }
##          })

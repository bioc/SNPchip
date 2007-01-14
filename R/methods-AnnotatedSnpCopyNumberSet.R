setMethod("chromosomeAnnotation", "AnnotatedSnpCopyNumberSet", function(object) object@chromosomeAnnotation)
setReplaceMethod("chromosomeAnnotation", c("AnnotatedSnpCopyNumberSet", "data.frame"), function(object, value){
  object@chromosomeAnnotation <- value
  object
})
setMethod("initialize", "AnnotatedSnpCopyNumberSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   featureData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   copyNumber = new("matrix"),
                   cnConfidence = new("matrix"),
                   chromosomeAnnotation = data.frame())
          {
            .Object@assayData <- assayDataNew(copyNumber = copyNumber, cnConfidence = cnConfidence)
            .Object@phenoData <- phenoData
            .Object@annotation <- annotation
            .Object@featureData <- featureData
            .Object@chromosomeAnnotation <- chromosomeAnnotation
            .Object
          })

#setMethod("alleleA", "AnnotatedSnpSet", function(object) alleleA(featureData(object)))
#setMethod("alleleB", "AnnotatedSnpSet", function(object) alleleB(featureData(object)))
#setMethod("chromosome", "AnnotatedSnpSet", function(object) chromosome(featureData(object)))
#setMethod("dbSnpId", "AnnotatedSnpSet", function(object) dbSnpId(featureData(object)))
#setMethod("enzyme", "AnnotatedSnpSet", function(object) enzyme(featureData(object)))
#setMethod("position", "AnnotatedSnpSet", function(object) position(featureData(object)))
#setMethod("probeSetId", "AnnotatedSnpSet", function(object) probeSetId(featureData(object)))
setAs("AnnotatedSnpCopyNumberSet", "AnnotatedSnpCopyNumberSetList",
      function(from){
        chr <- as.character(unique(chromosome(from)))
        snpSetList <- list()
        for(i in 1:length(chr)){
          snpSetList[[i]] <- from[chromosome(from) == chr[i], ]
        }
        new("AnnotatedSnpCopyNumberSetList", snpSetList=snpSetList)
      })

setMethod("plotSnp", "AnnotatedSnpCopyNumberSet",
          function(object, chromosomes, samples, ylim=NULL, xlim=NULL, col="black",
                   colCentromere="bisque", col.axis="brown", pch=".", cex=1, cex.axis=0.8,
                   cex.legend=1, cex.chr=0.8, oma=c(5, 3, 4, 0.5), mar=c(0, 0, 0, 0.2),
                   width.right=NULL, summaryPanel=FALSE, showLayout=TRUE, plotIt=TRUE,
                   digits=3, legend=TRUE, legend.stats="left", legend.pch="topleft",
                   legend.bty="o", legend.col="white", alternate.xaxis=TRUE,
                   xaxis=TRUE, jitter=TRUE, factor=0.1, yTicks=5, xTicks=2, bty="n",
                   bw=FALSE, ...){
            ##If black and white, change the default color scheme
            if(bw){
              col.axis <- gray(0)
              colCentromere <- gray(0.7)
            }
            chromosomes <- paste("chr", chromosomes, sep="")
            chromosomes[chromosomes == "chr23"] <- "chrX"
            chromosomes[chromosomes == "chr24"] <- "chrY"
            object <- object[chromosome(object) %in% chromosomes, samples]
            obj <- object
            
            if(is.null(ylim)) ylim <- c(ceiling(min(copyNumber(object), na.rm=TRUE)), floor(max(copyNumber(object), na.rm=TRUE)))
            samplenames <- sampleNames(object)
            chrAnn <- chromosomeAnnotation(object)[chromosomes,]

            object <- as(object, "AnnotatedSnpCopyNumberSetList")            
            object <- snpSetList(object)
            names(object) <- chromosomes
            k <- 1
            S <- dim(object[[1]])[2]
            N <- length(object)
            
            widths <- chrAnn$chromosomeSize
            widths <- widths/min(widths)
            if(summaryPanel){
              if(is.null(width.right)){
                width.right <- length(chromosomes)/1.5
              }
              N <- N+1
              widths <- c(widths, width.right)
            }
            nf <- layout(matrix(1:(S*N),
                                nc=N,
                                byrow=FALSE), widths=widths)

            ##Option to return just the layout
            if(!plotIt) {
              layout.show(nf)
              return()
            }
            par(mar=mar, oma=oma)
            for(chrom in chromosomes){
              for(i in 1:S){
                if(k == 1){ yaxis=TRUE; side.last <- 3} else yaxis=FALSE
                chromosomeSize <- chrAnn[chrom, 3]                
                xlim <- c(0, chromosomeSize)
                plotChromosome(object[[chrom]][, i], ylim=ylim, xlim=xlim,
                               yaxt="n", xaxt="n", yaxs="r", legendStats=FALSE,
                               legendPch=FALSE, centromereBorder=NA, xlab="", ylab="",
                               bty=bty, panel.xaxis=FALSE, panel.yaxis=yaxis,
                               yTicks=yTicks, mar=mar)
                if(side.last == 1) mtext(label, 3, line=2.5, cex=cex.chr)
                if(i == S) {
                  if(length(chromosomes) <= 6){
                    probs <- seq(0, 1, by=1/(xTicks+2))
                    probs <- probs[-c(1, length(probs))]
                    quants <- quantile(xlim, probs)                                                                                
                    tcl <- NULL
                    las <- 1
                  } else {
                    probs <- c(0, 0.5, 1)
                    quants <- quantile(xlim, probs)                                                            
                    labels <- c("", round(chromosomeSize/1e6, 0), "")
                    tcl <- 0
                    las <- 3
                  }
                  if(length(chromosomes) <= 6)
                    labels <- as.character(round(quants/1e6, 0))
                  if(side.last == 1) side <- 3 else side <- 1
                  axis(side, at=quants, outer=TRUE, labels=labels, tcl=tcl, cex.axis=cex.axis,
                       col=col.axis, col.axis=col.axis, las=las, line=0, lwd=1,
                       mgp=c(2, 0.5, 0))
                  if(side == 1){
                    if(chrom == chromosomes[1]){
                      label <- chrom
                    } else{
                      strsplit(chrom, "chr")[[1]][2]
                    }
                    mtext(label, 1, line=3.5, cex=cex.chr)
                  }
                  side.last <- side
                  if(!alternate.xaxis) side.last <- 3
                }
              }
              k <- k+1
            }
            mtext("Mb ", 1, at=0, line=0, outer=TRUE, cex=cex.axis,
                  col=col.axis, adj=1, las=las)
            if(alternate.xaxis & length(chromosomes) >= 2)
              mtext("Mb ", 3, at=0, line=0, outer=TRUE, cex=cex.axis,
                    col=col.axis, adj=0, las=las)            

            ###########################################################################
            ##Plot summary statistics
            object <- obj
            if(summaryPanel){
              if(sum(chromosome(object) != "chrX") > 0){
                obj <- object[chromosome(object) != "chrX", ]
              } else { obj <- object}
              cn <- colMeans(as.matrix(copyNumber(obj)))
              cn.sd <- apply(copyNumber(obj), 2, stats::sd)                
              stats <- cbind(cn, cn.sd)
              stats <- round(stats, digits)
            }            

            if(summaryPanel){
              showSummary <- function(x){
                par(mar=rep(0,4))
                plot(0:1, 0:1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
                legend(legend.stats, legend=c(
                                    paste(x["cn"], " avg CN"),
                                    paste(x["cn.sd"], " sd")), bty="n",
                       title=substr(x["samplenames"], 1, 10),
                       y.intersp=1.5,
                       cex=cex.legend,
                       text.col=c("black", "black", "black"))
              }
              stats <- data.frame(stats); stats$samplenames <- samplenames
              apply(stats, 1, showSummary)
            }
          })

setMethod("initialize", "AnnotatedSnpCopyNumberSetList",
          function(.Object,
                   snpSetList = list()){
            .Object@snpSetList <- snpSetList
            .Object
          })

setAs("AnnotatedSnpCopyNumberSetList", "AnnotatedSnpCopyNumberSet",
      function(from){
        chrAnn <- chromosomeAnnotation(snpSetList(from)[[1]])

        ##Combine assayData
        ad <- lapply(snpSetList(from), assayData)
        f.cn <- function(x) x$copyNumber
        f.cnConf <- function(x) x$cnConfidence
        cn.l <- lapply(ad, f.cn)
        cnConf.l <- lapply(ad, f.cnConf)
        cn <- do.call("rbind", cn.l)
        cnConf <- do.call("rbind", cnConf.l)

        ##Combine featureData
        fd <- lapply(snpSetList(from), featureData)
        pd <- lapply(fd, pData)
        fd.combined <- do.call("rbind", pd)
        fd <- new("AnnotatedDataFrame", data = fd.combined,
                  varMetadata = varMetadata(featureData(snpSetList(from)[[1]])))
        obj <- new("oligoSnpCopyNumberSet",
                   copyNumber = cn,
                   cnConfidence = cnConf,
                   featureData = fd,
                   phenoData = phenoData(snpSetList(from)[[1]]),
                   experimentData = experimentData(snpSetList(from)[[1]]))
        obj <- as(obj, "AnnotatedSnpCopyNumberSet")
        chromosomeAnnotation(obj) = chrAnn
        obj
      })

setMethod("snpSetList", "AnnotatedSnpCopyNumberSetList", function(object) object@snpSetList)
setMethod("plotCytoband", "AnnotatedSnpCopyNumberSet",
          function(object, ...){
          obj <- as(object, "AnnotatedSnpSet")
          plotCytoband(obj, ...)
        })

setMethod("plotChromosome", "AnnotatedSnpCopyNumberSet",
          function(object, col="black", colCentromere="bisque", centromereBorder=NA,
                   pch=".", cex=1, bty="o", legend.bty="n", xlim=NULL, ylim=NULL,                   
                   cex.axis=0.8, cex.legend=1, cex.main=1, xlab="Mb", ylab="copy number",
                   ps=16, digits=3, lwdCn=2, legendStats=TRUE, legendPch=TRUE,
                   legend.loc="topleft", xaxt="s", yaxt="s", xaxs="i", yaxs="r",
                   main="", mar=c(3, 3, 1, 0.2), panel.xaxis=FALSE, panel.yaxis=FALSE,
                   yTicks=5, xTicks=5, log="", bw=FALSE, ...){
            if(bw){
              col.axis <- "black"
              colCentromere <- gray(0.7)
            }
            chrom <- unique(SNPchip::chromosome(object))
            if(length(chrom) > 1) stop("object must contain only 1 chromosome")
            if(is.null(ylim)) ylim <- c(floor(min(copyNumber(object), na.rm=TRUE)), ceiling(max(copyNumber(object), na.rm=TRUE)))
            ##If more than 1 sample is present, it only plots the first
            if(dim(object)[2] > 1){
              warning("Only plotting the first sample")
              object <- object[,1]
            }

            y <- copyNumber(object)
            y[y < ylim[1]] <- ylim[1]
            y[y > ylim[2]] <- ylim[2]
            x <- position(object)

            chrAnn <- SNPchip::chromosomeAnnotation(object)[chrom, ]
            chromosomeSize <- chrAnn$chromosomeSize
            if(is.null(xlim)){
              xlim <- c(0, chromosomeSize)
              xlim[2] <- xlim[2]+2e6
            }
            if(log == "y") ylimit <- NULL else ylimit <- ylim
            par(las=1, ps=ps, mar=mar, "ylog")
            plot(x, y, pch=pch, cex.axis=cex.axis, cex.main=cex.main, xlab=xlab, ylab=ylab,
                 main=main, xaxs=xaxs, yaxs=yaxs, xlim=xlim, ylim=ylimit, xaxt="n",
                 yaxt=yaxt, bty=bty, log=log, col=col)
            if(panel.xaxis) axis(side=1, at=pretty(xlim, n=xTicks),
                                 labels=as.character(round(pretty(xlim, n=xTicks)/1e6,0)),
                                 cex.axis=cex.axis)
            if(panel.yaxis) axis(side=2, at=pretty(ylim, n=yTicks), outer=TRUE,
                                 cex.axis=cex.axis)

            ######################################################################
            ##Draw centromere
            ######################################################################
            centromere <- chrAnn[chrom, 1:2]
            rect(xleft=centromere[[1]], ybottom=ylim[1],
                 xright=centromere[[2]], ytop=ylim[2], col=colCentromere,
                 border=centromereBorder)
            cn <- mean(copyNumber(object))
            cn.sd <- sd(copyNumber(object))
            stats <- c(cn, cn.sd)
            stats <- round(stats, digits)
            names(stats) <- c("cn", "cn.sd")
            if(legendStats){
              legend(legend.loc, legend = c(substr(x["samplenames"], 1, min(nchar(x["samplenames"]),10)),
                                   paste(stats["cn"], " avg CN"),
                                   paste(stats["cn.sd"], " sd")), bty = legend.bty, cex = cex.legend,
                     text.col = c("black", "black", "black"))
            }
          })

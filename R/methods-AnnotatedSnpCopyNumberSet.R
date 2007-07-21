setMethod("chromosomeAnnotation", "AnnotatedSnpCopyNumberSet", function(object) object@chromosomeAnnotation)
setReplaceMethod("chromosomeAnnotation", c("AnnotatedSnpCopyNumberSet", "data.frame"), function(object, value){
  object@chromosomeAnnotation <- value
  object
})

setMethod("initialize", "AnnotatedSnpCopyNumberSet",
          function(.Object,
                   assayData = assayDataNew(
                     copyNumber=copyNumber,
                     cnConfidence = cnConfidence),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=character(),
                   copyNumber=new("matrix"),
                   cnConfidence=matrix(numeric(), nrow=nrow(copyNumber),
                     ncol=ncol(copyNumber),
                     dimnames=dimnames(copyNumber)),
                 chromosomeAnnotation = data.frame()){
            .Object@assayData <- assayData
            .Object@phenoData <- phenoData
            .Object@annotation <- annotation
            .Object@featureData <- featureData
            .Object@chromosomeAnnotation <- chromosomeAnnotation
            .Object
          })

setMethod("alleleA", "AnnotatedSnpCopyNumberSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "AnnotatedSnpCopyNumberSet", function(object) alleleB(featureData(object)))
setMethod("chromosome", "AnnotatedSnpCopyNumberSet", function(object) chromosome(featureData(object)))
setMethod("dbSnpId", "AnnotatedSnpCopyNumberSet", function(object) dbSnpId(featureData(object)))
setMethod("enzyme", "AnnotatedSnpCopyNumberSet", function(object) enzyme(featureData(object)))
setMethod("fragmentLength", "AnnotatedSnpCopyNumberSet", function(object) fragmentLength(fData(object)))
setMethod("position", "AnnotatedSnpCopyNumberSet", function(object) position(featureData(object)))

setMethod("plotSnp", "AnnotatedSnpCopyNumberSet",
          function(object,
                   chromosomes,
                   samples,
                   snpId=c(NA, 1e6),  ##SNP identifier, distance to plot on either side
                   oma=c(5, 4, 4, 0.5),
                   mar=c(0.5, 0, 0.5, 0.2),
                   width.right=NULL,
                   showLayout=FALSE,
                   plot=TRUE,
                   col="black",
                   bg="white",
                   cex=2,
                   pch=".", ##, length(unique(as.vector(calls(object))))),                                      
                   col.centromere=c("bisque", NA),
                   bw=FALSE,
                   lwd=2,
                   col.axis="brown",
                   cex.axis=1,
                   cex.main=1,
                   cex.lab=1,
                   xlim=NULL,                   
                   ylim=NULL,
                   log="",                   
                   xaxis.side=rep(c(1, 3), length.out=length(chromosomes)),
                   xaxs="i",                   
                   xaxt="s",
                   yaxs="i",
                   yaxt="n",
                   xlab="",
                   ylab="",
                   lab=c(2, 5, 7), ##see par
                   adj=0,
                   main="",
                   legend=c(TRUE, TRUE), ##legend for stats,  legend for plotting symbols
                   legend.panel=c(TRUE, FALSE),  ##plot legend on separate panel?
                   legend.location=c("topright", "bottomright"),
                   legend.bty=c("n", "n"),
                   legend.col=c("white", "white"),
                   ncol=2,
                   cex.legend=c(1, 1),
                   digits=3,
                   pt.cex=cex.legend*1.5,
                   bty="n",
                   addCytoband=FALSE,
                   height.cytoband=0.5, ##relative to plotting region for samples
                   ...){
            object <- as(object, "AnnotatedSnpSet")
            calls(object) <- matrix(1, nrow(object), ncol(object))
            plotSnp(object=object,
                    chromosomes=chromosomes,
                    samples=samples,
                    snpId=snpId,
                    oma=oma,
                    mar=mar,
                    width.right=width.right,
                    showLayout=showLayout,
                    plot=plot, 
                    col=col[1],
                    bg=bg[1],
                    cex=cex[1],
                    pch=pch[1],
                    col.centromere=col.centromere,
                    bw=bw,
                    lwd=lwd,
                    col.axis=col.axis,
                    cex.axis=cex.axis,
                    cex.main=cex.main,
                    cex.lab=cex.lab,
                    xlim=xlim,
                    ylim=ylim,
                    log=log,
                    xaxis.side=xaxis.side,
                    xaxs=xaxs,
                    xaxt=xaxt,
                    yaxs=yaxs,
                    yaxt=yaxt,
                    xlab=xlab,
                    ylab=ylab,
                    lab=lab,
                    adj=adj,
                    main=main,
                    legend=legend,
                    legend.panel=legend.panel,
                    legend.location=legend.location,
                    legend.bty=legend.bty,
                    legend.col=legend.col,
                    ncol=ncol,
                    cex.legend=cex.legend,
                    digits=digits,
                    pt.cex=pt.cex,
                    bty=bty,
                    addCytoband=addCytoband,
                    height.cytoband=height.cytoband, ...)
          })
          
setAs("AnnotatedSnpCopyNumberSet", "AnnotatedSnpSet",
      function(from){
        object <- new("AnnotatedSnpSet",
                      calls=matrix(NA, nrow(from), ncol(from)),
                      copyNumber=copyNumber(from),
                      callsConfidence=matrix(NA, nrow(from), ncol(from)),
                      cnConfidence=cnConfidence(from),
                      phenoData=phenoData(from),
                      featureData=featureData(from),
                      chromosomeAnnotation=chromosomeAnnotation(from),
                      annotation=annotation(from))
        object
      })

setMethod("show", "AnnotatedSnpCopyNumberSet", function(object) {
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


#setMethod("plotCytoband", "AnnotatedSnpCopyNumberSet",
#          function(object, ...){
#          obj <- as(object, "AnnotatedSnpSet")
#          plotCytoband(obj, ...)
#        })
#
##setMethod("plotChromosome", "AnnotatedSnpCopyNumberSet",
##          function(object, col="black", colCentromere="bisque", centromereBorder=NA,
##                   pch=".", cex=1, bty="o", legend.bty="n", xlim=NULL, ylim=NULL,                   
##                   cex.axis=0.8, cex.legend=1, cex.main=1, xlab="Mb", ylab="copy number",
##                   ps=16, digits=3, lwdCn=2, legendStats=TRUE, legendPch=TRUE,
##                   legend.loc="topleft", xaxt="s", yaxt="s", xaxs="i", yaxs="r",
##                   main="", mar=c(3, 3, 1, 0.2), panel.xaxis=FALSE, panel.yaxis=FALSE,
##                   yTicks=5, xTicks=5, log="", bw=FALSE, ...){
##            if(bw){
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
##            plot(x, y, pch=pch, cex.axis=cex.axis, cex.main=cex.main, xlab=xlab, ylab=ylab,
##                 main=main, xaxs=xaxs, yaxs=yaxs, xlim=xlim, ylim=ylimit, xaxt="n",
##                 yaxt=yaxt, bty=bty, log=log, col=col, cex=cex)
##            if(panel.xaxis) axis(side=1, at=pretty(xlim, n=xTicks),
##                                 labels=as.character(round(pretty(xlim, n=xTicks)/1e6,0)),
##                                 cex.axis=cex.axis)
##            if(panel.yaxis) axis(side=2, at=pretty(ylim, n=yTicks), outer=TRUE,
##                                 cex.axis=cex.axis)
##
##            ######################################################################
##            ##Draw centromere
##            ######################################################################
##            centromere <- chrAnn[chrom, 1:2]
##            rect(xleft=centromere[[1]], ybottom=ylim[1],
##                 xright=centromere[[2]], ytop=ylim[2], col=colCentromere,
##                 border=centromereBorder)
##            cn <- mean(copyNumber(object))
##            cn.sd <- sd(copyNumber(object))
##            stats <- c(cn, cn.sd)
##            stats <- round(stats, digits)
##            names(stats) <- c("cn", "cn.sd")
##            if(legendStats){
##              legend(legend.loc, legend = c(substr(x["samplenames"], 1, min(nchar(x["samplenames"]),10)),
##                                   paste(stats["cn"], " avg CN"),
##                                   paste(stats["cn.sd"], " sd")), bty = legend.bty, cex = cex.legend,
##                     text.col = c("black", "black", "black"))
##            }
##          })

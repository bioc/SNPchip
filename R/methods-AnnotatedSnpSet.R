setMethod("chromosomeAnnotation", "AnnotatedSnpSet", function(object) object@chromosomeAnnotation)
setReplaceMethod("chromosomeAnnotation", c("AnnotatedSnpSet", "data.frame"),
                 function(object, value){
                   object@chromosomeAnnotation <- value
                   object
                 })

setMethod("alleleA", "AnnotatedSnpSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "AnnotatedSnpSet", function(object) alleleB(featureData(object)))
##setMethod("chromosome", "AnnotatedSnpSet", function(object) chromosome(featureData(object)))
setMethod("dbSnpId", "AnnotatedSnpSet", function(object) dbSnpId(featureData(object)))
setMethod("enzyme", "AnnotatedSnpSet", function(object) enzyme(featureData(object)))
setMethod("fragmentLength", "AnnotatedSnpSet", function(object) fragmentLength(fData(object)))
setMethod("position", "AnnotatedSnpSet", function(object) position(featureData(object)))
##setMethod("probeSetId", "AnnotatedSnpSet", function(object) probeSetId(featureData(object)))

setMethod("initialize", "AnnotatedSnpSet",
          function(.Object,
                   assayData = assayDataNew(
                     calls=calls,
                     callsConfidence = callsConfidence,
                     copyNumber=copyNumber,
                     cnConfidence=cnConfidence, ...),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=character(),
                   chromosomeAnnotation=data.frame(),
                   calls=new("matrix"),
                   callsConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   copyNumber=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   cnConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   ...){
            data(chromosomeAnnotation)
            .Object@assayData <- assayData
            .Object@phenoData <- phenoData
            .Object@annotation <- annotation
            .Object@featureData <- featureData
            .Object@chromosomeAnnotation <- chromosomeAnnotation
            .Object@experimentData <- experimentData
            .Object
          })


setMethod("initialize", "SnpSet",
          function(.Object,
                   assayData = assayDataNew(
                     call = call,
                     callProbability = callProbability, ...),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   call = new("matrix"),
                   callProbability = matrix(
                     numeric(), nrow=nrow(call), ncol=ncol(call),
                     dimnames=dimnames(call)),
                   ... ) {
            callNextMethod(.Object,
                           assayData = assayData,
                           phenoData = phenoData,
                           featureData = featureData,
                           experimentData = experimentData,
                           annotation = annotation)
          })


setMethod("show", "AnnotatedSnpSet", function(object) {
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

setMethod("summary", "AnnotatedSnpSet", function(object, digits=3, ...){
  ##Calculate mean,sd copy number and prop no calls, prop het calls
  ##for each chromosome in each sample.  Return an S x 23 matrix for
  ##each.
  x <- list()
##  chrList <- as.list(unique(chromosome(object)))
##  if(sum(chromosome(object) != "X") > 0) chrList[[length(chrList) + 1]] <- "autosome"
  cn <- split(copyNumber(object), chromosome(object))
  cn <- lapply(cn, matrix, ncol=ncol(object), byrow=TRUE)
  x[[1]] <- round(sapply(cn, colMeans), digits)
##  chrStats <- list()
##  chromosomeMeanCn <- function(chrom, object){
##    if(chrom != "autosome"){
##      cn <- copyNumber(object)[chromosome(object) == chrom, ]      
##    } else  cn <- copyNumber(object)[chromosome(object) != "X",]
##    colMeans(as.matrix(cn))
##  }
##  chrStats[[1]] <- sapply(chrList, chromosomeMeanCn, object)
  colSds <- function(x) apply(x, 2, "sd")
  x[[2]] <- round(sapply(cn, colSds), digits)
  
  ##standard deviation of copy number for each chromosome
##  chromosomeSdCn <- function(chrom, object){
##    if(chrom != "autosome"){
##      cn <- copyNumber(object)[chromosome(object) == chrom, ]
##    } else  cn <- copyNumber(object)[chromosome(object) != "X",]
##    cn <- apply(as.matrix(cn), 2, stats::sd)
##  }
##  chrStats[[2]] <- sapply(chrList, chromosomeSdCn, object)

  ##Proportion of no calls for each chromosome
  calls <- split(calls(object), chromosome(object))
  calls <- lapply(calls, matrix, ncol=ncol(object), byrow=TRUE)
  calls.missing <- lapply(calls, "==", 4)
  x[[3]] <- round(sapply(calls.missing, colMeans), digits)
                              
#  propNoCalls <- function(chrom, object){
#    if(chrom != "autosome"){
#      calls <- calls(object)[chromosome(object) == chrom, ]
#    } else  calls <- calls(object)[chromosome(object) != "X", ]
#    calls <- as.matrix(calls)
#    colMeans(calls == 4)    
#  }
#  chrStats[[3]] <- sapply(chrList, propNoCalls, object)

  ##Proportion of homozygous calls given that a call was made
  homozygousCall <- function(x) x == 1 | x == 3
  calls.homozygous <- lapply(calls, homozygousCall)
  x[[4]] <- round(sapply(calls.homozygous, colMeans), digits)
  x[[5]] <- round(1-x[[4]]-x[[3]], digits)
  names(x) <- c("avgCN", "sdCN", "%NoCalls", "%Hom", "%Het")

  overall <- lapply(x, colMeans)
  x <- mapply("rbind", x, overall, SIMPLIFY=FALSE)

  labelRows <- function(x, sn){
    rownames(x) <- c(sn, "overall")
    return(x)
  }
  x <- lapply(x, labelRows, sn=sampleNames(object))

  reorderX <- function(x){
    y <- c(1:22, "X", "Y")
    y <- y[y %in% colnames(x)]
    x <- x[, y]
  }
  xx <- lapply(x, reorderX)
  
  ##Calculate grand average
##  if(dim(object)[2] > 1){
##    grand<- list()
##    grand[[1]] <- colMeans(chrStats$avgCopyNumber)
##    grand[[2]] <- apply(chrStats$sdCopyNumber, 2, stats::sd)
##    grand[[3]] <- colMeans(chrStats$propNoCalls)
##    grand[[4]] <- colMeans(chrStats$propHo)
##    grand <- do.call("rbind", grand)
##    rownames(grand) <- c("overall mean", "sd of means", "avg prop no calls", "avg prop AA/BB among calls")
##  } else grand <- NULL
##  stats <- list(chromosome=chrStats, overall=grand)
  return(xx)
})

setMethod("plotSnp", "AnnotatedSnpSet",
          function(object,
                   chromosomes,
                   samples,
#                   start.Mb=NULL, ##position in Mb
#                   stop.Mb=NULL,  ##position in Mb
                   snpId=c(NA, 1e6),  ##SNP identifier, distance to plot on either side
                   ##################################################
                   ##margins
                   oma=c(5, 4, 4, 0.5),
                   mar=c(0.5, 0, 0.5, 0.2),
##                   par=TRUE,
                   width.right=NULL,
                   showLayout=FALSE,
                   plot=TRUE,
                   ##################################################
                   ##plotting symbols and colors
                   col=c("royalblue", "red", "royalblue", "green3"), 
                   bg=rep("white", 4),
                   cex=c(2, 3, 2, 2),
                   pch=rep(".", 4), ##, length(unique(as.vector(calls(object))))),                                      
                   col.centromere=c("bisque", NA),
                   bw=FALSE,
                   lwd=2,
                   col.axis="brown",
                   cex.axis=1,
                   cex.main=1,
                   cex.lab=1,
                   ##################################################
                   ##Axis
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
                   ##################################################
                   ##legend
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
                   ##################################################
                   addCytoband=FALSE,
                   height.cytoband=0.5, ##relative to plotting region for samples
                   useLayout=TRUE,
                   ...){
            if(is.null(chromosomeAnnotation(object))){
              data(chromosomeAnnotation)
              chromosomeAnnotation(object) <- chromosomeAnnotation
            }

            chromosome.order <- as.character(chromosomes)
            if(any(legend.panel)) oma[4] <- 0
            
            ##If black and white, change the default color scheme
            if(bw){
              col <- c(gray(0.7, gray(0), gray(0.7), gray(0.9)))
              col.axis <- gray(0)
              col.centromere <- c(gray(0.7), NA)
            }
            chromosomeAnnotation <- chromosomeAnnotation(object)[chromosome.order, ]
            
            object <- object[chromosome(object) %in% chromosome.order, samples]
            S <- ncol(object)
            ##If plotting cytoband, need an extra row
            if(addCytoband){
              S <- S+1 
              xaxis.side <- rep(3, length.out=length(chromosome.order))
              oma <- c(6, 4, 4, 0.5)
            }
            obj <- object
            
            if(is.null(ylim)) ylim <- c(ceiling(min(copyNumber(object), na.rm=TRUE)), floor(max(copyNumber(object), na.rm=TRUE)))
            if(log == "y" & ylim[1] <= 0) {
              print("Must specify lower y-axis limit > 0 when plotting on log scale")
              ylim[1] <- 0.5
            }
            objList <- split(object, chromosome(object))

            k <- 1
            N <- length(objList)

            if(useLayout){
              par(oma=oma, mar=mar)
              widths <- chromosomeAnnotation$chromosomeSize
              widths <- widths/min(widths)
              if(any(legend.panel)){
                if(is.null(width.right)){
                  width.right <- length(chromosome.order)/1.5
                }
                N <- N+1
                widths <- c(widths, width.right)
              }
              if(addCytoband) heights <- c(rep(1, S-1), height.cytoband) else heights <- rep(1, S)
              nf <- layout(matrix(1:(S*N),
                                  nc=N,
                                  byrow=FALSE), widths=widths, heights=heights)
            }
            ##Option to return just the layout
            if(!plot) {
              layout.show(nf)
              return()
            }

            ##Make sure colors recycle correctly
            N <- sort(unique(as.vector(calls(object))))
            col <- col[N]
            pch <- pch[N]
            cex <- cex[N]
            bg <- bg[N]
            if(addCytoband){  data(cytoband); S <- S-1}
            names(xaxis.side) <- chromosome.order
            
            for(i in chromosome.order){
              for(j in 1:S){
                if(i == chromosome.order[1]) yaxis <- TRUE
                obj <- objList[[i]][, j]
                cn <- as.vector(copyNumber(obj))
                calls <- as.vector(calls(obj))
                chromosomeSize <- chromosomeAnnotation[i, 3]
##                if(is.null(xlim)) xlim <- c(0, chromosomeSize)
                if(is.null(xlim)) {
                  xlimit <- range(position(obj), na.rm=TRUE)
                } else xlimit <- xlim

                ##################################################
                ##Plot homozygous calls first
                homozygous <- calls == 1 | calls == 3
                plot(position(obj)[homozygous], cn[homozygous],
                     pch=pch[calls[homozygous]],
                     col=col[calls[homozygous]],
                     bg=bg[calls[homozygous]],
                     cex=cex[calls[homozygous]],
                     cex.axis=cex.axis,
                     cex.main=cex.main,
                     xlab=xlab,
                     ylab=ylab,
                     main=main,
                     xaxs=xaxs,                 
                     yaxs=yaxs,
                     xlim=xlimit,                 
                     ylim=ylim,
                     xaxt="n",
                     yaxt=yaxt,
                     bty=bty,
                     log=log)
                ##################################################
                ##Plot heterozygous calls last
                points(position(obj)[!homozygous], cn[!homozygous],
                       pch=pch[calls[!homozygous]],
                       col=col[calls[!homozygous]],
                       bg=bg[calls[!homozygous]],
                       cex=cex[calls[!homozygous]])                

                ######################################################################
                ##Draw centromere
                ######################################################################
                centromere <- chromosomeAnnotation[i, 1:2]
                if(log == "y"){
                  if(ylim[1] < 0.5) ylim[1] <- 0.5
                }
                rect(xleft=centromere[[1]], ybottom=ylim[1],
                     xright=centromere[[2]], ytop=ylim[2],
                     col=col.centromere[1],
                     border=col.centromere[2])
                
                ######################################################################
                ##Axes
                ######################################################################
                if(i == chromosome.order[1]){
                  axis(2, at=pretty(ylim, lab[2]), labels=pretty(ylim, lab[2]),
                       outer=FALSE, cex.axis=cex.axis, las=1)
                }
                if(xaxt != "n"){
                  if(j == 1 & xaxis.side[i] == 3) mtext(i, 3, line=2.5, cex=cex.lab, outer=FALSE)
                  if(j == S & xaxis.side[i] == 1) mtext(i, 1, line=2.5, cex=cex.lab, outer=FALSE)
                  axis(xaxis.side[i],
                       at=pretty(xlimit, lab[2]),
                       outer=TRUE,
                       labels=pretty(xlimit, lab[2])/1e6,
                       cex.axis=cex.axis,
                       col=col.axis,
                       col.axis=col.axis,
                       las=1, line=0,
                       lwd=1,
                       mgp=c(2, 0.5, 0))
                }
              }
              k <- k+1
              if(addCytoband){
                plotCytoband(obj, cytoBand=cytoband, cex.axis=cex.axis, xaxs=xaxs,
                             xlim=xlimit)
              }
            }
            if(i == chromosome.order[1])
              if(xaxt != "n"){
                mtext("Mb ", xaxis.side[i], line=2, outer=TRUE, cex=cex.lab,
                      col=col.axis, adj=0, las=1)
              }
            ###########################################################################
            ##Plot summary statistics
            ###########################################################################
            if(length(legend) == 1) legend <- rep(legend, 2)
            if(legend[1]){
              showSummary(object, where=legend.location[1], bty=legend.bty[1],
                          legend.panel=legend.panel[1], cex=cex.legend[1], col=col,
                          digits=digits)
            }
            ##plot legend
##          }
          if(legend[2]){
            if(legend.panel[2]) plot(0:1, 0:1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")              
            op <- par(bg=legend.col[2])
            if(any(pch == ".")) pt.cex[pch == "."] <- 6
            legend(legend.location[2],
                   pch=pch[1:2], 
                   col=col[1:2],
                   legend=c("AA/BB", "AB"),
                   cex=cex.legend[2],
                   y.intersp=1.5,
                   bty=legend.bty[2], pt.cex=pt.cex,
                   xjust=1, ncol=ncol)
          }
        })
            

##Courtesy of Jason Ting
setMethod("plotCytoband", "AnnotatedSnpSet",
          function(object, cytoBand=NULL, xlim=NULL, cex.axis=0.8, xaxs="r",
                   chromosome=NULL, main="", outer=TRUE, cytobandAxis=TRUE, ...){
            if(is.null(cytoBand)) {data(cytoband); cytoBand <- cytoband}
            if(is.null(chromosome))  chrom <- unique(chromosome(object))  else chrom <- as.character(chromosome)[1]

            cytoBand_chr <- cytoBand[as.character(cytoBand$chrom) == chrom,]
            cytoBand_chr_p <- cytoBand_chr[grep("^p",as.character(cytoBand_chr$name)),]
            cytoBand_chr_q <- cytoBand_chr[grep("^q",as.character(cytoBand_chr$name)),]
            p.bands <- length(cytoBand_chr_p$chromEnd)
            cut.left  <- c()
            cut.right <- c()
            ##  1st  band of arm or 1st  band after  "stalk"
            ##  last band of arm or last band before "stalk"
            for (band in 1:length(cytoBand_chr$chromEnd)) {
              if (band==1)                             { cut.left[band] <- TRUE; cut.right[band] <- FALSE} else
              if (band==p.bands)                       { cut.left[band] <- FALSE; cut.right[band] <- TRUE} else
              if (band==(p.bands+1))                   { cut.left[band] <- TRUE; cut.right[band] <- FALSE} else
              if (band==length(cytoBand_chr$chromEnd)) { cut.left[band] <- FALSE; cut.right[band] <- TRUE} else{
                cut.left[band] <- FALSE; cut.right[band] <- FALSE
              }
            }
            for (band in 1:length(cytoBand_chr$chromEnd)) {
              if (as.character(cytoBand_chr$gieStain[band])=="stalk") {
                cut.right[band-1] <- TRUE
                cut.left[band] <- NA
                cut.right[band] <- NA
                cut.left[band+1] <- TRUE
              }
            }
            plot(c(0,cytoBand_chr$chromEnd[length(cytoBand_chr$chromEnd)]),
                 c(0, 2), xlim=xlim, type="n", xlab="", ylab="",
                 axes=FALSE, xaxs=xaxs, main=main)
            for (i in 1:length(cytoBand_chr$chromEnd)) {
              start <- cytoBand_chr$chromStart[i]
              end   <- cytoBand_chr$chromEnd[i]
              delta = (end-start)/4
              stain <- as.character(cytoBand_chr$gieStain[i])
              if (stain=="gneg") {
                color <- "grey100"
              } else if (stain=="gpos25") {
                color <- "grey90"
              } else if (stain=="gpos50") {
                color <- "grey70"
              } else if (stain=="gpos75") {
                color <- "grey40"
              } else if (stain=="gpos100") {
                color <- "grey0"
              } else if (stain=="gvar") {
                color <- "grey100"
              } else if (stain=="acen") {
                color <- "brown4"
              } else if (stain=="stalk") {
                color <- "brown3"
              } else {
                color <- "white"
              }
              if (is.na(cut.left[i]) & is.na(cut.right[i])) {
                ## this is a "stalk", do not draw box. Draw two vertival lines instead
                delta <- (end-start)/3
                lines(c(start+delta,start+delta),c(0,2),col=color)
                lines(c(end-delta,end-delta),c(0,2),col=color)
              } else if (cut.left[i] & cut.right[i]) {      # cut both ends
                polygon(c(start,start+delta,end-delta,end,end,end-delta,start+delta,start),
                        c(0.3,0,0,0.3,1.7,2,2,1.7),col=color)
              } else if (cut.left[i]) {              # cut left end only
                polygon(c(start,start+delta,end,end,start+delta,start),
                        c(0.3,0,0,2,2,1.7),col=color)
              } else if (cut.right[i]) {             # cut right end only
                polygon(c(start,end-delta,end,end,end-delta,start),
                        c(0,0,0.3,1.7,2,2),col=color)
              } else {
                polygon(c(start, end, end, start),
                        c(0, 0, 2, 2), col=color)
              }
            }
            my.x <- (cytoBand_chr$chromStart+cytoBand_chr$chromEnd)/2
            if(cytobandAxis){
              axis(1, at=my.x, labels=as.character(cytoBand_chr$name), outer=outer, cex.axis=cex.axis,
                   line=1, las=3)
            }
          })

#create object of class snpscan with smoothed copynumbers and smoothed loh calls
#LOH calls are smoothed by setting homozygous = 0 and heterozygous = 1
setMethod("smoothSnp", "AnnotatedSnpSet",
          function(object, chromosomes,
                   samples,
                   span=1/10,
                   method="loess",
                   imputeNoCalls=TRUE,
                   verbose=TRUE){
            object <- object[chromosome(object) %in% chromosomes, samples]
            ##convert homozygous to 0 and heterozygous to 1
            calls(object)[calls(object) == 1 | calls(object) == 3] <- 0
            calls(object)[calls(object) == 2] <- 1
            smoothChromosome <- function(obj, span){
              loessX <- function(X, location, span){
                fit <- loess(X ~ location, span = span)$fitted
                return(fit)
              }
              ##Order by physical position before smoothing
              obj <- obj[order(position(obj)), ]
              cn.smooth <- apply(copyNumber(obj), 2, loessX, position(obj), span=span)
              call.smooth <- apply(calls(obj), 2, loessX, location=position(obj), span=span)
              rownames(cn.smooth) <- rownames(call.smooth) <- featureNames(obj) 
              copyNumber(obj) <- cn.smooth
              calls(obj) <- call.smooth
              obj
            }
            object.list <- split(object, chromosome(object))
            smooth.list <- lapply(object.list, smoothChromosome, span=span)
            smooth.set <- unsplitS4(smooth.list, featureData(object))
            return(smooth.set)
          })



.calculateYlim <- function(object, op){
  if("copyNumber" %in% ls(assayData(object))){
    ##only print this if there is more than one sample or more than 1 chromosome to plot
    if(length(unique(chromosome(object))) > 1 || ncol(object) > 1){
      print("one.ylim is FALSE. Calculating ylim based on the percentiles of the copy number distribution")
    }
    if(op$log == "y"){
      ylim <- range(copyNumber(object), na.rm=TRUE)
    } else{
      ##use 1 and 99 quantiles to avoid outliers
      ylim <- c(quantile(copyNumber(object), prob=0.001, na.rm=TRUE),
                quantile(copyNumber(object), prob=0.999, na.rm=TRUE))
    }
  } else{
    ##ylimits for genotypes??
    y <- .getY(object)
    ##jitter the genotype calls
    y <- jitter(y, amount=0.05)
    ylim <- range(y)
  }
  ylim
}

centromere <- function(chromosome){
  if(missing(chromosome) | !(chromosome %in% c(1:22, "X", "Y"))) stop("must specify chromosome 1-22, X or Y as a character string")
  data(chromosomeAnnotation, package="SNPchip", envir=environment())
  chromosomeAnnotation[chromosome, c("centromereStart", "centromereEnd")]
}

chromosomeSize <- function(chromosome){
  if(!is.character(chromosome)) stop("argument to chromosomeSize must be one of the following character strings: 1, ..., 22, X, or Y")
  if(any(!(chromosome %in% c(1:22, "X", "Y")))) stop("chromosome must be 1-22, X, or Y")  
  data(chromosomeAnnotation, package="SNPchip", envir=environment())
  chromosomeAnnotation[chromosome, "chromosomeSize"]
}

.drawCentromere <- function(object, op){
  data(chromosomeAnnotation, package="SNPchip", envir=environment())
  centromere <- chromosomeAnnotation[unique(chromosome(object)), ]
  xleft <- centromere["centromereStart"]
  xright <- centromere["centromereEnd"]
  rect(xleft=xleft, ybottom=op$ylim[1],
       xright=xright, ytop=op$ylim[2],
       col=op$col.centromere,
       border=op$border.centromere)  
}

.labelChromosome <- function(object, op, j){
  if(j == 1){
    if(op$label.chromosome)
      mtext(unique(chromosome(object)), side=3, outer=FALSE, line=op$line.label.chromosome, cex=op$cex.lab)
  }              
}

.drawYaxis <- function(object, op, j){
  if(unique(chromosome(object)) != op$firstChromosome) return()
  if(op$yaxt == "n") return()
  if("copyNumber" %in% ls(assayData(object))){
    ##Draw default y-axis
    at <- pretty(op$ylim)
    labels <- at
  }  else {
    at <- c(0, 1)
    labels <- c("AA/BB", "AB")
  }
  axis(side=2, at=at, labels=labels, las=1, cex.axis=op$cex.axis)  
}

.drawXaxis <- function(object, op, j){
  chromosomeName <- unique(chromosome(object))
  if(op$xaxt == "n") return()
  if(op$alternate.xaxis.side){
    side <- op$xaxis.side[[unique(chromosome(object))]]
  } else side <- op$xaxis.side
  if(side == 1 & j == ncol(object) | side == 3 & j == 1){
    axis(side,
         at=pretty(op$xlim[chromosomeName, ], op$lab[2]),
         outer=op$outer.axis,
         labels=pretty(op$xlim[chromosomeName, ], op$lab[2])/1e6,
         cex.axis=op$cex.axis,
         col=op$col.axis,
         col.axis=op$col.axis,
         las=1,
         line=op$line.axis,
         lwd=1,
         mgp=c(2, 0.5, 0))
    if(op$label.chromosome){
      mtext(unique(chromosome(object)), side=side, outer=FALSE, line=op$line.label.chromosome, cex=op$cex.lab)
    }
  }
}

.getCytoband <- function(object, op){
  if(op$add.cytoband){
    data(cytoband)              
    cytoband <- cytoband[cytoband[, "chrom"] == unique(chromosome(object)), ]
  }  else NULL
}

.drawCytobandWrapper <- function(S, cytoband, op, j, chromosomeName){
  if(j == S){
    if(op$add.cytoband){
      if(nrow(cytoband) > 0)
        plotCytoband(cytoband=cytoband,
                     xlim=op$xlim[chromosomeName, ],
                     xaxs=op$xaxs,
                     label.cytoband=op$label.cytoband,
                     cex.axis=op$cex.axis,
                     outer=op$outer.cytoband.axis)
    }
  }
}

.plot <- function(object, op, ...){
  ##could use a switch in the following statement to generate a
  ##class-specific par object
  if(missing(op)) op <- new("ParESet")
  chromosomeName <- unique(chromosome(object))
  
  if(!op$one.ylim){
    op$ylim <- .calculateYlim(object)
  }

  ##ensures homozygous genotypes are plotted first
  object <- .orderByGenotype(object)
  x <- .getX(object)##position
  y <- .getY(object, op)##calls or copy number

  if(!op$outer.ylab) ylab <- op$ylab else ylab <- ""

  ##Option to recycle graphical parameters by genotype call (when available)
  col <- .recycle(op$col, object)
  cex <- .recycle(op$cex, object)
  pch <- .recycle(op$pch, object)
  bg <- .recycle(op$bg, object)
  plot(x=x, y=y,
       xlim=op$xlim[chromosomeName, ],
       ylim=op$ylim,
       col=col,
       cex=cex,
       pch=pch,
       log=op$log,
       bg=bg,
       xaxt="n",
       xaxs=op$xaxs,
       main=op$main,
       ylab=ylab,
       yaxt="n",
       ...)
}

plotCytoband <- function(chromosome,
                         cytoband,
                         xlim,
                         xaxs="r",
                         new=TRUE,
                         label.cytoband=TRUE,
                         cex.axis=1,
                         outer=FALSE,
                         ...){
  def.par <- par(no.readonly=TRUE)
  on.exit(def.par)
  if(missing(cytoband)) data(cytoband, package="SNPchip", envir=environment())
  if(missing(chromosome)){
    if(length(unique(cytoband[, "chrom"])) > 1) stop("Must specify chromosome")
  }
  if(length(unique(cytoband$chrom)) > 1){
    cytoband <- cytoband[cytoband[, "chrom"] == chromosome, ]
  }
  rownames(cytoband) <- as.character(cytoband[, "name"])
  if(missing(xlim)) xlim <- c(0, chromosomeSize(unique(cytoband$chrom)))
  cytoband_p <- cytoband[grep("^p", rownames(cytoband), value=TRUE), ]
  cytoband_q <- cytoband[grep("^q", rownames(cytoband), value=TRUE), ]
  
  p.bands <- nrow(cytoband_p)
  cut.left  <- c()
  cut.right <- c()
  ##  1st  band of arm or 1st  band after  "stalk"
  ##  last band of arm or last band before "stalk"
  for (i in 1:nrow(cytoband)) {
    if (i == 1)                             { cut.left[i] <- TRUE; cut.right[i] <- FALSE} else
    if (i == p.bands)                       { cut.left[i] <- FALSE; cut.right[i] <- TRUE} else
    if (i == (p.bands+1))                   { cut.left[i] <- TRUE; cut.right[i] <- FALSE} else
    if (i == nrow(cytoband)) { cut.left[i] <- FALSE; cut.right[i] <- TRUE} else{
      cut.left[i] <- FALSE; cut.right[i] <- FALSE
    }
  }
  for (i in 1:nrow(cytoband)) {
    if (as.character(cytoband[i, "gieStain"]) == "stalk") {
      cut.right[i-1] <- TRUE
      cut.left[i] <- NA
      cut.right[i] <- NA
      cut.left[i+1] <- TRUE
    }
  }
  ##When plotting subregions of a chromosome, this prevents the cytobands from extending beyond the subsetted object
  ##exclude cytobands that end before the minimum plotting limits
  include <- cytoband[, "chromEnd"] > xlim[1] & cytoband[, "chromStart"] < xlim[2]            
  cytoband <- cytoband[include, ]
  cut.left <- cut.left[include]
  cut.right <- cut.right[include]
  if(new){
    plot(c(0, cytoband[nrow(cytoband), "chromEnd"]),
         c(0, 2),
         xlim=xlim,
         type="n",
         xlab="",
         ylab="",
         axes=FALSE,
         xaxs=xaxs,
         ...)
  }
  for (i in 1:nrow(cytoband)) {
    start <- cytoband[i, "chromStart"]
    last   <- cytoband[i, "chromEnd"]
    delta = (last-start)/4
    getStain <- function(stain){
      switch(stain,
             gneg="grey100",
             gpos25="grey90",
             gpos50="grey70",
             gpos75="grey40",
             gpos100="grey0",
             gvar="grey100",
             stalk="brown3",
             acen="brown4",
             "white")
    }
    color <- getStain(as.character(cytoband[i, "gieStain"]))
    if (is.na(cut.left[i]) & is.na(cut.right[i])) {
      ## this is a "stalk", do not draw box. Draw two vertical lines instead
      delta <- (last-start)/3
      lines(c(start+delta, start+delta), c(0,2), col=color)
      lines(c(last-delta, last-delta), c(0,2), col=color)
    } else if (cut.left[i] & cut.right[i]) {      # cut both lasts
      polygon(c(start, start+delta, last-delta, last, last, last-delta, start+delta, start),
              c(0.3, 0, 0, 0.3, 1.7, 2, 2, 1.7), col=color)
    } else if (cut.left[i]) {              # cut left last only
      polygon(c(start, start+delta, last, last, start+delta, start),
              c(0.3, 0, 0, 2, 2, 1.7), col=color)
    } else if (cut.right[i]) {             # cut right last only
      polygon(c(start, last-delta, last, last, last-delta, start),
              c(0, 0, 0.3, 1.7, 2, 2),col=color)
    } else {
      polygon(c(start, last, last, start),
              c(0, 0, 2, 2), col=color)
    }
  }
  my.x <- (cytoband$chromStart+cytoband$chromEnd)/2
  if(label.cytoband){
    axis(1, at=my.x,
         labels=rownames(cytoband),
         outer=outer,
         cex.axis=cex.axis,
         line=1, las=3)
  }
}

.isHomozygous <- function(object){
  calls(object) == 1 | calls(object) == 3
}

.orderByGenotype <- function(object){
  if(!("calls" %in% ls(assayData(object)))) return(object)
  gt <- as.vector(calls(object))
  gt[gt == 3] <- 1
  object[order(gt, decreasing=FALSE), ]
}
     
  

.recycle <- function(x, object){
  if(length(x) == nrow(object)) return(x)
  ##assume using 2 colors for homozygous and hets
  if(length(x) > 1){
    if("calls" %in% ls(assayData(object))){
      gt <- as.vector(calls(object))
      gt[is.na(gt)] <- 4
      ##assume that colors are to be recycled
      if(max(gt) > length(x)){
        warning("The length of the vector to be recycled is less than the number of genotype calls")
        print("Adding 'grey40' as a color for missing genotype.  See getPar for changing the defaults colors.")
        x <- c(x, "grey40")
      }
      if(max(gt) < length(x)){
        x <- x[sort(unique(gt))]
      }
      x <- x[gt]
    }
  }
  x
}

showSummary <- function(object, where, bty, legend.panel, cex, col, digits){
  f <- function(x, where, bty, legend.panel, cex, col){
    par(mar=rep(0,4))
    if(legend.panel) plot(0:1, 0:1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    legend(where, legend=c(
                    paste(x["ho"], " %AA/BB", sep=""),
                    paste(x["ht"], " %AB", sep=""),
                    paste(x["cn"], " avg CN"),
                    paste(x["cn.sd"], " sd")), bty=bty,
           title=substr(x["samplenames"], 1, 10),
           y.intersp=1.1,
           cex=cex,
           text.col=c("black", col[1], col[2], "black", "black"))
  }
  
  if(sum(chromosome(object) != "X") > 0){
    obj <- object[chromosome(object) != "X", ]
  } else { obj <- object}
  cn <- colMeans(as.matrix(copyNumber(obj)))
  cn.sd <- apply(copyNumber(obj), 2, sd)                
  ht <- colMeans(ifelse(calls(obj) == 2, 1, 0))
  ho <- colMeans(ifelse(calls(obj) == 1 | calls(obj) == 3, 1, 0))
  stats <- cbind(cn, cn.sd, ht, ho)
  stats <- round(stats, digits)
  stats <- data.frame(stats); stats$samplenames <- sampleNames(object)
  apply(stats, 1, f, where=where, bty=bty,
        legend.panel=legend.panel, cex=cex, col=col)
}








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

chromosome2numeric <- function(chromosome){
  chrom <- as.character(chromosome)
  chrom[chrom == "X"] <- 23
  chrom[chrom == "XY"] <- 24
  chrom[chrom == "Y"] <- 25
  chrom[chrom == "M"] <- 26
  chrom <- as.numeric(chrom)
  chrom
}

chromosomeSize <- function(chromosome){
  if(!is.character(chromosome)) stop("argument to chromosomeSize must be one of the following character strings: 1, ..., 22, X, or Y")
  if(any(!(chromosome %in% c(1:22, "X", "Y", "XY", "M")))) stop("chromosome must be 1-22, X, or Y")  
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



.isHomozygous <- function(object){
  calls(object) == 1 | calls(object) == 3
}

.orderByGenotype <- function(object){
  if(!("calls" %in% ls(assayData(object)))) return(object)
  gt <- as.vector(calls(object))
  gt[gt == 3] <- 1
  object[order(gt, decreasing=FALSE), ]
}
     
  

.recycle <- function(x, object, missing){
  if(length(x) == nrow(object)) return(x)
  ##assume using 2 colors for homozygous and hets
  if(length(x) > 1){
    if("calls" %in% ls(assayData(object))){
      gt <- as.vector(calls(object))
      gt[is.na(gt)] <- 4

      if(sum(!(gt %in% 1:4)) > 0){
        warning("Changing all genotypes that are not 1, 2, 3 to the value 4")
        gt[!(gt %in% 1:4)] <- 4
      }
      ##assume that colors are to be recycled
      if(max(gt) > length(x)){
        x <- c(x, missing)
      }
      if(max(gt) <= length(x)){
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








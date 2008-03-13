setMethod(".plotChromosome", "SnpLevelSet",
          function(object, op, hmmPredict, ...){

		  .getCytoband <- function(object, op){
			  if(op$add.cytoband){
				  data(cytoband)              
				  cytoband <- cytoband[cytoband[, "chrom"] == unique(chromosome(object)), ]
			  }  else NULL
		  }		  
		  cytoband <- .getCytoband(object, op)


		  .drawCytobandWrapper <- function(S, cytoband, op, j, chromosomeName){
			  if(!op$add.cytoband) return()
			  if(nrow(cytoband) > 0)
				  plotCytoband(cytoband=cytoband,
					       xlim=op$xlim[chromosomeName, ],
					       xaxs=op$xaxs,
					       label.cytoband=op$label.cytoband,
					       cex.axis=op$cex.axis,
					       outer=op$outer.cytoband.axis)
		  }

		  .drawYaxis <- function(object, op, j){
			  if(unique(chromosome(object)) != op$firstChromosome) return()
			  if(op$yaxt == "n") return()
			  if("copyNumber" %in% ls(assayData(object))){
				  at <- pretty(op$ylim)
				  at <- c(op$at, at)
				  labels <- at
			  }  else {
				  at <- c(0, 1)
				  labels <- c("AA/BB", "AB")
			  }
			  axis(side=2, at=at, labels=labels, las=1, cex.axis=op$cex.axis)  
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

		  if(!missing(hmmPredict)){
			  if(is.null(op$col.predict)){
				  require(RColorBrewer, quietly=TRUE) || stop("RColorBrewer package not available")
				  col.predict <- brewer.pal(length(states(hmmPredict)), "BrBG")
				  col.predict[states(hmmPredict) == "N"] <- "white"
				  print("col.predict not specified in list of graphical parameters. Using the following colors:")
				  print(col.predict)
				  op$col.predict <- col.predict
				  op$legend.fill.predict <- col.predict
			  }
			  if(is.null(op$height.predict)){
				  op$height.predict <- 0.2
			  }

			  if(!is.null(op$abline.v)){
				  if(length(sampleNames(hmmPredict)) == 1){
					  v <- breakpoints(hmmPredict)[, c("state", "start", "last")]
					  v <- v[v$state != "N", ]
					  v <- c(v$start, v$last)*1e6
					  op$abline.v <- v
				  }
			  }
		  }
		  
			  
		  if(op$cytoband.side == 3) cytobandOnTop <- TRUE else cytobandOnTop <- FALSE

##		  browser()
		  if(cytobandOnTop){
			  .drawCytobandWrapper(S=ncol(object),
					       cytoband=cytoband,
					       op=op,
					       j=j,
					       chromosomeName=unique(chromosome(object)))
		  }		  
		  for(j in 1:ncol(object)){
			  op$main <- op$main[j]
			  .plot(object[, j], op=op)
			  .drawYaxis(object=object, op=op)
			  if(!missing(hmmPredict)){
				  x <- split(hmmPredict[, j], chromosome(hmmPredict))
				  lapply(x, plotPredictions, op=op)
			  }
			  if(op$add.centromere)
				  .drawCentromere(object[, j], op)
			  if(cytobandOnTop)
				  .drawXaxis(object=object, op=op, j=j)
		  }
		  if(!cytobandOnTop){
			  .drawCytobandWrapper(S=ncol(object),
					       cytoband=cytoband,
					       op=op,
					       j=j,
					       chromosomeName=unique(chromosome(object)))
			  .drawXaxis(object=object, op=op, j=j)		  			  
		  }
		  return(op)
          })


.drawXaxis <- function(object, op, j){
	chromosomeName <- unique(chromosome(object))
	if(op$xaxt == "n") return()
	if(op$alternate.xaxis.side){
		side <- op$xaxis.side[[unique(chromosome(object))]]
	} else side <- op$xaxis.side
	if(side == 1 & j == ncol(object) | side == 3 & j == 1){
		##if labeling the cytoband, force the x-axis to be drawn on top
		if(op$cytoband.side==1 & op$label.cytoband){
##			op$outer.axis <- TRUE
##			op$line.axis <- 2
			side <- 3
		}
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
			if(!op$abbreviateChromosomeNames){
				chr <- paste("Chr", unique(chromosome(object)))
			} else{
				chr <- unique(chromosome(object))
			}
			mtext(chr,
			      side=side,
			      outer=FALSE,
			      line=op$line.label.chromosome,
			      cex=op$cex.lab)
		}
	}
}

setMethod("plotSnp", c("ParESet", "SnpLevelSet"),
          function(object, snpset, ...){

		  if(!missing(...)){
			  require(VanillaICE) || stop("VanillaICE package not available")
		  }
		  snpset <- snpset[!is.na(chromosome(snpset)), ]
		  if(object$useLayout){
			  layout(mat=object$mat,
				 widths=object$widths,
				 heights=object$heights,
				 respect=object$respect)
		  }
		  snpList <- split(snpset, chromosome(snpset))
		  names(snpList)[names(snpList) == "X"] <- "23"
		  names(snpList)[names(snpList) == "XY"] <- "24"            
		  names(snpList)[names(snpList) == "Y"] <- "25"
		  names(snpList)[names(snpList) == "M"] <- "26"            
		  snpList <- snpList[order(as.numeric(names(snpList)))]
		  names(snpList)[names(snpList) == "23"] <- "X"
		  names(snpList)[names(snpList) == "24"] <- "XY"            
		  names(snpList)[names(snpList) == "25"] <- "Y"
		  names(snpList)[names(snpList) == "26"] <- "M"

		  if(length(snpList) > 10) object$abbreviateChromosomeNames <- TRUE else object$abbreviateChromosomeNames <- FALSE

		  if(object$ylab == "copy number"){
			  if(any(apply(copyNumber(snpset), 2, "median") > 3) | any(apply(copyNumber(snpset), 2, "median") < 0)){
				  warning("The default ylabel 'copy number' may not be consistent with the quantity plotted on the vertical axes.  Typically, the median copy number is approximately 2 for autosomes or 1 for the male chromosome X")
			  }
		  }

		  par(allPlots(object))
		  for(i in 1:length(snpList)){
			  if(i == 1) par(yaxt="s") else par(yaxt="n")
			  object <- .plotChromosome(snpList[[i]], op=object, ...)
		  }
		  if(object$outer.ylab) mtext(object$ylab, side=object$side.ylab, outer=TRUE, las=3, cex=object$cex.ylab, line=object$line.ylab)
		  mtext(object$xlab, side=object$side.xlab, outer=object$outer.xlab, cex=object$cex.xlab, line=object$line.xlab, adj=0)
		  
		  if(length(sampleNames(snpset) == 1)){
			  mtext(object$main, side=3, outer=TRUE, cex=1.4)
		  }
		  return(object)
          })


setMethod("plotSnp", c("ParSnpCopyNumberSet", "SnpCopyNumberSet"),
          function(object, snpset, ...){
            old.par <- par(no.readonly=TRUE)
            on.exit(old.par)
            callNextMethod()
          })

setMethod("plotSnp", c("ParSnpCallSet", "SnpCallSet"),
          function(object, snpset, ...){
            old.par <- par(no.readonly=TRUE)
            on.exit(old.par)            
            callNextMethod()
          })

setMethod("plotSnp", c("ParSnpSet", "oligoSnpSet"),
          function(object, snpset, ...){
            old.par <- par(no.readonly=TRUE)
            on.exit(old.par)
            callNextMethod()
          })

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

.plot <- function(object, op, ...){
  ##could use a switch in the following statement to generate a
  ##class-specific par object
  if(missing(op)) op <- new("ParESet")
  chromosomeName <- unique(chromosome(object))

  if(!op$one.ylim){
    op$ylim <- .calculateYlim(object)
  }

  .orderByGenotype <- function(object){
	  if(!("calls" %in% ls(assayData(object)))) return(object)
	  gt <- as.vector(calls(object))
	  gt[gt == 3] <- 1
	  object[order(gt, decreasing=FALSE), ]
  }

  
  
  ##ensures homozygous genotypes are plotted first
  object <- .orderByGenotype(object)
  x <- .getX(object)##position
  y <- .getY(object, op)##calls or copy number

  if(!op$outer.ylab) ylab <- op$ylab else ylab <- ""


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
  
  ##Option to recycle graphical parameters by genotype call (when available)
  col <- .recycle(op$col, object, missing="grey40")
  cex <- .recycle(op$cex, object, missing=1)
  pch <- .recycle(op$pch, object, missing=".")
  bg <- .recycle(op$bg, object, missing="grey40")

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
  if(op$abline){
	  abline(h=op$abline.h, col=op$abline.col, lty=op$abline.lty, lwd=op$abline.lwd)
  }
  if(!is.null(op$abline.v)){
	  abline(v=op$abline.v, col=op$abline.v.col, lty=op$abline.v.lty, lwd=op$abline.v.lwd)
  }
  if(!is.null(op$legend)){
##	  ix <- order(unique(op$col))
	  legend(op$legend.location,
##		 fill=op$legend.fill,
		 col=op$legend.col,
		 pt.bg=op$legend.bg,
		 pch=op$legend.pch,
		 legend=op$legend,
		 bty=op$legend.bty)
  }
  return()
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
		xx <- c(0, cytoband[nrow(cytoband), "chromEnd"])
		yy <- c(0, 2)
		plot(xx,
		     yy,
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
	return()
}


plotPredictions <- function(object, op){
	breakpoints <- breakpoints(object)
	position <- position(object)
	chromosome <- chromosome(object)
	states <- states(object)
	predictions <- predictions(object)


	chr <- unique(chromosome)
	for(i in chr){
	  tmp <- breakpoints[breakpoints[, "chr"] == i, , drop=FALSE]
	  .drawRect <- function(x, position, op){
		  col <- op$col.predict
		  if(length(x) < 1) return()
		  start <- max(as.numeric(x["start"]) * 1e6, op$xlim[1])
		  last <- min(as.numeric(x["last"]) * 1e6, op$xlim[2])
		  predict <- predictions[position >= start & position <= last]
		  predict <- predict[!is.na(predict)]
		  if(length(unique(predict)) > 1) {
			  browser()
			  stop("predictions not unique")
		  }
		  col <- col[unique(predict)]
		  rect(xleft=start,
		       ybottom=op$ylim[1],
		       xright=last,
		       ytop=op$height.predict+op$ylim[1],
		       col=col,
		       border=col)
	  }
	  if(!is.null(tmp)){
		  ##best to draw the biggest regions first and the smallest regions last.
		  tmp <- tmp[order(tmp[, "size"], decreasing=TRUE), ]
		  apply(tmp, 1, .drawRect, position=position, op=op)
	  }
	  if(!is.null(op$legend.predict)){
		  legend(op$legend.location.predict,
			 fill=op$col.predict,
			 legend=op$legend.predict,
			 ncol=length(op$col.predict),
			 bty=op$legend.bty,
			 cex=op$cex.legend)
	  }
  }
	return(op)
}

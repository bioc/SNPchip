setMethod("initialize", "ParESet",
          function(.Object,
                   layout=TRUE,
                   col.axis="brown",
                   cex.main=1,
                   cex.axis=1,
                   cex.legend=1,
                   cex=2,
                   cex.lab=1,
                   pch=".",
                   col="black",
                   bg="white",
                   xaxs="r",
                   xaxt="s",
                   yaxs="r",
                   yaxt="s",
		   at=1:4,
                   lab=c(2, 5, 7), ##see par
                   adj=0,
                   bty="n",
                   ann=FALSE,
                   useLayout=TRUE,
                   mar=c(0.5, 0, 0.5, 0.2),
                   oma=c(4, 4, 4, 0.5),
                   las=1,
                   log="",
                   ylab="",
                   side.ylab=2,
                   outer.ylab=TRUE,
                   line.ylab=3,
                   cex.ylab=1,
                   xlab="Mb",
                   outer.xlab=TRUE,
                   side.xlab=1,
                   cex.xlab=1,
                   line.xlab=3,
                   outer.axis=TRUE,
                   line.axis=0,
                   main="",
                   col.centromere="bisque",
                   border.centromere="bisque",
                   xlim=NULL,
                   ylim=NULL,
                   one.ylim=TRUE,
                   add.cytoband=TRUE,
                   outer.cytoband=FALSE,
                   outer.cytoband.axis=FALSE,                   
                   label.cytoband=FALSE,
                   use.chromosome.size=FALSE, #for x-axis limits
                   label.chromosome=TRUE,
                   line.label.chromosome=2,
                   xaxis.side=1,
                   alternate.xaxis.side=FALSE,
                   mat=new("matrix", 1, 1),
                   heights=1,
                   widths=1,
                   respect=FALSE,
                   firstChromosome="1",
                   abline=FALSE,
                   abline.h=2,
                   abline.col="grey80",
                   abline.lty=1,
                   abline.lwd=1,
		   abline.v=NULL,
		   abline.v.col="grey20",
		   abline.v.lty=2,
		   abline.v.lwd=0.8,
                   ...){
            .Object@snpPar <- list(col.axis=col.axis,
                                   cex.main=cex.main,
                                   cex.axis=cex.axis,
                                   cex.legend=cex.legend,
                                   cex.lab=cex.lab,
                                   bty=bty,
                                   ann=ann,
                                   oma=oma,
                                   mar=mar,
                                   las=las,
                                   cex=cex,
                                   pch=pch,
                                   col=col,
                                   bg=bg,
                                   xaxs=xaxs,
                                   xaxt=xaxt,
                                   yaxt=yaxt,
				   at=at,
                                   yaxs=yaxs,
                                   lab=lab,
                                   adj=adj,
                                   log=log,
                                   xlab=xlab,
                                   side.xlab=side.xlab,
                                   outer.xlab=outer.xlab,
                                   cex.xlab=cex.xlab,
                                   line.xlab=line.xlab,
                                   ylab=ylab,
                                   side.ylab=side.ylab,
                                   outer.ylab=outer.ylab,
                                   cex.ylab=cex.ylab,
                                   line.ylab=line.ylab,
                                   main=main,
                                   xlim=xlim,
                                   ylim=ylim,
                                   col.centromere=col.centromere,
                                   border.centromere=border.centromere,
                                   one.ylim=one.ylim,
                                   add.cytoband=add.cytoband,
                                   outer.cytoband=outer.cytoband,
                                   outer.cytoband.axis=outer.cytoband.axis,
                                   label.cytoband=label.cytoband,
                                   use.chromosome.size=use.chromosome.size,
                                   label.chromosome=label.chromosome,
                                   line.label.chromosome=line.label.chromosome,
                                   xaxis.side=xaxis.side,
                                   alternate.xaxis.side=alternate.xaxis.side,
                                   ##Arguments to layout()
                                   useLayout=useLayout,
                                   mat=mat,
                                   heights=heights,
                                   widths=widths,
                                   respect=respect,
                                   firstChromosome=firstChromosome,
                                   abline=abline,
                                   abline.h=abline.h,
                                   abline.col=abline.col,
                                   abline.lty=abline.lty,
                                   abline.lwd=abline.lwd,				   
				   abline.v=abline.v,
				   abline.v.col=abline.v.col,
				   abline.v.lwd=abline.v.lwd,
				   abline.v.lty=abline.v.lty)
            .Object
          })

##updates graphical parameters with information from the data class
setMethod("getPar", c("ParESet", "SnpLevelSet"),
          function(object, snpset, add.cytoband, ...){
  snpset <- snpset[!is.na(chromosome(snpset)), ]
  ##layout
  chromosomeNames <- as.character(sort(chromosome2numeric(chromosome(snpset))))
  chromosomeNames[chromosomeNames == "23"] <- "X"
  chromosomeNames[chromosomeNames == "24"] <- "XY"
  chromosomeNames[chromosomeNames == "25"] <- "Y"
  chromosomeNames[chromosomeNames == "26"] <- "M"
  chromosomeNames <- unique(chromosomeNames)

  N <- length(chromosomeNames)
  S <- ncol(snpset)  
  data(chromosomeAnnotation, package="SNPchip", envir=environment())
  object$heights <- rep(1, ncol(snpset))
  if(!missing(add.cytoband)) object$add.cytoband <- add.cytoband
  if(object$add.cytoband){
    if(!object$outer.cytoband){
      S <- S+1
      if(S > 3){
        object$heights <- c(object$heights, 1/(ncol(snpset)*1.5))
      }
      if(S == 2 | S == 3){
        object$heights <- c(object$heights, 1/10)
      }
    } 
  }
  if(N > 10){
    object$alternate.xaxis.side <- TRUE
    side <- c(1, 3)[rep(1:2, N/2 + 1)]
    side <- side[1:N]
    options(warn=-1)
    names(side) <- chromosomeNames
    object$xaxis.side <- side
  }
  m <- matrix(1:(S*N), nc=N, byrow=FALSE)
  w <- chromosomeAnnotation[chromosomeNames, "chromosomeSize"]
  object$widths <- w/min(w)
  object$mat <- m

  ######################################################################
  ##graphical parameters
  if("copyNumber" %in% ls(assayData(snpset))){
    if(min(copyNumber(snpset), na.rm=TRUE) > 0) object$log <- "y" else object$log <- ""
    ##by default, we use the same ylimit on all the plots
    
    ##could make plot specific by adding an option one.ylim (or
    ##something to that effect) and calculating ylim in .plot()
    object$ylim <- .calculateYlim(snpset, object)
  } 
  def.op <- options(warn=-1)
  object$firstChromosome <- chromosomeNames[1]
  options(def.op)

  if(object$use.chromosome.size){
    object$xlim <- matrix(NA, nrow=length(chromosomeNames), ncol=2)    
    object$xlim[, 1] <- rep(0, nrow(object$xlim))
    object$xlim[, 2] <- chromosomeSize(chromosomeNames)
    rownames(object$xlim) <- chromosomeNames    
  } else{
    objList <- split(snpset, chromosome(snpset))
    objList <- objList[chromosomeNames]
    object$xlim <- t(sapply(objList, function(snpset) range(position(snpset))))
  }
  object
  ##set up defaults according to number of samples, chromosomes, position, etc.
})

setMethod("$", "ParESet", function(x, name){
  eval(substitute(snpPar(x)$NAME_ARG, list(NAME_ARG=name)))
})

setReplaceMethod("$", "ParESet",
                 function(x, name, value) {
                   snpPar(x)[[name]] = value
                   x
})

setMethod("snpPar", "ParESet", function(object) object@snpPar)

setReplaceMethod("snpPar", "ParESet", function(object, value) {
  object@snpPar <- value
  object
})

setMethod("allPlots", "ParESet", function(object){
  list(col.axis=object$col.axis,
       cex.main=object$cex.main,
       cex.lab=object$cex.lab,
       bty=object$bty,
       ann=object$ann,
       oma=object$oma,
       mar=object$mar,
       las=object$las,
       lab=object$lab)  
  })




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
                                   firstChromosome=firstChromosome)
            .Object
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


##setMethod("plotSpecific", "ParESet", function(object){
##  list(cex=object$cex,
##       pch=object$pch,
##       col=object$col,
##       bg=object$bg,
##       xaxs=object$xaxs,
##       xaxt=object$xaxt,
##       yaxt=object$yaxt,
##       lab=object$lab,
##       adj=object$adj)
##})

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


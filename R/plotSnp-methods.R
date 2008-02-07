setMethod(".plotChromosome", "SnpLevelSet",
          function(object, op, ...){
            cytoband <- .getCytoband(object, op)
            for(j in 1:ncol(object)){
              op$main <- op$main[j]
              .plot(object[, j], op=op)
              .drawYaxis(object=object, op=op)
              .drawCentromere(object[, j], op)              
              .drawCytobandWrapper(S=ncol(object), cytoband=cytoband, op=op, j=j, chromosomeName=unique(chromosome(object)))
              .drawXaxis(object=object, op=op, j=j)
            }
          })

setMethod("plotSnp", c("ParESet", "SnpLevelSet"),
          function(object, snpset){
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


            par(allPlots(object))
            for(i in 1:length(snpList)){
              if(i == 1) par(yaxt="s") else par(yaxt="n")
              .plotChromosome(snpList[[i]], op=object)
            }
            if(object$outer.ylab) mtext(object$ylab, side=object$side.ylab, outer=TRUE, las=3, cex=object$cex.ylab, line=object$line.ylab)
            mtext(object$xlab, side=object$side.xlab, outer=object$outer.xlab, cex=object$cex.xlab, line=object$line.xlab)
          })


setMethod("plotSnp", c("ParSnpCopyNumberSet", "SnpCopyNumberSet"),
          function(object, snpset){
            old.par <- par(no.readonly=TRUE)
            on.exit(old.par)
            callNextMethod()
          })

setMethod("plotSnp", c("ParSnpCallSet", "SnpCallSet"),
          function(object, snpset){
            old.par <- par(no.readonly=TRUE)
            on.exit(old.par)            
            callNextMethod()
          })

setMethod("plotSnp", c("ParSnpSet", "oligoSnpSet"),
          function(object, snpset){
            old.par <- par(no.readonly=TRUE)
            on.exit(old.par)            
            callNextMethod()
          })

setMethod("show", "ParESet", function(object) str(snpPar(object)))

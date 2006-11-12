setMethod("plotChromosome", "AnnotatedSnpSet",
          function(object, 
                   col="black",
                   colAA="blue",
                   colAB="red",
                   colNC="green3",
                   colCentromere="bisque",
                   centromereBorder=NA,
                   pch=".",
                   cex=1,
                   cexAA=1,
                   cexAB=1,
                   cexNC=1,
                   bty="o",
                   legend.bty="n",
                   xlim=NULL,
                   ylim=NULL,                   
                   cex.axis=0.8,
                   cex.legend=1,
                   cex.main=1,
                   xlab="Mb",
                   ylab="copy number",
                   ps=16,
                   digits=3,
                   lwdCn=2,
                   legendStats=TRUE,
                   legendPch=TRUE,
                   legend.loc="topleft",
                   xaxt="s",
                   yaxt="s",
                   xaxs="i",                   
                   yaxs="r",
                   main="",
                   mar=c(4, 4, 1, 0.2),
                   panel.xaxis=FALSE,
                   panel.yaxis=FALSE,
                   yTicks=5,
                   xTicks=5,
                   log="",
                   bw=FALSE, ...){
            if(bw){
              colAA <- gray(0.6)
              colAB <- "black"
              colNC <- gray(0.8)
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
            plot(x, y,
                 pch=pch,
                 cex.axis=cex.axis,
                 cex.main=cex.main,
                 xlab=xlab,
                 ylab=ylab,
                 main=main,
                 xaxs=xaxs,                 
                 yaxs=yaxs,
                 xlim=xlim,                 
                 ylim=ylimit,
                 xaxt="n",
                 yaxt=yaxt,
                 type="n",
                 bty=bty,
                 log=log)
            if(panel.xaxis) axis(side=1, at=pretty(xlim, n=xTicks),
                                 labels=as.character(round(pretty(xlim, n=xTicks)/1e6,0)),
                                 cex.axis=cex.axis)
            if(panel.yaxis) axis(side=2, at=pretty(ylim, n=yTicks), outer=TRUE,
                                 cex.axis=cex.axis)
            calls <- calls(object)
            hom <- calls == 1 | calls == 3
            het <- calls == 2
            points(jitter(x[hom], 1), jitter(y[hom], 1),
                   col = colAA,
                   cex = cexAA,
                   pch = pch)
            points(jitter(x[het],1), jitter(y[het], 1),
                   col = colAB,
                   cex = cexAB,
                   pch = pch)
            if(any(calls==4)){
              points(x[calls == 4], y[calls == 4, 1],
                     col = colNC,
                     cex = cexNC,
                     pch = pch)
            }
            ######################################################################
            ##Draw centromere
            ######################################################################
            centromere <- chrAnn[chrom, 1:2]
            rect(xleft=centromere[[1]], ybottom=ylim[1],
                 xright=centromere[[2]], ytop=ylim[2], col=colCentromere,
                 centromereBorder)
            cn <- mean(copyNumber(object))
            cn.sd <- sd(copyNumber(object))
            ht <- mean(ifelse(calls(object) == 2, 1, 0))
            ho <- mean(ifelse(calls(object) == 1 | calls(object) == 3, 1, 0))
            stats <- c(cn, cn.sd, ht, ho)
            stats <- round(stats, digits)
            names(stats) <- c("cn", "cn.sd", "ht", "ho")
            if(legendStats){
              par(bg="antiquewhite1")
              legend(legend.loc, legend = c(substr(x["samplenames"], 1, min(nchar(x["samplenames"]),10)),
                                  paste(stats["ho"], " %AA/BB", sep = ""),
                                  paste(stats["ht"], " %AB", sep = ""),
                                  paste(stats["cn"], " avg CN"),
                                  paste(stats["cn.sd"], " sd")), bty = legend.bty, cex = cex.legend,
                     text.col = c("black", colAA, colAB, "black", "black"))
            }
            if(legendPch){
              legend("topright",
                     pch=20,
                     col=c(colAA, colAB),
                     legend=c("AA/BB", "AB"),
                     cex=cex.legend,
                     bty="n",
                     pt.cex=cex.legend*1.5)
              par(bg="white")
            }
          })

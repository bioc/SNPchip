setMethod("plotSnp", "AnnotatedSnpSet",
          function(object, chromosomes, samples,
                   ylim=NULL,
                   xlim=NULL,
                   col="black",
                   colAA="blue",
                   colAB="red",
                   colNC="green3",
                   colCentromere="bisque",
                   col.axis="brown",
                   pch=".",
                   cex=1,
                   cexAA=1,
                   cexAB=1,
                   cexNC=1,
                   cex.axis=0.8,
                   cex.legend=1,
                   cex.chr=0.8,
                   oma=c(5, 3, 4, 0.5),
                   mar=c(0, 0, 0, 0.2),
                   width.right=NULL,
                   summaryPanel=FALSE,
                   showLayout=TRUE,
                   plotIt=TRUE,
                   digits=3,
                   legend=TRUE,
                   legend.location="topleft",
                   legend.bty="o",
                   legend.col="white",
                   alternate.xaxis=TRUE,
                   xaxis=TRUE,
                   jitter=TRUE,
                   factor=0.1,
                   yTicks=5,
                   xTicks=2,
                   bty="n",
                   bw=FALSE,
                   ...){

            ##If black and white, change the default color scheme
            if(bw){
              colAA <- gray(0.6)
              colAB <- "black"
              colNC <- gray(0.8)
              cex.axis <- "black"
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
            object <- as(object, "AnnotatedSnpSetList")            
            object <- SNPchip::snpSetList(object)
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
                plotChromosome(object[[chrom]][, i], colAA=colAA, colAB=colAB, colNC=colNC,
                               cexAA=cexAA, cexAB=cexAB, cexNC=cexNC,
                               ylim=ylim, xlim=xlim,
                               yaxt="n",
                               xaxt="n", 
                               yaxs="r",
                               legend=FALSE,
                               centromereBorder=NA,
                               xlab="", ylab="",
                               bty=bty,
                               panel.xaxis=FALSE,
                               panel.yaxis=yaxis,
                               yTicks=yTicks,
                               mar=mar)
                if(i == 1){
                  if(chrom == chromosomes[2]){
                    label <- chrom
                  } else {
                    label <- strsplit(chrom, "chr")[[1]][2]
                  }
                  if(side.last == 1) mtext(label, 3, line=2.5, cex=cex.chr)
                }
                if(i == S) {
                  if(length(chromosomes) <= 5){
                    probs <- seq(0, 1, by=1/(xTicks+2))
                    probs <- probs[-c(1, length(probs))]
                  } else {probs <- c(0, 1)}
                  quants <- quantile(xlim, probs)
                  if(side.last == 1) side <- 3 else side <- 1
                  axis(side, at=quants, outer=TRUE,
                       labels=as.character(round(quants/1e6, 0)),
                       cex.axis=cex.axis, col.axis=col.axis)
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
            mtext("Mb ", 1, at=0, line=1, outer=TRUE, cex=(cex.chr)*1.05, col=col.axis, adj=1)

            ###########################################################################
            ##Plot summary statistics
            object <- obj
            if(summaryPanel){
              if(sum(chromosome(object) != "chrX") > 0){
                obj <- object[chromosome(object) != "chrX", ]
              } else { obj <- object}
              cn <- colMeans(as.matrix(copyNumber(obj)))
              cn.sd <- apply(copyNumber(obj), 2, stats::sd)                
              ht <- colMeans(ifelse(calls(obj) == 2, 1, 0))
              ho <- colMeans(ifelse(calls(obj) == 1 | calls(obj) == 3, 1, 0))
              stats <- cbind(cn, cn.sd, ht, ho)
              stats <- round(stats, digits)
            }            

            if(summaryPanel){
              showSummary <- function(x){
                par(mar=rep(0,4))
                plot(0:1, 0:1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
                legend("topleft", legend=c(substr(x["samplenames"], 1, 10),
                                    paste(x["ho"], " %AA/BB", sep=""),
                                    paste(x["ht"], " %AB", sep=""),
                                    paste(x["cn"], " avg CN"),
                                    paste(x["cn.sd"], " sd")), bty="n",
                       cex=cex.legend,
                       text.col=c("black", colAA, colAB, "black", "black"))
              }
              stats <- data.frame(stats); stats$samplenames <- samplenames
              apply(stats, 1, showSummary)
            }
            ##plot legend
            if(legend){
              op <- par(bg=legend.col)
              if(summaryPanel){
                if(legend.location == "topleft") legend.location <- "bottomleft"
                legend.bty <- "n"
              }
              legend(legend.location, pch=20, 
                     col=c(colAA, colAB),
                     legend=c("AA/BB", "AB"), cex=cex.legend,
                     bty=legend.bty, pt.cex=cex.legend*1.5)
            }
          })

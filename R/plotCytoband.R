##Courtesy of Jason Ting
setMethod("plotCytoband", "AnnotatedSnpSet",
          function(object, cytoBand, mar=c(0, 0, 0, 0), xlim=NULL,
                   cex.axis=0.8, ...){
            chrom <- unique(SNPchip::chromosome(object))            
            cytoBand_chr <- cytoBand[as.character(cytoBand$chrom) == chrom,]
            cytoBand_chr_p <- cytoBand_chr[grep("^p",as.character(cytoBand_chr$name)),]
            cytoBand_chr_q <- cytoBand_chr[grep("^q",as.character(cytoBand_chr$name)),]
            p.bands <- length(cytoBand_chr_p$chromEnd)
            cut.left  <- c()
            cut.right <- c()
            ##  1st  band of arm or 1st  band after  "stalk"
            ##  last band of arm or last band before "stalk"
            for (band in 1:length(cytoBand_chr$chromEnd)) {
              if (band==1)                             { cut.left[band] = T; cut.right[band] = F } else
              if (band==p.bands)                       { cut.left[band] = F; cut.right[band] = T } else
              if (band==(p.bands+1))                   { cut.left[band] = T; cut.right[band] = F } else
              if (band==length(cytoBand_chr$chromEnd)) { cut.left[band] = F; cut.right[band] = T } else{
                cut.left[band] = F; cut.right[band] = F
              }
            }
            for (band in 1:length(cytoBand_chr$chromEnd)) {
              if (as.character(cytoBand_chr$gieStain[band])=="stalk") {
                cut.right[band-1] = T
                cut.left[band]    = NA
                cut.right[band]   = NA
                cut.left[band+1]  = T
              }
            }
            par(mar=mar)
            plot(c(0,cytoBand_chr$chromEnd[length(cytoBand_chr$chromEnd)]), c(0, 2),
                 xlim=xlim, type="n", xlab="", ylab="", axes=F, xaxs="i")
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
                ## this is a "stalk", do not drow box. Drow two vertival lines instead
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
            ##Figure margins too big in the next plot
            ##  par(mai=c(0.05,0.5,0.05,0.1), ps=6)
#            plot(c(0, cytoBand_chr$chromEnd[length(cytoBand_chr$chromEnd)]), c(0, 2),
#                 xlim=xlim, type="n", xlab="", ylab="", axes=F, xaxs="i")
#            browser()
            par(las=3)  # rotate text to vertical
            my.x <- (cytoBand_chr$chromStart+cytoBand_chr$chromEnd)/2            
            axis(1, at=my.x, labels=as.character(cytoBand_chr$name), outer=TRUE, cex.axis=cex.axis,
                 line=1)
#            for (i in 1:length(cytoBand_chr$name)) {
#             my.x <- (cytoBand_chr$chromStart[i]+cytoBand_chr$chromEnd[i])/2
##              text(x=my.x, y=1.7, labels=cytoBand_chr$name[i], adj=c(1,0.5))
#              axis(1, at=my.x, labels=cytoBand_chr$name[i], outer=TRUE)
##              lines(c(my.x, my.x), c(1.9, 2))
#            }
            par(ps=10, las=1)  # rotate text back to horizontal
          })


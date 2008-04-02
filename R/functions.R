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


.labelChromosome <- function(object, op, j){
  if(j == 1){
    if(op$label.chromosome)
      mtext(unique(chromosome(object)), side=3, outer=FALSE, line=op$line.label.chromosome, cex=op$cex.lab)
  }              
}

.isHomozygous <- function(object){
  calls(object) == 1 | calls(object) == 3
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








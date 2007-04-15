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

centromere <- function(annotation, chromosome){
  if(!all(as.character(chromosome) %in% rownames(annotation))){
     stop("Specified chromosome is not in the provided annotation")
   } else{
     annotation[as.character(chromosome), 1:2]
   }
}

chromosomeSize <- function(annotation, chromosome){
  if(!all(as.character(chromosome) %in% rownames(annotation))){
     stop("Specified chromosome is not in the provided annotation")
   } else{
     annotation[as.character(chromosome), 3]
   }
}

unsplitS4 <- function(value, featureData){
  ##Must be a better way
  obj <- new(class(value[[1]]),
             copyNumber=do.call("rbind", lapply(value, copyNumber)),
             cnConfidence=do.call("rbind", lapply(value, cnConfidence)),
             calls=do.call("rbind", lapply(value, calls)),
             callsConfidence=do.call("rbind", lapply(value, calls)),
             phenoData=phenoData(value[[1]]),
             annotation=annotation(value[[1]]),
             chromosomeAnnotation=chromosomeAnnotation(value[[1]]))
  featureData(obj) <- featureData[match(featureNames(obj), featureNames(featureData)), ]
  validObject(obj)
  obj
}

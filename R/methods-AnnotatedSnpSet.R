setAs("AnnotatedSnpSet", "AnnotatedSnpSetList",
      function(from){
        chr <- as.character(unique(chromosome(from)))
        snpSetList <- list()
        for(i in 1:length(chr)){
          snpSetList[[i]] <- from[chromosome(from) == chr[i], ]
        }
        new("AnnotatedSnpSetList", snpSetList=snpSetList)
      })

setMethod("chromosomeAnnotation", "AnnotatedSnpSet", function(object) object@chromosomeAnnotation)
setReplaceMethod("chromosomeAnnotation", c("AnnotatedSnpSet", "data.frame"),
                 function(object, value){
                   object@chromosomeAnnotation <- value
                   object
                 })

setMethod("alleleA", "AnnotatedSnpSet", function(object) alleleA(featureData(object)))
setMethod("alleleB", "AnnotatedSnpSet", function(object) alleleB(featureData(object)))
setMethod("chromosome", "AnnotatedSnpSet", function(object) chromosome(featureData(object)))
setMethod("dbSnpId", "AnnotatedSnpSet", function(object) dbSnpId(featureData(object)))
setMethod("enzyme", "AnnotatedSnpSet", function(object) enzyme(featureData(object)))
setMethod("position", "AnnotatedSnpSet", function(object) position(featureData(object)))
setMethod("probeSetId", "AnnotatedSnpSet", function(object) probeSetId(featureData(object)))

setMethod("initialize", "AnnotatedSnpSet",
          function(.Object,
                   assayData = assayDataNew(
                     calls=calls,
                     callsConfidence = callsConfidence,
                     copyNumber=copyNumber,
                     cnConfidence=cnConfidence, ...),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=character(),
                   chromosomeAnnotation=data.frame(),
                   calls=new("matrix"),
                   callsConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   copyNumber=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   cnConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   ...){
            data(chromosomeAnnotation)
            .Object@assayData <- assayData
            .Object@phenoData <- phenoData
            .Object@annotation <- annotation
            .Object@featureData <- featureData
            .Object@chromosomeAnnotation <- chromosomeAnnotation
            .Object@experimentData <- experimentData
            .Object
          })


setMethod("show", "AnnotatedSnpSet", function(object){
  tmp <- as(object, "eSet")
  show(tmp)

  cat("\nchromosomeAnnotation\n")
  N <- dim(chromosomeAnnotation(object))[1]
##  N <- as.character(unique(chromosome(object)))
##  if(length(N) > 2) N <- N[1:2]
  if(N > 2)  print(chromosomeAnnotation(object)[1:2, ])
  cat("...\n")
  print(chromosomeAnnotation(object)[N, ])
})

setMethod("summary", "AnnotatedSnpSet", function(object, digits=3, ...){
  
  ##Calculate mean,sd copy number and prop no calls, prop het calls
  ##for each chromosome in each sample.  Return an S x 23 matrix for
  ##each.
  chrList <- as.list(unique(chromosome(object)))
  if(sum(chromosome(object) != "chrX") > 0) chrList[[length(chrList) + 1]] <- "autosome"

  chrStats <- list()
  chromosomeMeanCn <- function(chrom, object){
    if(chrom != "autosome"){
      cn <- copyNumber(object)[chromosome(object) == chrom, ]      
    } else  cn <- copyNumber(object)[chromosome(object) != "chrX",]
    colMeans(as.matrix(cn))
  }
  chrStats[[1]] <- sapply(chrList, chromosomeMeanCn, object)

  ##standard deviation of copy number for each chromosome
  chromosomeSdCn <- function(chrom, object){
    if(chrom != "autosome"){
      cn <- copyNumber(object)[chromosome(object) == chrom, ]
    } else  cn <- copyNumber(object)[chromosome(object) != "chrX",]
    cn <- apply(as.matrix(cn), 2, stats::sd)
  }
  chrStats[[2]] <- sapply(chrList, chromosomeSdCn, object)

  ##Proportion of no calls for each chromosome
  propNoCalls <- function(chrom, object){
    if(chrom != "autosome"){
      calls <- calls(object)[chromosome(object) == chrom, ]
    } else  calls <- calls(object)[chromosome(object) != "chrX", ]
    calls <- as.matrix(calls)
    colMeans(calls == 4)    
  }
  chrStats[[3]] <- sapply(chrList, propNoCalls, object)

  ##Proportion of homozygous calls given that a call was made
  propHo <- function(chrom, calls, chromosomes){
    if(chrom != "autosome"){
      calls <- calls[chromosomes == chrom, ]
    } else calls <- calls[chromosomes != "chrX", ]

    perSample <- function(x){
      x <- x[x != 4]
      x[x == 2] <- 0
      x[x == 3] <- 1
      mean(x)
    }
    apply(as.matrix(calls), 2, perSample) 
  }
  chrStats[[4]] <- sapply(chrList, propHo, chromosomes=chromosome(object), calls=calls(object))
  
  colNames <- function(x, chrom){
    if(is.matrix(x)) colnames(x) <- chrom else{
      x <- t(as.matrix(x))
      colnames(x) <- chrom
    }
    x
  }
  chrStats <- lapply(chrStats, colNames, chrom=unlist(chrList))
  names(chrStats) <- c("avgCopyNumber", "sdCopyNumber", "propNoCalls", "propHo")

  ##Calculate grand average
  if(dim(object)[2] > 1){
    grand<- list()
    grand[[1]] <- colMeans(chrStats$avgCopyNumber)
    grand[[2]] <- apply(chrStats$sdCopyNumber, 2, stats::sd)
    grand[[3]] <- colMeans(chrStats$propNoCalls)
    grand[[4]] <- colMeans(chrStats$propHo)
    grand <- do.call("rbind", grand)
    rownames(grand) <- c("overall mean", "sd of means", "avg prop no calls", "avg prop AA/BB among calls")
  } else grand <- NULL

  stats <- list(chromosome=chrStats, overall=grand)
  return(stats)
})

            

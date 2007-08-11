setMethod("initialize", "oligoSnpSet",
          function(.Object,
                   assayData=assayDataNew(
                     calls=calls,
                     callsConfidence=callsConfidence,
                     copyNumber=copyNumber,
                     cnConfidence=cnConfidence, ...),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=character(),
                   calls=new("matrix"),
                   callsConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   copyNumber=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   cnConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)), ... ){
            .Object@assayData <- assayData
            .Object@phenoData <- phenoData
            .Object@annotation <- annotation
            .Object@featureData <- featureData
            .Object@experimentData <- experimentData
            .Object
          })

setValidity("oligoSnpSet", function(object) {
    msg <- validMsg(NULL, isValidVersion(object, "oligoSnpSet"))
    msg <- validMsg(msg, assayDataValidMembers(assayData(object), c("calls", "callsConfidence", "copyNumber", "cnConfidence")))
    if (is.null(msg)) TRUE else msg
})

setMethod(".getY", "oligoSnpSet", function(object, op, ...){
  y <- copyNumber(object)
  y[y < op$ylim[1]] <- op$ylim[1]
  y[y > op$ylim[2]] <- op$ylim[2]
  y
})

setMethod("show", "oligoSnpSet", function(object) {
  cat(class( object ), " (storageMode: ", storageMode(object), ")\n", sep="")
  adim <- dim(object)
  if (length(adim)>1)
    cat("assayData:",
        if (length(adim)>1)
        paste(adim[[1]], "features,",
              adim[[2]], "samples") else NULL,
        "\n")
  cat("  element names:",
      paste(assayDataElementNames(object), collapse=", "), "\n")
  cat("experimentData: use 'experimentData(object)'\n")
  pmids <- pubMedIds(object)
  if (length(pmids) > 0 && all(pmids != ""))
    cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
  cat("Annotation:", annotation(object), "\n")
  cat("phenoData\n")
  if(length(adim) > 1) show(phenoData(object))
  cat("featureData\n")      
  if(length(adim) > 1) show(featureData(object))
  cat("Annotation ")
  show(annotation(object))
})

setMethod("summary", "oligoSnpSet",
          function(object, digits=3, ...){
            ##Calculate mean,sd copy number and prop no calls, prop het calls
            ##for each chromosome in each sample.  Return an S x 23 matrix for
            ##each.
            x <- list()
            ##  chrList <- as.list(unique(chromosome(object)))
            ##  if(sum(chromosome(object) != "X") > 0) chrList[[length(chrList) + 1]] <- "autosome"
            cn <- split(copyNumber(object), chromosome(object))
            cn <- lapply(cn, matrix, ncol=ncol(object), byrow=TRUE)
            x[[1]] <- round(sapply(cn, colMeans), digits)
            ##  chrStats <- list()
            ##  chromosomeMeanCn <- function(chrom, object){
            ##    if(chrom != "autosome"){
            ##      cn <- copyNumber(object)[chromosome(object) == chrom, ]      
            ##    } else  cn <- copyNumber(object)[chromosome(object) != "X",]
            ##    colMeans(as.matrix(cn))
            ##  }
            ##  chrStats[[1]] <- sapply(chrList, chromosomeMeanCn, object)
            colSds <- function(x) apply(x, 2, "sd")
            x[[2]] <- round(sapply(cn, colSds), digits)
  
            ##standard deviation of copy number for each chromosome
            ##  chromosomeSdCn <- function(chrom, object){
            ##    if(chrom != "autosome"){
            ##      cn <- copyNumber(object)[chromosome(object) == chrom, ]
            ##    } else  cn <- copyNumber(object)[chromosome(object) != "X",]
            ##    cn <- apply(as.matrix(cn), 2, stats::sd)
            ##  }
            ##  chrStats[[2]] <- sapply(chrList, chromosomeSdCn, object)

            ##Proportion of no calls for each chromosome
            calls <- split(calls(object), chromosome(object))
            calls <- lapply(calls, matrix, ncol=ncol(object), byrow=TRUE)
            calls.missing <- lapply(calls, "==", 4)
            x[[3]] <- round(sapply(calls.missing, colMeans), digits)
                              
            ##  propNoCalls <- function(chrom, object){
            ##    if(chrom != "autosome"){
            ##      calls <- calls(object)[chromosome(object) == chrom, ]
            ##    } else  calls <- calls(object)[chromosome(object) != "X", ]
            ##    calls <- as.matrix(calls)
            ##    colMeans(calls == 4)    
            ##  }
            ##  chrStats[[3]] <- sapply(chrList, propNoCalls, object)

            ##Proportion of homozygous calls given that a call was made
            homozygousCall <- function(x) x == 1 | x == 3
            calls.homozygous <- lapply(calls, homozygousCall)
            x[[4]] <- round(sapply(calls.homozygous, colMeans), digits)
            x[[5]] <- round(1-x[[4]]-x[[3]], digits)
            names(x) <- c("avgCN", "sdCN", "%NoCalls", "%Hom", "%Het")

            overall <- lapply(x, colMeans)
            x <- mapply("rbind", x, overall, SIMPLIFY=FALSE)

            labelRows <- function(x, sn){
              rownames(x) <- c(sn, "overall")
              return(x)
            }
            x <- lapply(x, labelRows, sn=sampleNames(object))
            
            reorderX <- function(x){
              y <- c(1:22, "X", "Y")
              y <- y[y %in% colnames(x)]
              x <- x[, y]
            }
            xx <- lapply(x, reorderX)
  
            ##Calculate grand average
            ##  if(dim(object)[2] > 1){
            ##    grand<- list()
            ##    grand[[1]] <- colMeans(chrStats$avgCopyNumber)
            ##    grand[[2]] <- apply(chrStats$sdCopyNumber, 2, stats::sd)
            ##    grand[[3]] <- colMeans(chrStats$propNoCalls)
            ##    grand[[4]] <- colMeans(chrStats$propHo)
            ##    grand <- do.call("rbind", grand)
            ##    rownames(grand) <- c("overall mean", "sd of means", "avg prop no calls", "avg prop AA/BB among calls")
            ##  } else grand <- NULL
            ##  stats <- list(chromosome=chrStats, overall=grand)
            return(xx)
          })


#create object of class snpscan with smoothed copynumbers and smoothed loh calls
#LOH calls are smoothed by setting homozygous = 0 and heterozygous = 1
setMethod("smoothSnp", "oligoSnpSet",
          function(object, 
                   span=1/10,
                   method="loess",
                   imputeNoCalls=TRUE,
                   verbose=TRUE){
            ##convert homozygous to 0 and heterozygous to 1
            calls(object)[calls(object) == 1 | calls(object) == 3] <- 0
            calls(object)[calls(object) == 2] <- 1
            smoothChromosome <- function(obj, span){
              loessX <- function(X, location, span){
                fit <- loess(X ~ location, span = span)$fitted
                return(fit)
              }
              ##Order by physical position before smoothing
              obj <- obj[order(position(obj)), ]
              cn.smooth <- apply(copyNumber(obj), 2, loessX, position(obj), span=span)
              call.smooth <- apply(calls(obj), 2, loessX, location=position(obj), span=span)
              rownames(cn.smooth) <- rownames(call.smooth) <- featureNames(obj) 
              copyNumber(obj) <- cn.smooth
              calls(obj) <- call.smooth
              obj
            }
            object.list <- split(object, chromosome(object))
            smooth.list <- lapply(object.list, smoothChromosome, span=span)
            smooth.set <- unsplitS4(smooth.list, featureData(object))
            return(smooth.set)
          })




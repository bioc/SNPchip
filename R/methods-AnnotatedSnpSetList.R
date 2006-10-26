#setMethod("convert", "AnnotatedSnpSet",
#      function(from){
#        chr <- unique(chromosome(from))
#        snpSetList <- list()
#        for(i in 1:length(chr))
#          {
#            snpSetList[[i]] <- from[chromosome(from) == chr[i], ]
#          }
#        to <- new("AnnotatedSnpSetList", snpSetList = snpSetList)
#      })
setMethod("initialize", "AnnotatedSnpSetList",
          function(.Object,
                   snpSetList = list()){
            .Object@snpSetList <- snpSetList
            .Object
          })

setAs("AnnotatedSnpSetList", "AnnotatedSnpSet",
      function(from){
        chrAnn <- chromosomeAnnotation(snpSetList(from)[[1]])

        ##Combine assayData
        ad <- lapply(snpSetList(from), assayData)
        f.calls <- function(x) x$calls
        f.cn <- function(x) x$copyNumber
        f.callsConf <- function(x) x$callsConfidence
        f.cnConf <- function(x) x$cnConfidence
        calls.l <- lapply(ad, f.calls)
        cn.l <- lapply(ad, f.cn)
        callsConf.l <- lapply(ad, f.callsConf)
        cnConf.l <- lapply(ad, f.cnConf)
        calls <- do.call("rbind", calls.l)
        cn <- do.call("rbind", cn.l)
        callsConf <- do.call("rbind", callsConf.l)
        cnConf <- do.call("rbind", cnConf.l)

        ##Combine featureData
        fd <- lapply(snpSetList(from), featureData)
        pd <- lapply(fd, pData)
        fd.combined <- do.call("rbind", pd)
        fd <- new("AnnotatedDataFrame", data = fd.combined,
                  varMetadata = varMetadata(featureData(snpSetList(from)[[1]])))
        obj <- new("oligoSnpSet",
                   calls = calls,
                   copyNumber = cn,
                   callsConfidence = callsConf,
                   cnConfidence = cnConf,
                   featureData = fd,
                   phenoData = phenoData(snpSetList(from)[[1]]),
                   experimentData = experimentData(snpSetList(from)[[1]]))
        obj <- as(obj, "AnnotatedSnpSet")
        chromosomeAnnotation(obj) = chrAnn
        obj
      })

setMethod("snpSetList", "AnnotatedSnpSetList", function(object) object@snpSetList)



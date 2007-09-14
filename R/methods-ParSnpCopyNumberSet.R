setMethod("initialize", "ParSnpCopyNumberSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$ylab <- "copy number"
            .Object
          })

setMethod("plotSnp", c("ParSnpCopyNumberSet", "SnpCopyNumberSet"),
          function(object, snpset){
            old.par <- par(no.readonly=TRUE)
            on.exit(old.par)
            callNextMethod()
          })

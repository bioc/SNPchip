setMethod("initialize", "ParSnpCopyNumberSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$ylab <- "copy number"
            .Object
          })



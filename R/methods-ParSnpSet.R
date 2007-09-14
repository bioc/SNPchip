setMethod("initialize", "ParSnpSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$col <- c("lightblue", "red", "lightblue")
            .Object$bg <- c("lightblue", "red", "lightblue")            
            .Object$ylab <- "copy number"
            .Object
          })

setMethod("plotSnp", c("ParSnpSet", "oligoSnpSet"),
          function(object, snpset){
            old.par <- par(no.readonly=TRUE)
            on.exit(old.par)            
            callNextMethod()
          })

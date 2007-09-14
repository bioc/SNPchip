setMethod("initialize", "ParSnpCallSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$col <- c("lightblue", "red", "lightblue")
            .Object$ylab <- "genotype call"
            .Object
          })

setMethod("plotSnp", c("ParSnpCallSet", "SnpCallSet"),
          function(object, snpset){
            old.par <- par(no.readonly=TRUE)
            on.exit(old.par)            
            callNextMethod()
          })

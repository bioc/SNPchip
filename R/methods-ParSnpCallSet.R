setMethod("initialize", "ParSnpCallSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$col <- c("lightblue", "red", "lightblue")
            .Object$ylab <- "genotype call"
            .Object
          })

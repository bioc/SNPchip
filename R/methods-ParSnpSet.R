setMethod("initialize", "ParSnpSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$col <- c("lightblue", "red", "lightblue")
            .Object$bg <- c("lightblue", "red", "lightblue")            
            .Object$ylab <- "copy number"
            .Object
          })



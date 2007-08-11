setMethod(".combineChips", c("SnpCallSet", "SnpCallSet"),
          function(x, y, ...){
            fData <- rbind(fData(x), fData(y))
            featureData <- new("AnnotatedDataFrame",
                               data=fData,
                               varMetadata=fvarMetadata(x))
            warning("only using phenoData in first argument")
            new("SnpCallSet",
                featureData=featureData,
                phenoData=phenoData(x),
                calls=rbind(calls(x), calls(y)),
                callsConfidence=rbind(callsConfidence(x), callsConfidence(y)),
                experimentData=experimentData(x),
                annotation="mapping100k")
          })
setMethod("summary", signature(object = "SnpCallSet"),
          function(object, digits = 3, noCalls = FALSE, ...){
            het <- colMeans(ifelse(calls(object) == 2, 1, 0))
            hom <- colMeans(ifelse(calls(object) == 1 | calls(object) == 3, 1, 0))
            mat <- rbind(het, hom)
            mat
          })

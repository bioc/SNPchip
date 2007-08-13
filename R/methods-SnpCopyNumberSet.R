##setMethod(".getY", "SnpCopyNumberSet", function(object, op, ...){
##  y <- copyNumber(object)
##  y[y < op$ylim[1]] <- op$ylim[1]
##  y[y > op$ylim[2]] <- op$ylim[2]
##  y
##})

##setMethod("unsplitS4", c("SnpCopyNumberSet", "AnnotatedDataFrame"),
##          function(object, featureData){
##            obj <- new(class(object[[1]]),
##                       copyNumber=do.call("rbind", lapply(object, copyNumber)),
##                       cnConfidence=do.call("rbind", lapply(object, cnConfidence)),
##                       phenoData=phenoData(object[[1]]),
##                       annotation=annotation(object[[1]]),
##                       experimentData=experimentData(object[[1]]))
##            featureData(obj) <- featureData[match(featureNames(obj), featureNames(featureData)), ]
##            stopifnot(identical(rownames(copyNumber(obj)), rownames(fData(obj))))
##            obj
##          })


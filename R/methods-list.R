setMethod("unsplitSnpSet", c("list", "AnnotatedDataFrame"),
          function(from, annotatedDataFrame, ...){
            object <- switch(class(from[[1]]),
                             oligoSnpSet=new("oligoSnpSet", ...),
                             SnpCallSet=new("SnpCallSet", ...),
                             SnpCopyNumberSet=new("SnpCopyNumberSet", ...),
                             NULL)
            featureData(object) <- annotatedDataFrame[match(featureNames(object), featureNames(annotatedDataFrame)), ]
            stopifnot(identical(rownames(copyNumber(object)), rownames(fData(object))))
            object
          })

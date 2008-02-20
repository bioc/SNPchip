setMethod("chromosome", "SnpCallSet",
          function(object){
            if(!("chromosome" %in% fvarLabels(object))){

              fs <- featureNames(object)
              sql <- "SELECT man_fsetid, chrom FROM featureSet WHERE man_fsetid LIKE 'SNP%'"              
              ##Check if two objects have been combined
              pkgs <- strsplit(annotation(object), ",")[[1]]
              if(length(pkgs) > 1){
                object2 <- object1 <- object
                annotation(object1) <- pkgs[1]
                annotation(object2) <- pkgs[2]

                tmp1 <- dbGetQuery(db(object1), sql)
                tmp2 <- dbGetQuery(db(object2), sql)
                tmp <- rbind(tmp1, tmp2)
              } else {
                tmp <- dbGetQuery(db(object), sql)
              }
              idx <- match(fs, tmp[["man_fsetid"]])
              chr <- tmp[idx, "chrom"]
            } else {
              chr <- featureData(object)$chromosome
            }
            return(chr)
          })

setMethod("position", "SnpCallSet",
          function(object){
            if(!("position" %in% fvarLabels(object))){

              fs <- featureNames(object)
              sql <- "SELECT man_fsetid, physical_pos FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
              pkgs <- strsplit(annotation(object), ",")[[1]]
              if(length(pkgs) > 1){
                object2 <- object1 <- object
                annotation(object1) <- pkgs[1]
                annotation(object2) <- pkgs[2]                
              
                tmp1 <- dbGetQuery(db(object1), sql)
                tmp2 <- dbGetQuery(db(object2), sql)
                tmp <- rbind(tmp1, tmp2)
              } else{
                tmp <- dbGetQuery(db(object), sql)
              }
              idx <- match(fs, tmp[["man_fsetid"]])
              pos <- tmp[idx, "physical_pos"]
            } else {
              pos <- featureData(object)$position
            }
            return(pos)
          })





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

##setMethod("getSnpAnnotation", "SnpCallSet",
##          function(object){
##            callNextMethod()
##          })



setMethod("summary", signature(object = "SnpCallSet"),
          function(object, digits = 3, noCalls = FALSE, ...){
            het <- colMeans(ifelse(calls(object) == 2, 1, 0))
            hom <- colMeans(ifelse(calls(object) == 1 | calls(object) == 3, 1, 0))
            mat <- rbind(het, hom)
            mat
          })

##setMethod("unsplitS4", c("SnpCallSet", "AnnotatedDataFrame"),
##          function(object, featureData){
##            obj <- new(class(object[[1]]),
##                       calls=do.call("rbind", lapply(object, calls)),
##                       callsConfidence=do.call("rbind", lapply(object, calls)),
##                       experimentData=experimentData(object[[1]]),
##                       phenoData=phenoData(object[[1]]),
##                       annotation=annotation(object[[1]]))
##            featureData(obj) <- featureData[match(featureNames(obj), featureNames(featureData)), ]
##            stopifnot(identical(rownames(copyNumber(obj)), rownames(fData(obj))))
##            obj
##          })
